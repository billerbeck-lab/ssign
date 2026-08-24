#!/usr/bin/env python3
"""Cross-reference the ssign benchmark proteins against a tool's training set.

Three independent match routes, because training sets use different identifier
namespaces and sometimes carry a homolog rather than the exact protein:

  1. identifier  - benchmark accession universe (UniProt primary + secondary,
                   RefSeq, EMBL protein_id) intersected with the training-set IDs
  2. sequence    - identical amino-acid sequence, or the training record being a
                   prefix of the benchmark sequence (SignalP truncates to 70 aa)
  3. homology    - Smith-Waterman (pyopal, BLOSUM62, gap 11/1) of every benchmark
                   against every training record; best hit reported

Usage:
  match_training_set.py --db <training.fasta> --tool <name> --out <out.tsv>
"""

import argparse
import json
import os
import re
import sys
from collections import Counter

import pyopal
from _paths import (
    CLOSE_HOMOLOG,
    DISTANT_HOMOLOG,
    IN_TRAINING,
    NEAR_IDENTICAL_HOMOLOG,
    NO_MEANINGFUL_MATCH,
    WORK,
)

UNIPROT_RE = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
ACCLIKE_RE = re.compile(
    r"\b((?:WP|NP|YP|XP|AP|EAA|CAA|CAB|CAD|CAE|BAA|BAB|AAA|AAB|AAC|AAD|AAF|AAG|AAK|AAL|AAM|AAN|AAO|AAP|ABC|ACD|EDL|EFA)_?[0-9]{5,9}(?:\.[0-9]+)?)\b"
)


# SignalP distributes its datasets in a 3-line format: header, sequence, then a
# per-residue annotation string of equal length over this alphabet (S/T/L/P mark
# signal-peptide types, I cytoplasm, M transmembrane, O extracellular). Treating
# it as ordinary FASTA silently appends 70 junk characters to every sequence.
SP_ANNOTATION_ALPHABET = set("STLPIMO")

# A training record shorter than this is never treated as a truncated copy of a
# query; see the prefix_lens comment below.
MIN_PREFIX_LEN = 50


def _is_annotation(prev, line):
    return (
        prev is not None
        and len(line) == len(prev)
        and set(line) <= SP_ANNOTATION_ALPHABET
        and not set(prev) <= SP_ANNOTATION_ALPHABET
    )


def read_fasta(path):
    # Deliberately duplicated from ssign_lib.fasta_io: these audit scripts must run
    # standalone, without the ssign package installed.
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r").strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            elif line:
                if len(buf) == 1 and _is_annotation(buf[0], line):
                    continue
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def header_ids(header):
    """Every identifier-looking token in a FASTA header."""
    ids = set()
    first = header.split()[0] if header.split() else header
    for tok in re.split(r"[|\s,;]+", header):
        tok = tok.strip()
        if tok:
            ids.add(tok)
            ids.add(tok.split(".")[0])
    ids.add(first)
    ids.add(first.split(".")[0])
    for m in UNIPROT_RE.finditer(header):
        ids.add(m.group(1))
    for m in ACCLIKE_RE.finditer(header):
        ids.add(m.group(1))
        ids.add(m.group(1).split(".")[0])
    return {i for i in ids if len(i) >= 4}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--topn",
        type=int,
        default=150,
        help="candidates per ranking re-aligned in full mode (two rankings are unioned)",
    )
    args = ap.parse_args()

    bench = json.load(open(os.path.join(WORK, "benchmark_uniprot.json")))
    bench = [b for b in bench if b.get("sequence")]

    # ---- benchmark identifier universe, per protein
    bench_ids = []
    for b in bench:
        s = set()
        if b.get("primaryAccession"):
            s.add(b["primaryAccession"])
            s.add(b["primaryAccession"].split(".")[0])
        for k in ("secondaryAccessions", "xref_refseq", "embl_protein_ids"):
            for v in b.get(k) or []:
                s.add(v)
                s.add(v.split(".")[0])
        bench_ids.append(s)

    # ---- load training db
    db_names, db_seqs = [], []
    for h, s in read_fasta(args.db):
        db_names.append(h)
        db_seqs.append(re.sub(r"[^A-Z]", "", s.upper()))
    if not db_seqs:
        sys.exit(f"no sequences parsed from {args.db}")
    print(f"[{args.tool}] training records: {len(db_seqs)}", file=sys.stderr)

    db_id_index = {}
    for i, h in enumerate(db_names):
        for t in header_ids(h):
            db_id_index.setdefault(t, []).append(i)
    seq_index = {}
    for i, s in enumerate(db_seqs):
        seq_index.setdefault(s, []).append(i)
    # Lengths of the truncated records, for the truncated-training-set case:
    # SignalP keeps only the first 70 aa, so a full-length query matches by
    # prefix. The lower bound matters: without it any short peptide in the
    # database (DeepLocPro has records down to 8 residues) would match some
    # query's N-terminus by chance and be reported as a training hit.
    prefix_lens = sorted({len(s) for s in db_seqs if MIN_PREFIX_LEN <= len(s) <= 100})

    # pyopal defaults to BLOSUM50 with gap_open=3, which is ~3x cheaper than any
    # standard protein setting and lets local alignments stitch through unrelated
    # segments. Pin BLAST's protein defaults so the homolog tiers mean what they say.
    aligner = pyopal.Aligner("BLOSUM62", gap_open=11, gap_extend=1)
    # pyopal needs sequences over its alphabet; map anything odd to X
    ALPHA = set("ARNDCQEGHILKMFPSTWYVBZX*")
    clean = ["".join(c if c in ALPHA else "X" for c in s) for s in db_seqs]
    database = pyopal.Database(clean)

    rows = []
    for qi, b in enumerate(bench):
        qseq = "".join(c if c in ALPHA else "X" for c in b["sequence"].upper())

        # 1. identifier match
        id_hits = []
        for t in bench_ids[qi]:
            for j in db_id_index.get(t, []):
                id_hits.append((t, j))

        # 2. exact sequence / prefix
        seq_hits = list(seq_index.get(qseq, []))
        pref_hits = []
        if not seq_hits:
            for L in prefix_lens:
                cand = qseq[:L]
                for j in seq_index.get(cand, []):
                    pref_hits.append(j)

        # 3. homology: cheap SIMD score pass over the whole db, then one batched
        #    full-mode pass over the candidates to recover identity + coverage.
        #
        #    The candidate set is the union of two rankings, which is the whole
        #    point. Raw Smith-Waterman score favours long targets, but the best
        #    hit is chosen on `mfrac`, which is normalised by the SHORTER
        #    sequence. Those orders disagree badly on a length-heterogeneous
        #    database, so ranking by raw score alone pushes the true best hit
        #    out of the candidate set for most queries. Score divided by the
        #    shorter length is a close proxy for mfrac; taking both covers both.
        scores = [s.score for s in aligner.align(qseq, database, mode="score")]
        by_score = sorted(range(len(scores)), key=lambda k: scores[k], reverse=True)
        by_norm = sorted(
            range(len(scores)), key=lambda k: scores[k] / max(1, min(len(qseq), len(clean[k]))), reverse=True
        )
        order = list(dict.fromkeys(by_score[: args.topn] + by_norm[: args.topn]))
        subdb = pyopal.Database([clean[j] for j in order])
        full = aligner.align(qseq, subdb, mode="full")
        best = None
        for rank, r in enumerate(full):
            j = order[rank]
            # pyopal alignment string: M = match, X = mismatch, I/D = gaps.
            aln = r.alignment
            n_match = aln.count("M")
            ident = r.identity()  # matches / aligned columns, gaps excluded
            alen = len(aln)
            qcov = r.coverage("query")
            tcov = r.coverage("target")
            # Span coverage of whichever sequence is SHORTER. Training sets are
            # sometimes truncated (SignalP keeps only the first 70 aa), so query
            # coverage alone would understate a genuine full-length match.
            cov = qcov if len(qseq) <= len(clean[j]) else tcov
            # Gap-aware: fraction of the shorter sequence that is an identical
            # residue. Span coverage alone overstates gappy alignments between
            # long repeat-rich proteins, which is exactly what the T1SS RTX
            # toxins and T5bSS exoproteins in this list are.
            mfrac = n_match / max(1, min(len(qseq), len(clean[j])))
            cand = (ident, qcov, tcov, alen, j, r.score, cov, mfrac, n_match)
            if best is None or mfrac > best[7]:
                best = cand

        ident, qcov, tcov, alen, bj, bscore, cov, mfrac, n_match = best
        matched_by, db_rec, db_hdr = [], "", ""
        if id_hits:
            matched_by.append("identifier")
            db_rec = id_hits[0][0]
            db_hdr = db_names[id_hits[0][1]]
        if seq_hits:
            matched_by.append("exact_sequence")
            if not db_hdr:
                db_hdr = db_names[seq_hits[0]]
        if pref_hits:
            matched_by.append("sequence_prefix")
            if not db_hdr:
                db_hdr = db_names[pref_hits[0]]

        if matched_by:
            verdict = IN_TRAINING
        elif ident >= 0.90 and mfrac >= 0.80:
            verdict = NEAR_IDENTICAL_HOMOLOG
        elif ident >= 0.50 and mfrac >= 0.50:
            verdict = CLOSE_HOMOLOG
        elif ident >= 0.30 and mfrac >= 0.30:
            verdict = DISTANT_HOMOLOG
        else:
            verdict = NO_MEANINGFUL_MATCH

        rows.append(
            {
                "instance_id": b["instance_id"],
                "ss_type": b["ss_type"],
                "subtype": b["subtype"],
                "gene": b["gene"],
                "uniprot": b.get("primaryAccession", ""),
                "organism": (b.get("organism") or b.get("organism_csv") or "")[:70],
                "tool": args.tool,
                "verdict": verdict,
                "matched_by": ";".join(matched_by) or "-",
                "matched_identifier": db_rec or "-",
                "matched_record": (db_hdr or db_names[bj])[:120],
                "best_hit_header": db_names[bj][:120],
                "best_pct_identity": round(ident * 100, 1),
                "best_query_coverage": round(qcov * 100, 1),
                "best_target_coverage": round(tcov * 100, 1),
                "best_shorter_seq_coverage": round(cov * 100, 1),
                "best_identical_residue_fraction": round(mfrac * 100, 1),
                "best_identical_residues": n_match,
                "best_alignment_columns": alen,
                "best_sw_score": bscore,
            }
        )
        print(
            f"  {b['instance_id']:<10} {b['gene'][:18]:<18} {verdict:<24} "
            f"id={ident * 100:5.1f}% cov={cov * 100:5.1f}%  {db_names[bj][:55]}",
            file=sys.stderr,
        )

    if not rows:
        sys.exit(f"[{args.tool}] no benchmark proteins loaded; nothing written")
    cols = list(rows[0].keys())
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    print(f"\n[{args.tool}] verdicts: {dict(Counter(r['verdict'] for r in rows))}", file=sys.stderr)
    print(f"[{args.tool}] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
