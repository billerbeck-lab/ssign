#!/usr/bin/env python3
"""Run BLASTp locally against a user-provided database."""

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import TOOL_TIMEOUT_S
from ssign_lib.fasta_io import read_fasta
from ssign_lib.substrates import load_substrate_ids
from dedup_sequences import deduplicate_dict, expand_results_dict

# Terms to exclude from BLASTp hits
EXCLUDE_TERMS = [
    "hypothetical protein",
    "uncharacterized protein",
    "domain of unknown function",
    "duf",
    "unnamed protein product",
    "predicted protein",
]

# BLAST outfmt 6 column indices, mirroring the order in BLAST_OUTFMT below.
# Centralised so the parser doesn't drift if the format string changes.
BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen"
_COL_QSEQID = 0
_COL_SSEQID = 1
_COL_PIDENT = 2
_COL_ALN_LEN = 3
_COL_EVALUE = 10
_COL_BITSCORE = 11
_COL_STITLE = 12
_COL_QLEN = 13
_BLAST_MIN_FIELDS = 15  # all 15 outfmt fields must be present


def run_local_blastp(query_fasta, db_path, evalue, exclude_taxid, num_threads=4):
    """Run local BLAST+ blastp and return parsed hits."""
    cmd = [
        "blastp",
        "-query",
        query_fasta,
        "-db",
        db_path,
        "-outfmt",
        BLAST_OUTFMT,
        "-evalue",
        str(evalue),
        # -max_target_seqs is the documented BLAST+ footgun (Shah 2019,
        # bioinformatics 35:1613): does NOT return top-N best hits — returns
        # the first N to cross the e-value threshold during DB traversal.
        # 500 is the BLAST web default; widening the candidate pool plus
        # explicit bitscore-sort below recovers the actual best hit.
        "-max_target_seqs",
        "500",
        "-num_threads",
        str(num_threads),
    ]
    if exclude_taxid:
        cmd.extend(["-negative_taxids", str(exclude_taxid)])

    logger.info(f"Running local BLASTp against {db_path} ({num_threads} threads)")
    # FRAGILE: subprocess call requires BLAST+ (blastp) on PATH
    # If this breaks: install BLAST+ (sudo apt install ncbi-blast+)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TOOL_TIMEOUT_S)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"BLAST+ blastp binary not found: {e}\n"
            f"  Common causes:\n"
            f"    - NCBI BLAST+ is not installed or not on PATH\n"
            f"  How to fix:\n"
            f"    - Debian/Ubuntu: sudo apt install ncbi-blast+\n"
            f"    - macOS:         brew install blast\n"
            f"    - Conda:         conda install -c bioconda blast"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"BLASTp timed out after 4 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input sequences\n"
            f"    - Use a smaller database (e.g. Swiss-Prot instead of NR)"
        ) from e

    if result.returncode != 0:
        logger.error(f"BLASTp failed: {result.stderr[:500]}")
        raise RuntimeError(f"BLASTp exit code {result.returncode}")

    return parse_blast_tabular(result.stdout)


def parse_blast_tabular(output_text):
    """Parse BLAST tabular output (outfmt 6) and return best hit per query.

    Sorts every query's hits by bitscore (descending) and keeps the top
    one. Don't trust BLAST's row ordering: with -max_target_seqs widened
    to 500 the input contains many candidates per query, and BLAST does
    NOT guarantee they're sorted by score (the row order reflects DB
    traversal order). Picking by bitscore here ensures we keep the
    actually-best hit.
    """
    rows_by_query: dict[str, list[list[str]]] = {}
    for line in output_text.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < _BLAST_MIN_FIELDS:
            continue
        rows_by_query.setdefault(parts[_COL_QSEQID], []).append(parts)

    hits = {}
    for query_id, rows in rows_by_query.items():
        rows.sort(key=lambda r: float(r[_COL_BITSCORE]), reverse=True)
        best = rows[0]

        pident = float(best[_COL_PIDENT])
        qlen = int(best[_COL_QLEN])
        aln_len = int(best[_COL_ALN_LEN])
        qcov = (aln_len / qlen * 100) if qlen > 0 else 0

        # NCBI concatenates redundant subject titles with " >" in stitle —
        # take only the first segment so the description filter doesn't
        # see "hypothetical protein >RecName: hemolysin..." as hypothetical.
        primary_desc = best[_COL_STITLE].split(" >")[0]

        hits[query_id] = {
            "locus_tag": query_id,
            "blastp_hit_accession": best[_COL_SSEQID],
            "blastp_hit_description": primary_desc[:200],
            "blastp_pident": round(pident, 1),
            "blastp_qcov": round(qcov, 1),
            "blastp_evalue": float(best[_COL_EVALUE]),
            "blastp_bitscore": float(best[_COL_BITSCORE]),
        }

    return hits


def filter_hits(hits, min_pident, min_qcov):
    """Filter BLASTp hits by identity, coverage, and description."""
    filtered = {}
    for pid, hit in hits.items():
        if hit["blastp_pident"] < min_pident:
            continue
        if hit["blastp_qcov"] < min_qcov:
            continue

        desc_lower = hit["blastp_hit_description"].lower()
        if any(term in desc_lower for term in EXCLUDE_TERMS):
            continue

        filtered[pid] = hit

    return filtered


def main():
    parser = argparse.ArgumentParser(description="Run BLASTp locally")
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--db", required=True, help="Path to local BLAST+ database (e.g. NR)")
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min-pident", type=float, default=80)
    parser.add_argument("--min-qcov", type=float, default=80)
    parser.add_argument("--exclude-taxid", default="")
    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 4,
        help="Threads for blastp -num_threads (default: all CPUs)",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Extract substrate sequences into temp FASTA
    substrate_ids = load_substrate_ids(args.substrates)
    all_seqs = read_fasta(args.proteins)
    sub_seqs = {k: v for k, v in all_seqs.items() if k in substrate_ids}

    # Deduplicate to avoid redundant BLAST work
    unique_seqs, seq_groups = deduplicate_dict(sub_seqs)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        for pid, seq in unique_seqs.items():
            tmp.write(f">{pid}\n{seq}\n")
        tmp_path = tmp.name

    try:
        hits = run_local_blastp(tmp_path, args.db, args.evalue, args.exclude_taxid, args.threads)
    finally:
        os.unlink(tmp_path)

    filtered = filter_hits(hits, args.min_pident, args.min_qcov)
    # Expand results back to all duplicate proteins
    filtered = expand_results_dict(filtered, seq_groups)
    logger.info(f"{len(filtered)}/{len(hits)} hits pass filters for {args.sample}")

    fieldnames = [
        "locus_tag",
        "blastp_hit_accession",
        "blastp_hit_description",
        "blastp_pident",
        "blastp_qcov",
        "blastp_evalue",
        "blastp_bitscore",
    ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for hit in filtered.values():
            writer.writerow(hit)


if __name__ == "__main__":
    main()
