#!/usr/bin/env python3
"""Phase A: build the trusted 1-substrate-per-system GOLD LIST (the paper reference set).

One representative secreted protein per secretion-system instance, fully self-contained (host, gene,
accession, genome coordinates, the machinery it sits near, citation + DOI). This is the list we trust,
reference, and re-verify; it replaces the sprawling 337-effector table for downstream use.

Construction
------------
Machinery instances (90, from the systems audit): keep the verified substrate (from positives_all,
linked by sys_instance_id == instance_id) CLOSEST in gene-order to that instance's machinery -- the same
+/-N-gene proximity ssign uses, so the kept substrate is the one ssign is most likely to recover. The
26 instances with no verified substrate are DROPPED (untrusted). Cross-replicon instances (substrate on
a different molecule than the machinery, e.g. V. cholerae chr2 substrates vs chr1 Eps, Ralstonia
chromosomal Rip vs megaplasmid Hrp) are KEPT but flagged unreachable. An instance whose substrate
cannot be positioned on the genome at all is dropped.

T5SS (self-secreting; not in the machinery key) is added directly from the verified T5SS effectors:
each is its own instance, substrate = itself (T5a/T5c secrete their own passenger; T5b's secreted
protein is the TpsA effector, which is the corpus entry). reachable by construction (ssign emits T5SS
on MacSyFinder detection of the autotransporter/TpsB).

Distance metric uses the alias-normalised phase-1 gene-order index so locus-tag drift between tables is
moot.

Inputs : data/machinery/{instances,machinery_resolved}.tsv, data/dataset/positives_all.tsv,
         data/phase1/gene_order_index.tsv
Outputs: data/phase2/verification_phase_a/gold_list.tsv         (the trusted list)
         data/phase2/verification_phase_a/gold_list_dropped.tsv (instances dropped + why)
Run    : .venv/bin/python scripts/63_build_gold_list.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VDIR = BENCH / "data" / "phase2" / "verification_phase_a"
INF = 10**9

FIELDS = [
    "instance_id",
    "ss_type",
    "subtype",
    "organism",
    "gene",
    "uniprot",
    "effector_locus",
    "genome",
    "contig",
    "start",
    "stop",
    "strand",
    "stage_replicons",
    "nearest_machinery_gene",
    "nearest_machinery_locus",
    "distance_to_machinery_genes",
    "reachable_within_3",
    "distance_method",
    "n_candidate_substrates",
    "evidence_tier",
    "verification_status",
    "primary_ref",
    "citation_quote",
]


def _norm(lt: str) -> str:
    return (lt or "").replace("_", "").replace(" ", "").upper()


# PAO1 is staged under both its INSDC and RefSeq accessions in the gene-order index (the one known
# duplicate molecule in the panel), so machinery loci resolve to both; canonicalize to the RefSeq form
# so a single chromosome is never listed as two co-staging replicons. accession_base can't bridge these
# (different prefixes), hence the explicit map.
REPLICON_ALIAS = {"AE004091": "NC_002516.2"}


def _canon(acc: str) -> str:
    return REPLICON_ALIAS.get(acc, REPLICON_ALIAS.get(acc.split(".")[0], acc))


def gene_order():
    """normalized locus -> [(record_acc, ordinal, start, end, strand)] incl. aliases (old tags)."""
    idx = defaultdict(list)
    for r in read_tsv(BENCH / "data" / "phase1" / "gene_order_index.tsv"):
        entry = (
            r["record_acc"],
            int(r["ordinal"]),
            int(r["start"]) if r.get("start") else None,
            int(r["end"]) if r.get("end") else None,
            r.get("strand", ""),
        )
        for k in {_norm(r["locus_tag"])} | {_norm(a) for a in (r.get("aliases", "") or "").split(";") if a.strip()}:
            idx[k].append(entry)
    return idx


def main() -> int:
    inst_meta = {r["instance_id"]: r for r in read_tsv(BENCH / "data" / "machinery" / "instances.tsv")}
    pos = read_tsv(BENCH / "data" / "dataset" / "positives_all.tsv")
    mach = read_tsv(BENCH / "data" / "machinery" / "machinery_resolved.tsv")
    gidx = gene_order()

    mach_loci = defaultdict(list)  # instance -> [normalized machinery locus]
    for m in mach:
        if m.get("locus_tag", "").strip():
            mach_loci[m["instance_id"]].append(_norm(m["locus_tag"]))
    # instance -> {replicon carrying machinery}, taken from the INDEX (same accession namespace as the
    # substrate contig) so INSDC/RefSeq aliases of one chromosome (PAO1 AE004091 vs NC_002516) collapse
    # to a single replicon instead of looking like two.
    mach_acc = defaultdict(set)
    for inst_id, loci in mach_loci.items():
        for mlt in loci:
            for rec, *_ in gidx.get(mlt, []):
                mach_acc[inst_id].add(rec)

    sub_by_inst = defaultdict(list)
    for r in pos:
        sid = r.get("sys_instance_id", "").strip()
        if sid in inst_meta:
            sub_by_inst[sid].append(r)

    def locus_of(r):
        return (r.get("placement_effector_locus") or r.get("locus_tag") or "").strip()

    def position(r):
        """Best (record, ordinal, start, end, strand) for a substrate, or None if unplaceable."""
        for rec, ordn, s, e, strand in gidx.get(_norm(locus_of(r)), []):
            return (rec, ordn, s, e, strand)
        # fallback to the corpus's own coords (T5SS / coord-only rows)
        if r.get("contig") and r.get("start") and r.get("stop"):
            return (r["contig"], None, int(r["start"]), int(r["stop"]), r.get("strand", ""))
        if r.get("refseq_genome") and r.get("placement_start") and r.get("placement_stop"):
            return (
                r["refseq_genome"],
                None,
                int(r["placement_start"]),
                int(r["placement_stop"]),
                r.get("placement_strand", ""),
            )
        return None

    def nearest_machinery(inst_id, srec, sordn):
        """(min same-replicon gene distance, machinery gene, machinery locus) or (INF, '', '')."""
        best, bg, bl = INF, "", ""
        macc_by_lt = {}
        for m in mach:
            if m["instance_id"] == inst_id and m.get("locus_tag", "").strip():
                macc_by_lt[_norm(m["locus_tag"])] = (m.get("gene", ""), m["locus_tag"])
        if sordn is None:
            return best, bg, bl
        for mlt in mach_loci.get(inst_id, []):
            for mrec, mordn, *_ in gidx.get(mlt, []):
                if mrec == srec and abs(sordn - mordn) < best:
                    best = abs(sordn - mordn)
                    bg, bl = macc_by_lt.get(mlt, ("", mlt))
        return best, bg, bl

    rows, dropped = [], []

    def emit(inst_id, ss, subtype, organism, sub, pos_tuple, dist, dmethod, nmg, nml, ncand):
        rec, ordn, s, e, strand = pos_tuple
        # Replicons that MUST be staged together for ssign: the substrate's replicon + every replicon
        # carrying this instance's machinery. For cross-replicon systems (machinery and substrate on
        # different chromosomes/plasmids) this lists both, so the assembly is run whole and the pieces
        # are never split into separate ssign runs.
        replicons = ";".join(sorted({_canon(rec)} | {_canon(a) for a in mach_acc.get(inst_id, set())}))
        rows.append(
            {
                "instance_id": inst_id,
                "ss_type": ss,
                "subtype": subtype,
                "organism": organism,
                "gene": sub.get("gene", ""),
                "uniprot": sub.get("uniprot", ""),
                "effector_locus": locus_of(sub),
                "genome": sub.get("refseq_genome", ""),
                "contig": rec,
                "start": s if s is not None else "",
                "stop": e if e is not None else "",
                "strand": strand,
                "stage_replicons": replicons,
                "nearest_machinery_gene": nmg,
                "nearest_machinery_locus": nml,
                "distance_to_machinery_genes": "" if dist >= INF else dist,
                "reachable_within_3": "yes" if (dist != INF and dist <= 3) or dmethod == "self" else "no",
                "distance_method": dmethod,
                "n_candidate_substrates": ncand,
                "evidence_tier": sub.get("evidence_tier", ""),
                "verification_status": sub.get("verification_status", ""),
                "primary_ref": sub.get("primary_ref", ""),
                "citation_quote": sub.get("citation_quote") or sub.get("quote", ""),
            }
        )

    # --- machinery instances (90): collapse to the closest verified substrate ---
    for inst_id, meta in inst_meta.items():
        subs = sub_by_inst.get(inst_id, [])
        if not subs:
            dropped.append({"instance_id": inst_id, "ss_type": meta["ss_type"], "reason": "no_verified_substrate"})
            continue
        scored = []
        for r in subs:
            p = position(r)
            if p is None:
                continue
            dist, nmg, nml = nearest_machinery(inst_id, p[0], p[1])
            scored.append((dist, locus_of(r), r, p, nmg, nml))
        if not scored:
            dropped.append({"instance_id": inst_id, "ss_type": meta["ss_type"], "reason": "substrate_unpositionable"})
            continue
        scored.sort(key=lambda x: (x[0], x[1]))  # closest, then deterministic by locus
        dist, _, sub, p, nmg, nml = scored[0]
        dmethod = "gene" if dist < INF else "cross_replicon"
        emit(
            inst_id,
            meta["ss_type"],
            sub.get("subtype", ""),
            meta.get("organism", ""),
            sub,
            p,
            dist,
            dmethod,
            nmg,
            nml,
            len(subs),
        )

    # --- T5SS (self-secreting): each verified effector is its own instance, substrate = itself ---
    t5 = [r for r in pos if r["ss_type"] == "T5SS"]
    for i, r in enumerate(sorted(t5, key=lambda x: (x.get("subtype", ""), x.get("gene", ""))), 1):
        p = position(r)
        if p is None:
            dropped.append({"instance_id": f"T5SS_{i:02d}", "ss_type": "T5SS", "reason": "substrate_unpositionable"})
            continue
        emit(
            f"T5SS_{i:02d}",
            "T5SS",
            r.get("subtype", ""),
            r.get("organism", ""),
            r,
            p,
            0,
            "self",
            "(self-secreting)",
            locus_of(r),
            1,
        )

    rows.sort(key=lambda r: (r["ss_type"], r["instance_id"]))
    write_tsv(VDIR / "gold_list.tsv", FIELDS, rows)
    write_tsv(VDIR / "gold_list_dropped.tsv", ["instance_id", "ss_type", "reason"], dropped)

    from collections import Counter

    print(f"gold_list.tsv: {len(rows)} trusted system<->substrate pairs")
    print(f"  by ss_type: {dict(Counter(r['ss_type'] for r in rows))}")
    print(f"  reachable within +/-3 genes: {sum(1 for r in rows if r['reachable_within_3'] == 'yes')} / {len(rows)}")
    print(f"  distance method: {dict(Counter(r['distance_method'] for r in rows))}")
    print(f"dropped: {len(dropped)} -> {dict(Counter(d['reason'] for d in dropped))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
