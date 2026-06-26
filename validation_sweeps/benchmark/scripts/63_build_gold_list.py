#!/usr/bin/env python3
"""Phase A: build the trusted 1-substrate-per-system GOLD LIST (the paper reference set).

One representative secreted protein per secretion-system instance, fully self-contained (host, gene,
accession, genome coordinates, the machinery it sits near + gene-distance, reachability, found-by-ssign,
DOI, citation). This is the list we trust, reference, and re-verify; it replaces the 337-effector table.

Built from the SAME answer key the recall figure (52_system_recall) uses, so the gold list and the
figure share one instance basis:
  * instance set + per-effector distance to machinery come from `ceiling_per_effector.tsv`
    (`instance_id` covers the machinery-audit instances AND the T1SS rescue instances that carry the
    RTX toxins / proteases; `nearest_dist` + `reachable_n3` are the trusted proximity calls).
  * the cleaned, citation-verified effector set + ssign_call (found) come from `actual_per_effector`
    via clean_dataset (T3SS from the detection-on tag, others default -- mirrors fig 06).
  * provenance (DOI, quote, organism, subtype) + coordinates come from positives_all / the gene-order
    index.

Per instance we keep the effector CLOSEST to the machinery (min nearest_dist) -- the one ssign's
+/-N-gene window is most likely to reach -- making system<->substrate 1:1. T5SS (self-secreting; not in
the proximity tables) is added directly from the verified T5SS effectors, substrate = itself.

`stage_replicons` records every replicon that must be staged together for ssign (substrate replicon +
machinery replicons), so multi-replicon hosts (V. cholerae chr2, Ralstonia megaplasmid, ...) are run
whole, never split.

Inputs : data/phase1/ceiling_per_effector.tsv, data/phase2/actual_per_effector.{panel_genbank_default,
         panel_genbank_t3ss}.tsv, data/dataset/positives_all.tsv, data/machinery/machinery_resolved.tsv,
         data/phase1/gene_order_index.tsv
Outputs: data/phase2/verification_phase_a/gold_list.tsv
Run    : .venv/bin/python scripts/63_build_gold_list.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import clean_dataset  # noqa: E402
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
VDIR = P2 / "verification_phase_a"
INF = 10**9
REPLICON_ALIAS = {"AE004091": "NC_002516.2"}  # PAO1 INSDC/RefSeq duplicate -> one chromosome

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
    "found_by_ssign",
    "n_candidate_substrates",
    "evidence_tier",
    "verification_status",
    "primary_ref",
    "citation_quote",
]


def _norm(lt: str) -> str:
    return (lt or "").replace("_", "").replace(" ", "").upper()


def _canon(acc: str) -> str:
    return REPLICON_ALIAS.get(acc, REPLICON_ALIAS.get(acc.split(".")[0], acc))


def gene_order():
    """normalized locus -> [(record_acc, start, end, strand)] incl. aliases (old tags)."""
    idx = defaultdict(list)
    for r in read_tsv(BENCH / "data" / "phase1" / "gene_order_index.tsv"):
        idx[_norm(r["locus_tag"])].append(
            (
                r["record_acc"],
                int(r["start"]) if r.get("start") else None,
                int(r["end"]) if r.get("end") else None,
                r.get("strand", ""),
            )
        )
        for a in (r.get("aliases", "") or "").split(";"):
            if a.strip():
                idx[_norm(a)].append((r["record_acc"], None, None, ""))
    return idx


def main() -> int:
    gidx = gene_order()
    ceiling = {}  # (ss_type, effector_locus) and (ss_type, uniprot) -> ceiling row
    for r in read_tsv(BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"):
        if r.get("effector_locus", "").strip():
            ceiling[(r["ss_type"], r["effector_locus"])] = r
        if r.get("uniprot", "").strip() and r["uniprot"] != "-":
            ceiling.setdefault((r["ss_type"], r["uniprot"]), r)
    pos_by_uni, pos_by_loc, pos_by_gene = {}, {}, {}
    for r in read_tsv(BENCH / "data" / "dataset" / "positives_all.tsv"):
        if r.get("uniprot", "").strip() and r["uniprot"] != "-":
            pos_by_uni.setdefault(r["uniprot"], r)
        for lc in (r.get("placement_effector_locus"), r.get("locus_tag")):
            if lc and lc.strip():
                pos_by_loc.setdefault(lc.strip(), r)
        pos_by_gene.setdefault((r.get("gene", "").lower(), r["ss_type"]), r)
    mach_acc = defaultdict(set)  # instance -> {machinery replicon}, from the index (alias-safe)
    for m in read_tsv(BENCH / "data" / "machinery" / "machinery_resolved.tsv"):
        for rec, *_ in gidx.get(_norm(m.get("locus_tag", "")), []):
            mach_acc[m["instance_id"]].add(_canon(rec))

    dropped = clean_dataset.dropped_id()

    def cei(r):
        return ceiling.get((r["ss_type"], r.get("effector_locus", ""))) or ceiling.get(
            (r["ss_type"], r.get("uniprot", ""))
        )

    def coords(eff_locus):
        for rec, s, e, strand in gidx.get(_norm(eff_locus), []):
            if s is not None:
                return rec, s, e, strand
        return "", None, None, ""

    # --- spine: cleaned, verified effectors grouped into instances (T3SS detection-on tag, else default) ---
    inst_effs = defaultdict(list)  # (ss_type, instance_id) -> [effector dicts enriched with ceiling]
    for tag, only in (
        ("panel_genbank_default", lambda ss: ss != "T3SS"),
        ("panel_genbank_t3ss", lambda ss: ss == "T3SS"),
    ):
        for r in clean_dataset.load_clean_actual(P2 / f"actual_per_effector.{tag}.tsv"):
            if not only(r["ss_type"]):
                continue
            if (r["gene"], r["uniprot"]) in dropped or (r["gene"], r.get("effector_locus", "")) in dropped:
                continue
            c = cei(r)
            if not c or not c.get("instance_id", "").strip():
                continue
            inst_effs[(r["ss_type"], c["instance_id"])].append({**r, "_cei": c})

    rows, dropped_inst = [], []
    for (ss, inst_id), effs in inst_effs.items():
        if not any(e["testable"] == "yes" for e in effs):
            dropped_inst.append({"instance_id": inst_id, "ss_type": ss, "reason": "no_testable_effector"})
            continue

        def dist(e):
            d = e["_cei"].get("nearest_dist", "").strip()
            return int(d) if d.isdigit() else INF

        # Representative = the CLOSEST-to-machinery effector, but preferring one that preserves the
        # instance's status: a found effector over an unfound one, then a reachable over an unreachable.
        # So the 1:1 collapse keeps "instance found/reachable" intact (matches fig 06's any-effector
        # logic) while still picking the nearest qualifying substrate ("if ssign did detect it, it
        # still will").
        def rank(e):
            found = e["ssign_call"] == "emitted_secreted"
            reach = e["reachable_n3"] == "true"  # actual_per_effector value (post answer-key anchor fixes)
            return (0 if found else 1, 0 if reach else 1, dist(e), e.get("effector_locus", ""))

        rep = min(effs, key=rank)
        c = rep["_cei"]
        p = (
            pos_by_uni.get(rep["uniprot"])
            or pos_by_loc.get(rep.get("effector_locus", ""))
            or pos_by_gene.get((rep["gene"].lower(), ss))
            or {}
        )
        rec, s, e, strand = coords(rep.get("effector_locus", ""))
        rec = _canon(rec) if rec else _canon(rep.get("unit_id", ""))
        replicons = ";".join(sorted({rec} | mach_acc.get(inst_id, set()) - {""}))
        d = c.get("nearest_dist", "").strip()
        rows.append(
            {
                "instance_id": inst_id,
                "ss_type": ss,
                "subtype": p.get("subtype", ""),
                "organism": p.get("organism", "") or c.get("refseq_genome", ""),
                "gene": rep["gene"],
                "uniprot": rep["uniprot"],
                "effector_locus": rep.get("effector_locus", ""),
                "genome": c.get("refseq_genome", ""),
                "contig": rec,
                "start": s if s is not None else "",
                "stop": e if e is not None else "",
                "strand": strand,
                "stage_replicons": replicons,
                "nearest_machinery_gene": c.get("nearest_gene", ""),
                "nearest_machinery_locus": c.get("nearest_locus", ""),
                "distance_to_machinery_genes": d if d.isdigit() else "",
                "reachable_within_3": "yes" if rep["reachable_n3"] == "true" else "no",
                "found_by_ssign": "yes" if rep["ssign_call"] == "emitted_secreted" else "no",
                "n_candidate_substrates": len(effs),
                "evidence_tier": p.get("evidence_tier", ""),
                "verification_status": p.get("verification_status", ""),
                "primary_ref": p.get("primary_ref", ""),
                "citation_quote": p.get("citation_quote") or p.get("quote", ""),
            }
        )

    # --- T5SS (self-secreting): each verified effector is its own instance, substrate = itself ---
    # found comes from the tier-2 rerun's real emission (rerun_coords), not 53_t5ss_recall (which infers
    # via local MacSyFinder and over-counts, see 4.2). reachable = subtype has a TXSScan model: T5a/b/c
    # yes, T5d/T5e (plpD/eae) no -- the paper builds no profile for those, so ssign cannot detect them.
    try:
        from rerun_coords import RerunIndex  # noqa: PLC0415

        ridx = RerunIndex()
    except Exception:
        ridx = None
    t5 = [r for r in read_tsv(BENCH / "data" / "dataset" / "positives_all.tsv") if r["ss_type"] == "T5SS"]
    for i, r in enumerate(sorted(t5, key=lambda x: (x.get("subtype", ""), x.get("gene", ""))), 1):
        rec, s, e, strand = coords(r.get("placement_effector_locus") or r.get("locus_tag", ""))
        if not rec and r.get("contig"):
            rec, s, e, strand = (
                r["contig"],
                int(r["start"]) if r.get("start") else None,
                int(r["stop"]) if r.get("stop") else None,
                r.get("strand", ""),
            )
        has_model = r.get("subtype", "") not in ("T5dSS", "T5eSS")
        found = ""
        if ridx and r.get("contig") and r.get("start") and r.get("stop"):
            j = ridx.join(r["contig"], int(r["start"]), int(r["stop"]))
            found = "yes" if (j and j["emitted"]) else "no"
        rows.append(
            {
                "instance_id": f"T5SS_{i:02d}",
                "ss_type": "T5SS",
                "subtype": r.get("subtype", ""),
                "organism": r.get("organism", ""),
                "gene": r.get("gene", ""),
                "uniprot": r.get("uniprot", ""),
                "effector_locus": r.get("locus_tag", ""),
                "genome": r.get("refseq_genome", ""),
                "contig": rec or r.get("contig", ""),
                "start": s if s is not None else "",
                "stop": e if e is not None else "",
                "strand": strand,
                "stage_replicons": _canon(rec or r.get("contig", "")),
                "nearest_machinery_gene": "(self-secreting)" if has_model else "(no TXSScan model)",
                "nearest_machinery_locus": r.get("locus_tag", ""),
                "distance_to_machinery_genes": "0" if has_model else "",
                "reachable_within_3": "yes" if has_model else "no",
                "found_by_ssign": found,
                "n_candidate_substrates": 1,
                "evidence_tier": r.get("evidence_tier", ""),
                "verification_status": r.get("verification_status", ""),
                "primary_ref": r.get("primary_ref", ""),
                "citation_quote": r.get("citation_quote") or r.get("quote", ""),
            }
        )

    rows.sort(key=lambda r: (r["ss_type"], r["instance_id"]))
    write_tsv(VDIR / "gold_list.tsv", FIELDS, rows)

    from collections import Counter

    print(f"gold_list.tsv: {len(rows)} trusted system<->substrate pairs")
    print(f"  by ss_type: {dict(Counter(r['ss_type'] for r in rows))}")
    print(f"  reachable within +/-3 genes: {sum(1 for r in rows if r['reachable_within_3'] == 'yes')}")
    print(f"  found by ssign (non-T5): {sum(1 for r in rows if r['found_by_ssign'] == 'yes')}")
    print(
        f"  T1SS rescue (RTX toxins/proteases) included: {sum(1 for r in rows if r['instance_id'].startswith('T1SS_R'))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
