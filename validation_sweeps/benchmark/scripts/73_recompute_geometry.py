#!/usr/bin/env python3
"""Recompute gold-list geometry (nearest machinery / distance / reachable) for ALL 90 rows from the
secretion-system components ssign actually DETECTED in the tier-2 rerun.

WHY (replaces machinery_resolved.tsv for these columns): the curated `data/machinery/machinery_resolved.tsv`
is a hand-resolved, PARTIAL component list. For several rows it places the effector hundreds-to-thousands of
genes from "its machinery" while the rerun shows the effector sits 1-3 genes from the detected apparatus
(e.g. CopN curated 520 vs detected 1; BipC 11 vs 3; TplE 1454 vs 1; Tle1 265 vs 1). Those produced
`reachable_within_3=no` for effectors ssign clearly emitted via proximity. The honest, internally-consistent
source is the machinery TXSScan/MacSyFinder detected in the same run that produced `found_by_ssign`.

Distance is measured to the nearest detected component of the effector's COGNATE system type (gene-order
distance on the same contig), exactly the notion ssign's per-component proximity window uses:
  T1SS->T1SS, T2SS->T2SS, T3SS->T3SS, T4SS->pT4SSt, T6SS->T6SSi*, T5bSS->T5bSS translocator.
  T5aSS/T5cSS/T5dSS/T5eSS are single-polypeptide self-secretors: no separate machinery, distance 0, reach="self".
  (T4aP / MSH detected labels are type-IV / MSHA pili, NOT secretion systems -- never matched.)
Cognate-only is deliberate: VirA stays reach=no (24 genes from its T3SS apparatus) even though it emitted via a
T5aSS autotransporter 1 gene away -- a real "recovered, but not by cognate proximity" case, not a near-miss.

Statuses: recomputed | self_secreting | self_no_model (T5d/T5e, no TXSScan profile) | no_cognate_machinery
(cognate system not detected on the effector's contig -> reach=no) | contig_absent | effector_not_in_order
(can't map the effector to a rerun ORF -> stored value kept by scripts/65).

This changes ONLY the three geometry columns; it never touches found_by_ssign, so recall is unaffected.
Read-only: writes geometry_recompute.tsv; scripts/65 applies it. Run: .venv/bin/python scripts/73_recompute_geometry.py
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402
from rerun_coords import CONTIG_ALIAS, RerunIndex, _emitted_path, _raw_path  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
OUT = VD / "geometry_recompute.tsv"

SELF_TYPES = {"T5aSS", "T5cSS", "T5dSS", "T5eSS"}  # single-polypeptide autotransporters / inverse-AT
# cognate detected-component matcher per effector type (matched on the component's `ss_type` column).
# pili (T4aP, MSH) are intentionally absent -- they are not secretion machinery.
MACHINERY_FOR = {
    "T1SS": lambda t: t.startswith("T1SS"),
    "T2SS": lambda t: t.startswith("T2SS"),
    "T3SS": lambda t: t.startswith("T3SS"),
    "T4SS": lambda t: t == "pT4SSt",
    "T5bSS": lambda t: t == "T5bSS",
    "T6SS": lambda t: t.startswith("T6SS"),
}
GEOM_KEYS = ("nearest_machinery_gene", "distance_to_machinery_genes", "reachable_within_3")


def gene_order_index(unit: Path) -> dict[str, dict[str, int]]:
    """{contig: {locus: idx}} -- ORFs ranked by start within each contig of the unit's raw results."""
    by_contig: dict[str, list[dict]] = {}
    with open(_raw_path(unit), newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("contig") and r.get("start"):
                by_contig.setdefault(r["contig"], []).append(r)
    idx: dict[str, dict[str, int]] = {}
    for contig, rows in by_contig.items():
        rows.sort(key=lambda r: int(r["start"]))
        idx[contig] = {r["locus_tag"]: i for i, r in enumerate(rows)}
    return idx


def detected_components(unit: Path) -> list[dict]:
    """[{ss_type, locus, gene_name}] from the `# Secretion Systems` section (component rows with a real locus)."""
    path = _emitted_path(unit)
    if not path:
        return []
    lines = Path(path).read_text().split("\n")
    start = next((i for i, ln in enumerate(lines) if ln.lstrip("# ").lower().startswith("secretion systems")), None)
    if start is None:
        return []
    # NB: do NOT stop at the first blank line. results.csv splits secretion systems into TWO subsections --
    # "(with secreted proteins)" and "(other)" -- divided by a blank line, and a system with no emitted
    # substrate (frequently the T6SS) lives only in "(other)". Read to EOF and keep just `component` rows;
    # the blank, the second "# Secretion Systems (other)" header, and its repeated column header all fail the
    # record_type=="component" filter. (This is why detected_components canNOT share _read_emitted_loci's
    # single-section reader -- that one must terminate, this one must not.)
    block = [ln for ln in lines[start + 1 :] if ln.strip()]
    return [
        {"ss_type": r["ss_type"], "locus": r["locus_tag"], "gene_name": r.get("gene_name", "")}
        for r in csv.DictReader(block)
        if r.get("record_type") == "component" and r.get("locus_tag")
    ]


def main() -> int:
    gold = read_tsv(GOLD)
    ridx = RerunIndex()
    order_cache: dict[str, dict] = {}
    comp_cache: dict[str, list] = {}

    out = []
    for r in gold:
        ss, subtype = r["ss_type"], r["subtype"]
        eff_type = subtype or ss
        contig, start, stop = r["contig"].strip(), int(r["start"]), int(r["stop"])
        row = {"instance_id": r["instance_id"], "ss_type": ss, "subtype": subtype, "gene": r["gene"]}
        for k in GEOM_KEYS:
            row[f"stored_{k}"] = r.get(k, "")
        blank = {"rc_nearest": "", "rc_locus": "", "rc_dist": "", "rc_reach": ""}

        rc = contig if contig in ridx._contig2unit else CONTIG_ALIAS.get(contig)
        unit = ridx._contig2unit.get(rc)
        if unit is None:
            out.append({**row, **blank, "status": "contig_absent"})
            continue
        order = order_cache.setdefault(str(unit), gene_order_index(unit))
        comps = comp_cache.setdefault(str(unit), detected_components(unit))
        contig_order = order.get(rc, {})

        eo = ridx.emitted_overlap(contig, start, stop, strand=r.get("strand"))
        eff_locus = eo["best_locus"] if eo else ""
        if eff_locus not in contig_order:
            # the effector overlaps no rerun ORF -> Bakta missed the gene (e.g. NleD in an IS-dense region).
            # An un-annotated gene can't be near machinery or emitted, so reach="no" is the honest answer.
            out.append(
                {
                    **row,
                    **blank,
                    "rc_nearest": "(effector ORF absent from rerun)",
                    "rc_reach": "no",
                    "status": "effector_not_in_order",
                }
            )
            continue
        ei = contig_order[eff_locus]

        if eff_type in SELF_TYPES:
            is_self = any(c["locus"] == eff_locus and c["ss_type"] == eff_type for c in comps)
            out.append(
                {
                    **row,
                    "rc_nearest": f"self ({eff_type} autotransporter)"
                    if is_self
                    else f"self ({eff_type}, no TXSScan model)",
                    "rc_locus": eff_locus,
                    "rc_dist": "0",
                    "rc_reach": "self",
                    "status": "self_secreting" if is_self else "self_no_model",
                }
            )
            continue

        match = MACHINERY_FOR.get(eff_type, lambda t: False)
        mach = [c for c in comps if match(c["ss_type"]) and c["locus"] in contig_order and c["locus"] != eff_locus]
        if not mach:
            out.append(
                {
                    **row,
                    **blank,
                    "rc_nearest": "(no cognate machinery detected)",
                    "rc_reach": "no",
                    "status": "no_cognate_machinery",
                }
            )
            continue
        nearest = min(mach, key=lambda c: abs(contig_order[c["locus"]] - ei))
        dist = abs(contig_order[nearest["locus"]] - ei)
        out.append(
            {
                **row,
                "rc_nearest": nearest["gene_name"],
                "rc_locus": nearest["locus"],
                "rc_dist": str(dist),
                "rc_reach": "yes" if dist <= 3 else "no",
                "status": "recomputed",
            }
        )

    header = [
        "instance_id",
        "ss_type",
        "subtype",
        "gene",
        *[f"stored_{k}" for k in GEOM_KEYS],
        "rc_nearest",
        "rc_locus",
        "rc_dist",
        "rc_reach",
        "status",
    ]
    out.sort(key=lambda x: (x["ss_type"], x["subtype"], x["instance_id"]))
    write_tsv(OUT, header, out)

    print(f"geometry_recompute.tsv: {len(out)} rows (all gold rows)")
    print(f"  status: {dict(Counter(x['status'] for x in out))}")
    changed = [x for x in out if x["stored_reachable_within_3"] != x["rc_reach"] and x["status"] != "contig_absent"]
    print(f"  reachable changed vs stored: {len(changed)}")
    for x in changed:
        print(
            f"    {x['instance_id']:9} {x['subtype'] or x['ss_type']:7} {x['gene']:14} "
            f"reach {x['stored_reachable_within_3']}->{x['rc_reach']:4} "
            f"dist {x['stored_distance_to_machinery_genes']}->{x['rc_dist']} ({x['status']})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
