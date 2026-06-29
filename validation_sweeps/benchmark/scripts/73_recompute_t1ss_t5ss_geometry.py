#!/usr/bin/env python3
"""Recompute geometry for the 32 RTX-T1SS + T5SS gold rows that machinery_resolved.tsv could not resolve.

`scripts/65 Reanchor.geometry()` measures gene-order distance against `data/machinery/machinery_resolved.tsv`,
which only holds the proximity-panel (T1/T2/T3/T4/T6) instances -- it has NO RTX-T1SS or T5SS instance, so
those 32 rows carry whatever scripts/63 wrote and were never re-derivable by the audit. This script recomputes
them from the SAME source ssign's own proximity step uses: the secretion-system components MacSyFinder/TXSScan
actually DETECTED in the tier-2 rerun (the `# Secretion Systems` section of each `<unit>_results.csv`).

Method per row (gene-order distance, identical notion to proximity_analysis):
  1. map the gold (contig,start,stop) to its rerun ORF locus via the max-overlap join (RerunIndex), so effector
     and machinery loci are in the SAME Bakta namespace (the gold RefSeq locus is renamed on reannotation);
  2. build a gene-order index for that rerun contig (ORFs sorted by start);
  3. read the detected components and, for the row's ss_type, pick the machinery loci:
       T1SS                      -> all T1SS components (mfp/abc/omf); distance = nearest, on the same contig;
       T5bSS                     -> the T5bSS_translocator (TpsB pore) -- TPS is TWO-partner, the TpsA effector
                                    is NOT its own machinery, so this REPLACES the wrong "(self-secreting)" tag;
       T5aSS/T5cSS/T5dSS/T5eSS   -> autotransporter / inverse-AT: the effector IS the secretion apparatus.
                                    If its own locus is the detected component, that's genuinely self-secreting
                                    (distance 0, reachable N/A -- there is no separate machinery to be near).
  4. distance = |idx(effector) - idx(nearest machinery locus)|; reachable_within_3 = distance <= 3.

Self-secretors keep distance 0 but get nearest_machinery_gene = "self (T5xSS autotransporter)" so the column is
HONEST about why (the prior "(self-secreting)" with reachable=yes conflated self-secretion with proximity).
A genome with no detected component of the needed type (or a cross-contig translocator) is reported with the
reason, not a fabricated number.

Read-only diff: prints stored-vs-recomputed and writes t1ss_t5ss_geometry_recompute.tsv. Apply (if adopted) is
a separate step in scripts/65. Run: .venv/bin/python scripts/73_recompute_t1ss_t5ss_geometry.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402
from rerun_coords import CONTIG_ALIAS, RerunIndex, _emitted_path, _raw_path  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
AUDIT = VD / "geometry_overlap_audit.tsv"
OUT = VD / "t1ss_t5ss_geometry_recompute.tsv"

SELF_TYPES = {"T5aSS", "T5cSS", "T5dSS", "T5eSS"}  # autotransporters / inverse-AT: effector == apparatus
# which detected components count as "machinery" for a gold row of a given (ss_type, subtype). Matched on the
# system-type column (`ss_type`); within a T5bSS system the lone mandatory component IS the TpsB translocator.
MACHINERY_FOR = {
    "T1SS": lambda comp_ss_type: comp_ss_type.startswith("T1SS"),
    "T5bSS": lambda comp_ss_type: comp_ss_type == "T5bSS",
    "T6SS": lambda comp_ss_type: comp_ss_type.startswith("T6SS"),  # T6SSi/T6SSii components
}

# Rows outside the no-machinery set that still need rerun-detected geometry. T6SS_17/18 were re-anchored
# (2026-06-29 id re-review) to the correct effector; machinery_resolved.tsv holds only a partial T6SS
# component set for their instances (so its gene-distance came out 5-7 = "not reachable"), which contradicts
# the rerun where both effectors sit 1 gene from a detected T6SS component and emit. Recompute from the
# detected components, consistent with the T1SS/T5SS rows above.
EXTRA_TARGETS = ("T6SS_17", "T6SS_18")


def gene_order_index(unit: Path) -> dict[str, dict[str, dict]]:
    """{contig: {locus: {idx, start, end, strand}}} from the unit's raw results, ORFs ranked by start per contig."""
    by_contig: dict[str, list[dict]] = {}
    with open(_raw_path(unit), newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("contig") and r.get("start") and r.get("end"):
                by_contig.setdefault(r["contig"], []).append(r)
    idx: dict[str, dict[str, dict]] = {}
    for contig, rows in by_contig.items():
        rows.sort(key=lambda r: int(r["start"]))
        idx[contig] = {
            r["locus_tag"]: {"idx": i, "start": int(r["start"]), "end": int(r["end"]), "strand": r.get("strand", "")}
            for i, r in enumerate(rows)
        }
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
    # stop at the first blank line, the section boundary -- same as rerun_coords._read_emitted_loci /
    # bench_runout._read_secreted_chunk, so a future trailing `#`-section can't leak rows into this DictReader.
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == ""), len(lines))
    block = [ln for ln in lines[start + 1 : end] if ln.strip()]
    comps = []
    for r in csv.DictReader(block):
        if r.get("record_type") == "component" and r.get("locus_tag"):
            comps.append({"ss_type": r["ss_type"], "locus": r["locus_tag"], "gene_name": r.get("gene_name", "")})
    return comps


def main() -> int:
    gold = {r["instance_id"]: r for r in read_tsv(GOLD)}
    audit = {r["instance_id"]: r for r in read_tsv(AUDIT)}
    targets = [iid for iid, a in audit.items() if a["geom_status"] in ("no_machinery_instance", "unindexed")]
    targets += [iid for iid in EXTRA_TARGETS if iid not in targets]

    ridx = RerunIndex()
    order_cache: dict[str, dict] = {}
    comp_cache: dict[str, list] = {}

    out = []
    for iid in targets:
        r = gold[iid]
        ss, subtype = r["ss_type"], r["subtype"]
        eff_type = subtype or ss  # T5SS rows carry the real subtype; T1SS rows use ss_type
        contig, start, stop = r["contig"].strip(), int(r["start"]), int(r["stop"])
        row = {
            "instance_id": iid,
            "ss_type": ss,
            "subtype": subtype,
            "gene": r["gene"],
            "stored_nearest": r["nearest_machinery_gene"],
            "stored_dist": r["distance_to_machinery_genes"],
            "stored_reach": r["reachable_within_3"],
        }

        rc = contig if contig in ridx._contig2unit else CONTIG_ALIAS.get(contig)
        unit = ridx._contig2unit.get(rc)
        if unit is None:
            row.update(rc_nearest="", rc_locus="", rc_dist="", rc_reach="", status="contig_absent")
            out.append(row)
            continue

        ukey = str(unit)
        order = order_cache.setdefault(ukey, gene_order_index(unit))
        comps = comp_cache.setdefault(ukey, detected_components(unit))
        contig_order = order.get(rc, {})

        # effector's rerun locus = the max-overlap ORF at the gold coords on this contig
        eo = ridx.emitted_overlap(contig, start, stop, strand=r.get("strand"))
        eff_locus = eo["best_locus"] if eo else ""
        eff = contig_order.get(eff_locus)
        if eff is None:
            row.update(rc_nearest="", rc_locus="", rc_dist="", rc_reach="", status="effector_not_in_order")
            out.append(row)
            continue

        if eff_type in SELF_TYPES:
            # Single-polypeptide self-secretor (classical AT T5aSS/T5cSS, patatin-like T5dSS, inverse-AT T5eSS):
            # the effector IS the secretion apparatus, so there is no separate machinery to be "near" -> distance 0,
            # reachable not applicable. T5dSS/T5eSS have no TXSScan profile, so TXSScan won't detect the effector as
            # a component; that is a model gap, not evidence of a partner, so they stay self-secreting + a note.
            is_self = any(c["locus"] == eff_locus and c["ss_type"] == eff_type for c in comps)
            row.update(
                rc_nearest=f"self ({eff_type} autotransporter)" if is_self else f"self ({eff_type}, no TXSScan model)",
                rc_locus=eff_locus,
                rc_dist="0",
                rc_reach="self",
                status="self_secreting" if is_self else "self_no_model",
            )
            out.append(row)
            continue

        match = MACHINERY_FOR.get(eff_type, lambda t: False)
        mach = [c for c in comps if match(c["ss_type"]) and c["locus"] in contig_order and c["locus"] != eff_locus]
        if not mach:
            row.update(
                rc_nearest="", rc_locus="", rc_dist="", rc_reach="no", status=f"no_{eff_type}_component_on_contig"
            )
            out.append(row)
            continue

        nearest = min(mach, key=lambda c: abs(contig_order[c["locus"]]["idx"] - eff["idx"]))
        dist = abs(contig_order[nearest["locus"]]["idx"] - eff["idx"])
        row.update(
            rc_nearest=nearest["gene_name"],
            rc_locus=nearest["locus"],
            rc_dist=str(dist),
            rc_reach="yes" if dist <= 3 else "no",
            status="recomputed",
        )
        out.append(row)

    header = [
        "instance_id",
        "ss_type",
        "subtype",
        "gene",
        "stored_nearest",
        "stored_dist",
        "stored_reach",
        "rc_nearest",
        "rc_locus",
        "rc_dist",
        "rc_reach",
        "status",
    ]
    out.sort(key=lambda x: (x["ss_type"], x["subtype"], x["instance_id"]))
    write_tsv(OUT, header, out)

    from collections import Counter

    print(f"t1ss_t5ss_geometry_recompute.tsv: {len(out)} rows")
    print(f"  status: {dict(Counter(x['status'] for x in out))}")
    print(f"\n  {'id':9} {'sub':7} {'gene':10} stored->recomputed")
    for x in out:
        changed = (x["stored_dist"] != x["rc_dist"]) or (x["stored_reach"] != x["rc_reach"])
        flag = "  <-- CHANGED" if changed else ""
        print(
            f"  {x['instance_id']:9} {x['subtype'] or x['ss_type']:7} {x['gene']:10} "
            f"[{x['stored_nearest'][:18]} d{x['stored_dist']} {x['stored_reach']}] -> "
            f"[{x['rc_nearest'][:22]} d{x['rc_dist']} {x['rc_reach']}] ({x['status']}){flag}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
