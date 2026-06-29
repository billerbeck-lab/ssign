#!/usr/bin/env python3
"""Phase A final audit: recompute geometry + found_by_ssign for ALL 90 gold rows and diff vs stored.

The blind-review sweeps (scripts 67-70) checked the *labels* on each row (organism, accession, locus
identity, citation). They did NOT re-check the two columns recall is actually computed from:

  geometry      - nearest_machinery_gene / distance_to_machinery_genes / reachable_within_3
  found_by_ssign- did ssign emit a secreted protein for this effector in the tier-2 rerun?

Only the 12 re-anchored rows had those recomputed (scripts/65 Reanchor); the other 78 carry whatever
scripts/63 wrote at build time. This audit recomputes BOTH uniformly and flags every disagreement.

(a) GEOMETRY: reuses scripts/65 `Reanchor.geometry()` verbatim (same gene-order index + machinery table),
    so a corrected row stays internally consistent and a build-time error surfaces as a diff. Rows whose
    effector_locus is absent from the phase-1 index (the 12 T5SS/espP genomes it never covered) can't be
    recomputed -> reported as `unindexed`, stored value kept.

(b) FOUND_BY_SSIGN: the rule Teo fixed on 2026-06-29 -- *any* base-pair overlap between the gold span and
    an EMITTED secreted protein on the molecule = found; no overlap with any emitted protein = ssign failed
    (including the case where Bakta missed the ORF entirely, which is a real ssign failure mode the
    benchmark should capture). Molecule identity is reconciled first: gold contigs are a RefSeq/INSDC mix,
    and the rerun stages whichever the input GenBank used; `rerun_coords.CONTIG_ALIAS` maps the 3 verified
    RefSeq<->INSDC pairs (same DNA, lenratio 1.000). A contig that can't be reconciled would be a DATA FIX,
    not a fail -- but after aliasing all 90 reconcile, so none occur.

For transparency the audit also reports, per row, what the two *existing* matchers say, to show where they
diverge from the agreed rule:
  join_found       - rerun_coords.RerunIndex.join (max-overlap; can mis-report a 0-overlap protein).
  three_prime_found- scripts/24 bridge (3'-stop within 3 bp, same strand; stricter than any-overlap).

Inputs : data/phase2/verification_phase_a/gold_list_final.tsv
         data/phase1/gene_order_index.tsv, data/machinery/machinery_resolved.tsv  (geometry)
         rerun/ + rerun_fullasm/                                                  (found_by_ssign)
Outputs: data/phase2/verification_phase_a/geometry_overlap_audit.tsv  (+ printed summary)
Run    : .venv/bin/python scripts/71_geometry_overlap_audit.py
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402
from rerun_coords import RerunIndex  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
OUT = VD / "geometry_overlap_audit.tsv"

# Reuse scripts/65 Reanchor.geometry() verbatim (module name starts with a digit -> import by path).
_spec = importlib.util.spec_from_file_location("apply_corrections65", Path(__file__).parent / "65_apply_corrections.py")
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)
Reanchor = _m.Reanchor

GEOM_KEYS = ("nearest_machinery_gene", "distance_to_machinery_genes", "reachable_within_3")


def main() -> int:
    gold = read_tsv(GOLD)
    engine = Reanchor()
    ridx = RerunIndex()

    out = []
    for r in gold:
        iid, eff = r["instance_id"], r["effector_locus"].strip()
        row = {"instance_id": iid, "ss_type": r["ss_type"], "gene": r["gene"], "effector_locus": eff}

        # (a) geometry. Only rows whose instance has resolved machinery (machinery_resolved.tsv) can
        # be recomputed; the RTX (T1SS_R*) and T5SS rows were enriched from curation/convention and have
        # no resolved instance, so a recompute "disagreement" there is expected, not a mis-anchor. The
        # only alarming bucket is `DIFF`: a row WITH resolved machinery whose recompute != stored.
        try:
            g = engine.geometry(iid, eff)
            for k in GEOM_KEYS:
                row[f"stored_{k}"] = r.get(k, "")
                row[f"recomputed_{k}"] = g[k]
            if not engine.mach.get(iid):
                row["geom_status"] = "no_machinery_instance"
            else:
                row["geom_status"] = "match" if all(str(r.get(k, "")) == str(g[k]) for k in GEOM_KEYS) else "DIFF"
        except ValueError:
            row["geom_status"] = "unindexed"
            for k in GEOM_KEYS:
                row[f"stored_{k}"] = r.get(k, "")
                row[f"recomputed_{k}"] = ""

        # (b) found_by_ssign -- molecule-reconciled any-overlap with an emitted protein
        contig, start, stop = r["contig"].strip(), int(r["start"]), int(r["stop"])
        eo = ridx.emitted_overlap(contig, start, stop, strand=r.get("strand"))
        stored_found = r.get("found_by_ssign", "")
        if eo is None:
            row.update(
                molecule="UNRECONCILABLE",
                stored_found=stored_found,
                recomputed_found="",
                found_status="data_fix",
                overlap_reason="contig_absent_from_rerun",
            )
        else:
            recomputed = eo["found"]
            j = ridx.join(contig, start, stop)
            row.update(
                molecule=eo["contig"],
                stored_found=stored_found,
                recomputed_found=recomputed,
                found_status=("match" if stored_found == recomputed else "DIFF"),
                overlap_reason=eo["reason"],
                n_overlap=eo["n_overlap"],
                overlap_bp=eo["best_overlap_bp"],
                overlap_locus=eo["best_locus"],
                emitted_locus=eo["emitted_locus"],
                join_found="yes" if (j and j["emitted"]) else "no",
                three_prime_found=eo["three_prime_found"],
            )
        out.append(row)

    header = [
        "instance_id",
        "ss_type",
        "gene",
        "effector_locus",
        "geom_status",
        *[f"stored_{k}" for k in GEOM_KEYS],
        *[f"recomputed_{k}" for k in GEOM_KEYS],
        "molecule",
        "stored_found",
        "recomputed_found",
        "found_status",
        "overlap_reason",
        "n_overlap",
        "overlap_bp",
        "overlap_locus",
        "emitted_locus",
        "join_found",
        "three_prime_found",
    ]
    # geometry DIFFs first, then found DIFFs, then the rest -- the residual to eyeball.
    rank = {"DIFF": 0, "data_fix": 1, "no_machinery_instance": 2, "unindexed": 3, "match": 4}
    out.sort(key=lambda x: (rank.get(x["geom_status"], 5), rank.get(x["found_status"], 5), x["instance_id"]))
    write_tsv(OUT, header, out)

    geom = Counter(x["geom_status"] for x in out)
    found = Counter(x["found_status"] for x in out)
    print(f"geometry_overlap_audit.tsv: {len(out)} rows")
    print(f"  geometry      : {dict(geom)}")
    print(f"  found_by_ssign: {dict(found)}")

    gd = [x for x in out if x["geom_status"] == "DIFF"]
    fd = [x for x in out if x["found_status"] == "DIFF"]
    if gd:
        print(f"\n  GEOMETRY DIFFS ({len(gd)}):")
        for x in gd:
            print(
                f"    {x['instance_id']:9} {x['gene']:14} "
                f"reach {x['stored_reachable_within_3']}->{x['recomputed_reachable_within_3']} "
                f"dist {x['stored_distance_to_machinery_genes']}->{x['recomputed_distance_to_machinery_genes']} "
                f"near {x['stored_nearest_machinery_gene']}->{x['recomputed_nearest_machinery_gene']}"
            )
    if fd:
        print(f"\n  FOUND_BY_SSIGN DIFFS ({len(fd)}):")
        for x in fd:
            print(
                f"    {x['instance_id']:9} {x['gene']:14} {x['stored_found']}->{x['recomputed_found']} "
                f"({x['overlap_reason']}, n_overlap={x['n_overlap']}, bp={x['overlap_bp']}, "
                f"emit={x['emitted_locus'] or '-'}; join={x['join_found']} 3'={x['three_prime_found']})"
            )

    # divergence of the existing matchers from the agreed any-overlap rule (informational)
    rec = [x for x in out if x["found_status"] in ("match", "DIFF")]
    jdiv = [x for x in rec if x["join_found"] != x["recomputed_found"]]
    tdiv = [x for x in rec if x["three_prime_found"] != x["recomputed_found"]]
    print(f"\n  matcher divergence from any-overlap rule: join={len(jdiv)}  three_prime(scripts/24)={len(tdiv)}")
    for x in tdiv:
        print(
            f"    3' diverges: {x['instance_id']:9} {x['gene']:14} any-overlap={x['recomputed_found']} 3'={x['three_prime_found']}"
        )

    testable = [x for x in out if x["found_status"] != "data_fix"]
    yes = sum(x["recomputed_found"] == "yes" for x in testable)
    print(
        f"\n  recall under any-overlap rule: {yes}/{len(testable)} = {100 * yes / len(testable):.0f}% (all 90, no SS-type exclusion)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
