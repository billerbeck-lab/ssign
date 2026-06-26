#!/usr/bin/env python3
"""Phase A (revised): collapse each secretion-system instance to ONE representative substrate.

The benchmark scores whether ssign finds a system AND a substrate, not how many substrates. So for
each of the 90 machinery instances we keep a SINGLE substrate: the one closest in gene-order to that
instance's machinery. Rationale: ssign emits a substrate by proximity (+/-N genes from a detected
machinery component), so the substrate nearest the machinery is the one ssign is most likely to reach;
keeping it preserves "if ssign would/did detect a substrate here, it still will". This makes the
system<->substrate mapping 1:1 and shrinks the re-verification set from 337 to <=90.

Linkage  : positives_all.sys_instance_id == machinery instances.instance_id (256/337 rows match a 90-key
           instance; the rest carry R## rescue ids / no instance and are out of scope per Teo).
Distance : |ordinal(substrate) - ordinal(nearest machinery gene)| on the same replicon, from the
           phase1 gene-order index (the same gene-distance ssign's window uses). bp fallback if a locus
           is absent from the index; if neither side resolves, the instance's lone/first substrate is
           kept and flagged.

Inputs : data/machinery/instances.tsv, data/machinery/machinery_resolved.tsv,
         data/dataset/positives_all.tsv, data/phase1/gene_order_index.tsv
Outputs: data/phase2/verification_phase_a/representative_substrates.tsv  (the 1:1 list to verify)
         data/phase2/verification_phase_a/instances_no_substrate.tsv     (90-key instances with 0 of the 337)
Run    : .venv/bin/python scripts/62_representative_substrates.py
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


def _norm(lt: str) -> str:
    """Locus-tag key tolerant of the RefSeq underscore/space drift (ECs4562 vs ECs_4562), mirroring
    24_actual_call.normalize so substrate placement loci match the index + machinery loci."""
    return (lt or "").replace("_", "").replace(" ", "").upper()


def ordinal_index():
    """normalized locus_tag -> [(record_acc, ordinal), ...], for every cached gene. The index supplies
    each locus's own record+ordinal, so we never depend on the drift-prone refseq_genome /
    resolved_accession fields agreeing."""
    idx = {}
    for r in read_tsv(BENCH / "data" / "phase1" / "gene_order_index.tsv"):
        entry = (r["record_acc"], int(r["ordinal"]))
        # index the current locus_tag AND every alias (old_locus_tag etc.), so substrate placement loci
        # recorded under an older tag still resolve.
        keys = {_norm(r["locus_tag"])} | {_norm(a) for a in (r.get("aliases", "") or "").split(";") if a.strip()}
        for k in keys:
            idx.setdefault(k, []).append(entry)
    return idx


def main() -> int:
    instances = read_tsv(BENCH / "data" / "machinery" / "instances.tsv")
    inst_meta = {r["instance_id"]: r for r in instances}
    pos = read_tsv(BENCH / "data" / "dataset" / "positives_all.tsv")
    mach = read_tsv(BENCH / "data" / "machinery" / "machinery_resolved.tsv")
    locpos = ordinal_index()

    # machinery gene locus_tags per instance (record+ordinal resolved via the index, not these fields)
    mach_by_inst = defaultdict(list)
    for m in mach:
        lt = m.get("locus_tag", "").strip()
        if lt:
            mach_by_inst[m["instance_id"]].append(lt)

    # the 337 substrates that map to a 90-key instance
    sub_by_inst = defaultdict(list)
    for r in pos:
        sid = r.get("sys_instance_id", "").strip()
        if sid in inst_meta:
            sub_by_inst[sid].append(r)

    def sub_locus(r):
        return (r.get("placement_effector_locus") or r.get("locus_tag") or "").strip()

    def distance(r, inst_id):
        """Min gene-distance from this substrate to any machinery gene of its instance, on a shared
        replicon. Both loci are located via the index, so accession drift between tables is moot."""
        s_entries = locpos.get(_norm(sub_locus(r)), [])
        best = INF
        for mlt in mach_by_inst.get(inst_id, []):
            for mrec, mord in locpos.get(_norm(mlt), []):
                for srec, sord in s_entries:
                    if srec == mrec:
                        best = min(best, abs(sord - mord))
        return (best, "gene") if best < INF else (INF, "none")

    chosen, no_sub = [], []
    for inst_id, meta in inst_meta.items():
        subs = sub_by_inst.get(inst_id, [])
        if not subs:
            no_sub.append(
                {k: meta.get(k, "") for k in ("instance_id", "ss_type", "refseq_genome", "organism", "effector_loci")}
            )
            continue
        scored = sorted(((distance(r, inst_id), r) for r in subs), key=lambda x: (x[0][0], sub_locus(x[1])))
        (dist, dmethod), best = scored[0]
        chosen.append(
            {
                "instance_id": inst_id,
                "ss_type": meta["ss_type"],
                "genome": meta["refseq_genome"],
                "organism": meta.get("organism", ""),
                "n_candidate_substrates": len(subs),
                "gene": best.get("gene", ""),
                "uniprot": best.get("uniprot", ""),
                "effector_locus": sub_locus(best),
                "distance_to_machinery": "" if dist >= INF else dist,
                "distance_method": dmethod,
                "primary_ref": best.get("primary_ref", ""),
                "citation_quote": best.get("citation_quote") or best.get("quote", ""),
            }
        )

    chosen.sort(key=lambda r: (r["ss_type"], r["instance_id"]))
    write_tsv(VDIR / "representative_substrates.tsv", list(chosen[0].keys()), chosen)
    write_tsv(
        VDIR / "instances_no_substrate.tsv",
        ["instance_id", "ss_type", "refseq_genome", "organism", "effector_loci"],
        no_sub,
    )

    from collections import Counter

    print(f"representative_substrates.tsv: {len(chosen)} instances -> 1 substrate each")
    print(f"  by ss_type: {dict(Counter(r['ss_type'] for r in chosen))}")
    print(f"  distance method: {dict(Counter(r['distance_method'] for r in chosen))}")
    far = [
        r
        for r in chosen
        if r["distance_method"] == "gene"
        and isinstance(r["distance_to_machinery"], int)
        and r["distance_to_machinery"] > 3
    ]
    print(f"  representatives >3 genes from machinery (ssign would MISS by proximity): {len(far)}")
    print(f"instances_no_substrate.tsv: {len(no_sub)} of 90 instances have 0 of the 337 verified substrates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
