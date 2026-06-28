#!/usr/bin/env python3
"""Phase A second-pass review: build the two agent-sweep batch sets from the 90-row gold list.

Two independent sweeps, each a separate dir under gold_review2/ with its own batches/ + verdicts/ + SCHEMA:

  sweep1_identity  - the rows the deterministic precompute (scripts/66) could not clear: every row whose
                     identity_signals flag contains `no_acc` or `len_mismatch`. Agents find the correct
                     genome-matched UniProt accession (preferring a reviewed Swiss-Prot entry) or confirm
                     none exists, and adjudicate the two length mismatches. Batches carry the genome
                     coordinates AND the precomputed signal so agents adjudicate, not re-derive.
  sweep2_citation  - ALL 90 rows. Agents check the primary reference: organism match, documents SECRETION
                     (strict: supernatant/translocation/assay, not mere virulence), quote verbatim.

Per-sweep batches are chunked small (<=8 identity, <=10 citation) so each agent does careful lookups, and
the citation sweep is grouped by ss_type so an agent sees one system family at a time.

Inputs : data/phase2/verification_phase_a/gold_list_final.tsv
         data/phase2/verification_phase_a/identity_signals.tsv
Outputs: data/phase2/verification_phase_a/gold_review2/sweep1_identity/batches/*.tsv (+ batch_plan.json)
         data/phase2/verification_phase_a/gold_review2/sweep2_citation/batches/*.tsv (+ batch_plan.json)
Run    : .venv/bin/python scripts/67_build_review2_batches.py
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
SIGNALS = VD / "identity_signals.tsv"
GR2 = VD / "gold_review2"

IDENTITY_CHUNK = 8
CITATION_CHUNK = 10

# identity_signals flags (scripts/66) that the precompute could not auto-clear -> need an agent.
# Kept as a named constant so a flag rename in 66 fails the guard below instead of silently shrinking the sweep.
TARGET_FLAGS = ("no_acc", "len_mismatch")
SS_TYPE_RE = re.compile(r"^[A-Za-z0-9]+$")  # ss_type is used as a batch filename stem; keep it path-safe

IDENTITY_COLS = [
    "row_id",
    "ss_type",
    "gene",
    "organism",
    "genome",
    "contig",
    "start",
    "stop",
    "strand",
    "effector_locus",
    "uniprot_in",
    "uniprot_len",
    "expected_aa",
    "len_ratio",
    "signal_flags",
    "uniprot_protein",
]
CITATION_COLS = [
    "row_id",
    "ss_type",
    "subtype",
    "gene",
    "uniprot",
    "organism",
    "primary_ref",
    "citation_quote",
]


def _chunk(rows: list, n: int) -> list[list]:
    return [rows[i : i + n] for i in range(0, len(rows), n)]


def _write_sweep(name: str, cols: list[str], batches: dict[str, list[dict]]) -> list[dict]:
    bdir = GR2 / name / "batches"
    bdir.mkdir(parents=True, exist_ok=True)
    (GR2 / name / "verdicts").mkdir(parents=True, exist_ok=True)
    plan = []
    for bid, brows in batches.items():
        write_tsv(bdir / f"{bid}.tsv", cols, brows)
        plan.append({"id": bid, "n": len(brows)})
    (GR2 / name / "batch_plan.json").write_text(json.dumps(plan, indent=2))
    return plan


def main() -> int:
    gold = {r["instance_id"]: r for r in read_tsv(GOLD)}
    signals = {r["instance_id"]: r for r in read_tsv(SIGNALS)}

    # ---- sweep 1: identity (only rows the precompute could not clear) ----
    flagged = [iid for iid, s in signals.items() if set(TARGET_FLAGS) & set(s["flags"].split(";"))]
    if not flagged:
        raise SystemExit(f"no rows carry any of {TARGET_FLAGS} -- did a flag get renamed in scripts/66?")
    flagged.sort(key=lambda i: (gold[i]["ss_type"], i))
    id_rows = []
    for iid in flagged:
        g, s = gold[iid], signals[iid]
        id_rows.append(
            {
                "row_id": iid,
                "ss_type": g["ss_type"],
                "gene": g["gene"],
                "organism": g["organism"],
                "genome": g["genome"],
                "contig": g["contig"],
                "start": g["start"],
                "stop": g["stop"],
                "strand": g["strand"],
                "effector_locus": g["effector_locus"],
                "uniprot_in": g["uniprot"],
                "uniprot_len": s["uniprot_len"],
                "expected_aa": s["expected_aa"],
                "len_ratio": s["len_ratio"],
                "signal_flags": s["flags"],
                "uniprot_protein": s["uniprot_protein"],
            }
        )
    id_batches = {f"s1_{i + 1}": c for i, c in enumerate(_chunk(id_rows, IDENTITY_CHUNK))}
    p1 = _write_sweep("sweep1_identity", IDENTITY_COLS, id_batches)

    # ---- sweep 2: citation (all 90, grouped by ss_type) ----
    by_type: dict[str, list[dict]] = defaultdict(list)
    for iid, g in gold.items():
        by_type[g["ss_type"]].append(g)
    cit_batches: dict[str, list[dict]] = {}
    for st in sorted(by_type):
        if not SS_TYPE_RE.match(st):
            raise SystemExit(f"ss_type {st!r} is not filename-safe; refusing to build batch ids from it")
        rows = sorted(by_type[st], key=lambda r: r["instance_id"])
        for i, chunk in enumerate(_chunk(rows, CITATION_CHUNK), 1):
            cit_batches[f"{st}_{i}"] = [
                {
                    "row_id": r["instance_id"],
                    "ss_type": r["ss_type"],
                    "subtype": r["subtype"],
                    "gene": r["gene"],
                    "uniprot": r["uniprot"],
                    "organism": r["organism"],
                    "primary_ref": r["primary_ref"],
                    "citation_quote": r["citation_quote"],
                }
                for r in chunk
            ]
    p2 = _write_sweep("sweep2_citation", CITATION_COLS, cit_batches)

    print(f"sweep1_identity: {len(id_rows)} flagged rows -> {len(p1)} batches {[b['id'] for b in p1]}")
    print(f"  flagged: {flagged}")
    print(f"sweep2_citation: {sum(b['n'] for b in p2)} rows -> {len(p2)} batches {[b['id'] for b in p2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
