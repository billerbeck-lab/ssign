#!/usr/bin/env python3
"""Split gold_list_final.tsv into N blind review batches for the exhaustive holistic re-read (Teo, 2026-06-30).

Blind = the agent sees only the row's CLAIM (organism/gene/uniprot/locus/coords/ref/quote), NOT the `correction`
column or any of my prior verdict columns (found/reachable/distance/verification_status). The agent verifies the
claim from scratch (UniProt entry + the cited paper) so its verdict isn't anchored on what I already concluded.

Writes blind_batch_{1..N}.tsv to the scratchpad. Read-only on the gold list.
Run: .venv/bin/python scripts/76_make_blind_batches.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

VD = Path(__file__).resolve().parents[1] / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
OUTDIR = Path(
    "/tmp/claude-1000/-home-teo-Desktop-billerbeck-lab-ssign-package/9cd99a1f-5b90-4b24-956f-27d4a979dd71/scratchpad"
)
N_BATCHES = 6
BLIND_COLS = [
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
    "primary_ref",
    "citation_quote",
]


def main() -> int:
    gold = read_tsv(GOLD)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    # round-robin so each batch spans SS types (no agent gets only T1 or only T3)
    batches: list[list[dict]] = [[] for _ in range(N_BATCHES)]
    for i, r in enumerate(gold):
        batches[i % N_BATCHES].append({c: r.get(c, "") for c in BLIND_COLS})
    for b, rows in enumerate(batches, 1):
        out = OUTDIR / f"blind_batch_{b}.tsv"
        write_tsv(out, BLIND_COLS, rows)
        print(f"{out}  ({len(rows)} rows: {', '.join(r['instance_id'] for r in rows)})")
    print(f"\n{len(gold)} rows -> {N_BATCHES} batches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
