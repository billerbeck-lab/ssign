#!/usr/bin/env python3
"""Split the 57-genome rerun panel into CX3 batch-jobs for the established run_benchmark_batch.pbs.

The benchmark runs a BATCH of genomes serially in one job (run_benchmark_batch.pbs), and you submit
several batch-jobs at once -- PBS schedules up to your concurrent-job cap and queues the rest. This
writes:
  rerun_batches/rerun_NN.txt : the 57 rerun unit-ids, split into N_BATCHES files
  SUBMIT_rerun.sh            : qsubs every batch with the rerun config (clean-slate Bakta + T3SS on
                               + whole-genome enrichment + annotation on)

Unit-ids come from data/phase2/rerun_inputs.txt (one .gbff per distinct genome). N_BATCHES defaults to
12 so you get 12 jobs in flight; PBS runs as many concurrently as the cap allows.

Run (locally, before git push):  <repo>/.venv/bin/python scripts/cx3/make_rerun_batches.py
"""

from __future__ import annotations

from pathlib import Path

BENCH = Path(__file__).resolve().parents[2]
INPUTS = BENCH / "data" / "phase2" / "rerun_inputs.txt"
CXDIR = BENCH / "scripts" / "cx3"
BATCHDIR = CXDIR / "rerun_batches"
N_BATCHES = 12


def main() -> int:
    units = [line[:-5] if line.endswith(".gbff") else line for line in INPUTS.read_text().split()]
    units = [u for u in units if u]
    BATCHDIR.mkdir(parents=True, exist_ok=True)
    for old in BATCHDIR.glob("rerun_*.txt"):
        old.unlink()

    # even split into N_BATCHES (round-robin keeps batch sizes within 1 of each other)
    batches: list[list[str]] = [[] for _ in range(N_BATCHES)]
    for i, u in enumerate(units):
        batches[i % N_BATCHES].append(u)
    batches = [b for b in batches if b]

    names = []
    for n, batch in enumerate(batches, 1):
        name = f"rerun_{n:02d}.txt"
        (BATCHDIR / name).write_text("\n".join(batch) + "\n")
        names.append(name)

    submit = """#!/bin/bash
# benchmark-final-validation rerun: submit every batch-job at once. PBS runs up to your concurrent-job
# cap and queues the rest, so you do NOT trickle or babysit. Config (clean slate, all genomes identical):
#   genbank input + REANNOTATE=1 -> Bakta re-annotates from scratch (NO --use-input-annotations)
#   INCLUDE_T3SS=1  -> T3SS detected   |   ENRICH=1 + WHOLE_GENOME=1 -> whole-genome DLP/DSE + stats
#   ANNOT=1         -> InterProScan + EggNOG + pLM-BLAST + ProtParam on
# Results land in ~/runs/benchmark_phase2/rerun/<unit>/.  Override REPO/WALLTIME via env if needed.
set -eu
REPO="${REPO:-$HOME/blastp_t5a/ssign}"
BENCH="$REPO/validation_sweeps/benchmark"
PBS="$BENCH/scripts/cx3/run_benchmark_batch.pbs"
BATCHES="$BENCH/scripts/cx3/rerun_batches"
GB="$BENCH/inputs_gb"
# each batch = 4-5 genomes x ~30-48 min ~= 3-4 h; 8 h leaves margin. Over-requesting walltime only
# lengthens the queue wait (PBS backfills shorter jobs first). Override via env if a batch runs long.
WALLTIME="${WALLTIME:-08:00:00}"
for b in "$BATCHES"/rerun_*.txt; do
    qsub -l select=1:ncpus=32:mem=64gb:ngpus=1:gpu_type=RTX6000 -l walltime="$WALLTIME" \\
        -v "BATCH_FILE=$b,INPUT_MODE=genbank,REANNOTATE=1,INCLUDE_T3SS=1,WHOLE_GENOME=1,ENRICH=1,ANNOT=1,INPUT_DIR_GB=$GB,SSIGN_VENV=$REPO/.venv,RUN_TAG=rerun" \\
        "$PBS"
done
echo "submitted $(ls "$BATCHES"/rerun_*.txt | wc -l) batch jobs. Watch: qstat -u \\$USER"
"""
    (CXDIR / "SUBMIT_rerun.sh").write_text(submit)

    print(f"wrote {len(names)} batches to {BATCHDIR.relative_to(BENCH)}/ + SUBMIT_rerun.sh")
    print(f"  {len(units)} genomes across {len(names)} batches (sizes: {sorted(len(b) for b in batches)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
