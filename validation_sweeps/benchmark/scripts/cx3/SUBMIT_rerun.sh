#!/bin/bash
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
WALLTIME="${WALLTIME:-72:00:00}"
for b in "$BATCHES"/rerun_*.txt; do
    qsub -l select=1:ncpus=32:mem=64gb:ngpus=1:gpu_type=RTX6000 -l walltime="$WALLTIME" \
        -v "BATCH_FILE=$b,INPUT_MODE=genbank,REANNOTATE=1,INCLUDE_T3SS=1,WHOLE_GENOME=1,ENRICH=1,ANNOT=1,INPUT_DIR_GB=$GB,SSIGN_VENV=$REPO/.venv,RUN_TAG=rerun" \
        "$PBS"
done
echo "submitted $(ls "$BATCHES"/rerun_*.txt | wc -l) batch jobs. Watch: qstat -u \$USER"
