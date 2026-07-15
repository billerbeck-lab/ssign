#!/bin/bash
# Submit ONE PBS job that runs ssign from the CONTAINER IMAGE (the self-contained
# .sif) at a given tier. Validation counterpart to run_container_extended.pbs.
#
# Mirrors submit_batched_overnight.sh but for the container: it builds the qsub in
# an array and runs it FROM THIS SCRIPT, so the long line is never pasted into a
# terminal (paste-mangling is the #1 cause of stuck/garbage container jobs).
#
# Usage (on CX3, from the repo dir ~/blastp_t5a/ssign):
#   bash scripts/cx3/submit_container.sh $HOME/xantho_gbff/GCF_000.gbff
#   bash scripts/cx3/submit_container.sh --sif $EPHEMERAL/ssign_v7.sif --tier extended GENOME.gbff
#   bash scripts/cx3/submit_container.sh --tier full --walltime 24:00:00 GENOME.gbff   # full panel
#   bash scripts/cx3/submit_container.sh --dry-run GENOME.gbff     # print the qsub, do not submit
#
# The genome path MUST be absolute (the PBS input check runs before any cd on the
# compute node). gpu_type=RTX6000 is pinned in the -l select below (and also baked
# into the PBS directive), so the job is placeable however it is submitted.

set -eu

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PBS_SCRIPT="$SCRIPT_DIR/run_container_extended.pbs"
test -f "$PBS_SCRIPT" || { echo "FATAL: $PBS_SCRIPT not found"; exit 1; }

GPU="RTX6000"
WALLTIME="06:00:00"
NCPUS="32"
MEM="64gb"
TIER="extended"
SIF="${SIF:-${EPHEMERAL:-$HOME/ephemeral}/ssign_v7.sif}"
DRY="0"

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --sif) SIF="$2"; shift 2 ;;
        --gpu) GPU="$2"; shift 2 ;;
        --walltime) WALLTIME="$2"; shift 2 ;;
        --tier) TIER="$2"; shift 2 ;;
        --dry-run|--print) DRY="1"; shift ;;
        -h|--help) sed -n '2,18p' "$0" | sed 's/^# \?//'; exit 0 ;;
        --) shift; break ;;
        -*) echo "Unknown flag: $1"; exit 1 ;;
        *) break ;;
    esac
done

[ "$#" -ge 1 ] || { echo "FATAL: pass a genome path (absolute). See --help."; exit 1; }
GENOME="$1"
case "$GENOME" in
    /*) ;;
    *) echo "FATAL: genome path must be absolute (got: $GENOME)"; exit 1 ;;
esac
test -f "$GENOME" || { echo "FATAL: genome not found: $GENOME"; exit 1; }
test -f "$SIF"    || { echo "FATAL: image not found: $SIF (pass --sif)"; exit 1; }

echo "Submitting container ssign job:"
echo "  image:   $SIF"
echo "  genome:  $GENOME"
echo "  tier=$TIER  gpu=$GPU  ncpus=$NCPUS mem=$MEM  walltime=$WALLTIME"
echo

# gpu_type in -l select is MANDATORY on v1_gpu72 (omitting it queues the job
# forever). Building the command in an array keeps the -v string in one token so a
# terminal can't mangle it on paste, and lets --dry-run print exactly what runs.
QSUB_ARGS=(
    -l "select=1:ncpus=${NCPUS}:mem=${MEM}:ngpus=1:gpu_type=$GPU"
    -l "walltime=$WALLTIME"
    -N "ssign_container_${TIER}"
    -v "GENOME=${GENOME},SIF=${SIF},TIER=${TIER}"
    "$PBS_SCRIPT"
)
if [ "$DRY" = "1" ]; then
    printf 'DRY RUN: would submit:\n  qsub'
    printf ' %q' "${QSUB_ARGS[@]}"
    printf '\n'
    exit 0
fi
jid=$(qsub "${QSUB_ARGS[@]}")
echo "Submitted: $jid"
echo
echo "Confirm it is placeable (must show gpu_type=RTX6000):"
echo "  qstat -f $jid | grep Resource_List.select"
echo "Watch:  qstat -u \$USER      (Q becomes R when it starts)"
