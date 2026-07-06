#!/bin/bash
# Submit ONE PBS job that runs ssign in batched multi-genome mode
# (pooled predictions + annotations across N genomes via MultiGenomeRunner).
#
# Usage:
#   # Default — 2-genome smoke test (K-12 + PAO1) on RTX6000:
#   bash scripts/cx3/submit_batched_overnight.sh
#
#   # All 4 tutorial genomes (K-12 + PAO1 + Vc + Salm):
#   bash scripts/cx3/submit_batched_overnight.sh --tutorial-all
#
#   # Specific genomes:
#   bash scripts/cx3/submit_batched_overnight.sh /path/g1.gbff /path/g2.gbff
#
#   # Override GPU type:
#   bash scripts/cx3/submit_batched_overnight.sh --gpu L40S
#
# The 2-genome default is sized for fast verification runs — K-12 and
# PAO1 differ enough (4314 vs 5680 proteins, different SS profiles) to
# exercise per-genome / pool / split routing without paying the
# annotation cost of 4 genomes.
#
# RTX6000 is the default because it places reliably on Imperial CX3
# v1_gpu72 at the 64-core / 120-GB spec — L40S and others have failed
# the placement-set check at that size as of 2026-06-04.

set -eu

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PBS_SCRIPT="$SCRIPT_DIR/run_batched_multi.pbs"
test -f "$PBS_SCRIPT" || { echo "FATAL: $PBS_SCRIPT not found"; exit 1; }

GPU="RTX6000"
# Per-tool timeouts now scale with input size (openspec size-aware-tool-timeouts): a large
# multi-genome --enrichment-stats run pools whole-genome predictions (e.g. ~160k proteins ->
# ~11h DLP cap), so the PBS walltime is the binding limit. Pass --walltime 24:00:00 (or more)
# for a full panel; the 8h default only fits a handful of genomes.
WALLTIME="8:00:00"
USE_TUTORIAL_ALL="0"
NCPUS="64"
MEM="120gb"
DRY="0"
EXTRA="${SSIGN_EXTRA_ARGS:-}"
# Install tier the run targets. full adds BLASTp-vs-NR + HH-suite-vs-UniRef30
# (the PBS wrapper exports SSIGN_BLAST_NR + SSIGN_HHSUITE_UNICLUST for those).
# Default extended keeps every existing submit unchanged.
TIER="extended"

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --gpu) GPU="$2"; shift 2 ;;
        --walltime) WALLTIME="$2"; shift 2 ;;
        --tier) TIER="$2"; shift 2 ;;
        --tutorial-all) USE_TUTORIAL_ALL="1"; shift ;;
        # 32c/64gb places more easily than the 64c/120gb default when the big
        # nodes are busy; enough for the extended tier on one genome at a time.
        --small) NCPUS="32"; MEM="64gb"; shift ;;
        # Turn on the circular-shift enrichment test (+ per-type/pooled figures).
        --enrichment-stats) EXTRA="${EXTRA:+$EXTRA }--enrichment-stats"; shift ;;
        # Print the qsub command instead of submitting (verify before you fire).
        --dry-run|--print) DRY="1"; shift ;;
        -h|--help) sed -n '2,25p' "$0" | sed 's/^# \?//'; exit 0 ;;
        --) shift; break ;;
        -*) echo "Unknown flag: $1"; exit 1 ;;
        *) break ;;
    esac
done

if [ "$#" -gt 0 ]; then
    GENOMES=("$@")
elif [ "$USE_TUTORIAL_ALL" = "1" ]; then
    GENOMES=(
        "$HOME/ssign-tutorial/ecoli_k12.gbff"
        "$HOME/ssign-tutorial/pseudomonas_pao1.gbff"
        "$HOME/ssign-tutorial/salmonella_typhimurium_lt2.gbff"
        "$HOME/ssign-tutorial/vibrio_cholerae_n16961.gbff"
    )
else
    GENOMES=(
        "$HOME/ssign-tutorial/ecoli_k12.gbff"
        "$HOME/ssign-tutorial/pseudomonas_pao1.gbff"
    )
fi

for g in "${GENOMES[@]}"; do
    test -f "$g" || { echo "FATAL: genome not found: $g"; exit 1; }
done

# Build colon-separated INPUT_GBFFS in one safe assignment — no terminal
# wrap during the qsub line gets to mangle a long quoted string.
GBFFS=""
for g in "${GENOMES[@]}"; do
    GBFFS="${GBFFS}${GBFFS:+:}${g}"
done

echo "Submitting batched ssign job:"
echo "  GPU: $GPU"
echo "  walltime: $WALLTIME"
echo "  ${#GENOMES[@]} genomes:"
for g in "${GENOMES[@]}"; do echo "    $g"; done
echo

# gpu_type is MANDATORY on v1_gpu72 — omitting it queues the job forever
# ("Placement set is too small"). Building the whole command in an array keeps
# the long -v string in one token so a terminal can't mangle it on paste, and
# lets --dry-run print exactly what would be submitted.
QSUB_ARGS=(
    -l "select=1:ncpus=${NCPUS}:mem=${MEM}:ngpus=1:gpu_type=$GPU"
    -l "walltime=$WALLTIME"
    -N "ssign_batched_${#GENOMES[@]}genomes"
    -v "INPUT_GBFFS=${GBFFS},GPU_TYPE=${GPU},TIER=${TIER},SSIGN_EXTRA_ARGS=${EXTRA}"
    "$PBS_SCRIPT"
)
echo "  ncpus=$NCPUS mem=$MEM  tier=$TIER  extra-args: ${EXTRA:-(none)}"
echo
if [ "$DRY" = "1" ]; then
    printf 'DRY RUN — would submit:\n  qsub'
    printf ' %q' "${QSUB_ARGS[@]}"
    printf '\n'
    exit 0
fi
jid=$(qsub "${QSUB_ARGS[@]}")
echo "Submitted: $jid"
echo
echo "Watch with:"
echo "  qstat -u \$USER"
echo "Live log (after it starts):"
echo "  tail -f \"\$(ls -td \$HOME/runs/batched_* | head -1)/ssign.run.log\""
