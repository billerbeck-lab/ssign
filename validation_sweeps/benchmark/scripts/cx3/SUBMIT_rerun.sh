#!/bin/bash
# benchmark-final-validation: submit the 57-genome tier-2 rerun on CX3.
#
# One job per genome (the cx3-submit skill's single-genome pattern), trickled 2 at a time to respect
# the gpu72 ~2-concurrent-job cap. Each job: extended tier, Bakta RE-ANNOTATION (clean slate -- NO
# --use-input-annotations), annotation ON, --enrichment-stats (forces whole-genome DLP/DSE). Per-genome
# wallclock ceiling ~4.4h (pLM-BLAST + EggNOG dominate), so this trickles over ~2-5 days.
#
# Idempotent: records each submitted genome in .rerun_submitted and skips it on re-run, so you can stop
# and restart safely. Run it under tmux so it survives a dropped login:
#   tmux new -s rerun 'bash validation_sweeps/benchmark/scripts/cx3/SUBMIT_rerun.sh'
#
# Optional: export SSIGN_EXTRA_ARGS="--skip-plmblast" before running to drop pLM-BLAST (the 4h pole and
# the weakest annotator); cuts per-genome time ~4x. WALLTIME and REPO are overridable via env.
set -eu

REPO="${REPO:-$HOME/blastp_t5a/ssign}"
GB="$REPO/validation_sweeps/benchmark/inputs_gb"
LIST="$REPO/validation_sweeps/benchmark/data/phase2/rerun_inputs.txt"
SUB="$REPO/scripts/cx3/submit_batched_overnight.sh"
DONE="$REPO/validation_sweeps/benchmark/scripts/cx3/.rerun_submitted"
WALLTIME="${WALLTIME:-12:00:00}"
CAP="${CAP:-2}"

test -f "$LIST" || { echo "FATAL: $LIST not found (rsync inputs + git pull first)"; exit 1; }
test -f "$SUB"  || { echo "FATAL: $SUB not found (git pull the repo on CX3)"; exit 1; }
touch "$DONE"

n=0
total=$(grep -cve '^[[:space:]]*$' "$LIST")
while read -r g; do
    [ -z "$g" ] && continue
    if grep -qxF "$g" "$DONE"; then echo "skip (already submitted): $g"; continue; fi
    test -f "$GB/$g" || { echo "FATAL: input not found: $GB/$g"; exit 1; }
    # block until a job slot frees up (gpu72 hidden ~2-job cap)
    while [ "$(qstat -u "$USER" 2>/dev/null | grep -c ssign_batched)" -ge "$CAP" ]; do sleep 300; done
    n=$((n + 1))
    echo "[$n/$total] submitting $g"
    bash "$SUB" --small --enrichment-stats --walltime "$WALLTIME" "$GB/$g"
    echo "$g" >>"$DONE"
    sleep 30
done <"$LIST"
echo "all $total genomes submitted. Watch: qstat -u \$USER"
