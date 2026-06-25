#!/bin/bash
# benchmark-final-validation: re-run the 4 genomes the tier-2 rerun ran on FRAGMENT inputs.
#
# Why: scripts/58 (now fixed) used to read inputs_gb/<acc>.gbff for every unit. For 4 genomes the corpus
# accession is only a plasmid/WGS contig (hlyA 175-CDS plasmid, apxIA 33-CDS contig, ltxA 242-CDS contig,
# Serralysin 254-CDS contig), so MacSyFinder never saw the chromosomal T1SS machinery and the toxins
# could not emit. script 50 had already staged the PARENT whole assemblies in inputs_gb_fullasm/; this
# job runs THOSE, with the identical rerun config, into a separate RUN_TAG so it never clobbers rerun/.
#
# PREREQ (one-time): the 4 full-assembly inputs were NOT in the original rsync. From the LAPTOP repo root:
#   rsync -av \
#     validation_sweeps/benchmark/inputs_gb_fullasm/{NZ_CP031766.1,NZ_CBDBTK010000022.1,NZ_JABJZG010000001.1,NZ_JBCGCZ010000007.1}.gbff \
#     ttr25@login.cx3.hpc.ic.ac.uk:blastp_t5a/ssign/validation_sweeps/benchmark/inputs_gb_fullasm/
#   then on CX3:  cd ~/blastp_t5a/ssign && git pull   (gets this script + batch + the script-58 fix)
#
# Config matches SUBMIT_rerun.sh exactly (genbank + REANNOTATE=1 + T3SS on + whole-genome enrichment +
# annotation on) so the results reconcile 1:1 with rerun/. Results land in
#   ~/runs/benchmark_phase2/rerun_fullasm/<unit>/   (retrieve to rerun_fullasm/<unit>/ on the laptop).
set -eu
REPO="${REPO:-$HOME/blastp_t5a/ssign}"
BENCH="$REPO/validation_sweeps/benchmark"
PBS="$BENCH/scripts/cx3/run_benchmark_batch.pbs"
BATCH="$BENCH/scripts/cx3/rerun_fullasm_batch.txt"
GB="$BENCH/inputs_gb_fullasm"
# 4 genomes x ~30-48 min each ~= 3 h; 8 h leaves margin. Over-requesting walltime only lengthens the
# queue wait (PBS backfills shorter jobs first), it does not buy anything. Override via env if a genome runs long.
WALLTIME="${WALLTIME:-08:00:00}"
qsub -l select=1:ncpus=32:mem=64gb:ngpus=1:gpu_type=RTX6000 -l walltime="$WALLTIME" \
    -v "BATCH_FILE=$BATCH,INPUT_MODE=genbank,REANNOTATE=1,INCLUDE_T3SS=1,WHOLE_GENOME=1,ENRICH=1,ANNOT=1,INPUT_DIR_GB=$GB,SSIGN_VENV=$REPO/.venv,RUN_TAG=rerun_fullasm" \
    "$PBS"
echo "submitted 1 batch job (4 full-assembly genomes). Watch: qstat -u \$USER"
