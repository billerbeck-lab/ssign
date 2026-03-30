#!/bin/bash
# Full multi-genome test: 2 genomes with HHpred + orthologs + BLAST+
set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"
source ~/ssign_test/bin/activate

OUTDIR="/tmp/ssign_multi_full"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

GENOME1="/mnt/d/ssign_package/tests/data/Xanthobacter_tagetidis_TagT2C_genomic.gbff"
GENOME2="/mnt/d/ssign_package/tests/data/Roseixanthobacter_finlandensis_VTT_E-85241_genomic.gbff"

echo "=== MULTI-GENOME TEST: $(date) ==="
echo "Genome 1: Xanthobacter_tagetidis"
echo "Genome 2: Roseixanthobacter_finlandensis"
echo ""

# Verify BLAST+
echo "=== BLAST+ check ==="
blastp -version | head -1

cd /mnt/d/ssign_package

# Run genome 1
echo ""
echo "========================================="
echo "=== GENOME 1: Xanthobacter ($(date +%H:%M:%S)) ==="
echo "========================================="
START1=$(date +%s)

python -u src/ssign_app/core/runner.py \
    --input "$GENOME1" \
    --sample xanthobacter \
    --outdir "$OUTDIR/xanthobacter" \
    --run-blastp \
    --run-hhsuite \
    --run-interproscan \
    --run-protparam \
    2>&1

END1=$(date +%s)
echo "Genome 1 time: $((END1 - START1))s"

# Run genome 2
echo ""
echo "========================================="
echo "=== GENOME 2: Roseixanthobacter ($(date +%H:%M:%S)) ==="
echo "========================================="
START2=$(date +%s)

python -u src/ssign_app/core/runner.py \
    --input "$GENOME2" \
    --sample roseixanthobacter \
    --outdir "$OUTDIR/roseixanthobacter" \
    --run-blastp \
    --run-hhsuite \
    --run-interproscan \
    --run-protparam \
    2>&1

END2=$(date +%s)
echo "Genome 2 time: $((END2 - START2))s"

# Cross-genome ortholog grouping
echo ""
echo "========================================="
echo "=== CROSS-GENOME ORTHOLOGS ($(date +%H:%M:%S)) ==="
echo "========================================="
START3=$(date +%s)

python -u -c "
import sys
sys.path.insert(0, 'src')
from ssign_app.core.runner import run_cross_genome_orthologs

result = run_cross_genome_orthologs(
    genome_outdirs=['$OUTDIR/xanthobacter', '$OUTDIR/roseixanthobacter'],
    output_dir='$OUTDIR',
    min_pident=40.0,
    min_qcov=70.0,
    progress_callback=lambda s, p, m: print(f'  [{p}%] {s}: {m}'),
)
print()
print('=== CROSS-GENOME RESULTS ===')
for k, v in result.items():
    print(f'  {k}: {v}')
"

END3=$(date +%s)
echo "Cross-genome time: $((END3 - START3))s"

# Final summary
echo ""
echo "========================================="
echo "=== FINAL SUMMARY ==="
echo "========================================="
echo "Genome 1 output:"
ls -la "$OUTDIR/xanthobacter/"*.csv 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
echo ""
echo "Genome 2 output:"
ls -la "$OUTDIR/roseixanthobacter/"*.csv 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
echo ""
echo "Cross-genome output:"
ls -la "$OUTDIR/"*.csv "$OUTDIR/"*.faa 2>/dev/null | awk '{print "  " $NF " (" $5 ")"}'
echo ""

# Check for cross-genome ortholog columns in integrated CSVs
echo "Cross-genome columns in integrated CSVs:"
for d in xanthobacter roseixanthobacter; do
    CSV=$(ls "$OUTDIR/$d/"*integrated*.csv 2>/dev/null | head -1)
    if [ -n "$CSV" ]; then
        COLS=$(head -1 "$CSV" | tr ',' '\n' | grep xg_ || echo "none")
        echo "  $d: $COLS"
    fi
done

TOTAL=$((END3 - START1))
echo ""
echo "Total wall time: ${TOTAL}s ($((TOTAL/60))m)"
echo "=== TEST COMPLETE: $(date) ==="
