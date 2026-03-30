#!/bin/bash
cd /mnt/d/ssign_package
source ~/ssign_test/bin/activate

echo "=== DLP TIMING TEST: 50 proteins ==="
START=$(date +%s)

python src/ssign_app/scripts/run_deeplocpro.py \
  --input /tmp/dlp_test50.faa \
  --sample TEST50 \
  --output /tmp/dlp_test50_results.tsv \
  --mode remote

RC=$?
END=$(date +%s)
ELAPSED=$((END - START))

echo ""
echo "Exit code: $RC"
echo "TIMING: ${ELAPSED} seconds for 50 proteins"
echo "Per protein: $(echo "scale=2; $ELAPSED / 50" | bc) seconds"
echo ""
NLINES=$(wc -l < /tmp/dlp_test50_results.tsv 2>/dev/null || echo 0)
echo "Result lines: $NLINES"
head -3 /tmp/dlp_test50_results.tsv 2>/dev/null
