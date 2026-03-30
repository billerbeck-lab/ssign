#!/bin/bash
cd /mnt/d/ssign_package
source ~/ssign_test/bin/activate

echo "=== DLP TIMING TEST: 5 proteins ==="
START=$(date +%s)

python src/ssign_app/scripts/run_deeplocpro.py \
  --input /tmp/dlp_test5.faa \
  --sample TEST \
  --output /tmp/dlp_test5_results.tsv \
  --mode remote

RC=$?
END=$(date +%s)
ELAPSED=$((END - START))

echo ""
echo "Exit code: $RC"
echo "TIMING: ${ELAPSED} seconds for 5 proteins"
echo "Per protein: $(echo "scale=1; $ELAPSED / 5" | bc) seconds"
echo ""
echo "=== Results ==="
cat /tmp/dlp_test5_results.tsv 2>/dev/null || echo "No results file"
