#!/bin/bash
# Run multi-genome test from local filesystem to avoid WSL D: mount drops
export PATH="$HOME/.local/bin:$PATH"
source ~/ssign_test/bin/activate

# Use local copy to avoid D: drive disconnect
cd /tmp/ssign_local
pip install -e . -q 2>&1 | tail -2

# Verify import
python -c "from ssign_app.core.runner import run_cross_genome_orthologs; print('Import OK')"

# Run test
export SSIGN_TEST_OUTDIR="/tmp/ssign_multi_full"
python -u tests/test_multi_genome_full.py 2>&1
echo "EXIT CODE: $?"
