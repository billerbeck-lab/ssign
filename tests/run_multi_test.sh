#!/bin/bash
export PATH="$HOME/.local/bin:$PATH"
source ~/ssign_test/bin/activate
cd /mnt/d/ssign_package
python -u tests/test_multi_genome_full.py 2>&1
echo "EXIT CODE: $?"
