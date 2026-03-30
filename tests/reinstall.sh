#!/bin/bash
export PATH="$HOME/.local/bin:$PATH"
source ~/ssign_test/bin/activate
cd /mnt/d/ssign_package
pip install -e . -q 2>&1 | tail -5
python -c "from ssign_app.core.runner import run_cross_genome_orthologs; print('Import OK')"
