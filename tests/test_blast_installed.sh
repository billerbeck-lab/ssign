#!/bin/bash
export PATH="$HOME/.local/bin:$PATH"
source ~/ssign_test/bin/activate
which blastp
blastp -version | head -1
which makeblastdb
echo "BLAST+ ready"
