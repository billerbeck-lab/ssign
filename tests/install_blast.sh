#!/bin/bash
# Download and install BLAST+ without sudo
set -e

BLAST_VER="2.16.0"
BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz"
INSTALL_DIR="$HOME/.local"

cd /tmp
echo "Downloading BLAST+ ${BLAST_VER}..."
wget -q "$BLAST_URL" -O blast.tar.gz || curl -sL "$BLAST_URL" -o blast.tar.gz
echo "Download complete: $(ls -lh blast.tar.gz | awk '{print $5}')"

echo "Extracting..."
tar xzf blast.tar.gz

echo "Installing to ${INSTALL_DIR}/bin/..."
mkdir -p "${INSTALL_DIR}/bin"
cp ncbi-blast-${BLAST_VER}+/bin/blastp "${INSTALL_DIR}/bin/"
cp ncbi-blast-${BLAST_VER}+/bin/makeblastdb "${INSTALL_DIR}/bin/"

# Add to PATH if not already there
if ! echo "$PATH" | grep -q "${INSTALL_DIR}/bin"; then
    echo "export PATH=\"${INSTALL_DIR}/bin:\$PATH\"" >> ~/.bashrc
    echo "Added ${INSTALL_DIR}/bin to PATH in ~/.bashrc"
fi

# Verify
export PATH="${INSTALL_DIR}/bin:$PATH"
echo ""
echo "=== Verification ==="
which blastp
blastp -version | head -1
which makeblastdb
makeblastdb -version | head -1

# Cleanup
rm -rf blast.tar.gz ncbi-blast-${BLAST_VER}+
echo "Done!"
