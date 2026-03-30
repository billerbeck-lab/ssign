FROM python:3.11-slim

LABEL maintainer="Billerbeck Lab"
LABEL description="ssign pLM-BLAST image — protein language model BLAST against ECOD70"

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

# PyTorch CPU (pLM-BLAST doesn't need GPU)
RUN pip install --no-cache-dir \
    torch==2.1.2 --index-url https://download.pytorch.org/whl/cpu

# pLM-BLAST from GitHub
RUN pip install --no-cache-dir \
    git+https://github.com/labstructbioinf/pLM-BLAST.git \
    pandas==2.2.3 \
    biopython==1.84 \
    requests==2.32.3

# ECOD70 database mount point (~5-10GB):
#   -v /path/to/ecod70:/db/ecod70
# NOTE: pLM-BLAST ECOD70 on MPI server is BROKEN. Use local mode only.

RUN groupadd -r ssign && useradd -r -g ssign ssign
USER ssign

WORKDIR /work
