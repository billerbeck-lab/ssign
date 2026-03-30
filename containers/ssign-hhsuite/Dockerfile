FROM ubuntu:22.04

LABEL maintainer="Billerbeck Lab"
LABEL description="ssign HH-suite image — hhblits, hhsearch with database mount points"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    hhsuite \
    python3 \
    python3-pip \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    requests==2.32.3 \
    biopython==1.84 \
    pandas==2.2.3

# Database mount points — users download and mount these:
# Pfam-A (~3GB):      -v /path/to/pfam:/db/pfam
# PDB70 (~20GB):      -v /path/to/pdb70:/db/pdb70
# UniClust30 (~25GB): -v /path/to/uniclust30:/db/uniclust30
#
# Download instructions:
#   Pfam-A:      https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
#   PDB70:       https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
#   UniClust30:  https://wwwuser.gwdg.de/~compbiol/uniclust/2023_02/

RUN groupadd -r ssign && useradd -r -g ssign ssign
USER ssign

WORKDIR /work
