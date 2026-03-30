FROM ubuntu:22.04

LABEL maintainer="Billerbeck Lab"
LABEL description="ssign annotation image — Prodigal, MacSyFinder v2, HMMER, BLAST+, InterProScan"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    prodigal \
    hmmer \
    ncbi-blast+ \
    wget \
    curl \
    procps \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# MacSyFinder v2
RUN pip3 install --no-cache-dir \
    macsyfinder==2.1.4 \
    pandas==2.2.3 \
    biopython==1.84

# Install TXSScan models for MacSyFinder
RUN macsydata install TXSScan || true

# InterProScan (user can bind-mount a local install instead)
# The full install is ~80GB; we provide the wrapper and let users mount their install.
# To use: -v /path/to/interproscan-5.x:/opt/interproscan
ENV INTERPROSCAN_HOME=/opt/interproscan
ENV PATH="${INTERPROSCAN_HOME}:${PATH}"

# Bakta is optional — too large for default image (~30GB database)
# Users who need it should: pip3 install bakta && bakta_db download --output /db/bakta
# Then mount: -v /path/to/bakta_db:/db/bakta

RUN groupadd -r ssign && useradd -r -g ssign ssign
USER ssign

WORKDIR /work
