# ssign v1.0.0 bundle image — skeleton.
#
# This Dockerfile freezes the ssign Python stack against a SHA-pinned CUDA
# base for ≥5-year reproducibility (per the longevity-commitment plan).
# Model weights and reference databases are NOT baked in — they are fetched
# on the host with scripts/fetch_weights.sh and scripts/fetch_databases.sh,
# then bind-mounted at run time. This keeps the image at ~3-5 GB instead of
# 600+ GB and avoids licence redistribution issues with DTU and EggNOG
# (covered conservatively below; see plan addendum E.6 for the toggles
# that flip if their replies arrive positive).
#
# Build:
#   docker build -f containers/Dockerfile -t ssign:1.0.0 .
#
# Run:
#   docker run --gpus all --rm \
#     -v $HOME/.ssign/databases:/home/ssign/.ssign/databases:ro \
#     -v $HOME/.ssign/models:/home/ssign/.ssign/models:ro \
#     -v $PWD:/work \
#     ssign:1.0.0 run /work/genome.gbff --outdir /work/output
#
# See containers/README.md for full usage notes including Singularity.

# ---------------------------------------------------------------------------
# Base image — SHA-pinned for reproducibility
# ---------------------------------------------------------------------------
# nvidia/cuda runtime image: CUDA 12.4 + cuDNN, no compiler toolchain.
# Runtime variant is sufficient because PyTorch ships its own CUDA kernels.
#
# FRAGILE: SHA digest below is a PLACEHOLDER. Resolve and lock at build
# time with `docker pull nvidia/cuda:12.4.1-runtime-ubuntu22.04` followed
# by `docker inspect --format='{{index .RepoDigests 0}}' ...`. Locking
# the digest is what makes the image deterministic across years; without
# it, the same `docker build` command can produce different stacks.
ARG CUDA_BASE=nvidia/cuda:12.4.1-runtime-ubuntu22.04@sha256:0000000000000000000000000000000000000000000000000000000000000000
FROM ${CUDA_BASE}

# SSIGN_VERSION is the *image tag* — the Python package version comes from
# pyproject.toml and is set by the release process. Until pyproject is
# bumped to 1.0.0, `pip show ssign` inside the image will report 0.9.0
# even when this ARG (and the OCI label below) say 1.0.0.
ARG SSIGN_VERSION=1.0.0
ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="ssign" \
      org.opencontainers.image.description="Secretion-system Identification for Gram-negative bacterial genomes" \
      org.opencontainers.image.source="https://github.com/billerbeck-lab/ssign" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.version="${SSIGN_VERSION}"

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# ---------------------------------------------------------------------------
# System packages
# ---------------------------------------------------------------------------
# python3.10 is Ubuntu 22.04's default — also the minimum supported by
# pyproject.toml. build-essential + python3-dev cover source-build fallbacks
# for any pip dep without a manylinux wheel; revisit at build time if no
# such fallback fires (drop ~400 MB). ncbi-blast+ provides the `blastp`
# binary used at full tier; everything else (bakta, eggnog-mapper, …)
# arrives via the [extended] pip extra below.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.10 \
        python3.10-venv \
        python3-pip \
        python3-dev \
        build-essential \
        wget \
        tar \
        unzip \
        git \
        ca-certificates \
        libgomp1 \
        ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Non-root user
# ---------------------------------------------------------------------------
# uid 1000 matches the typical first-user uid on Ubuntu hosts so bind-mounted
# output directories carry the host user's ownership without needing
# `--user $(id -u):$(id -g)` for the common case. Override at run time if
# the host user has a different uid.
RUN groupadd --gid 1000 ssign \
    && useradd --uid 1000 --gid ssign --create-home --shell /bin/bash ssign

# ---------------------------------------------------------------------------
# ssign install
# ---------------------------------------------------------------------------
# Copy the source tree, install with the [extended] extra (matches the
# v1.0.0 default offering: prediction tools + Bakta + EggNOG-mapper +
# pLM-BLAST + PLM-Effector). Install non-editable so the image is frozen
# at build time.
#
# FRAGILE: pyproject.toml dependency versions are NOT yet locked via
# uv.lock. Step 9 of the publication roadmap pins every transitive
# dep; until then, two builds days apart can resolve different versions.
# The skeleton accepts this — locking is the build-time job.
WORKDIR /opt/ssign
COPY --chown=ssign:ssign pyproject.toml README.md LICENSE ./
COPY --chown=ssign:ssign src ./src
COPY --chown=ssign:ssign scripts ./scripts

RUN pip install --upgrade pip setuptools wheel \
    && pip install '.[extended]'

# ---------------------------------------------------------------------------
# Bind-mount points + TXSScan profiles
# ---------------------------------------------------------------------------
# /home/ssign/.ssign/{databases,models} match the default paths in ssign,
# so the user mounts the host-side directories produced by fetch_*.sh
# straight onto them — no env-var or CLI overrides needed. /work is the
# convention for input/output.
RUN mkdir -p /home/ssign/.ssign/databases /home/ssign/.ssign/models /work \
    && chown -R ssign:ssign /home/ssign /work

# Conservative DTU + EggNOG handling (plan addendum E.4 + addendum 10 J.1/J.2):
# - DTU SignalP: NOT bundled. DTU confirmed 2026-05-07 the licence does not
#   permit redistribution (final). Users either bind-mount a local install
#   onto /opt/signalp6/ and pass --signalp-path, or override the new
#   offline-first default with --signalp-mode remote.
# - DTU DeepLocPro: NOT bundled. Awaiting Ole's redistribution reply; if
#   permitted, add `COPY deeplocpro/ /opt/deeplocpro/` here per plan E.6.
#   Same opt-in pattern as SignalP for now.
# - EggNOG database: NOT bundled. The user runs `download_eggnog_data.py`
#   inside `~/.ssign/databases/eggnog/` on the host before launching the
#   container. fetch_databases.sh --tier extended performs that step.

USER ssign
WORKDIR /home/ssign

# Pre-install TXSScan HMM profiles into ~/.macsyfinder/. runner.py would
# otherwise trigger this on first invocation via `macsydata install --user
# TXSScan==1.1.4`, which fails offline. Baking them in keeps the image
# self-contained. Pin matches runner.py:823.
RUN macsydata install --user TXSScan==1.1.4

WORKDIR /work

ENTRYPOINT ["ssign"]
CMD ["--help"]
