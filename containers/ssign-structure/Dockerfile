FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

LABEL maintainer="Billerbeck Lab"
LABEL description="ssign structure image — ESMFold (GPU), Foldseek"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    wget \
    procps \
    && rm -rf /var/lib/apt/lists/*

# PyTorch with CUDA 11.8 support
RUN pip3 install --no-cache-dir \
    torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118

# ESMFold dependencies
RUN pip3 install --no-cache-dir \
    fair-esm==2.0.0 \
    transformers==4.36.2 \
    biopython==1.84 \
    pandas==2.2.3 \
    scipy==1.14.1

# Foldseek — download static Linux binary
RUN wget -q https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz \
    && tar xzf foldseek-linux-avx2.tar.gz \
    && mv foldseek/bin/foldseek /usr/local/bin/ \
    && rm -rf foldseek foldseek-linux-avx2.tar.gz

# ESMFold model weights are downloaded on first run (~2GB)
# To pre-cache: mount a volume at /root/.cache/huggingface

RUN groupadd -r ssign && useradd -r -g ssign ssign

WORKDIR /work
