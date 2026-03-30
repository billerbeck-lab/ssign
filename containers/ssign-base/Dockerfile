FROM python:3.11-slim

LABEL maintainer="Billerbeck Lab"
LABEL description="ssign base image — Python dependencies for data processing, figures, and reporting"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    pandas==2.2.3 \
    numpy==1.26.4 \
    scipy==1.14.1 \
    matplotlib==3.10.0 \
    seaborn==0.13.2 \
    biopython==1.84 \
    pyyaml==6.0.2 \
    goatools==1.4.12 \
    obonet==0.3.1 \
    sankeyflow==1.3.0 \
    networkx==3.4.2 \
    requests==2.32.3 \
    jinja2==3.1.4

# Create non-root user
RUN groupadd -r ssign && useradd -r -g ssign ssign
USER ssign

WORKDIR /work
