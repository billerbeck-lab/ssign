# Installation

## Easy Mode (Recommended for most users)

The simplest way to run ssign. Uses a Streamlit GUI and remote web APIs for
annotation tools. No Nextflow, Docker, or large database downloads needed.

### Requirements

- Python >= 3.10

### Install

**Ubuntu / Debian** (22.04+, Debian 12+):

Modern Ubuntu/Debian prevents direct `pip install` into the system Python
(PEP 668). Use pipx, which auto-manages an isolated environment:

```bash
sudo apt install pipx       # one-time setup
pipx install ssign
ssign
```

**macOS** (Homebrew Python):

```bash
pip install ssign
ssign
```

**Any system with conda**:

```bash
conda create -n ssign python=3.12
conda activate ssign
pip install ssign
ssign
```

**Developer install** (editable mode — source changes take effect immediately):

```bash
sudo apt install python3-venv    # if on Ubuntu/Debian
python3 -m venv ~/.ssign-env
source ~/.ssign-env/bin/activate
pip install -e /path/to/ssign_package
ssign
```

No system packages needed beyond Python — ssign ships its own HMMER implementation
via pyhmmer. The GUI walks you through uploading genomes and configuring the pipeline.

> **Note:** If you encounter issues with the bundled hmmsearch shim, you can
> install the original HMMER as a fallback: `sudo apt install hmmer`

### What runs where (Easy Mode)

| Tool | How it runs | License |
|------|------------|---------|
| MacSyFinder | Locally (pip-installed) | GPL |
| HMMER (hmmsearch) | Locally via pyhmmer shim (pip-installed) | BSD |
| Prodigal | Locally via pyrodigal (pip-installed) | GPL |
| DeepLocPro | BioLib remote API | Free (no DTU license) |
| SignalP 6.0 | BioLib remote API | Free (no DTU license) |
| DeepSecE | Locally (optional, 7.3 GB model) | MIT |
| BLASTp | NCBI web API (3 req/sec) | Free |
| InterProScan | EBI web API (30 req/sec) | Free |
| HH-suite | MPI Toolkit API (200 jobs/hr) | Free |
| ProtParam | Locally via BioPython | BSD |

### Optional: DeepSecE

DeepSecE predicts which secretion system type a protein is secreted by.
It requires downloading a ~7.3 GB ESM language model on first run.
Enable it in the GUI settings — it is off by default.

```bash
pip install deepsece    # if not already installed with ssign[full]
```

---

## Power Mode (HPC / batch processing)

For processing many genomes with local databases. Requires Nextflow + containers.

### Requirements

- [Nextflow](https://nextflow.io/) >= 23.04
- [Java](https://adoptium.net/) >= 11 (required by Nextflow)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/  # or /usr/local/bin/

# Verify
nextflow -version
```

### Pipeline Installation

```bash
git clone https://github.com/billerbeck-lab/ssign.git
cd ssign
```

### Container Images

Docker images are pulled automatically on first run. To pre-pull:

```bash
docker pull ghcr.io/billerbeck-lab/ssign-base:latest
docker pull ghcr.io/billerbeck-lab/ssign-annotation:latest
docker pull ghcr.io/billerbeck-lab/ssign-structure:latest
docker pull ghcr.io/billerbeck-lab/ssign-hhsuite:latest
docker pull ghcr.io/billerbeck-lab/ssign-plmblast:latest
```

For Singularity (HPC):

```bash
singularity pull ssign-base.sif docker://ghcr.io/billerbeck-lab/ssign-base:latest
singularity pull ssign-annotation.sif docker://ghcr.io/billerbeck-lab/ssign-annotation:latest
# ... etc
```

### DTU Tools (local mode only)

For Power Mode with local DeepLocPro or SignalP installs (optional — remote mode
works without these):

**DeepLocPro**
1. Apply for license: https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/
2. Download and install following DTU instructions
3. Provide path: `--deeplocpro_path /path/to/deeplocpro`

**SignalP 6.0**
1. Apply for license: https://services.healthtech.dtu.dk/services/SignalP-6.0/
2. Download and install following DTU instructions
3. Provide path: `--signalp_path /path/to/signalp6`

### Optional Databases (local annotation)

| Database | Size | Required for |
|----------|------|-------------|
| Swiss-Prot (BLAST) | ~1.5 GB | BLASTp local (recommended) |
| nr (BLAST) | ~300 GB | BLASTp local (comprehensive) |
| Pfam-A (HH-suite) | ~3 GB | HH-suite local |
| PDB70 (HH-suite) | ~20 GB | HH-suite local |
| UniClust30 | ~25 GB | HH-suite local |
| InterProScan | ~80 GB | InterProScan local |
| Foldseek AF/SP | ~10 GB | Foldseek |
| ECOD70 | ~5-10 GB | pLM-BLAST |
| Bakta DB | ~30 GB | Bakta (only for raw contigs) |

**Minimum install (remote annotation):** ~2 GB (pipeline + containers, no databases)
**Recommended install:** ~15 GB (+ Swiss-Prot + Foldseek)
**Full install:** ~450+ GB (all databases)
