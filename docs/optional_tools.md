# Optional Tool Installation Guide

ssign works out of the box with `pip install ssign`. The tools below extend its capabilities and can be installed at any time.

---

## DeepSecE — Additional Secreted Protein Prediction

DeepSecE predicts which secretion system type a protein is secreted by, providing cross-validation against DeepLocPro. Download is approximately **7.3 GB** (PyTorch + ESM protein language model).

**Important:** Install inside your ssign virtual environment.

### All platforms

```bash
# Activate your ssign environment first
source ~/.ssign-env/bin/activate    # Linux/macOS
# or on Windows WSL: source ~/.ssign-env/bin/activate

pip install ssign[deepsece]
```

ssign will automatically detect DeepSecE on next launch.

---

## Bakta — Higher-Quality Gene Prediction

If your input is raw FASTA contigs (not GenBank), ssign uses Pyrodigal by default. Bakta provides richer functional annotation. Requires a **~2 GB** database download.

### All platforms

```bash
source ~/.ssign-env/bin/activate

pip install ssign[bakta]
bakta_db download --output ~/bakta_db --type light
```

Then enable Bakta in the ssign GUI and enter the database path (`~/bakta_db`).

---

## BLAST+ — Ortholog Grouping Across Genomes

When analysing multiple genomes, BLAST+ enables ortholog grouping to identify shared secreted proteins. BLAST+ is **not pip-installable** and must be installed through your system package manager. Download is approximately **200 MB**.

### Linux (Ubuntu / Debian)

```bash
sudo apt install ncbi-blast+
```

### macOS

```bash
brew install blast
```

### Conda (cross-platform)

```bash
conda install -c bioconda blast
```

ssign will automatically detect BLAST+ on next launch.

---

## DeepLocPro (Local Mode) — Faster Localization Prediction

By default, ssign runs DeepLocPro via the free DTU cloud API. Local mode is faster but requires a DTU academic license and approximately **5 GB** of model files. GPU recommended.

### All platforms

1. Register for a free academic license at [DTU DeepLocPro](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)
2. Download and install following DTU's instructions
3. In ssign, select "Local install" for DeepLocPro and enter the install path

---

## SignalP 6.0 (Local Mode) — Faster Signal Peptide Prediction

By default, ssign runs SignalP via the free DTU cloud API. Local mode is faster but requires a DTU academic license and approximately **1 GB** download.

### All platforms

1. Register for a free academic license at [DTU SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
2. Download and install following DTU's instructions
3. In ssign, select "Local install" for SignalP and enter the install path

---

## Install Everything At Once

To install all pip-installable optional dependencies:

```bash
source ~/.ssign-env/bin/activate
pip install ssign[full]
```

This installs DeepSecE + Bakta. BLAST+, DeepLocPro (local), and SignalP (local) still need separate installation as described above.

You can also combine specific extras: `pip install ssign[deepsece,bakta]`
