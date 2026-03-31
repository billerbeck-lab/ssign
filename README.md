# ssign — Secretion-System Identification for Gram Negatives

**ssign** links established bioinformatic tools into a single pipeline that identifies secretion systems and their associated secreted proteins in gram-negative bacterial genomes. Given one or more genome files, ssign detects which secretion systems are present, predicts which proteins are secreted through them, and optionally annotates those proteins with functional and structural information.

Built for the [Billerbeck Lab](https://www.billerbecklab.com/) at Imperial College London. Version 0.1.0 (Beta).

---

## Installation

ssign is a Python package with a browser-based GUI. Installation takes a few minutes and requires no bioinformatics experience.

### Linux (Ubuntu / Debian)

1. **Open a terminal.** Search for "Terminal" in your app launcher, or press `Ctrl+Alt+T`.

2. **Install Python** (if not already installed). Most Linux distributions include Python, but you may need pip and venv:

   ```bash
   sudo apt install python3 python3-pip python3-venv
   ```

3. **Create a virtual environment.** This creates an isolated Python environment so ssign's dependencies don't conflict with other software on your system:

   ```bash
   python3 -m venv ~/.ssign-env && source ~/.ssign-env/bin/activate
   ```

4. **Install ssign:**

   ```bash
   pip install ssign
   ```

5. **Launch ssign:**

   ```bash
   ssign
   ```

   This opens a browser-based interface where you upload genome files and configure the pipeline.

> **Tip:** Each time you open a new terminal, reactivate the environment before running ssign:
> ```bash
> source ~/.ssign-env/bin/activate && ssign
> ```

### macOS

1. **Open Terminal.** Search for "Terminal" in Spotlight (`Cmd+Space`), or find it in Applications > Utilities.

2. **Install Python** (if not already installed). The recommended way is via [Homebrew](https://brew.sh/):

   ```bash
   brew install python
   ```

3. **Create a virtual environment.** This creates an isolated Python environment so ssign's dependencies don't conflict with other software on your system:

   ```bash
   python3 -m venv ~/.ssign-env && source ~/.ssign-env/bin/activate
   ```

4. **Install ssign:**

   ```bash
   pip install ssign
   ```

5. **Launch ssign:**

   ```bash
   ssign
   ```

   This opens a browser-based interface where you upload genome files and configure the pipeline.

> **Tip:** Each time you open a new terminal, reactivate the environment before running ssign:
> ```bash
> source ~/.ssign-env/bin/activate && ssign
> ```

### Windows (via WSL)

ssign requires the Windows Subsystem for Linux (WSL), which provides a full Linux environment inside Windows.

1. **Open PowerShell as Administrator.** Right-click the Start button and select "Terminal (Admin)" or "PowerShell (Admin)".

2. **Install WSL:**

   ```powershell
   wsl --install
   ```

   This installs Ubuntu by default. The download is approximately 1-2 GB.

3. **Restart your computer** when prompted.

4. **Open the Ubuntu app** from the Start menu. On first launch, it will ask you to create a Linux username and password.

5. **Follow the Linux instructions above** (starting from step 2) inside the Ubuntu terminal.

---

## Optional Dependencies

ssign works out of the box with `pip install ssign`. The following optional tools extend its capabilities.

### DeepSecE — additional secreted protein prediction

DeepSecE is a deep learning model that predicts secretion system effectors. It provides an additional layer of evidence alongside DeepLocPro and SignalP. The download is approximately 7.3 GB (ESM protein language model).

```bash
pip install ssign[deepsece]
```

Same command on all platforms. Install inside your ssign virtual environment.

### Bakta — higher-quality ORF prediction

If your input is raw FASTA contigs (not GenBank), ssign uses Pyrodigal for gene prediction by default. Bakta provides more comprehensive annotation. After installing the Python package, you also need to download its database (~2 GB).

```bash
pip install ssign[bakta]
bakta_db download --output ~/bakta_db --type light
```

Same command on all platforms. Install inside your ssign virtual environment.

### BLAST+ — ortholog grouping across genomes

When analysing multiple genomes, BLAST+ enables ortholog grouping to identify shared secreted proteins. BLAST+ is not pip-installable and must be installed through your system package manager (~200 MB).

**Linux (Ubuntu / Debian):**
```bash
sudo apt install ncbi-blast+
```

**macOS:**
```bash
brew install blast
```

**Conda (cross-platform):**
```bash
conda install -c bioconda blast
```

### pLM-BLAST — protein language model structural search

pLM-BLAST enables sequence-structure similarity searches using protein language model embeddings. It is not pip-installable and requires local setup. See the [pLM-BLAST GitHub repository](https://github.com/labstructbioinf/pLM-BLAST) for installation instructions.

> **Note:** pLM-BLAST and Foldseek are **ssign-power only** features. They require local databases and/or 3D structure files and are not available through the cloud-based ssign GUI.

### Install all pip-installable optional dependencies

```bash
pip install ssign[full]
```

Or combine specific extras: `pip install ssign[deepsece,bakta]`

---

## Pipeline Overview

ssign runs in four stages:

### Stage 1 — Secretion System Identification

Detects secretion systems (T1SS, T2SS, T3SS, T4SS, T5SS, T6SS, T9SS, etc.) using **MacSyFinder v2** with the **TXSScan** models. Only systems meeting a completeness threshold (default 80%) are retained. Flagellum, Tad, and T3SS are excluded by default (configurable).

### Stage 2 — Secreted Protein Identification

Predicts which proteins are secreted through the detected systems using a combination of:

- **DeepLocPro** — deep learning subcellular localisation predictor for prokaryotes (required)
- **DeepSecE** — secretion system effector predictor (optional, ~7.3 GB model)
- **SignalP 6.0** — signal peptide predictor (optional)

Secreted proteins are identified through guilt-by-association: proteins encoded near secretion system components that are predicted to be extracellular are flagged as candidate secreted proteins. Cross-validation between tools reduces false positives.

### Stage 3 — Secreted Protein Annotation (Optional)

Each detected secreted protein can be annotated using any combination of five tools:

| Tool | What it provides | API provider |
|------|-----------------|--------------|
| **BLASTp** | Sequence homology hits | NCBI |
| **HHpred** | Remote homology / domain detection | MPI Bioinformatics Toolkit |
| **InterProScan** | Domain and family annotations | EBI |
| **Foldseek** | Structural similarity search | ssign-power only |
| **ProtParam** | Physicochemical properties | Local (Biopython) |

Each tool can be independently enabled or disabled in the GUI.

> **Note on well-annotated genomes:** If your genome already has thorough functional annotations, you can skip the optional annotation tools to save significant time. The existing annotations from your input file are already included in the results. You can always re-run with annotation tools enabled later using the Resume feature.

### Stage 4 — Generate Data and Figures

Merges all results into output tables and generates publication-ready summary figures.

---

## Supported Input Formats

| Format | Extensions | Notes |
|--------|-----------|-------|
| **GenBank** (recommended) | `.gbff`, `.gbk`, `.gb` | Includes gene annotations — no ORF prediction needed |
| **GFF3** | `.gff`, `.gff3` | Requires associated FASTA sequence |
| **FASTA contigs** | `.fasta`, `.fna`, `.fa` | Requires ORF prediction (Pyrodigal, or optionally Bakta) |

GenBank format is recommended because it contains both the nucleotide sequence and gene annotations, allowing ssign to skip the ORF prediction step entirely.

---

## Output

ssign produces the following output files:

| File | Description |
|------|-------------|
| `ssign_results.csv` | Main results table, chunked into three sections: (1) secreted proteins with annotations, (2) secretion systems with their associated proteins, (3) other detected secretion systems |
| `ssign_results_raw.csv` | Complete unfiltered results for all proteins |
| `ssign_summary.txt` | Plain-text summary of detected systems and secreted proteins |
| `figures/` | Publication-ready figures (system diagrams, annotation summaries) |

---

## Important: API Fragility

Most annotation and prediction tools in ssign run via free cloud APIs provided by DTU, NCBI, MPI Bioinformatics Toolkit, and EBI. While this means no local database setup is needed, there are trade-offs:

- **Speed**: A full run with all annotation tools typically takes **1-3 hours** depending on genome size, number of secreted proteins, and server load.
- **Rate limits**: APIs enforce request limits (e.g., NCBI allows 3 requests/second, MPI Toolkit allows 200/hour). ssign respects these limits automatically.
- **Reliability**: Cloud services occasionally experience downtime. ssign includes retry logic with exponential backoff, but extended outages may require reruns.
- **Reproducibility**: API databases are updated periodically, so results may vary slightly between runs months apart.

If speed and reliability are critical, consider **ssign-power** (see below).

---

## ssign-power

For local or HPC execution with full local databases (no API dependency), install **ssign-power**. This uses Nextflow DSL2 with Docker or Singularity containers and runs all tools locally, including Foldseek and pLM-BLAST which are not available in the standard ssign GUI.

```bash
nextflow run ssign --input genome.gbff --outdir results -profile docker
```

ssign-power requires Nextflow, Java 11+, and Docker or Singularity. It is not needed for normal use — the standard `pip install ssign` workflow covers most use cases.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Excluded systems | `Flagellum, Tad, T3SS` | Secretion system types to exclude from analysis |
| Confidence threshold | `0.8` | Minimum DeepLocPro extracellular probability |
| Proximity window | `3` | Number of genes (+/-) around each SS component to search |
| System completeness | `0.8` | Minimum MacSyFinder wholeness score to accept a system |

All parameters are configurable through the GUI settings panel.

---

## Citations

If you use ssign in your research, please cite the underlying tools:

- **MacSyFinder v2**: Neron B, Denise R, Coluzzi C, Touchon M, Rocha EPC, Abby SS. MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes. *Peer Community Journal*. 2023;3:e28. [doi:10.24072/pcjournal.250](https://doi.org/10.24072/pcjournal.250)

- **TXSScan**: Abby SS, Cury J, Guglielmini J, Neron B, Touchon M, Rocha EPC. Identification of protein secretion systems in bacterial genomes. *Scientific Reports*. 2016;6:23080. [doi:10.1038/srep23080](https://doi.org/10.1038/srep23080)

- **DeepLocPro**: Moreno J, Nielsen H, Winther O, Teufel F. Predicting the subcellular location of prokaryotic proteins with DeepLocPro. *Bioinformatics*. 2024;40(12):btae677. [doi:10.1093/bioinformatics/btae677](https://doi.org/10.1093/bioinformatics/btae677)

- **DeepSecE**: Zhang Y, Guan J, Li C, Wang Z, Deng Z, Gasser RB, Song J, Ou HY. DeepSecE: a deep-learning-based framework for multiclass prediction of secreted proteins in gram-negative bacteria. *Research*. 2023;6:0258. [doi:10.34133/research.0258](https://doi.org/10.34133/research.0258)

- **SignalP 6.0**: Teufel F, Almagro Armenteros JJ, Johansen AR, Gislason MH, Piber SI, Tsirigos KD, Winther O, Brunak S, von Heijne G, Nielsen H. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nature Biotechnology*. 2022;40(7):1023-1025. [doi:10.1038/s41587-021-01156-3](https://doi.org/10.1038/s41587-021-01156-3)

- **BLAST+**: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

- **HH-suite3**: Steinegger M, Meier M, Mirdita M, Vohringer H, Haunsberger SJ, Soding J. HH-suite3 for fast remote homology detection and deep protein annotation. *BMC Bioinformatics*. 2019;20:473. [doi:10.1186/s12859-019-3019-7](https://doi.org/10.1186/s12859-019-3019-7)

- **InterProScan 5**: Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, McWilliam H, Maslen J, Mitchell A, Nuka G, Pesseat S, Quinn AF, Sangrador-Vegas A, Scheremetjew M, Yong SY, Lopez R, Hunter S. InterProScan 5: genome-scale protein function classification. *Bioinformatics*. 2014;30(9):1236-1240. [doi:10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031)

- **Foldseek**: van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Soding J, Steinegger M. Fast and accurate protein structure search with Foldseek. *Nature Biotechnology*. 2024;42(2):243-246. [doi:10.1038/s41587-023-01773-0](https://doi.org/10.1038/s41587-023-01773-0)

- **Pyrodigal**: Larralde M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. *Journal of Open Source Software*. 2022;7(72):4296. [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)

- **Bakta**: Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*. 2021;7(11):000685. [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685)

- **Biopython**: Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422-1423. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

---

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).

Copyright (C) 2026 Billerbeck Lab, Imperial College London.
