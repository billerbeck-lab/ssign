# ssign — Secretion-system Identification for Gram-Negative Bacteria

[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Status: Beta](https://img.shields.io/badge/status-beta-orange)](#roadmap-to-v100)

ssign detects secretion systems in Gram-negative bacterial genomes, identifies
the proteins they secrete, and annotates those proteins with functional and
structural information from the major bioinformatics databases. Built for the
[Billerbeck Lab](https://www.billerbecklab.com/) at Imperial College London.

**Version `0.9.0` — pre-publication baseline.** The next release (`1.0.0`)
will be the publication version: fully offline-capable, SHA-pinned Docker
image, Zenodo-DOI'd. See [Roadmap to v1.0.0](#roadmap-to-v100).

---

## 🌐 Hosted web service — coming soon

A public web service that lets you submit a genome in the browser and receive
the full ssign report without installing anything locally is **planned for
release alongside v1.0.0**. If you don't have the hardware or CLI
experience to run ssign yourself, that will be the easiest entry point. This
section will be updated with the URL once the service is live.

In the meantime, you can run ssign locally (below) or via Google Colab (see
[`colab/`](colab/), notebook shipping with v1.0.0).

---

## Quickstart

```bash
# Create an isolated environment
python3 -m venv .venv && source .venv/bin/activate

# Install (from PyPI once v1.0.0 ships; currently from source)
pip install git+https://github.com/billerbeck-lab/ssign.git@v0.9.0-prerefactor

# Launch the GUI
ssign
```

Opens a browser-based interface for uploading genomes and configuring the
pipeline. Command-line mode is also supported (see `ssign --help`).

**System requirements:** Linux or macOS, Python ≥ 3.10. CUDA-capable GPU
recommended for DeepSecE and (in v1.0.0) PLM-Effector.

**Full install instructions** — including WSL on Windows, optional tool
extras (Bakta, DeepSecE, BLAST+), and dependency management — are in
[`docs/optional_tools.md`](docs/optional_tools.md). Will be restructured into
`docs/how-to/install.md` at v1.0.0.

---

## What ssign does

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Stage 1: Secretion-System Detection                                    │
│    MacSyFinder v2 + TXSScan models → validated secretion systems        │
│                                                                         │
│  Stage 2: Secreted-Protein Prediction                                   │
│    DeepLocPro + DeepSecE + SignalP (+ PLM-Effector in v1.0.0)           │
│    → candidate proteins + guilt-by-association with SS component        │
│       neighborhoods                                                     │
│                                                                         │
│  Stage 3: Cross-Validation + Proximity Analysis                         │
│    Per-SS-component ±3-gene window, same contig only                    │
│    → ranked list of candidate secreted proteins                         │
│                                                                         │
│  Stage 4: Optional Functional Annotation                                │
│    BLASTp, HH-suite, InterProScan, Bakta, EggNOG, pLM-BLAST, ProtParam  │
│    → integrated annotations + consensus voting across tools             │
│                                                                         │
│  Stage 5: Integration + Reporting                                       │
│    → HTML report, result tables, summary figures                        │
└─────────────────────────────────────────────────────────────────────────┘
```

Detailed description of each stage, threshold choices, and tool selection
rationale lives in [`docs/`](docs/) (will be restructured per Diataxis at
v1.0.0).

---

## Supported inputs

| Format        | Extensions              |
| ------------- | ----------------------- |
| Genbank       | `.gbff`, `.gbk`, `.gb`  |
| FASTA contigs | `.fasta`, `.fna`, `.fa` |
| Protein FASTA | `.faa`                  |

---

## Output

- `ssign_results.csv` — main results, three sections: (1) secreted proteins
  with annotations, (2) secretion systems with associated proteins, (3)
  other systems detected but without high-confidence substrates.
- `ssign_results_raw.csv` — complete unfiltered per-protein results.
- `ssign_summary.txt` — plain-text summary.
- `figures/` — summary figures produced by the pipeline (system diagrams,
  annotation heatmaps, enrichment plots). These are summary-quality, not
  paper-quality; publication figures are regenerated separately from scripts
  in the top-level `figures/` directory.
- HTML report with embedded interactive tables.

---

## Key parameters

| Parameter             | Default                | Meaning                                                                     |
| --------------------- | ---------------------- | --------------------------------------------------------------------------- |
| `excluded_systems`    | `Flagellum, Tad, T3SS` | System types to skip (T3SS excluded by default — DeepSecE unreliable on it) |
| `conf_threshold`      | `0.8`                  | DeepLocPro minimum extracellular probability                                |
| `proximity_window`    | `3`                    | +/- N genes around each SS component (same contig only)                     |
| `wholeness_threshold` | `0.8`                  | Minimum MacSyFinder completeness to accept a system                         |

All configurable in the GUI or via CLI flags. Full parameter reference in
[`docs/configuration.md`](docs/configuration.md).

---

## Install tiers

ssign ships in three tiers. Pick the one that matches your storage capacity
and use case. You can always upgrade later by running the tier-aware database
fetcher with a new `--tier`.

| Tier         | Disk    | What's included                                                                                   | Install                       |
| ------------ | ------- | ------------------------------------------------------------------------------------------------- | ----------------------------- |
| **base**     | ~17 GB  | Secretion-system detection + secreted-protein prediction (DLP, DSE, SignalP, PLM-E) + Bakta light | `pip install ssign`           |
| **extended** | ~130 GB | base + EggNOG + HH-suite (Pfam + PDB70) + InterProScan + pLM-BLAST                                | `pip install ssign[extended]` |
| **full**     | ~630 GB | extended + BLAST NR + Bakta full DB + HH-suite UniRef30                                           | `pip install ssign[full]`     |

After pip install, download the matching database bundle:

```bash
bash scripts/fetch_databases.sh --tier base       # or: extended / full
```

The fetcher pulls from **pinned Zenodo DOIs** — identical bytes on every
run, forever. See [`data/README.md`](data/README.md) for what each tier
downloads.

### pip vs Docker

Pick one of:

- **pip** (the commands above) — installs into your Python environment. Flexible
  but depends on your system's Python, CUDA, and libraries staying compatible.
- **Docker** (`docker pull billerbeck-lab/ssign:1.0.0`) — frozen SHA-pinned
  environment, guaranteed reproducible for 5+ years. Recommended for paper-
  reproducibility and webserver deployments. Available from v1.0.0 onwards.

Both pip and Docker work with any install tier — the tier is controlled by
which database bundle you fetch, not by which environment you choose.

### Cherry-picking individual tools

If none of the three tiers matches your situation, pick individual extras:

```bash
pip install ssign[deepsece]          # just DeepSecE on top of base
pip install ssign[bakta,deepsece]    # combine any pip extras
pip install ssign[dev]               # test + lint dependencies (for contributors)
```

System binaries (BLAST+, HH-suite, InterProScan) are installed separately
per your platform:

```bash
# BLAST+
sudo apt install ncbi-blast+      # Debian/Ubuntu
brew install blast                 # macOS
conda install -c bioconda blast    # Conda (cross-platform)

# HH-suite (extended + full)
sudo apt install hhsuite
conda install -c bioconda hhsuite

# InterProScan (extended + full) — Java, manual install
# See docs/optional_tools.md for step-by-step instructions.
```

Full platform-specific install guide: [`docs/optional_tools.md`](docs/optional_tools.md)
(will migrate to `docs/how-to/install.md` at v1.0.0).

---

## Roadmap to v1.0.0

v1.0.0 is the publication release and will include:

- **Fully offline operation** — no external API dependencies. All tools run
  from local binaries and databases. The current baseline still uses BioLib
  (DeepLocPro, SignalP), NCBI remote BLAST, EBI InterProScan, and the MPI
  Toolkit HHpred web service — all replaced or made local for v1.0.0.
- **New tools** — Bakta + EggNOG (whole-genome), PLM-Effector (first-class
  prediction), pLM-BLAST / ECOD70 (substrate annotation).
- **Docker bundle image** — SHA-pinned, reproducible for 5+ years, published
  to Docker Hub / GHCR.
- **Zenodo deposits** — separate DOIs for source code, model weights, and
  database bundle. Paper cites all three.
- **FAIR-compliant repository layout** — per the
  [FAIR4RS principles](https://doi.org/10.1038/s41597-022-01710-x) (Barker
  et al. 2022, _Scientific Data_).
- **Public hosted web service** (post-publication) — Flask-based, BLAST-style
  submission form, job queue, results page.

Track progress in [`CHANGELOG.md`](CHANGELOG.md).

---

## Citing ssign

If you use ssign in your research, please cite the software itself in addition
to the underlying tools. At v1.0.0 we will have a Zenodo DOI and, once the
manuscript is published, a paper DOI — both will be listed here. For the
pre-publication baseline, cite via [`CITATION.cff`](CITATION.cff) or the
GitHub tag `v0.9.0-prerefactor`.

---

## Citing the underlying tools

ssign integrates many excellent open-source tools. If your analysis uses a
given tool, please cite it alongside ssign.

<details>
<summary>Full list (click to expand)</summary>

- **MacSyFinder v2**: Neron B, Denise R, Coluzzi C, Touchon M, Rocha EPC,
  Abby SS. _Peer Community Journal_. 2023;3:e28. [doi:10.24072/pcjournal.250](https://doi.org/10.24072/pcjournal.250)
- **TXSScan**: Abby SS, Cury J, Guglielmini J, Neron B, Touchon M, Rocha EPC.
  _Scientific Reports_. 2016;6:23080. [doi:10.1038/srep23080](https://doi.org/10.1038/srep23080)
- **DeepLocPro**: Moreno J, Nielsen H, Winther O, Teufel F. _Bioinformatics_.
  2024;40(12):btae677. [doi:10.1093/bioinformatics/btae677](https://doi.org/10.1093/bioinformatics/btae677)
- **DeepSecE**: Zhang Y, Guan J, Li C, Wang Z, Deng Z, Gasser RB, Song J,
  Ou HY. _Research_. 2023;6:0258. [doi:10.34133/research.0258](https://doi.org/10.34133/research.0258)
- **SignalP 6.0**: Teufel F, Almagro Armenteros JJ, Johansen AR, et al.
  _Nature Biotechnology_. 2022;40(7):1023-1025. [doi:10.1038/s41587-021-01156-3](https://doi.org/10.1038/s41587-021-01156-3)
- **PLM-Effector**: _(v1.0.0 — citation on integration)_
- **BLAST+**: Camacho C, Coulouris G, Avagyan V, et al. _BMC Bioinformatics_.
  2009;10:421. [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
- **HH-suite3**: Steinegger M, Meier M, Mirdita M, et al. _BMC Bioinformatics_.
  2019;20:473. [doi:10.1186/s12859-019-3019-7](https://doi.org/10.1186/s12859-019-3019-7)
- **InterProScan 5**: Jones P, Binns D, Chang HY, et al. _Bioinformatics_.
  2014;30(9):1236-1240. [doi:10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031)
- **Bakta**: Schwengers O, Jelonek L, Dieckmann MA, et al. _Microbial Genomics_.
  2021;7(11):000685. [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685)
- **EggNOG-mapper**: _(v1.0.0 — citation on integration)_
- **pLM-BLAST**: _(v1.0.0 — citation on integration)_
- **Pyrodigal**: Larralde M. _JOSS_. 2022;7(72):4296. [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
- **Biopython**: Cock PJA, Antao T, Chang JT, et al. _Bioinformatics_.
  2009;25(11):1422-1423. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

</details>

---

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to file issues, propose
features, and submit pull requests. Contributions welcome, especially for
documentation and new tool integrations.

---

## License

ssign is distributed under the **GNU General Public License v3.0 or later**
(GPL-3.0-or-later). See [`LICENSE`](LICENSE).

---

## Authors

- **M. Teo Reid** — Department of Bioengineering, Imperial
  College London. ORCID: [0009-0009-9239-5743](https://orcid.org/0009-0009-9239-5743)
- **Owen Terpstra** — Molecular Microbiology, University of Groningen.
  ORCID: [0000-0002-8767-4061](https://orcid.org/0000-0002-8767-4061)
- **Karan Kumar** — Industrial Systems Biotechnology Research Group, iAMB,
  RWTH Aachen University. ORCID: [0000-0003-0012-8314](https://orcid.org/0000-0003-0012-8314)
- **Dr. Sonja Billerbeck** — Department of Bioengineering,
  Imperial College London. ORCID: [0000-0002-3092-578X](https://orcid.org/0000-0002-3092-578X)

Correspondence: t.reid25@imperial.ac.uk, s.billerbeck@imperial.ac.uk
