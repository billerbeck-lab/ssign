# ssign: Secretion-system Identification for Gram-Negative Bacteria

[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Status: Beta](https://img.shields.io/badge/status-beta-orange)](#roadmap-to-v100)

ssign detects secretion systems in Gram-negative bacterial genomes, identifies
the proteins they secrete, and annotates those proteins with functional and
structural information from the major bioinformatics databases. Built for the
[Billerbeck Lab](https://www.billerbecklab.com/) at Imperial College London.

**Version `0.9.0` (pre-publication baseline).** The next release (`1.0.0`)
is the publication version: fully offline-capable, SHA-pinned Docker image,
Zenodo-DOI'd. See [Roadmap to v1.0.0](#roadmap-to-v100).

---

## Hosted web service: coming soon

A public web service for browser-based genome submission is planned alongside
v1.0.0. Until then, run ssign locally (below) or via Google Colab (see
[`colab/`](colab/), notebook ships with v1.0.0).

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

For an end-to-end walkthrough on a real genome, see
[`docs/tutorials/first_run.md`](docs/tutorials/first_run.md). Full install
instructions (WSL, optional tool extras, dependency management) in
[`docs/how-to/install.md`](docs/how-to/install.md).

---

## What ssign does

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Stage 1: Secretion-System Detection                                    │
│    MacSyFinder v2 + TXSScan models → validated secretion systems        │
│                                                                         │
│  Stage 2: Secreted-Protein Prediction                                   │
│    DeepLocPro + DeepSecE + SignalP + PLM-Effector                       │
│    → candidate proteins, voted on by independent predictors             │
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

Per-stage detail in [`docs/explanation/pipeline_overview.md`](docs/explanation/pipeline_overview.md).
Per-decision rationale in [`docs/explanation/design_decisions.md`](docs/explanation/design_decisions.md).

---

## Supported inputs

| Format        | Extensions              |
| ------------- | ----------------------- |
| Genbank       | `.gbff`, `.gbk`, `.gb`  |
| FASTA contigs | `.fasta`, `.fna`, `.fa` |
| Protein FASTA | `.faa`                  |

---

## Output

Per genome, ssign writes to your `--outdir`:

- `<sample-id>_results.csv`: main results (chunked: secreted proteins, then
  associated systems, then other systems).
- `<sample-id>_results_raw.csv`: every column with no filtering.
- `<sample-id>_summary.txt`: plain-text report and enrichment summary.
- `figures/<sample-id>/*.png`: five summary figures.

In multi-genome batches the GUI also writes combined `ssign_results.csv`,
`ssign_results_raw.csv`, and `ssign_summary.txt` at the outdir root.

Full column reference: [`docs/reference/output_files.md`](docs/reference/output_files.md).

---

## Key parameters

| Parameter             | Default              | Meaning |
| --------------------- | -------------------- | ------- |
| `excluded_systems`    | `Flagellum Tad T3SS` | System types to skip (T3SS excluded by default; DeepSecE unreliable on it). |
| `conf_threshold`      | `0.8`                | DeepLocPro minimum extracellular probability. |
| `proximity_window`    | `3`                  | +/-N genes around each SS component (same contig only). |
| `wholeness_threshold` | `0.8`                | Minimum MacSyFinder completeness to accept a system. |

All configurable in the GUI or via CLI flags. Full parameter reference in
[`docs/reference/cli.md`](docs/reference/cli.md).

---

## Install tiers

ssign ships in three tiers; pick the one matching your storage budget.
Upgrade later by re-running the database fetcher with a new `--tier`.

| Tier         | Disk    | What's included                                                                                   | Install                       |
| ------------ | ------- | ------------------------------------------------------------------------------------------------- | ----------------------------- |
| **base**     | ~17 GB  | Secretion-system detection + secreted-protein prediction (DLP, DSE, SignalP, PLM-E) + Bakta light | `pip install ssign`           |
| **extended** | ~130 GB | base + EggNOG + HH-suite (Pfam + PDB70) + InterProScan + pLM-BLAST                                | `pip install ssign[extended]` |
| **full**     | ~630 GB | extended + BLAST NR + Bakta full DB + HH-suite UniRef30                                           | `pip install ssign[full]`     |

After pip install, fetch the matching database bundle (pulled from pinned
Zenodo DOIs for long-term reproducibility):

```bash
bash scripts/fetch_databases.sh --tier base       # or: extended / full
```

See [`data/README.md`](data/README.md) for per-tier contents.

### pip vs Docker

- **pip** installs into your Python environment; depends on local Python,
  CUDA, and system libraries staying compatible.
- **Docker** (`docker pull billerbeck-lab/ssign:1.0.0`, available from v1.0.0)
  is SHA-pinned and reproducible for 5+ years. Recommended for paper-
  reproducibility and webserver deployments.

Tier is chosen by the database bundle you fetch, not by pip vs Docker; both
work with any tier.

### Cherry-picking individual tools

If none of the three tiers matches your situation, pick individual extras:

```bash
pip install ssign[deepsece]          # just DeepSecE on top of base
pip install ssign[bakta,deepsece]    # combine any pip extras
pip install ssign[dev]               # test + lint dependencies (for contributors)
```

System binaries (BLAST+, HH-suite, InterProScan) are installed separately per
your platform:

```bash
# BLAST+
sudo apt install ncbi-blast+      # Debian/Ubuntu
brew install blast                 # macOS
conda install -c bioconda blast    # cross-platform

# HH-suite (extended + full)
sudo apt install hhsuite
conda install -c bioconda hhsuite

# InterProScan (extended + full): Java, manual install.
# See docs/how-to/install.md for step-by-step instructions.
```

Full platform-specific install guide: [`docs/how-to/install.md`](docs/how-to/install.md).

---

## Roadmap to v1.0.0

v1.0.0 is the publication release. Most pipeline-side work has already landed
on `main`; remaining items are packaging, documentation, and the release
itself.

**Pipeline (mostly landed on main)**

- ✅ Offline-first: DeepLocPro and SignalP run from local installs by
  default. Users without a DTU academic licence can opt into the DTU
  webserver as a fallback (no licence needed on the user's part). NCBI
  remote BLAST, EBI InterProScan webserver, and MPI Toolkit HHpred remote
  modes removed in favour of local binaries.
- ✅ Bakta + EggNOG whole-genome annotation, PLM-Effector first-class
  secretion prediction, pLM-BLAST / ECOD70 substrate annotation.
- ✅ Re-annotate inputs with Bakta by default; `--use-input-annotations`
  preserves curated GenBank annotations.
- ✅ Cross-validation rule: DLP, DSE, and PLM-Effector treated as equal
  secretion predictors (any one flagging means a candidate). SignalP is
  evidence-only. `n_prediction_tools_agreeing` column carried through.
- ✅ Pipeline order: `enrichment_testing` runs before substrate filtering;
  stats filter default ON for ≥10 genomes.

**Packaging and distribution (in flight)**

- **Docker bundle image**: SHA-pinned, reproducible for 5+ years, published
  to Docker Hub / GHCR.
- **Zenodo deposits**: separate DOIs for source code, model weights, and
  database bundle; paper cites all three.
- ✅ Tier-aware database fetcher (`scripts/fetch_databases.sh --tier
  {base,extended,full}`; Zenodo DOIs land at release time).
- **FAIR-compliant repository layout** per
  [FAIR4RS](https://doi.org/10.1038/s41597-022-01710-x) (Barker et al. 2022,
  _Scientific Data_).
- **Diataxis documentation** (tutorials / how-to / reference / explanation).
- **`bio.tools` registration** for FAIR findability.

**Hosted web service (post-publication):** Flask-based submission form, job
queue, results page.

Track progress in [`CHANGELOG.md`](CHANGELOG.md).

---

## Citing ssign

Cite via [`CITATION.cff`](CITATION.cff) or the GitHub tag
`v0.9.0-prerefactor`. Zenodo and paper DOIs will be added here at v1.0.0
release.

---

## Citing the underlying tools

ssign integrates many open-source tools. Please cite any tool your analysis
uses alongside ssign.

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
- **PLM-Effector**: Zheng D. _Briefings in Bioinformatics_. 2026;27(2):bbag143.
  [doi:10.1093/bib/bbag143](https://doi.org/10.1093/bib/bbag143)
- **BLAST+**: Camacho C, Coulouris G, Avagyan V, et al. _BMC Bioinformatics_.
  2009;10:421. [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
- **HH-suite3**: Steinegger M, Meier M, Mirdita M, et al. _BMC Bioinformatics_.
  2019;20:473. [doi:10.1186/s12859-019-3019-7](https://doi.org/10.1186/s12859-019-3019-7)
- **InterProScan 5**: Jones P, Binns D, Chang HY, et al. _Bioinformatics_.
  2014;30(9):1236-1240. [doi:10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031)
- **Bakta**: Schwengers O, Jelonek L, Dieckmann MA, et al. _Microbial Genomics_.
  2021;7(11):000685. [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685)
- **EggNOG-mapper v2**: Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P,
  Huerta-Cepas J. _Molecular Biology and Evolution_. 2021;38(12):5825-5829.
  [doi:10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293)
- **pLM-BLAST**: _(v1.0.0: citation on integration)_
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

- **M. Teo Reid**, Department of Bioengineering, Imperial College London.
  ORCID: [0009-0009-9239-5743](https://orcid.org/0009-0009-9239-5743)
- **Owen Terpstra**, Molecular Microbiology, University of Groningen.
  ORCID: [0000-0002-8767-4061](https://orcid.org/0000-0002-8767-4061)
- **Karan Kumar**, Industrial Systems Biotechnology Research Group, iAMB,
  RWTH Aachen University. ORCID: [0000-0003-0012-8314](https://orcid.org/0000-0003-0012-8314)
- **Dr. Sonja Billerbeck**, Department of Bioengineering, Imperial College
  London. ORCID: [0000-0002-3092-578X](https://orcid.org/0000-0002-3092-578X)

Correspondence: t.reid25@imperial.ac.uk, s.billerbeck@imperial.ac.uk
