# ssign: Secretion-system Identification for Gram-Negative Bacteria

[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Status: v1.0.0 released](https://img.shields.io/badge/status-v1.0.0%20released-brightgreen)](CHANGELOG.md)
[![Tests](https://github.com/billerbeck-lab/ssign/actions/workflows/test.yml/badge.svg)](https://github.com/billerbeck-lab/ssign/actions/workflows/test.yml)
[![Lint](https://github.com/billerbeck-lab/ssign/actions/workflows/lint.yml/badge.svg)](https://github.com/billerbeck-lab/ssign/actions/workflows/lint.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21441317.svg)](https://doi.org/10.5281/zenodo.21441317)

ssign predicts secretion systems in Gram-negative bacterial genomes, 
the proteins they secrete, and annotates those proteins.

**Version `1.0.0`**

## What ssign does

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Stage 1: Secretion-System Prediction                                   │
│    MacSyFinder v2 + TXSScan models                                      │
│                                                                         │
│  Stage 2: Secreted Protein Prediction                                   │
│    DeepLocPro + DeepSecE + SignalP                                      │
│                                                                         │
│  Stage 3: Proximity Analysis & Optional Statistical Analysis            │
│    Per-SS-component ± gene window                                       │
│                                                                         │
│  Stage 4: Optional Functional Annotation                                │
│    BLASTp, HH-suite, InterProScan, Bakta, EggNOG, pLM-BLAST, ProtParam  │
│                                                                         │
│  Stage 5: Integration + Reporting                                       │
│      HTML report, result tables, summary figures                        │
└─────────────────────────────────────────────────────────────────────────┘
```

Per-stage detail in [pipeline_overview.md](docs/explanation/pipeline_overview.md).
Per-decision rationale in [design_decisions.md](docs/explanation/design_decisions.md).

---

## Supported inputs

`.gbff`, `.gbk`, `.gb`, `.fasta`, `.fna`, `.fa`, `.faa`, `.gff`

## Output

Per genome, ssign writes to your `--outdir`:

- `<sample-id>_results.csv`: main results (chunked: secreted proteins, then
  associated systems, then other systems).
- `<sample-id>_results_raw.csv`: raw outputs of all tools run.
- `<sample-id>_summary.txt`: plain-text report and enrichment summary.
- `figures/<sample-id>/*.png`: summary set of overview figures & enrichment fold/significance figures when `--enrichment-stats` is on.

Multi-genome batches write each genome into its own `<outdir>/<sample-id>/`
subdirectory, plus a combined `combined_results.csv` and `combined_summary.txt`
at the outdir root.

Full column reference: [output_files.md](docs/reference/output_files.md).

---

## Key parameters

| Parameter             | Default              | Meaning |
| --------------------- | -------------------- | ------- |
| `--excluded_systems`    | `Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P` | Surface/uptake appendages skipped by default. |
| `--conf_threshold`      | `0.8`                | DeepLocPro minimum extracellular probability to be considered a secreted protein. |
| `--proximity_window`    | `3`                  | +/-N genes around each SS component. |
| `--wholeness_threshold` | `0.8`                | Minimum MacSyFinder completeness to accept a system. |
| `--enrichment-stats`  | `off`                | Increases compute time to have on but tells you secreted protein enrichment. |

All configurable via CLI flags. Full reference:
[cli.md](docs/reference/cli.md) · [run.md](docs/how-to/run.md).

`ssign` is a command-line tool. Run `ssign run --help` for every option and `ssign doctor --tier (your tier)` for verification and troubleshooting.

---

## Install tiers

ssign's image itself is ~20gb and then you add one of three tiers, differing in which databases you fetch and which tools are automatically run. Pick the one matching your storage budget and compute access; upgrade later by re-running the fetcher
with a new `--tier`.

| Tier         | DB disk | Adds |
| ------------ | ------- | ---- |
| **base**     | ~4 GB   | Secretion-system detection + secreted-protein prediction (DeepLocPro, DeepSecE, SignalP) + Bakta light |
| **extended** | ~100 GB | base + EggNOG + InterProScan + pLM-BLAST |
| **full**     | ~500 GB | extended + Bakta full DB + HH-suite (Pfam + PDB70 + UniRef30) + BLASTp-vs-Swiss-Prot |

**System requirements:** Linux (base, extended, full) or macOS (base only), Python >= 3.10. A CUDA-capable GPU is
recommended for the neural predictors (DeepLocPro, DeepSecE) and pLM-BLAST.

Full install instructions are in the
[install guide](docs/how-to/install.md).

Maintainer build and publish steps: [containers/README.md](containers/README.md).

---

## Citing

Cite via [CITATION.cff](CITATION.cff). The container image is archived on Zenodo:
[doi:10.5281/zenodo.21441317](https://doi.org/10.5281/zenodo.21441317).

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
- **pLM-BLAST**: Kaminski K, Ludwiczak J, Pawlicki K, et al. _Bioinformatics_.
  2023;39(10):btad579. [doi:10.1093/bioinformatics/btad579](https://doi.org/10.1093/bioinformatics/btad579)
- **Pyrodigal**: Larralde M. _JOSS_. 2022;7(72):4296. [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
- **Biopython**: Cock PJA, Antao T, Chang JT, et al. _Bioinformatics_.
  2009;25(11):1422-1423. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

</details>

---

## AI usage disclosure

Parts of the ssign codebase, documentation, and tests were drafted with the
assistance of large language models. All AI-assisted output was reviewed, edited,
and verified by human authors. Code changes were
validated against the test (`pytest tests/unit/`) and, for pipeline-affecting changes, against real-genome validation
runs. Scientific claims, default thresholds, and biological rationale were
verified against the primary literature, see
[design_decisions.md](docs/explanation/design_decisions.md).

---

## Contributing

See [CONTRIBUTING.md](.github/CONTRIBUTING.md) for how to file issues, propose
features, and submit pull requests. A [Code of Conduct](.github/CODE_OF_CONDUCT.md)
applies to all project spaces; security issues should be reported privately per
[SECURITY.md](.github/SECURITY.md).

---

## License

ssign is distributed under the **GNU General Public License v3.0 or later**
(GPL-3.0-or-later). See [LICENSE](LICENSE). Note that the container image is for
non-commercial research use (it bundles DeepLocPro under CC BY-NC-SA 4.0); see
[licensing.md](docs/explanation/licensing.md) for the full per-component
breakdown.

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
