# ssign: Secretion-system Identification for Gram-Negative Bacteria

[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Status: v1.0.0 pre-release](https://img.shields.io/badge/status-v1.0.0%20pre--release-orange)](#roadmap-to-v100)
[![Tests](https://github.com/billerbeck-lab/ssign/actions/workflows/test.yml/badge.svg)](https://github.com/billerbeck-lab/ssign/actions/workflows/test.yml)
[![Lint](https://github.com/billerbeck-lab/ssign/actions/workflows/lint.yml/badge.svg)](https://github.com/billerbeck-lab/ssign/actions/workflows/lint.yml)

ssign predicted secretion systems in Gram-negative bacterial genomes, 
the proteins they secrete, and annotates those proteins.

**Version `1.0.0`**

## What ssign does

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Stage 1: Secretion-**System** Prediction                               │
│    MacSyFinder v2 + TXSScan models                                      │
│                                                                         │
│  Stage 2: Secreted **Protein** Prediction                               │
│    DeepLocPro + DeepSecE + SignalP                                      │
│                                                                         │
│  Stage 3: Proximity Analysis                                            │
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

| Format        | Extensions              |
| ------------- | ----------------------- |
| GenBank       | `.gbff`, `.gbk`, `.gb`  |
| FASTA contigs | `.fasta`, `.fna`, `.fa` |
| Protein FASTA | `.faa`                  |
| GFF3          | `.gff` (with a same-stem nucleotide FASTA alongside) |

---

## Key parameters

| Parameter             | Default              | Meaning |
| --------------------- | -------------------- | ------- |
| `excluded_systems`    | `Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P` | Surface/uptake appendages skipped by default (not protein secretion systems). T3SS is detected by default; DeepSecE is never trusted for T3SS (it over-calls flagellar proteins), so T3SS relies on MacSyFinder + DeepLocPro + proximity. Add `T3SS` to exclude it entirely. |
| `conf_threshold`      | `0.8`                | DeepLocPro minimum extracellular probability. |
| `proximity_window`    | `3`                  | +/-N genes around each SS component (same contig only). |
| `wholeness_threshold` | `0.8`                | Minimum MacSyFinder completeness to accept a system. |

All configurable via CLI flags. Full reference:
[cli.md](docs/reference/cli.md) · [run.md](docs/how-to/run.md).

---

## Install tiers

ssign ships in three tiers, differing only in which databases you fetch. Pick
the one matching your storage budget; upgrade later by re-running the fetcher
with a new `--tier`.

| Tier         | DB disk | Adds |
| ------------ | ------- | ---- |
| **base**     | ~4 GB   | Secretion-system detection + secreted-protein prediction (DeepLocPro, DeepSecE, SignalP) + Bakta light |
| **extended** | ~100 GB | base + EggNOG + InterProScan + pLM-BLAST |
| **full**     | ~500 GB | extended + Bakta full DB + HH-suite (Pfam + PDB70 + UniRef30) + BLASTp-vs-Swiss-Prot (NR opt-in) |

Add ~18 GB of model weights (ESM2, DeepSecE checkpoint + ESM-1b backbone,
DeepLocPro, ProtT5), downloaded automatically on first run and shared across
tiers. Container users get them pre-baked, so nothing downloads at run time.

After installing, fetch the matching database bundle and verify:

```bash
ssign fetch-databases --tier extended    # or: base / full
ssign doctor --tier extended             # reports any missing DB / binary / weight, with the fix
```
---

`ssign` is a command-line tool. Run `ssign run --help` for every option, or
`ssign doctor` to check your install.

**System requirements:** Linux or macOS, Python >= 3.10. A CUDA-capable GPU is
recommended for the neural predictors (DeepLocPro, DeepSecE) and pLM-BLAST.

Full install instructions (tiers, optional tool extras, HPC) are in the
[install guide](docs/how-to/install.md).

Full container guide (hardware, HPC/PBS, verification, troubleshooting):
[install.md](docs/how-to/install.md). Maintainer build and
publish steps: [containers/README.md](containers/README.md).

---

## Output

Per genome, ssign writes to your `--outdir`:

- `<sample-id>_results.csv`: main results (chunked: secreted proteins, then
  associated systems, then other systems).
- `<sample-id>_results_raw.csv`: every column, no filtering.
- `<sample-id>_summary.txt`: plain-text report and enrichment summary.
- `figures/<sample-id>/*.png`: the curated summary set, numbered `01`-`06`
  (secreted proteins by SS type; size and physicochemical properties; and four
  functional-category figures, COG / KEGG / EggNOG / curated consensus, when
  annotation tools have run), plus the enrichment fold/significance figures
  (per-tool and combined) when `--enrichment-stats` is on.

Multi-genome batches write each genome into its own `<outdir>/<sample-id>/`
subdirectory, plus a combined `combined_results.csv` and `combined_summary.txt`
at the outdir root and the curated figure set over all genomes as `0N_pooled_*`
(figure `01` is the cross-genome overview, one stacked bar per genome). With
BLAST+ installed they also get cross-genome ortholog groups
(`cross_genome_ortholog_groups.csv`) and a conservation figure
(`figures/07_cross_genome_orthologs.png`).

Full column reference: [output_files.md](docs/reference/output_files.md).

---

## Citing ssign

Cite via [CITATION.cff](CITATION.cff). The Zenodo and paper DOIs will be added
here at the v1.0.0 release.

### Citing the underlying tools

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
and verified by the human authors: code changes were
validated against the test (`pytest tests/unit/`) and, for pipeline-affecting changes, against real-genome validation
runs. Scientific claims, default thresholds, and biological rationale were
verified against the primary literature cited in this README and in
[design_decisions.md](docs/explanation/design_decisions.md).

---

## Contributing

See [CONTRIBUTING.md](.github/CONTRIBUTING.md) for how to file issues, propose
features, and submit pull requests. Contributions are welcome, especially for
documentation and new tool integrations. A [Code of Conduct](.github/CODE_OF_CONDUCT.md)
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
