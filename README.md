# ssign — Secretion-system Identification for Gram Negatives

A pipeline for identifying secretion system substrates in gram-negative bacterial genomes.

## Overview

```
Input Genome(s)          ssign Pipeline                      Output
─────────────── ──────────────────────────────── ──────────────────────
                 ┌─────────────────────────────┐
 FASTA contigs ──┤ 1. Input Processing         │
 GenBank/GBFF ──┤    (Prodigal / Bakta / parse)│
 GFF3 + FASTA ──┤                              │
                 ├─────────────────────────────┤
                 │ 2. SS Detection             │ → valid_systems.tsv
                 │    (MacSyFinder v2 + TXSScan)│
                 ├─────────────────────────────┤
                 │ 3. Secreted Protein Predict  │
                 │    (DeepLocPro + DeepSecE    │
                 │     + SignalP 6.0)           │
                 ├─────────────────────────────┤
                 │ 4. Substrate Identification  │ → substrates.tsv
                 │    (proximity + T5SS + stats)│
                 ├─────────────────────────────┤
                 │ 5. Annotation (optional)     │
                 │    BLASTp, HH-suite,         │
                 │    Foldseek, InterProScan,   │
                 │    ProtParam                 │
                 ├─────────────────────────────┤
                 │ 6. Integration & Reporting   │ → master_substrates.csv
                 │    (consensus + figures)     │   ssign_report.html
                 └─────────────────────────────┘   figures/
```

## Two Ways to Run ssign

### Easy Mode (Streamlit GUI + web APIs)

No Nextflow, Docker, or system packages needed. Annotation tools run via
remote web APIs (NCBI, EBI, BioLib, MPI Toolkit). Best for most users.

```bash
# On Ubuntu/Debian (one-time setup)
sudo apt install pipx
pipx install ssign

# On macOS / other systems
pip install ssign
```

Then just run:
```bash
ssign
```

This opens a Streamlit GUI where you upload genomes and configure the pipeline.
DeepSecE is optional (7.3 GB model download, can be enabled in the GUI).

### Power Mode (Nextflow + local databases)

For HPC environments or batch processing of many genomes with local databases.
Requires Nextflow, Java 11+, and Docker or Singularity.

```bash
# Single genome
nextflow run ssign --input genome.gbff --outdir results -profile docker

# Multiple genomes via samplesheet
nextflow run ssign --input samplesheet.csv --outdir results -profile docker

# Skip annotation (faster, identification only)
nextflow run ssign --input genome.gbff --skip_annotation -profile docker

# HPC with PBS scheduler
nextflow run ssign --input samplesheet.csv -profile hpc_pbs
```

## Samplesheet Format

```csv
sample,input_1,input_2
genome_a,/path/to/genome_a.gbff,
genome_b,/path/to/contigs.fasta,
genome_c,/path/to/annotation.gff3,/path/to/contigs.fasta
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Genome file or samplesheet CSV |
| `--outdir` | `./results` | Output directory |
| `--excluded_systems` | `Flagellum,Tad,T3SS` | SS types to exclude |
| `--conf_threshold` | `0.8` | DeepLocPro extracellular probability minimum |
| `--proximity_window` | `3` | +/- N genes per SS component |
| `--skip_annotation` | `false` | Skip all optional annotation tools |
| `--blastp_mode` | `remote` | `local` or `remote` for BLASTp |

See `nextflow.config` for all parameters (Power Mode) or the GUI settings panel (Easy Mode).

## Annotation Tools

Each tool can be independently skipped and supports local/remote modes:

| Tool | Local DB Size | Remote API | Default Mode |
|------|--------------|------------|--------------|
| BLASTp | nr: ~300GB / swissprot: ~1.5GB | NCBI (3 req/sec) | remote |
| HH-suite | Pfam-A ~3GB + PDB70 ~20GB | MPI Toolkit (200/hr) | off |
| InterProScan | ~80GB | EBI (30 req/sec) | remote |
| Foldseek | ~10GB | Web API | local |
| pLM-BLAST | ECOD70 ~5-10GB | broken | off |
| ProtParam | none (BioPython) | N/A | on |

## Prediction Tools

| Tool | Install | License | Easy Mode |
|------|---------|---------|-----------|
| DeepLocPro | pip (pybiolib) or DTU download | Free via BioLib / DTU academic | BioLib remote (no license) |
| DeepSecE | pip (deepsece) | MIT | Optional (7.3 GB ESM model) |
| SignalP 6.0 | pip (pybiolib) or DTU download | Free via BioLib / DTU academic | BioLib remote (no license) |

## License

GPL-3.0. See [LICENSE](LICENSE).

## Citation

If you use ssign, please cite the tools it wraps:
- MacSyFinder: Abby et al., 2014
- DeepLocPro: Thumuluri et al., 2022
- DeepSecE: Zhang et al., 2023
- SignalP 6.0: Teufel et al., 2022
