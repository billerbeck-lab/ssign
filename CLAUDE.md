# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ssign** (Secretion-system Identification for Gram Negatives) — identifies secretion system substrates in gram-negative bacterial genomes. Built for the Billerbeck Lab at Imperial College London. Version 0.1.0 (Beta).

Two execution modes:
- **Easy Mode**: Pure Python via Streamlit GUI. Uses web APIs for annotation tools. Install: `pip install ssign` + `sudo apt install hmmer`. No Nextflow/Docker needed.
- **Power Mode**: Nextflow DSL2 pipeline with Docker/Singularity containers and local databases. For HPC batch processing.

## Running (Easy Mode — primary focus)

```bash
pip install ssign
ssign                     # launches Streamlit GUI
```

## Running (Power Mode)

```bash
nextflow run main.nf --input genome.gbff --outdir results -profile docker
```

## Tests

```bash
pytest tests/unit/ -v                           # unit tests
pytest tests/unit/test_proximity_analysis.py -v  # single test
pytest tests/integration/test_pipeline_e2e.py -v # integration (needs network)
```

## Architecture

### Two parallel orchestration layers (same scripts, different runners)

**Easy Mode** (primary development focus):
```
ssign CLI → Streamlit GUI (Home.py) → PipelineRunner (core/runner.py) → Python scripts → web APIs
```

**Power Mode**:
```
nextflow run main.nf → workflows/*.nf → modules/local/*.nf → bin/ scripts → containers
```

Both modes call the same underlying Python scripts. The scripts live in two mirrored locations:
- `bin/` — used by Nextflow (auto-added to PATH)
- `src/ssign_app/scripts/` — packaged with pip install

### Pipeline Phases (6 steps)

1. **Input Processing** → detect format (GenBank/GFF3/FASTA), extract proteins + gene order
2. **SS Detection** → MacSyFinder v2 + TXSScan → validated secretion systems + components
3. **Secreted Protein Prediction** → DeepLocPro + DeepSecE + SignalP → cross-validated predictions
4. **Substrate Identification** → per-component proximity analysis + T5SS handling + enrichment stats → filtered substrates
5. **Optional Annotation** → BLASTp, HH-suite, Foldseek, InterProScan, ProtParam (each independently skippable)
6. **Integration & Reporting** → merge annotations + HTML report + publication figures

## Key Code Locations

- `src/ssign_app/Home.py` — Streamlit GUI entry point
- `src/ssign_app/core/runner.py` — **PipelineRunner**: pure-Python orchestrator (Easy Mode). Handles step sequencing, resume/checkpoint, progress callbacks, core vs optional failure handling.
- `src/ssign_app/cli.py` — CLI entry point (`ssign` command)
- `src/ssign_app/scripts/` — all pipeline scripts (packaged copy of `bin/`)
- `bin/` — all Python scripts (used by Nextflow Power Mode)
- `bin/ssign_lib/` — shared library (constants, FASTA I/O, PDB utils, GO utils, retry, manifest)
- `bin/ssign_lib/constants.py` — **single source of truth** for thresholds
- `modules/local/` — Nextflow process definitions (Power Mode only)
- `workflows/` — Nextflow subworkflows (Power Mode only)
- `containers/` — Dockerfiles for 5 images (Power Mode only)

## Dependencies (Easy Mode)

- **No system packages required** — everything is pip-installable
- **Python**: streamlit, pandas, numpy, biopython, matplotlib, seaborn, scipy, macsyfinder, pybiolib, pyrodigal, pyhmmer
- **Optional**: `deepsece` (7.3 GB ESM model, user-enabled in GUI)
- DeepLocPro and SignalP run via BioLib remote API (no DTU license needed)

## Critical Bug Fixes (preserve these)

1. **Per-component proximity** (proximity_analysis.py): Uses +/-N genes from each individual SS component, NOT the full system boundary span. The old approach caused ~26 false positives.
2. **DSE cross-genome leakage** (proximity_analysis.py, system_filtering.py): DeepSecE may predict an SS type that doesn't exist in that genome. Must validate via `dse_type_in_genome()`.
3. **pLDDT scale mismatch** (pdb_utils.py): ESMFold outputs 0-1, AlphaFold DB uses 0-100. Auto-detection + normalization prevents ESMFold structures from being filtered out.
4. **BLASTp hit splitting** (run_blastp.py): NCBI concatenates descriptions with " >"; must split before filtering.
5. **DeepSecE T3SS unreliability**: MacSyFinder found 0 T3SS across 74 genomes; DeepSecE predicted 1,808 (mostly flagellar misclassification). T3SS excluded by default.
6. **HH-suite remote mode**: Uses "alignment" parameter, NOT "sequence".
7. **Foldseek metric**: Uses qtmscore (query-normalized), NOT alntmscore.

## Key Parameters

- `conf_threshold = 0.8` — DeepLocPro extracellular probability minimum
- `wholeness_threshold = 0.8` — MacSyFinder system completeness minimum
- `proximity_window = 3` — +/-N genes per SS component
- `required_fraction_correct = 0.8` — fraction of SS components correctly localized
- `excluded_systems = 'Flagellum,Tad,T3SS'` — excluded by default
- Each annotation tool has `skip_*` and `*_mode` (local/remote) config options

## Shims (src/ssign_app/shims/)

`hmmsearch` — pyhmmer-based drop-in for the HMMER3 binary. Registered as a console_script so `pip install ssign` puts `hmmsearch` on PATH. MacSyFinder finds it automatically. Only implements the subset of hmmsearch flags that MacSyFinder uses (`--cpu`, `-o`, `--cut_ga`, `-E`). If it breaks due to MacSyFinder or pyhmmer updates, fallback is `sudo apt install hmmer`.

## Error Handling Pattern

All fragile external dependencies (subprocess calls, web APIs, optional imports) use this pattern:
```python
# FRAGILE: <what can fail and why>
# If this breaks: <how to fix>
```
PipelineRunner.check_dependencies() runs pre-flight checks before the pipeline starts and warns about missing tools. Individual steps produce definitive errors with specific fix instructions.

## Conventions

- Tool wrappers follow `run_<tool>.py` naming; output standardized TSV/CSV
- Scripts in `src/ssign_app/scripts/` and `bin/` should stay in sync
- PipelineRunner (runner.py) is the primary orchestrator for Easy Mode development
- DeepLocPro and SignalP default to remote mode (BioLib) — no DTU license needed
- DeepSecE is optional and off by default in Easy Mode (large model download)
