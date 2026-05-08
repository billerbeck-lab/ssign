# `scripts/`: install fetchers and one-off analysis scripts

Scripts here either install ssign's external assets (databases, weights) or
perform single-purpose analyses that produced results or figures in the
ssign paper. They are **not** part of the ssign pipeline itself; the
pipeline code lives in `src/ssign_app/` and is importable from the
installed package.

## When to put a script here vs in `src/ssign_app/`

- **Here (`scripts/`):** one-off paper-specific analysis, figure-regeneration
  wrappers, dataset-curation helpers. Not imported by other code. Runnable
  standalone.
- **In `src/ssign_app/`:** reusable modules, CLI entry points, code that any
  ssign user would invoke as part of the pipeline.

If a script becomes generally useful to multiple analyses, promote it into
`src/ssign_app/` as a module.

## Already shipped

- `fetch_weights.sh`: downloads DeepSecE checkpoint, ProtT5, ESM, and
  PLM-Effector weights from their canonical sources (Zenodo at v1.0.0
  release; upstream mirrors as fallback).
- `fetch_databases.sh`: tier-aware (`--tier base|extended|full`) downloader
  for Bakta, EggNOG, HH-suite, InterProScan, BLAST NR, ECOD70, taxdump.

## To be populated in Phase 8

- Per-figure regeneration scripts referenced from `figures/`.
- Validation-set curation scripts.
- Benchmark driver scripts.

## Running

Every script should be runnable from the repo root:

```bash
cd /path/to/ssign
python scripts/<script>.py [args]
# or
bash scripts/<script>.sh
```

and should document its inputs, outputs, and any required environment at the
top of the file.

## Style

- Python scripts: use the same deps as the main package; no separate venv.
- Shell scripts: POSIX-compatible, `set -euo pipefail` at the top.
- Document non-obvious assumptions in comments.
- Record random seeds for any stochastic analysis.

## Status

Install fetchers (`fetch_weights.sh`, `fetch_databases.sh`) shipped.
Paper-analysis scripts populated as the corresponding analyses are
finalised before publication.
