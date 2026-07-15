# `scripts/`: helper scripts for installing and running ssign

Small user-facing helpers that sit alongside the pip package. They are **not** part
of the pipeline itself (that lives in `src/ssign_app/` and is importable from the
installed package); these wrap external-asset setup and the container launcher.

| Script | What it does |
|---|---|
| `fetch_databases.sh` | Tier-aware (`--tier base\|extended\|full`) downloader for the reference databases (Bakta, EggNOG, HH-suite, InterProScan, ECOD30, taxdump; BLAST NR opt-in). Model weights auto-download on first run or are baked into the container. |
| `ssign-run` | One-line launcher for the ssign Apptainer/Singularity container: binds databases + SignalP, sets the tier and RAM budget, adds `--nv`, stages the image. |
| `ssign-setup-dtu` | Installs the DTU predictors (DeepLocPro + SignalP 6) locally so ssign runs them on your machine, not the DTU webserver. |

Usage is documented in the how-to guides:
[`install.md`](../docs/how-to/install.md) (native/pip),
[`install_container.md`](../docs/how-to/install_container.md) (container),
[`run_on_hpc.md`](../docs/how-to/run_on_hpc.md) (clusters).

## Style

- Shell: `set -euo pipefail` at the top; document inputs, outputs, and required
  environment in the header comment.
- Python (in the package, not here): same deps as the main package; no separate venv.
