# `scripts/`: helper scripts for installing and running ssign

Small user-facing helpers that sit alongside the pip package.

| Script | What it does |
|---|---|
| `fetch_databases.sh` | Tier-aware (`--tier base\|extended\|full`) downloader for the reference databases. Model weights auto-download on first run or are baked into the container. |
| `ssign-run` | One-line launcher for the ssign Apptainer/Singularity container. |
| `ssign-setup-dtu` | Installs the DTU predictors (DeepLocPro + SignalP 6). |

Usage is documented in the how-to guides:
[`install.md`](../docs/how-to/install.md) (container + native install) and
[`run.md`](../docs/how-to/run.md) (running ssign).
