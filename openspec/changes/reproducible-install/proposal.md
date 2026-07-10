## Why

ssign orchestrates ~a dozen external tools (MacSyFinder, the DTU predictors, DeepSecE, Bakta, EggNOG, InterProScan, pLM-BLAST, HH-suite, BLAST, plus the torch/transformers stack). Any of them updating can break a run, so a reproducible, version-locked install must be the **norm**, not an afterthought. Today it is not: `pyproject.toml` pins everything with `>=` (0 exact pins), and `pip install ssign` ignores the committed `uv.lock`, so a fresh install drifts to whatever is latest. The earlier Docker skeleton is the wrong default (Linux-only to build, banned on most HPCs, cannot run on Mac). We evaluated the options and chose **Singularity/Apptainer** as the default reproducible artifact and a **base-only conda environment** for Mac.

## What Changes

- **Default = Singularity/Apptainer `.sif`** (Linux + HPC, all tiers). Built from the SHA-pinned base + the frozen `uv.lock`; runs with no root; one archivable file (Zenodo DOI for the paper's long-term-reproducibility promise). This replaces the Docker image as the shipped artifact.
- **Mac = a base-tier-only conda `environment.yml`.** Base enables only `deeplocpro / deepsece / signalp / protparam` (+ core MacSyFinder/pyrodigal/pyhmmer), all Mac-compatible. The Linux-only tools (InterProScan, HH-suite, BLAST, Bakta) are extended/full only, so they never appear on Mac.
- **Windows: unsupported.** Several core tools have no Windows build; Windows users use WSL2, which is just the Linux path.
- **Docker: retired.** The pinning work already done (refreshed `uv.lock`, SHA-pinned CUDA base, "export the lock then pip install" recipe) ports into the Singularity build file; the `Dockerfile` + `build_and_test.sh` are removed.
- The Python stack stays pinned via `uv.lock` (shared basis for both the `.sif` and the conda env).

## Capabilities

### New Capabilities
- `reproducible-install`: version-locked, reproducible install artifacts for ssign — a Singularity/Apptainer image (default; Linux/HPC; all tiers; Zenodo-archivable) and a base-only conda environment (Mac), both built from the pinned dependency set, with documented per-tier database mounting.

### Modified Capabilities
<!-- none: supersedes the (unmerged) docker-release-image change; no existing spec's requirements change. -->

## Impact

- **Files**: `containers/ssign.def` (new Singularity definition), `containers/environment.yml` (new, base-only conda env for Mac), remove `containers/Dockerfile` + `containers/build_and_test.sh`, `containers/README.md` (rewrite for `.sif` + conda), `docs/development/install_test_runbook.md` Part B (Singularity build test), a build/test helper. `uv.lock` already refreshed.
- **No pipeline/runtime code changes**; base already excludes InterProScan and the other Linux-only tools (verified in `TIER_TOOL_DEFAULTS`), so Mac base needs no gating change.
- **Supersedes** the `docker-release-image` change (removed).
- **Out of scope**: Windows-native support; full tier on Mac (InterProScan is Linux-only); bundling DTU weights or databases (licence/size — bind-mounted); publishing the artifact to a registry / minting the Zenodo DOI (release-time).
