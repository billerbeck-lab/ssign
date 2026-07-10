## Why

The longevity commitment (a "zero-maintenance" pledge to the supervisor, defensible only if the offline + Docker + Zenodo mitigation stack is built before Oct 2026) needs a **reproducible container image** as one of its legs. A skeleton exists (`containers/Dockerfile` + `containers/README.md`) but has **never been built or run**, and it cannot yet deliver reproducibility: the CUDA base image SHA is a placeholder (all zeros) and the Python dependency versions are unpinned, so two `docker build`s days apart can resolve different stacks. An unpinned image defeats the entire purpose. The full-tier code path is already smoke-validated (run 3256040, 16/16), and the image is independent of the in-flight 74/75-genome scale runs, so the build can proceed now in parallel.

## What Changes

- **Pin for reproducibility**: resolve and lock the CUDA base image digest (replace the placeholder `ARG CUDA_BASE` SHA) and pin the transitive Python dependency tree (a lockfile / constraints file wired into the image's `pip install`), so a rebuild years later resolves the identical stack.
- **First real build + offline smoke-run**: build `ssign:1.0.0` and run the minimal T5aSS golden fixture through it with databases/models bind-mounted, confirming the known substrate call (BIMENO_04457) offline. This is Part B of `docs/development/install_test_runbook.md`.
- **Per-tier mount story**: document exactly which host directories to bind-mount for base vs extended vs full (one image; `extended` and `full` share identical pip deps, so the tier is chosen entirely by the DB/weight volume, not a different image).
- **In-container DTU / EggNOG paths**: verify SignalP/DeepLocPro via bind-mount (`--signalp-path` / `--deeplocpro-path`) and the `--*-mode remote` webserver fallback work inside the container; EggNOG stays host-installed (biopython pin conflict, unchanged).
- **Singularity/Apptainer**: verify the `docker build` → `.sif` conversion path the README describes, for HPC users who cannot run Docker.

## Capabilities

### New Capabilities
- `release-container`: a reproducible, SHA/dependency-pinned container image for ssign that runs the pipeline offline with host-provided databases/models bind-mounted, selectable per install tier, with a documented build + in-container smoke-test procedure and a Singularity conversion path.

### Modified Capabilities
<!-- none: no existing spec's requirements change; this adds the container-release capability and validates the existing skeleton. -->

## Impact

- **Files**: `containers/Dockerfile` (lock CUDA SHA, wire the lockfile into `pip install`), `containers/README.md` (per-tier mount matrix, verified build/run/sif commands), a new dependency lockfile (`uv.lock` or `containers/constraints.txt`), an optional `containers/build_and_test.sh` helper, and `docs/development/install_test_runbook.md` Part B (fill in the verified steps).
- **No pipeline/runtime code changes** and no new Python runtime deps: this packages and validates the existing stack. The pin may surface version conflicts that require adjusting `pyproject.toml` constraints, but not application logic.
- **Build/host requirement**: a CUDA-capable machine to build + GPU-run the image (e.g. the RTX 4070 box); the laptop can prepare the lockfile and Dockerfile edits but not GPU-run the image.
- **Ties to**: the longevity-commitment plan (Docker leg of the Oct 2026 stack) and publication-roadmap step 9 (dependency pinning).
- **Out of scope** (separate work): the CX3 full-tier scale confirmation (in-flight; a runner bug it surfaces would rebuild the now-pinned image cheaply); registry publishing / Zenodo DOI deposit (release-time); bundling DTU weights (licence-blocked, DTU confirmed 2026-05-07).
