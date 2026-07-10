## 1. Pin dependencies (reproducibility)

- [x] 1.1 Lock format = `uv.lock` (the repo already standardises on it — `research_software_audit.md` cites `uv.lock` as the pin mechanism, and a committed `uv.lock` exists).
- [x] 1.2 Refreshed `uv.lock` (was **stale** vs pyproject: `uv lock --check` failed; `uv lock` added accelerate/fairscale/h5py/sentencepiece — a build from the stale lock would have shipped a **broken extended image** missing those). uv resolution is cross-platform/deterministic, so refreshed on the laptop; the CUDA-host build (2.1) validates the resolved set actually installs+runs.
- [x] 1.3 Wired the lock into `containers/Dockerfile`: `uv export --frozen --extra extended --no-hashes --no-emit-project` → `pip install -r`, then `pip install --no-deps .` for ssign itself.
- [x] 1.4 CUDA base digest locked: `sha256:517da2300c184c9999ec203c2665244bdebd3578d12fcc7065e83667932643d9` (`nvidia/cuda:12.4.1-runtime-ubuntu22.04`, resolved 2026-07-10; re-pin procedure in the Dockerfile comment).

## 2. First build + offline smoke run

- [ ] 2.1 **[CUDA host]** Build the image: `containers/build_and_test.sh` (or `docker build -f containers/Dockerfile -t ssign:1.0.0 .`); fix any dependency-conflict fallout from the pin (adjust `pyproject.toml` constraints only, not app logic — then `uv lock` + rebuild).
- [ ] 2.2 **[CUDA host]** Run the minimal T5aSS golden fixture in-container offline (DeepLocPro bind-mounted, all other tools off, `--network none`); confirm it reports substrate BIMENO_04457. `build_and_test.sh <deeplocpro-dir>` does this.
- [x] 2.3 Added `containers/build_and_test.sh` (uv-lock check → `docker build` → offline golden run asserting BIMENO_04457).

## 3. Per-tier mount matrix

- [ ] 3.1 Document in `containers/README.md` the exact host directories to bind-mount for base / extended / full (one image; tier = which `.ssign/databases` + `.ssign/models` subset is mounted). Include a worked `docker run` per tier.
- [ ] 3.2 Verify a base-tier run (mount only base DBs) and an extended-tier run (larger DB mount) both work against the single image.

## 4. Licence-restricted tools in-container

- [ ] 4.1 Verify SignalP/DeepLocPro via bind-mounted host install (`--signalp-path` / `--deeplocpro-path`) run offline in-container.
- [ ] 4.2 Verify the DTU webserver fallback (`--signalp-mode remote` / `--deeplocpro-mode remote`) works from inside the container (network on).
- [ ] 4.3 Confirm the EggNOG host-install note is accurate (biopython pin conflict) and the `.ssign/databases/eggnog/` mount is picked up when present.

## 5. HPC (Singularity/Apptainer)

- [ ] 5.1 Convert the image to `.sif` and run the golden fixture on the target HPC (CX3) with DBs bind-mounted; confirm parity with the Docker run. Record the Apptainer version.
- [ ] 5.2 Update `containers/README.md` with the verified `.sif` build + run commands.

## 6. Runbook + close-out

- [ ] 6.1 Fill in `docs/development/install_test_runbook.md` Part B with the verified build/run/sif steps (replace runbook-only placeholders with the commands that worked; log any deviation).
- [ ] 6.2 Record in the change: resolved CUDA digest, lockfile hash, image size, and the golden-run result. Note remaining release-time steps (registry push, Zenodo DOI) as out of scope.
