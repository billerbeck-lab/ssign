## 1. Shared dependency pin (carried over)

- [x] 1.1 Lock format = `uv.lock` (repo standard); refreshed the stale lock (added the 4 missing `[extended]` deps). Cross-platform, deterministic.
- [x] 1.2 CUDA base digest resolved + locked: `nvidia/cuda:12.4.1-runtime-ubuntu22.04@sha256:517da2300c184c9999ec203c2665244bdebd3578d12fcc7065e83667932643d9` (reused as the Singularity bootstrap base).

## 2. Singularity/Apptainer image (default)

- [ ] 2.1 Write `containers/ssign.def`: `Bootstrap: docker` from the pinned CUDA base; `%post` mirrors the retired Dockerfile recipe (`python3 -m pip install ... uv` → `uv export --locked --extra extended --no-hashes --no-emit-project` → `pip install -r` → `pip install --no-deps .`), bake TXSScan profiles, create the `.ssign/{databases,models}` + `/work` mount points.
- [ ] 2.2 **[this laptop, Linux]** Build unprivileged: `apptainer build ssign.sif containers/ssign.def` (install apptainer first: `sudo pacman -S apptainer`); fix any pin/build fallout (adjust pyproject constraints only, then `uv lock` + rebuild).
- [ ] 2.3 **[this laptop]** Offline golden run: `apptainer run --no-home --containall <mounts> ssign.sif run <fixture> ...` with DeepLocPro bind-mounted, all other tools off; confirm substrate BIMENO_04457. (DeepLocPro isn't on this laptop → if unavailable here, do the container-starts + `ssign --version` smoke here and defer the golden run to a DLP box.)
- [ ] 2.4 Add a build/test helper capturing the `uv lock` check → `apptainer build` → golden run.

## 3. Base-only conda environment (macOS)

- [ ] 3.1 Write `containers/environment.yml`: base-tier deps only (DeepLocPro/SignalP host-provided or the DTU channel, DeepSecE/ProtParam/MacSyFinder/pyrodigal/pyhmmer/torch), pinned to the `uv.lock` versions. Explicitly exclude InterProScan / HH-suite / BLAST+ / Bakta.
- [ ] 3.2 Verify the env resolves for osx-64 + osx-arm64 and contains no Linux-only tool (dry-run solve).
- [ ] 3.3 (User, on a Mac) create the env and run the golden fixture at base tier → BIMENO_04457.

## 4. Retire Docker

- [ ] 4.1 Remove `containers/Dockerfile` and `containers/build_and_test.sh`; remove the `docker-release-image` openspec change (superseded).
- [ ] 4.2 Rewrite `containers/README.md` for the `.sif` (build + run + per-tier mount matrix) and the Mac conda env; drop Docker instructions.

## 5. Databases, DTU tools, HPC

- [ ] 5.1 Document the per-tier host directories to bind-mount (base/extended/full) against the single `.sif`; verify a base and an extended run.
- [ ] 5.2 Verify the DTU bind-mount (`--signalp-path`/`--deeplocpro-path`) and `--*-mode remote` fallback inside the `.sif`.
- [ ] 5.3 **[CX3]** Run the `.sif` on the target HPC with `--nv` + DBs bind-mounted; confirm parity with a native run. Record the Apptainer version.

## 6. Docs + close-out

- [ ] 6.1 Fill `docs/development/install_test_runbook.md` Part B with the verified Singularity build/run steps; add the Windows→WSL2 note.
- [ ] 6.2 Record the resolved base digest, lock state, `.sif` size, and golden-run result. Note release-time steps (Zenodo DOI deposit of the `.sif`) as out of scope.
