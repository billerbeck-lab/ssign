## 1. Shared dependency pin (carried over)

- [x] 1.1 Lock format = `uv.lock` (repo standard); refreshed the stale lock (added the 4 missing `[extended]` deps). Cross-platform, deterministic.
- [x] 1.2 CUDA base digest resolved + locked: `nvidia/cuda:12.4.1-runtime-ubuntu22.04@sha256:517da2300c184c9999ec203c2665244bdebd3578d12fcc7065e83667932643d9` (reused as the Singularity bootstrap base).

## 2. Singularity/Apptainer image (default)

- [x] 2.1 Wrote `containers/ssign.def`: `Bootstrap: docker` from the **digest-only** pinned CUDA base (Apptainer rejects `tag@digest`); `%post` mirrors the Dockerfile recipe (`python3 -m pip install ... uv` → `uv export --locked --extra extended --no-hashes --no-emit-project` → `pip install -r` → `pip install --no-deps .`); bakes TXSScan under a fixed `HOME=/opt/ssign_home`; creates `.ssign/{databases,models}` + `/work`.
- [x] 2.2 **[this laptop, 2026-07-10]** Built unprivileged (`apptainer build --fakeroot`, apptainer 1.5.2), **5.0 GB `.sif`**. Gotcha found + worked around: the build sandbox + torch overflow a RAM-backed tmpfs `/tmp` (`[Errno 122] disk quota exceeded`) — must point `APPTAINER_TMPDIR`/`CACHEDIR`/`TMPDIR` at real disk (the `build_sif.sh` helper does this). No pyproject change needed; the pinned stack built clean.
- [x] 2.3 In-container smoke (laptop, offline, `--containall`): `ssign --version`, `doctor --tier base` = **18/18 Python packages OK**, MacSyFinder 2.1.6 + hmmsearch shim on PATH, **TXSScan 1.1.4 baked + discoverable**, HOME resolves, DTU tools correctly optional. **End-to-end on CX3 (2026-07-10): 17/17 steps PASSED offline** on the T5aSS fixture — MacSyFinder/TXSScan + DeepLocPro (from the mounted conda env) + all downstream through report/figures. Runtime findings baked into the run recipe / to-doc:
  - **DTU tools**: mount the host conda env at its **real absolute path** (`-B $ENV:$ENV`) so the tool's bundled Python/libs resolve; point `--deeplocpro-path $ENV/bin`. Confirmed working.
  - **ESM weights**: DeepLocPro/DeepSecE download the ~2.5 GB ESM2 model at runtime; air-gapped nodes can't fetch it. Mount the host torch cache to the container's cache path (`-B $HOME/.cache/torch:/opt/ssign_home/.cache/torch`). For the archival image, **bake the ESM weights in** (self-contained).
  - **Tier marker**: the image has none, so it defaults to `extended` and expects host-provided tools (EggNOG etc.) → **container runs must pass `--tier` explicitly**, or bake a base marker / change the default.
  - The BIMENO_04457 substrate call differs laptop-vs-CX3 (borderline DeepLocPro localization, 0.60 extracellular) — a DeepLocPro-version difference, not a container issue; the container reproduces CX3-native output (parity check pending).
- [x] 2.4 Added `containers/build_sif.sh` (uv-lock check → real-disk-tmpdir `apptainer build` → `--version`/`doctor` smoke → optional golden run).

## 3. macOS base install

- [ ] 3.1 **Reconsidered — conda may be unnecessary.** Every base-tier tool is pip-installable (MacSyFinder, pyhmmer, pyrodigal, DeepSecE, torch, biopython; DeepLocPro/SignalP are DTU, host-provided either way). Base has **no** Linux-only binary that conda would lock, so macOS base = a **pinned `pip install`** (lock-derived), simpler than a conda env. Decision pending Teo's OK; if yes, ship a base requirements lock instead of `environment.yml`.
- [ ] 3.2 Provide + verify the chosen macOS base install (pinned pip requirements, or conda env) resolves on osx-64 + osx-arm64 with no Linux-only tool.
- [ ] 3.3 (User, on a Mac) install + run the golden fixture at base tier.

## 4. Retire Docker

- [x] 4.1 Removed `containers/Dockerfile` + `build_and_test.sh` and the `docker-release-image` change (superseded).
- [x] 4.2 Rewrote `containers/README.md` for the `.sif`: build-off-cluster + real-disk tmpdir, the CX3-validated run recipe (DTU conda-env mount at real path, ESM cache mount, `--tier`, `--nv/--writable-tmpfs/--containall`), per-tier DB mount matrix, DTU webserver fallback, HPC build-elsewhere-run-here, macOS base = pinned pip, Windows→WSL2.

## 5. Databases, DTU tools, HPC

- [ ] 5.1 Document the per-tier host directories to bind-mount (base/extended/full) against the single `.sif`; verify a base and an extended run.
- [ ] 5.2 Verify the DTU bind-mount (`--signalp-path`/`--deeplocpro-path`) and `--*-mode remote` fallback inside the `.sif`.
- [ ] 5.3 **[CX3]** Run the `.sif` on the target HPC with `--nv` + DBs bind-mounted; confirm parity with a native run. Record the Apptainer version.

## 6. Docs + close-out

- [ ] 6.1 Fill `docs/development/install_test_runbook.md` Part B with the verified Singularity build/run steps; add the Windows→WSL2 note.
- [ ] 6.2 Record the resolved base digest, lock state, `.sif` size, and golden-run result. Note release-time steps (Zenodo DOI deposit of the `.sif`) as out of scope.
