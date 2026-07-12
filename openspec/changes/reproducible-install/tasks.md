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

- [x] 3.1 **Decided: pinned pip, no conda** (Teo OK'd). Base is all-pip, so macOS base = `pip install -c containers/requirements-base.lock.txt ssign`. Shipped `containers/requirements-base.lock.txt` (108 deps, `uv export` core-only from `uv.lock`; no extended/Linux-only tools). README updated.
- [ ] 3.2 (Needs a Mac) verify the base lock resolves on osx-64 + osx-arm64 (uv export carries cross-platform markers, so it should; confirm on Amine's Mac).
- [ ] 3.3 (User, on a Mac) install + run the golden fixture at base tier.

## 4. Retire Docker

- [x] 4.1 Removed `containers/Dockerfile` + `build_and_test.sh` and the `docker-release-image` change (superseded).
- [x] 4.2 Rewrote `containers/README.md` for the `.sif`: build-off-cluster + real-disk tmpdir, the CX3-validated run recipe (DTU conda-env mount at real path, ESM cache mount, `--tier`, `--nv/--writable-tmpfs/--containall`), per-tier DB mount matrix, DTU webserver fallback, HPC build-elsewhere-run-here, macOS base = pinned pip, Windows→WSL2.

## 5. Databases, DTU tools, HPC

- [x] 5.1 **base** per-tier mounts documented + verified on CX3 (17/17). **extended/full image gaps closed in `ssign.def` (2026-07-10 rebuild):** (a) **Java 11** baked (`openjdk-11-jre-headless`, matches native `Java/11.0.25`) so InterProScan's `interproscan.sh` finds `java` under `--containall`; (b) **pLM-BLAST baked** — pinned clone `4dddea3` at `/opt/pLM-BLAST`, `SSIGN_PLMBLAST_SCRIPT` defaulted in `%environment`. Key finding: `run_plm_blast.py` runs it as `sys.executable plmblast.py`, so it uses ssign's own Python; its full runtime stack (torch/transformers/**numba/scikit-learn**, all already in `uv.lock`) needs no extra pip install, and the ECOD DB stays host-mounted. Host-provided at run time (documented in README extended example): EggNOG-mapper conda env + `EGGNOG_DATA_DIR` (prepend its bin via `APPTAINERENV_PREPEND_PATH` so `emapper.py` resolves under `--containall`), InterProScan install via `SSIGN_INTERPROSCAN_PATH`, ECOD via `SSIGN_ECOD_DB`. **Two container-only read-only traps found + fixed while smoke-testing the rebuilt image (2026-07-10):** (i) pLM-BLAST's `alntools` uses numba `@njit(cache=True)`, which writes its cache next to the read-only source → baked `NUMBA_CACHE_DIR=/tmp/ssign_numba_cache` in `%environment` (verified `import alntools` now succeeds under `--containall`); (ii) InterProScan writes scratch under its install dir → README mounts `$DBROOT` **rw** (matches native). Image is 5.1 GB (Java + pLM-BLAST added ~0.1 GB), Java 11.0.31 + `plmblast.py`/`embeddings.py`/`alntools` all present. **Remaining: validate the extended run on one real Xanthobacter genome on CX3 (task 5.3) — the minimal T5aSS fixture yields 0 substrates on CX3's DeepLocPro, so annotation tools would no-op; use a substrate-yielding genome.**
- [ ] 5.2 Verify the DTU bind-mount (`--signalp-path`/`--deeplocpro-path`) and `--*-mode remote` fallback inside the `.sif`.
- [~] 5.3 **[CX3]** Run the `.sif` on the target HPC with `--nv` + DBs bind-mounted; confirm parity with a native run. Record the Apptainer version. **First extended run (job 3268527, 2026-07-10) FAILED at step 2 (extract_proteins): Bakta shells out to external binaries (tRNAscan-SE, aragorn, INFERNAL, diamond, real HMMER, amrfinder) that are NOT pip-installable and NOT baked** (base smoke used `--use-input-annotations`, so Bakta was never exercised in-container before). Downstream tools (DLP/IPS/pLM-BLAST/EggNOG) were skipped-after-core-fail, so still unproven. Fix (no image rebuild): `run_container_extended.pbs` now mounts the host `~/.conda/envs/bakta-deps` and PREPENDs it (real toolchain, same as native; its real HMMER also covers MacSyFinder). **Third container-only gap** after numba + IPS-temp. Rerun pending. **For the self-contained archival image, Bakta's toolchain must be baked** (real HMMER + ~6 bioconda binaries + amrfinder via micromamba) — deferred alongside the ESM bake; both are "mount now, bake for archival".

## 6. Docs + close-out

- [ ] 6.1 Fill `docs/development/install_test_runbook.md` Part B with the verified Singularity build/run steps; add the Windows→WSL2 note.
- [ ] 6.2 Record the resolved base digest, lock state, `.sif` size, and golden-run result. Note release-time steps (Zenodo DOI deposit of the `.sif`) as out of scope.
