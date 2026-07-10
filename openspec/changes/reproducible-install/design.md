## Context

`pyproject.toml` declares every dependency with `>=` (0 exact pins, 27 ranged, 10 unpinned), and `pip install ssign` does not read the committed `uv.lock`, so the default install drifts. A Docker skeleton existed but Docker is the wrong default: Linux-only to build, root to build, banned on most HPCs (needs root to run), and it cannot run on macOS at all. We compared pip-lockfile / conda / Singularity across install friction, lock completeness, platform reach, and whether they lock the non-Python binaries. Decision (with the user): Singularity default + base-only conda for Mac, Docker dropped.

Confirmed in code: `TIER_TOOL_DEFAULTS` (constants.py) enables `deeplocpro/deepsece/signalp/protparam` at base; `eggnog/interproscan/plmblast` are extended, `blastp/hhsuite` full. So base carries no Linux-only tool and is Mac-safe.

## Goals / Non-Goals

**Goals:**
- A version-locked default: a Singularity/Apptainer `.sif` built from a SHA-pinned base + frozen `uv.lock`, running the pipeline offline with host-mounted databases; rebuild resolves the identical stack.
- A base-only conda `environment.yml` so Mac users get a reproducible base install without any Linux-only tool.
- One archivable artifact (the `.sif`) for the paper's long-term reproducibility (Zenodo DOI) and locked-down HPCs.

**Non-Goals:**
- Windows-native support (WSL2 = the Linux path).
- Full tier on Mac (InterProScan is Linux-only; full tier is an HPC/Linux job by design).
- Bundling DTU weights or databases (licence / 500 GB — bind-mounted).
- Registry publishing / Zenodo DOI minting (release-time).

## Decisions

1. **Singularity/Apptainer is the default, Docker is dropped.** Apptainer runs with no root, is preinstalled on essentially every HPC, and produces a single archivable `.sif`. It builds from the same SHA-pinned CUDA base via `Bootstrap: docker` (no local Docker daemon needed) and runs the same "export `uv.lock` → pip install" recipe. *Alternative:* keep Docker as the build source and convert. Rejected: the target users (HPC) can't run Docker, and building it needs Docker installed.
2. **Conda is Mac-only and base-only.** Base's four tools are all Mac-compatible; the Linux-only tools never appear at base. *Alternative:* extended on Mac via conda. Rejected: InterProScan has no Mac/Windows build and extended/full are HPC-scale jobs anyway.
3. **`uv.lock` is the shared pin.** Both artifacts consume the same frozen Python versions (the `.sif` via `uv export`, the conda env by pinning the same versions). The lock was stale and is refreshed.
4. **Windows via WSL2 only**, documented, no separate artifact.
5. **Unprivileged `.sif` build.** Apptainer 1.x builds without root (user namespaces / `--fakeroot`), so the image can be built + tested on an ordinary Linux box (including this laptop).

## Risks / Trade-offs

- **Apptainer version drift on target HPCs** (older Singularity) → verify the `.sif` on the real HPC (CX3) once; record the version.
- **Conda reproducibility is version- not build-exact** unless pinned to explicit builds → pin exact versions; accept host-glibc variance (still far better than unpinned pip).
- **DTU tools + databases stay outside both artifacts** (licence/size) → no container fully "locks everything"; document the expected DTU-tool version + DB versions separately.
- **`.sif` is large + Linux-only** → that's why Mac gets conda; the `.sif` is for Linux/HPC + archival.
- **GPU**: the `.sif` uses the host NVIDIA driver via `--nv`; the Mac conda env is CPU/MPS torch (no CUDA on Mac).

## Migration Plan

- Additive + a removal: add `containers/ssign.def` and `containers/environment.yml`, remove `containers/Dockerfile` + `containers/build_and_test.sh`, rewrite `containers/README.md`. Supersede/remove the `docker-release-image` change.
- Validation gate: the `.sif` builds from the pinned base + lock, and the in-container golden run reports the known T5aSS substrate (BIMENO_04457) offline; the Mac conda env resolves and runs the base pipeline on the same fixture.

## Open Questions

- **Conda env pin granularity**: pin exact versions only, or exact builds (`=ver=build`) for stricter reproducibility at the cost of solve fragility? Default to exact versions.
- **`.sif` base**: reuse the CUDA runtime base (GPU-capable, larger) or a slim CPU base for a smaller image, with a separate GPU variant? Default to the single CUDA base for now (matches full-tier GPU tools); revisit if image size is a problem for archival.
