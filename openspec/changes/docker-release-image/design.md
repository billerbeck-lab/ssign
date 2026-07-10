## Context

`containers/Dockerfile` (May) freezes the ssign stack on a SHA-pinned CUDA base and installs `.[extended]` non-editably; databases/weights are bind-mounted, not baked (keeps the image ~3–5 GB vs 600+ GB and dodges DTU/EggNOG redistribution). It is a **skeleton**: never built, CUDA base SHA is a placeholder of zeros, and Python deps are unpinned (no lockfile). `docs/development/install_test_runbook.md` Part B (Docker build test) is runbook-only, gated on Part A + the full-tier smoke test. The full-tier smoke test passed (run 3256040), so image correctness no longer waits on anything; the 74/75-genome CX3 reruns are a scheduling/scale concern orthogonal to the image.

Confirmed: `extended` and `full` optional-dependency sets are byte-identical (18 packages each); tiers differ only by which databases/weights exist on the host. So there is one code image, and the tier is a mount-time choice.

## Goals / Non-Goals

**Goals:**
- A `docker build` that resolves the **identical** stack years apart (base image digest + Python deps both pinned).
- The image actually built once and proven to run the pipeline **offline** on the golden T5aSS fixture with bind-mounted DBs/models, producing the known substrate call.
- A documented per-tier mount matrix (base / extended / full) against the single image.
- Verified in-container paths for the not-bundled tools (DTU bind-mount + webserver fallback; EggNOG host-installed) and a verified Singularity/Apptainer conversion.

**Non-Goals:**
- Confirming full tier at 74/75-genome scale (in-flight CX3 reruns; a bug there is a runner fix + cheap rebuild, not an image-design change).
- Publishing to a registry or minting a Zenodo DOI (release-time).
- Bundling DTU SignalP/DeepLocPro weights (licence-blocked; DTU confirmed non-redistributable 2026-05-07).
- Changing pipeline/runtime code or which tools run.

## Decisions

1. **Dependency lock mechanism.** Generate a fully-resolved lockfile and wire it into the image build so `pip install` is deterministic. *Option A: `uv.lock`* (fast, hash-pinned, but adds `uv` to the build) vs *Option B: a `containers/constraints.txt`* from `pip freeze` of a known-good venv, passed as `pip install -c`. Leaning **B** (constraints.txt): no new build tool, transparent, and it pins the exact versions the smoke test validated. The lockfile is committed so the pinned set is auditable. (Confirm against publication-roadmap step 9, which may standardise on `uv.lock` repo-wide — if so, use A.)
2. **One image, per-tier by mount.** Do not build three images. The Dockerfile installs `.[extended]` (== full code); base users simply mount only the base DB set. Document the mount matrix rather than branching the Dockerfile. *Alternative:* a slim base image without the extended pip extras. Rejected for now: it doubles the build/test surface for a marginal image-size win, and the extras are pip-only (no heavy system deps).
3. **CUDA base SHA locked at pin-time, not floating.** Resolve `nvidia/cuda:12.4.1-runtime-ubuntu22.04` to its current digest and hard-code it. Accept that a future CVE in the base means a deliberate re-pin (documented), which is the correct trade-off for reproducibility.
4. **Smoke test = the golden fixture, offline.** Reuse `tests/fixtures/golden/t5ass_minimal` (DeepLocPro-local, every other tool off) as the in-container acceptance run: it needs only a small DLP mount, finishes in ~2 min, and has a known-correct expected substrate. This makes the Docker smoke test a real scientific check, not just "container starts."
5. **Build happens on a CUDA host.** The lockfile + Dockerfile edits are prepared on the laptop; the actual `docker build` + `--gpus all` run happen on the RTX 4070 box (or CX3 via Singularity). The change delivers the exact build/run commands, not a laptop-run image.

## Risks / Trade-offs

- **Pinning surfaces a dependency conflict** (e.g. torch/CUDA vs a transitive pin) → caught at the first real build; may require a `pyproject.toml` constraint tweak (not app logic). Better found now than at release.
- **Base image digest ages** (CVEs) → documented re-pin procedure; reproducibility is worth the manual bump.
- **GPU-host dependency** means the smoke test can't run in CI on the laptop → mitigated by the golden fixture being CPU-friendly (DeepLocPro CPU) so the *pipeline* part runs anywhere Docker+DLP exist; only GPU-accelerated tools need the GPU.
- **Singularity drift** on HPC (older Apptainer) → verify the `.sif` build on the actual target (CX3) once, document the version.

## Migration Plan

- Additive: no existing behaviour changes. The lockfile and Dockerfile edits are new/localised; rollback = revert the change (the skeleton returns).
- Validation gate: the image builds cleanly and the in-container golden run matches the expected T5aSS substrate call. Record the resolved CUDA digest + lockfile hash in the change.

## Open Questions

- ~~**Lockfile format**~~ **Resolved: `uv.lock`.** The repo already standardises on it (`research_software_audit.md` cites `uv.lock`; a committed `uv.lock` exists). The Dockerfile exports the frozen lock to a pip requirements file at build time (keeps the install pip-based, lock authoritative). The committed lock was stale and has been refreshed.
- **Base-image slimming**: is a separate minimal base-tier image wanted for the webserver/GUI-only user, or is "one image, mount less" sufficient for v1.0.0? Default to the latter unless the supervisor wants a smaller base artifact.
