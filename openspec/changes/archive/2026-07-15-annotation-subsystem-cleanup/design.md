## Context

The annotation phase runs up to six tools (BLASTp, HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam) concurrently inside `runner.parallel_group(N)`, which sets `SSIGN_PARALLEL_GROUP_SIZE=N`. `ssign_lib/resources.py` exposes `effective_ram_gb()` (min of SSIGN_MAX_RAM_GB / scheduler alloc / cgroup / physical) and `parallel_share_ram_gb()` (= effective / N). Tools inside the group must size memory to the *share*, not the whole allocation. A 2026-07-15 audit found EggNOG force-loads its 38.5 GB DB via a CLI default that bypasses its own RAM autodetect, and three siblings size to the whole node. Separately, PLM-Effector has been off at every tier since 2026-06-17 (over-calls ~25% of a proteome) but is still fully vendored and wired in.

## Goals / Non-Goals

**Goals:**
- Extended/full runs fit their RAM allocation on ≤32 GB nodes without manual flags.
- Every annotation tool sizes memory/CPU to its parallel-group share.
- PLM-Effector fully removed from the codebase (Teo confirmed full deletion incl. vendored package, 2026-07-15).
- A single `--skip-annotation` flag replaces enumerating six `--skip-*`.

**Non-Goals:**
- Changing prediction tools (DLP/DSE/SignalP) — audited CLEAN.
- Serializing the annotation group (the share-sizing fix makes concurrency fit; no need to run tools one-at-a-time).
- Touching the DBs or the image bake (HH-suite bake + build-cache are separate, in reproducible-install).

## Decisions

- **Size to the share, not the node.** EggNOG already uses `parallel_share_ram_gb()`; make InterProScan (Java `-Xmx`), HH-suite (worker fan-out), and pLM-BLAST (CPU budget when GPU-less) match. When a tool runs standalone (group size 1) `parallel_share_ram_gb()` returns the full amount, so single-tool behavior is unchanged, only concurrent behavior is corrected.
- **EggNOG dbmem default `True` → `None`.** The RAM autodetect (`_autodetect_dbmem`, share ≥ 50 GB) already exists; the CLI default was defeating it. `None` lets it govern: dbmem on when RAM allows, off (mmap the local-cached DB) otherwise. Alternative (keep `True`, document `--no-eggnog-dbmem`) rejected: it leaves a footgun that every small-node user must know to disarm.
- **Full PLM-E deletion.** Alternative (keep the package importable for the classifier) was offered and declined; the classifier will source that feature elsewhere. Remove wiring from `cli.py`, `runner.py` (config + `_step_plm_effector` + effort model), `constants.py`, `doctor.py`, `dependency_manifest.py`, `Home.py`, `system_filtering.py`, `proximity_analysis.py`, `cross_validate_predictions.py`, `multi_runner.py`, and delete `run_plm_effector.py`, `merge_plm_effector_outputs.py`, and the vendored `plm_effector/` dir + its `coefficients.json`/`effort_model.py` entries.
- **`--skip-annotation` is a convenience OR over the six skip flags.** It skips only the optional annotation phase; predictions and the T5aSS whole-protein pass that depend on annotation degrade gracefully (empty annotation columns), matching how each individual `--skip-*` already behaves.

## Risks / Trade-offs

- [Removing PLM-E from the effort model / `coefficients.json` leaves a dangling key lookup] → grep every `_PREDICT`/`_WHOLE`/coefficient reference; add a test that the effort model builds with no PLM-E entry.
- [HH-suite RAM clamp reduces worker parallelism on big-RAM many-core nodes] → it only clamps *down* when RAM-limited; a 64 GB+ node keeps the CPU-derived count.
- [`--skip-annotation` masking a tier that expected annotations] → it is an explicit user opt-out; document that it empties annotation columns, same as the per-tool flags.
- [Classifier feature-source loss] → accepted by Teo; note it in `design_decisions.md` so the classifier project is aware.

## Migration Plan

- PLM-E was off by default, so no default-output change; anyone who passed `--no-skip-plm-effector` now gets an unknown-flag error (intended — it's gone). Note in the changelog.
- Ships in the next image rebuild (v7); no runtime migration.
