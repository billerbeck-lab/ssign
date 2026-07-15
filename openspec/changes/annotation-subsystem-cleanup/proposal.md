## Why

Two problems in the annotation/prediction subsystem surfaced during the v6 container validation. (1) A 32 GB run OOM-killed EggNOG because `--eggnog-dbmem` is force-enabled by a CLI default, bypassing the RAM-aware autodetect that already exists; an audit found three sibling tools that also size memory/CPU to the whole node instead of their parallel-group share. (2) PLM-Effector remains vendored and wired in despite being off at every tier since 2026-06-17 (it over-calls ~25% of a proteome); carrying dead-but-loadable code cuts against the "works for a long time, low upkeep" release goal. Fixing both now lets them ride the same image rebuild.

## What Changes

- **BREAKING (CLI): Remove PLM-Effector entirely.** Delete `run_plm_effector.py`, the vendored `plm_effector/`, the `--skip-plm-effector`/`--plme-batch-size` flags, its `PipelineConfig` fields, its `constants` wiring (`_OPT_IN` + `max_stacking` gate), the `_step_plm_effector` runner step, its effort-model coefficient, and its enrichment/figure references. Leave a one-paragraph "tested, over-called ~25% of a proteome, removed" note in `design_decisions.md`.
- **Fix four annotation resource-sizing hazards** so extended/full stop OOMing on <64 GB nodes:
  - EggNOG: `--eggnog-dbmem` CLI default `True` → `None` (auto), so `run_eggnog._autodetect_dbmem()` (RAM-share ≥ 50 GB) governs instead of forcing the 38.5 GB in-RAM load.
  - InterProScan: size the Java `-Xmx` from `parallel_share_ram_gb()` not `effective_ram_gb()` (keep the `[4,64]` clamp).
  - HH-suite: clamp worker fan-out by `parallel_share_ram_gb()`, not by CPU count alone.
  - pLM-BLAST: count it in the annotation CPU budget's `n_heavy` when no GPU is present (avoid CPU oversubscription in the shared group).
- **Add `--skip-annotation`**: one flag that skips all six annotation tools (BLASTp, HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam) at once, instead of six separate `--skip-*` flags.
- **Docs:** drop the now-unnecessary `--no-eggnog-dbmem` workaround guidance from `install_container.md` and `install.md` (auto-detection handles small nodes).

## Capabilities

### New Capabilities
- `annotation-resource-sizing`: annotation tools size RAM/CPU to their parallel-group share (not the whole node), so a concurrent annotation group fits the job's allocation instead of OOM-killing on small nodes.
- `skip-annotation`: a single CLI flag to skip the entire optional annotation phase.

### Modified Capabilities
- `plme-prediction`: **removed** — PLM-Effector is deleted from the pipeline; the capability and its spec are retired.

## Impact

- Code: `cli.py`, `core/runner.py` (`PipelineConfig`, `_step_plm_effector`, `_annotation_cpu_budget`, step list), `scripts/run_eggnog.py`, `run_interproscan.py`, `run_hhsuite.py`, `run_plm_blast.py`, `scripts/run_plm_effector.py` (deleted), vendored `plm_effector/` (deleted), `ssign_lib/constants.py`, `runtime/coefficients.json` + effort-model (drop PLM-E entry).
- Tests: delete PLM-E tests; add unit tests for dbmem-off-at-low-RAM and IPS-heap-uses-share and `--skip-annotation`.
- Docs: `design_decisions.md` (PLM-E removal note), `install_container.md`, `install.md`, CLI reference.
- Image: no new bake; PLM-E drops out of the source on the next rebuild (v7). Behavior change is opt-out flag removal (PLM-E was already off by default, so no default-output change).
