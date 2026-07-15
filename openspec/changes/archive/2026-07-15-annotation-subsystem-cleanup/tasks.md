# Tasks

## 1. PLM-Effector removal — voting/consensus (do this carefully, per Teo's "clean removal")
- [x] 1.1 `cross_validate_predictions.py`: drop PLM-Effector as the third vote. Remove `_plm_effector_flag`, the `plm_e_data` param + `--plm-effector` arg + its loader, the `if plm_e_secreted: evidence.append("PLM-Effector")` branch, and the `plm_effector_secreted`/`plm_effector_type`/`plm_effector_max_prob` output columns. `n_prediction_tools_agreeing` now ranges **0-2** (DLP + DSE); `is_secreted = n_agreeing >= 1` stays; SignalP stays evidence-only. Update the module docstring (lines ~18-20: "0-3" → "0-2", "three predictors" → "DeepLocPro + DeepSecE").
- [x] 1.2 `proximity_analysis.py`: remove the `plm_effector_secreted`/`plm_effector_type`/`plm_effector_max_prob` reads (lines ~177-181, 219-221) and the three column names from the substrate output header (~272-274) + the module docstring mention (~52). This is a substrate-TSV schema change.
- [x] 1.3 `system_filtering.py`, `multi_runner.py`, `doctor.py`, `dependency_manifest.py`, `Home.py` (GUI): remove every PLM-Effector reference/branch. For each, confirm the surrounding logic still holds with PLM-E gone (no dangling "3 predictors" copy, no empty GUI section).
- [x] 1.4 Grep sweep: `grep -rniE 'plm_?effector|plme|max_stacking' src/` returns nothing outside deleted files. Fix any stragglers.

## 2. PLM-Effector removal — wiring, config, effort model
- [x] 2.1 `cli.py`: remove `--skip-plm-effector`/`--no-skip-plm-effector`, `--plme-batch-size`, `--plm-effector-types`, and the `plm_effector_types`/`skip_plm_effector` passthrough (~line 605).
- [x] 2.2 `core/runner.py`: remove `_step_plm_effector`, the PLM-E step from the DAG/step list, the `plm_effector_types` PipelineConfig field (~345), `skip_plm_effector`, `plme_batch_size`, and the per-type-paths/merge logic (~2518, 2539). Keep the `secretion_evidence` output column (~3297) — it just no longer contains "PLM-Effector".
- [x] 2.3 `ssign_lib/constants.py`: remove `_OPT_IN = {"plm_effector"}` (~314) and the surrounding skip-resolution special-casing; remove any `max_stacking`/PLM-E threshold constants. Verify the tier resolver still resolves cleanly without the `_OPT_IN` member.
- [x] 2.4 `runtime/effort_model.py`: remove `"plm_effector": "gpu"` (~34) and `plm_effector` from `_PREDICT` (~45). `runtime/coefficients.json`: delete the `plm_effector` entry (~169). Add a guard/test that the effort model builds with no PLM-E key.

## 3. Delete PLM-Effector files
- [x] 3.1 Delete `src/ssign_app/scripts/run_plm_effector.py`, `src/ssign_app/scripts/merge_plm_effector_outputs.py`, and the vendored `src/ssign_app/scripts/plm_effector/` directory (ensemble.py, feature_extraction.py, inference.py, …).
- [x] 3.2 `design_decisions.md`: condense the existing PLM-E section (§3.1) to a short "tested, found to over-call ~25% of a proteome, removed 2026-07-15" record; note the secretion-classifier must source that feature independently. (Also removed §3.4 "Why PLM-Effector is vendored", the tier-size + citation mentions; renumbered §3.5→§3.4 and its cross-ref.)

## 4. Annotation resource-sizing fixes (size to the group share)
- [x] 4.1 EggNOG: `cli.py` `--eggnog-dbmem` default `True` → `None` (auto). Confirm `None` → `PipelineConfig.eggnog_dbmem=None` → runner passes neither flag → `run_eggnog._autodetect_dbmem()` governs.
- [x] 4.2 InterProScan `run_interproscan.py:168`: source the Java `-Xmx` from `parallel_share_ram_gb()` not `effective_ram_gb()`; keep the `[4,64]` clamp.
- [x] 4.3 HH-suite `run_hhsuite.py:278` + `runner.py:2298`: clamp worker fan-out by `parallel_share_ram_gb()` (e.g. `min(cpu_derived, int(share/~3.5))`), not by CPU alone.
- [x] 4.4 pLM-BLAST: count it in `_annotation_cpu_budget` `n_heavy` (runner.py ~2162) when no GPU is visible, so it takes its CPU share instead of the full `cpu_per_genome`.

## 5. `--skip-annotation` flag
- [x] 5.1 `cli.py`: add `--skip-annotation` (BooleanOptionalAction, default off). When set, resolve all six `skip_*` (blastp/hhsuite/interproscan/plmblast/eggnog/protparam) to True unless a per-tool inverse (`--no-skip-<tool>`) was explicitly passed. Resolved in `_config_from_args` via an `ann_skip` dict so per-tool `--no-skip-<tool>` (False, not None) survives.
- [x] 5.2 `core/runner.py`: verified — the flag just sets the six existing `skip_*` config fields, so the annotation steps take the same no-op path the individual flags already exercise. Confirmed end-to-end (all-six-off alone, `--no-skip-eggnog` override, tier defaults intact when off).

## 6. Docs
- [x] 6.1 `docs/how-to/install_container.md` (RAM row + Troubleshooting) and `docs/how-to/install.md`: drop the `--no-eggnog-dbmem`-on-small-nodes guidance (autodetect now handles it); state ≥64 GB is only for the in-RAM speedup. Also stripped every PLM-E mention from install.md (intro, tier table, vendored section, break-list, skip-tools list, weights export).
- [x] 6.2 `docs/reference/cli.md`: removed the PLM-Effector flag section, added `--skip-annotation` (Misc annotation table) + `--eggnog-dbmem` (EggNOG table, auto default).
- [x] 6.3 `docs/reference/output_files.md`: dropped the `plm_effector_*` example mention; documented `n_prediction_tools_agreeing` as 0-2 (DLP+DSE; SignalP evidence-only).

## 7. Tests + verification
- [x] 7.1 Deleted PLM-E unit/integration test files (all `test_plm_effector_*` + `test_run_plm_effector*` + `test_merge_plm_effector*`, already git-rm'd). Updated consumers: `test_cross_validate_predictions.py`, `test_proximity_analysis.py` (+ shared `_helpers.py` PREDICTIONS_FIELDS / make_prediction_row), `test_runner.py`, `test_multi_runner.py`, `test_doctor.py`, `test_wrapper_main_smoke.py`, `test_cli_run_subcommand.py`, `test_audit_database_sizes.py`, `test_pipeline_e2e_golden.py`, `test_imports.py`, `test_scaled_timeout.py`, `test_resources.py` — all PLM-E expectations dropped.
- [x] 7.2 Golden fixture `t5ass_minimal_predictions.tsv`: surgically dropped the three `plm_effector_*` columns (deterministic column deletion, byte format preserved). No downstream golden carried them. NOTE: the fixture is also stale for a *separate* pre-existing reason (missing `cytoplasmic_membrane_prob`); a full DeepLocPro-driven regen is still owed and can only run on a licensed-DLP host — tracked in NOTES.md.
- [x] 7.3 Unit tests: (a) already covered by `test_run_eggnog.py::TestAutodetectParallelGroupAware`; (b) added `test_run_interproscan.py::test_heap_sizes_to_group_share_not_whole_node`; (c) added `test_cli_run_subcommand.py::TestSkipAnnotation` (3 cases).
- [x] 7.4 `pytest tests/unit -q` green (1365 passed). `import ssign_app` OK; §1.4 grep sweep clean (only pLM-BLAST + the intentional design_decisions "tested, removed" note remain).
