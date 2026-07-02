## Why

Every external tool runs under a flat `TOOL_TIMEOUT_S = 14400` (4h) cap, sized for one genome. On the 52-genome benchmark single-job run (`MultiGenomeRunner`, `--enrichment-stats`), segment B pools every genome's proteome into **one** whole-genome prediction call (160,831 proteins). DeepLocPro/DeepSecE/SignalP each ran past 4h and were killed, cross-validate got no localization output, and **all 52 genomes emitted 0 secreted proteins**. Verified it is a *time* wall, not memory (all three stream in batches). The failure is a fixed cap that doesn't scale with input size, not a structural flaw: the bundled effort model predicts DLP whole-genome for 160,831 proteins ≈ 5.4h, DSE ≈ 5.0h, so 4h was ~10-30% too short. ssign already ships a fitted per-tool runtime model (`runtime/effort_model.py` + `coefficients.json`); the cap should be derived from it instead of hardcoded.

## What Changes

- Add a `scaled_timeout(tool, size, regime)` helper that returns `max(floor, margin × predicted_seconds)`, where `predicted_seconds` comes from the existing effort model and `floor` is the current `TOOL_TIMEOUT_S`. Unmodelled tools return `floor` unchanged.
- The **runner** computes each tool's timeout from its actual input size (`count_sequences`) + regime and passes it to `run_script(..., timeout=...)`, replacing the ~10 hardcoded `14400`/`7200` values. The runner's `_run_script_streaming` layer is what killed the failed run, so this is the layer that must scale.
- Each tool **wrapper** accepts `--timeout N` so its inner `subprocess.run(timeout=...)` uses the same scaled cap (defaulting to `TOOL_TIMEOUT_S` for standalone use). Otherwise the inner 4h fires first and defeats the fix.
- Applies to all 11 modelled tools (bakta, deeplocpro, deepsece, eggnog, extract_proteins, interproscan, macsyfinder, plm_blast, plm_effector, protparam, signalp), not just the whole-genome predictors. Only-substrate tools (EggNOG especially, predicted ~4.4h at 819 pooled substrates) are the next most likely to hit the flat cap at panel scale, so the guard is general.
- **Margin policy:** `margin = 2` for well-calibrated fits (`low_confidence = False`); a wider margin for `low_confidence = True` fits so a thin/noisy fit can't false-kill. The floor covers the degenerate cases (e.g. SignalP whole-genome is a garbage fit but predicts tiny, and SignalP is fast, so the floor is safe).
- **Fail-loud:** a prediction-step timeout surfaces as a clear error rather than silently cascading to "No DeepLocPro output" + 0 substrates for every genome.
- Small/neighborhood runs are unchanged: predicted time is tiny, so `max(floor, …)` keeps the current 4h cap (zero regression).

## Capabilities

### New Capabilities
- `size-aware-tool-timeouts`: derive each external tool's subprocess timeout from the fitted effort model (scaled by input size, with a floor and a confidence-aware margin), enforced consistently at the runner and wrapper layers.

### Modified Capabilities
<!-- None. This capability consumes `runtime-effort-model` unchanged; it does not alter that model's requirements. -->

## Impact

- **Code**: new helper in `src/ssign_app/runtime/` (consumes `effort_model.effort`/`resolve_regime`); `core/runner.py` (~10 `run_script(..., timeout=)` call sites compute the scaled value); tool wrappers with an inner `subprocess.run` timeout (`run_deeplocpro`, `run_signalp`, `run_interproscan`, `run_eggnog`, `run_blastp`, `run_bakta`, `run_ortholog_grouping`) gain a `--timeout` arg. DeepSecE needs only the runner-side change (no inner timeout).
- **Behavior**: large-input tool calls get a proportionally larger cap and complete instead of being killed; small runs unchanged. `TOOL_TIMEOUT_S` becomes a floor, not the value. `plm_blast` keeps its 24h floor.
- **Constants/data**: `TOOL_TIMEOUT_S` semantics change from "the cap" to "the floor"; `coefficients.json` is now read on the hot path (runner) as well as the GUI.
- **Tests**: unit tests for `scaled_timeout` (floor, margin, unmodelled→floor, low-confidence widening); wrapper `--timeout` plumbing; a runner test that a modelled step receives a size-scaled timeout.
- **Out of scope**: PBS `--walltime` (set manually at CX3 submit; docs note a full panel needs `--walltime 24:00:00`); the `MultiGenomeRunner` per-genome / chunking rework (explicitly rejected as unnecessary); DeepSecE `batch_size=1` → VRAM-auto batching (a real speedup, tracked as a separable stretch task).
