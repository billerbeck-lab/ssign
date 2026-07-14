## Why

The effort model already projects per-tool wall-clock, and the `runtime-eta-online` capability already **specifies** a machine-adaptive ETA (a prior at run start, online rate inference per limiting factor, a replay harness), but it was never implemented: no ETA is shown to the user, and `scaled_timeout` still uses the raw reference-machine estimate. As we saw live (a legit tool on a slow network filesystem added ~77 min; a weak machine could be killed by a timeout sized for the reference box), the machine actually running matters. This change implements that capability and feeds its machine-adjusted projections into both a live ETA and the kill-timeout.

## What Changes

- Implement the spec-only `runtime-eta-online` estimator: tool-set-agnostic ETA composition (per-tool effort sum, `max` over the concurrent predictor block), a prior ETA at run start, and online machine-rate inference (`rate = effort(size) / observed_wallclock`) applied to not-yet-run tools sharing the same limiting factor, refining as each tool completes.
- Refine the "prior at t=0" requirement: derive the initial prior from a **coarse resource ratio** (the current machine's GPU present/type + CPU core count vs the calibration reference machine's profile), in addition to the spec's historical `_pipeline_total` rows, so a first estimate exists even with no history and on a never-seen machine. Requires encoding the reference machine profile (cores/GPU) alongside the coefficients.
- Wire the runner step loop to feed each completed tool's observed duration into the estimator and print a live-refining `estimated remaining ~Y min` after each step.
- **Modify `size-aware-tool-timeouts`**: `scaled_timeout` consumes the estimator's machine-adjusted per-tool prediction (once a limiting-factor rate is inferred) instead of the raw reference estimate, so a slow machine gets proportionally longer kill-timeouts. The floor, the confidence-aware margins, the no-fit fallback, and the wrapper/runner agreement are all preserved.
- Ship the replay harness the spec requires (offline convergence check from a finished run's per-tool wall-clocks).
- RAM stays out of scope for runtime (treated as noise); the spec's `io_factor` already captures the storage-speed effect that dominated the 77-min startup.

## Capabilities

### New Capabilities

(none: the `runtime-eta-online` capability already exists as a spec; this change implements and refines it.)

### Modified Capabilities

- `runtime-eta-online`: implement the estimator; refine the t=0 prior to use a resource ratio (GPU/cores) with the reference machine profile encoded; add the runner-side live display.
- `size-aware-tool-timeouts`: `scaled_timeout` consumes the machine-adjusted per-tool prediction (online-inferred limiting-factor rate) rather than the raw reference-machine estimate, keeping floor / margin / no-fit behaviour intact.

## Impact

- `src/ssign_app/runtime/`: a new online-estimator module (prior + per-limiting-factor rate inference + ETA composition); reference machine profile + resource detection; `timeouts.scaled_timeout` gains an optional machine-rate adjustment; reuse the existing `limiting_factor` / `effort` / `resolve_regime`.
- `src/ssign_app/core/runner.py`: feed observed per-step durations to the estimator, print the live ETA, pass machine-adjusted timeouts. Resource detection already exists for `doctor`.
- `calibration/`: record the reference machine profile (cores/GPU) so ratios are possible; add the replay-convergence harness.
- Tests: deterministic unit tests (fixed coefficients + a fake machine profile + simulated observed durations, asserting the prior ETA, the converged ETA, and the machine-adjusted `scaled_timeout`), plus a replay-convergence test.
