## 1. Reference profile + resource detection

- [ ] 1.1 Record the calibration reference machine's profile (CPU core count + GPU model) in `coefficients.json` `_meta`, sourced from the machine that produced the current fit.
- [ ] 1.2 Add a `machine_profile()` helper in `runtime/` that detects the current machine's cores (`os.cpu_count`/cgroup-effective) and GPU (present/model), reusing the detection `doctor` already does.
- [ ] 1.3 Add `resource_ratio(reference, current)` returning a per-limiting-factor prior multiplier (cpu ~ core-count ratio with diminishing returns; gpu ~ present/model factor; io ~ 1.0).

## 2. Online estimator (pure functions in `runtime/`)

- [ ] 2.1 New `runtime/eta.py` with a stateful `RuntimeEstimator`: holds per-limiting-factor rate estimates (init from the resource prior), the tool set to run, and their regimes/sizes.
- [ ] 2.2 `prior_eta(tools, sizes, profile)`: compose the per-tool `effort()` sum with `max` over the concurrent predictor block, scaled by `resource_ratio`; return `(eta_seconds, ci)`.
- [ ] 2.3 `observe(tool, observed_seconds)`: infer `rate = effort(size)/observed` for the tool's limiting factor, update that class's rate (clamped to a sane range), and recompute the remaining ETA.
- [ ] 2.4 `remaining_eta()`: machine-adjusted projection over not-yet-run tools; tighten CI as more classes are observed.
- [ ] 2.5 Honour tier + `skip_*`: only tools that will actually run contribute; unmodelled tools (no `effort` fit) contribute a floor-based guess flagged in the CI.

## 3. Timeout integration

- [ ] 3.1 Extend `scaled_timeout` to accept an optional machine-rate adjustment (the inferred class rate) and apply it to `predicted_seconds` before the `max(floor, margin × ...)` formula; unchanged when no rate is available.
- [ ] 3.2 Thread the estimator's current class rate into the runner's `_scaled_tool_timeout` so each tool's kill-timeout tracks the observed machine speed. Preserve floor / margin / no-fit / low-confidence behaviour.

## 4. Runner wiring + live display

- [ ] 4.1 Build the `RuntimeEstimator` once `n_proteins` and the tool set are known (after input parsing); print the prior ETA.
- [ ] 4.2 After each completed step, call `observe()` with the measured duration and print a one-line `estimated remaining ~Y min`.
- [ ] 4.3 Ensure the parallel predictor block reports each tool's own duration (not the block's) so per-class rates are inferred correctly.

## 5. Replay harness (spec requirement)

- [ ] 5.1 Offline harness that replays a finished run's per-tool wall-clocks through the estimator in completion order, recording ETA + CI after each step, using only data available up to that step.
- [ ] 5.2 Report ETA error per step to show convergence toward the observed total.

## 6. Tests + docs

- [ ] 6.1 Deterministic unit tests: fixed coefficients + a fake reference/current profile + simulated observed durations, asserting (a) the prior ETA, (b) the converged ETA after class observations, (c) the machine-adjusted `scaled_timeout` (slow-machine widens, no-rate falls back).
- [ ] 6.2 Replay-convergence test on a recorded run from `calibration/runs.jsonl`.
- [ ] 6.3 Update `calibration/README.md` (reference profile field) and note the live ETA in the run docs.
