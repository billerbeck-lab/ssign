## 1. Reference profile + resource detection

- [x] 1.1 Record the calibration reference machine's profile (CPU core count + GPU model) in `coefficients.json` `_meta`, sourced from the machine that produced the current fit. (`_meta.reference_profile` = {cpu_cores:24, gpu:"A40", gpu_vram_gb:48}.)
- [x] 1.2 Add a `machine_profile()` helper in `runtime/machine.py` that detects the current machine's cores + GPU, reusing `resources.effective_cpu_count` + `probe_cuda_device` (doctor's probes).
- [x] 1.3 Add `resource_ratio(reference, current)` returning a per-limiting-factor prior multiplier (cpu = core-count ratio ^0.5 diminishing returns; gpu = CPU-fallback penalty 0.02 when GPU absent else parity; io = 1.0). 9 unit tests in `test_machine_profile.py`.

## 2. Online estimator (pure functions in `runtime/`)

- [x] 2.1 DEVIATION: the stateful estimator already existed as `runtime/estimator.py` `Estimator` (observe/rate/eta over a plan of concurrent stages). Extended it rather than adding a duplicate `eta.py`: added a `prior_rates` seed + `Estimator.from_profile(coeffs, current)` that seeds t=0 rates from `resource_ratio`.
- [x] 2.2 Prior ETA = `Estimator.from_profile(...).eta(plan)` before any observation (plan-based API supersedes the spec's `(tools, sizes)` tuple form: sum serial stages, `max` over each concurrent stage; `EtaResult` carries the CI). Prior projections keep the wide `PRIOR_MACHINE_U` uncertainty.
- [x] 2.3 `observe(stage_idx, tool, regime, size, wallclock)`: infers `rate = effort/observed` for the tool's limiting factor; the inferred rate then overrides the prior for remaining same-factor tools.
- [x] 2.4 `eta(plan)`: machine-adjusted projection over not-yet-run steps (observed steps contribute their real wall-clock); CI tightens as more factors are observed (inferred rate -> `u_machine=0`).
- [x] 2.5 Only tools in the plan contribute (the runner builds the plan from the actual tier/`skip_*` stage list, §4); unmodelled tools are surfaced via `EtaResult.n_unmodeled` so the ETA reads as a lower bound. Covered by `test_unmodelled_tool_omitted_not_crash`.

## 3. Timeout integration

- [x] 3.1 Extend `scaled_timeout` to accept an optional machine-rate adjustment (the inferred class rate) and apply it to `predicted_seconds` before the `max(floor, margin × ...)` formula; unchanged when no rate is available.
- [x] 3.2 `_scaled_tool_timeout` now passes `machine_rate=estimator.rate(limiting_factor(tool))` (the *observed* rate, not the prior) into `scaled_timeout`; None until a same-factor tool completes, so pre-estimator behaviour (floor/margin/no-fit) is preserved until then. Defensive: any failure -> machine_rate=None.

## 4. Runner wiring + live display

- [x] 4.1 `Estimator.from_profile()` built at the top of `_execute_stages`; the first ETA (prior) prints as `estimated total ~Y min` once extract_proteins yields the protein count. Sizes computed lazily in `_eta_sizes()` (proteins/neighborhood via cached FASTA counts, substrates via measured count or a 2% prior before filtering).
- [x] 4.2 `_eta_step(stage_idx, step_id, duration)` after each successful step observes it + prints `estimated remaining ~Y min (range lo-hi)`. Plan rebuilt each call from current sizes (refines once filtering lands the real substrate count). Wrapped so a bug never breaks a run.
- [x] 4.3 Already satisfied: the parallel path keys each future to its own `t_step_start`, so `result.duration_s` is per-tool. The hook fires per-tool inside the `as_completed` loop with the shared stage_idx (distinct `(stage_idx, tool)` keys). New logic in pure `runtime/run_eta.py`; tests in `test_run_eta.py` (12).

## 5. Replay harness (spec requirement)

- [x] 5.1 `run_eta.replay(coeffs, plan, completions)` replays ordered `(stage_idx, tool, regime, size, wallclock)` completions through a fresh bare `Estimator`, recording the projected total + CI after each step (only data up to that step).
- [x] 5.2 `ReplayResult.steps[k].error_pct` = |projected total - true total| / true total; true total is the all-observed ETA. `test_replay_converges_to_true_total` asserts error shrinks to 0 and the CI tightens.

## 6. Tests + docs

- [x] 6.1 Done across the suite (synthetic coeffs, no mocks): (a) prior ETA -> `test_runtime_estimator::test_from_profile_slow_machine_widens_prior_eta` + `test_prior_rate_scales_projection_before_any_observation`; (b) converged ETA -> `test_observation_overrides_prior` + `test_replay_*`; (c) machine-adjusted timeout -> `test_scaled_timeout::test_machine_rate_slow_widens_timeout` / `_fast_never_below_floor` / `_none_is_unchanged`.
- [x] 6.2 DEVIATION: `runs.jsonl` is in the private memory calibration dir (never in git), so the in-repo convergence test uses a deterministic synthetic run (`test_run_eta::test_replay_converges_to_true_total`, 0.5x machine -> error -> 0). The memory-dir `replay.py` still covers real recorded runs.
- [x] 6.3 Memory-dir `calibration/README.md` notes the reference_profile encoding + resource_ratio prior; `docs/tutorials/first_run.md` describes the live `estimated remaining ~Y min` line (informational, wide-then-sharpening, CPU-only caveat). EggNOG refit VERIFIED (targeted, per Teo): OLS over the 74 fittable runs.jsonl rows reproduces the shipped eggnog a/b/LOO exactly (137% LOO is intrinsic cache variance, mitigated by the scratch-staging fix); recorded in coefficients.json _meta, no coefficient change.
