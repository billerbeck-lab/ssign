## Context

The `runtime-eta-online` spec defines a machine-adaptive ETA (prior at t=0, online per-limiting-factor rate inference, replay harness) but was never implemented. `src/ssign_app/runtime/` already provides the reference-machine math: `effort(tool, size, regime)` (returns projected seconds or None), `resolve_regime()`, `limiting_factor(tool)` (cpu / gpu / io), and `scaled_timeout()` (used only for kill-timeouts, off the raw reference estimate). Cross-machine wall-clocks are recorded in `calibration/runs.jsonl`, but the coefficients do not encode the reference machine's own profile (cores / GPU), so no current-vs-reference ratio can be computed yet.

## Goals / Non-Goals

**Goals:**
- Implement the online estimator: prior ETA, online per-limiting-factor rate inference, tool-set-agnostic composition.
- Show a live-refining "estimated remaining ~Y min" in the runner.
- Feed the machine-adjusted per-tool prediction into `scaled_timeout` (kill-timeouts track the real machine).
- Encode the reference machine profile and add a coarse resource prior at t=0.
- Deterministic unit tests + the replay-convergence harness the spec requires.

**Non-Goals:**
- RAM-based runtime modelling (noise; the `io_factor` already captures storage-speed effects).
- Refitting the effort coefficients (separate calibration work).
- Predicting scheduler queue wait (out of scope; this estimates compute time only).

## Decisions

1. **Limiting-factor classes, reuse existing `limiting_factor()`.** Each tool is cpu-, gpu-, or io-bound. As a tool finishes, infer that class's machine rate `r = effort(size) / observed_wallclock` and apply it to the remaining tools of the same class. This is the per-class factor idea, and the third (`io`) class captures network-filesystem slowness (the 77-min startup case).
2. **Prior at t=0, reconcile the resource ratio with historical rows.** Use a coarse resource ratio (current cores + GPU present/type vs the reference profile) as the portable default so a first ETA exists on a never-seen machine, and tighten with historical `_pipeline_total` rows when present. Wide CI until observations arrive.
3. **Reference machine profile** (core count + GPU model) is recorded alongside the coefficients (in `_meta`), sourced from the calibration machine that produced the fit. Ratio per class: cpu ~ core-count ratio with diminishing returns, gpu ~ present/model factor, io ~ 1.0 prior then learned.
4. **Composition:** sum serial tools; take the `max` over the concurrent predictor block (DeepLocPro / DeepSecE / SignalP). Only tools that will actually run (tier + `skip_*`) contribute.
5. **Timeout integration:** `scaled_timeout` gains an optional machine-rate adjustment. Before any observation for a class, use the resource prior; after, use the inferred rate. Floor, confidence margins, and no-fit fallback are unchanged.
6. **Display cadence:** print the prior once after input parsing (`n_proteins` known), then a one-line refined "estimated remaining ~Y min" after each completed step. Not per-log-line spam.
7. **Separation:** the estimator is a small stateful object owned by the runner; the math (prior, rate inference, composition) is pure functions in `runtime/`, testable without the runner.

## Risks / Trade-offs

- **Thin early observations.** One completion gives a noisy per-class rate; mitigate with a wide CI, clamp the rate to a sane range, and tighten as more tools of the class finish.
- **`io_factor` is confounded** (a tool can be slow for non-IO reasons); treat io as a conservatively-learned residual, not a first-class prediction.
- **Single reference machine** makes the t=0 ratio approximate; that is acceptable because the online correction supersedes the prior within the first couple of tools.
- **Un-fit tools** (blastp / hh_suite priors) contribute floor-based guesses; the ETA is explicitly approximate and carries a CI, never presented as exact.
