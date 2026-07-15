# runtime-eta-online Specification

## Purpose
TBD - created by archiving change runtime-eta-estimator. Update Purpose after archive.
## Requirements
### Requirement: Tool-set-agnostic ETA composition

The estimator SHALL compute a run's total ETA by summing the per-tool effort projections for exactly the tools that will run, accounting for tools that execute in parallel (the predictor block) versus serially. It MUST handle any tier or cherry-picked tool subset, including combinations not present in the calibration data, by composing the available per-tool models.

#### Scenario: Cherry-picked subset never run before
- **WHEN** a run enables a tool combination with no matching historical `_pipeline_total` row
- **THEN** the estimator still returns an ETA built from the individual per-tool effort models

#### Scenario: Parallel predictor block
- **WHEN** DeepLocPro, SignalP, and DeepSecE run concurrently
- **THEN** their contribution to the ETA is the max of their projected times, not the sum

### Requirement: Prior ETA at run start

At t=0, once `n_proteins` and the tool set are known but no tool has finished, the estimator SHALL emit a prior ETA with a wide confidence interval. The prior SHALL be derived from a **coarse resource ratio** between the current machine (GPU present/type + CPU core count) and the calibration reference machine's recorded profile, applied to the per-tool effort sum, and MAY additionally be tightened with historical `_pipeline_total` rows at similar size/tier when available. The machine's true per-limiting-factor rates are unknown at this point and MUST be treated as latent, refined online as tools complete.

#### Scenario: First estimate before any tool completes on a never-seen machine
- **WHEN** input parsing finishes, `n_proteins` is known, and there is no historical `_pipeline_total` row for this machine
- **THEN** the estimator still returns a prior ETA, scaled by the current-vs-reference resource ratio, with a confidence interval wide enough to cover machine-speed uncertainty

#### Scenario: Reference machine profile is available for the ratio
- **WHEN** the effort coefficients are loaded
- **THEN** they carry the reference machine's core count and GPU model, so the current-vs-reference resource ratio can be computed

### Requirement: Online machine-rate inference

As each tool completes, the estimator SHALL infer this machine's rate for that tool's limiting factor (`cpu_rate`, `gpu_rate`, or `io_factor`) as `effort_T(size) / observed_wallclock`, and apply the inferred rate to all not-yet-run tools sharing that limiting factor. Subsequent completions SHALL refine the rate estimate.

#### Scenario: First CPU tool tightens the estimate
- **WHEN** the first CPU-bound tool (e.g. MacSyFinder or Bakta) finishes
- **THEN** the estimator updates `cpu_rate` and narrows the ETA for all remaining CPU-bound tools

#### Scenario: First GPU tool tightens the estimate
- **WHEN** the first GPU-bound tool finishes
- **THEN** the estimator updates `gpu_rate` and narrows the ETA for all remaining GPU-bound tools

#### Scenario: IO factor inferred from EggNOG
- **WHEN** EggNOG finishes far slower than its CPU effort predicts
- **THEN** the estimator raises `io_factor` (inferring NFS-bound storage) and widens IO-sensitive tool estimates accordingly

### Requirement: Replay harness for convergence validation

The system SHALL provide an offline harness that replays a finished run's per-tool wallclocks through the estimator in completion order, recording the ETA and CI after each step, so convergence toward the true total can be measured. The harness MUST use only the per-tool wallclocks available up to each step, never the final total.

#### Scenario: Replay reports convergence
- **WHEN** a finished run's tool completions are replayed in order
- **THEN** the harness reports the ETA error after each completion, showing the estimate converging to the observed total

### Requirement: Live ETA displayed during a run

The runner SHALL surface the estimator's ETA to the user: once as a prior after input parsing, and then a refined "estimated remaining" after each completed step, on a single log line, using the machine-adjusted projection rather than the raw reference estimate.

#### Scenario: ETA refines after the first tools of each class
- **WHEN** the first CPU-bound tool and the first GPU-bound tool have completed
- **THEN** the displayed "estimated remaining" reflects this machine's inferred `cpu_rate` and `gpu_rate`, not the reference-machine prior

#### Scenario: Estimate is presented as approximate
- **WHEN** the ETA is shown
- **THEN** it is labelled as an estimate (e.g. "~Y min") and never as an exact time, consistent with the wide early confidence interval

