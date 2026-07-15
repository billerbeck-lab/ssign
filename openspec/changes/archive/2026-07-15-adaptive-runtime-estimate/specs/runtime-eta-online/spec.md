## MODIFIED Requirements

### Requirement: Prior ETA at run start

At t=0, once `n_proteins` and the tool set are known but no tool has finished, the estimator SHALL emit a prior ETA with a wide confidence interval. The prior SHALL be derived from a **coarse resource ratio** between the current machine (GPU present/type + CPU core count) and the calibration reference machine's recorded profile, applied to the per-tool effort sum, and MAY additionally be tightened with historical `_pipeline_total` rows at similar size/tier when available. The machine's true per-limiting-factor rates are unknown at this point and MUST be treated as latent, refined online as tools complete.

#### Scenario: First estimate before any tool completes on a never-seen machine
- **WHEN** input parsing finishes, `n_proteins` is known, and there is no historical `_pipeline_total` row for this machine
- **THEN** the estimator still returns a prior ETA, scaled by the current-vs-reference resource ratio, with a confidence interval wide enough to cover machine-speed uncertainty

#### Scenario: Reference machine profile is available for the ratio
- **WHEN** the effort coefficients are loaded
- **THEN** they carry the reference machine's core count and GPU model, so the current-vs-reference resource ratio can be computed

## ADDED Requirements

### Requirement: Live ETA displayed during a run

The runner SHALL surface the estimator's ETA to the user: once as a prior after input parsing, and then a refined "estimated remaining" after each completed step, on a single log line, using the machine-adjusted projection rather than the raw reference estimate.

#### Scenario: ETA refines after the first tools of each class
- **WHEN** the first CPU-bound tool and the first GPU-bound tool have completed
- **THEN** the displayed "estimated remaining" reflects this machine's inferred `cpu_rate` and `gpu_rate`, not the reference-machine prior

#### Scenario: Estimate is presented as approximate
- **WHEN** the ETA is shown
- **THEN** it is labelled as an estimate (e.g. "~Y min") and never as an exact time, consistent with the wide early confidence interval
