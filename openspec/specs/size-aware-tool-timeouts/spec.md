# size-aware-tool-timeouts Specification

## Purpose
TBD - created by archiving change size-aware-tool-timeouts. Update Purpose after archive.
## Requirements
### Requirement: Size-aware timeout derivation

The system SHALL provide a `scaled_timeout(tool, size, regime)` helper that returns an integer subprocess timeout in seconds equal to `max(floor, ceil(margin × predicted_seconds))`, where `predicted_seconds` is the effort model's reference-machine estimate for `(tool, size, regime)`, `floor` defaults to `TOOL_TIMEOUT_S`, and `margin` defaults to 2.

#### Scenario: Large input scales the cap above the floor
- **WHEN** `scaled_timeout` is called for a tool whose predicted time exceeds `floor / margin` (e.g. DeepLocPro whole_genome at 160,831 proteins, predicted ≈ 5.4h)
- **THEN** the returned timeout is `ceil(margin × predicted_seconds)` (≈ 11h), strictly greater than the floor

#### Scenario: Small input stays at the floor
- **WHEN** `scaled_timeout` is called for an input whose `margin × predicted_seconds` is below `floor` (e.g. a neighborhood prediction of a few hundred proteins)
- **THEN** the returned timeout equals `floor`, leaving current behavior unchanged

### Requirement: Unmodelled tools fall back to the floor

When the effort model has no fit for `(tool, regime)` (e.g. blastp, hh_suite, orthologs-blast), `scaled_timeout` SHALL return `floor` unchanged, so no unmodelled tool loses its existing timeout.

#### Scenario: Tool with no coefficient block
- **WHEN** `scaled_timeout` is called for a tool/regime with no entry in `coefficients.json`
- **THEN** the returned timeout equals `floor` (the tool's current fixed cap)

### Requirement: Confidence-aware margin

For fits flagged `low_confidence`, `scaled_timeout` SHALL apply a wider margin than the default so a thin or noisy fit cannot produce a timeout that is likely to kill a still-running tool.

#### Scenario: Low-confidence fit widened
- **WHEN** `scaled_timeout` is called for a `(tool, regime)` whose effort fit has `low_confidence = True`
- **THEN** the effective margin is greater than the default confident-fit margin (2), producing a more generous cap

#### Scenario: Degenerate low fit still protected by the floor
- **WHEN** a low-confidence fit predicts a value far below the true runtime (e.g. SignalP whole_genome, a flat ≈13 min fit)
- **THEN** the returned timeout is at least `floor`, and since that tool's true runtime is below the floor, it is not killed prematurely

### Requirement: Runner derives each tool timeout from actual input size

For every modelled tool it invokes, the runner SHALL compute the subprocess timeout via `scaled_timeout` using the tool's actual input size (sequence count of the FASTA being processed) and its resolved regime, instead of a hardcoded constant.

#### Scenario: Whole-genome prediction on a pooled proteome
- **WHEN** the runner invokes DeepLocPro/DeepSecE/SignalP on a pooled whole-genome FASTA of N proteins with enrichment enabled
- **THEN** each tool's timeout is `scaled_timeout(tool, N, "whole_genome")`, large enough that a run whose predicted time is under the cap completes rather than being killed at 4h

### Requirement: Wrapper and runner timeout layers agree

A tool wrapper that enforces its own inner `subprocess.run` timeout SHALL accept the runner-computed value (via a `--timeout` argument) and use it, defaulting to `TOOL_TIMEOUT_S` when invoked standalone, so the inner cap never fires before the runner's outer cap.

#### Scenario: Runner passes a scaled timeout to a wrapper
- **WHEN** the runner computes a scaled timeout of T seconds for a tool and invokes its wrapper
- **THEN** the wrapper's inner `subprocess.run` uses T (not the hardcoded `TOOL_TIMEOUT_S`)

#### Scenario: Wrapper run standalone
- **WHEN** a wrapper is invoked directly without `--timeout`
- **THEN** it uses `TOOL_TIMEOUT_S` as before

### Requirement: Timeouts fail loud

When a tool exceeds its (scaled) timeout, the pipeline SHALL surface the timeout as an explicit, attributable error identifying the tool and the cap, rather than allowing it to be masked as a downstream "missing input" failure.

#### Scenario: Prediction timeout is reported directly
- **WHEN** a whole-genome prediction tool exceeds its scaled timeout
- **THEN** the run reports that tool timing out at its cap, not only the downstream "No DeepLocPro output" cross-validate failure

