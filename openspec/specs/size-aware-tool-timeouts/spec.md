# size-aware-tool-timeouts Specification

## Purpose
TBD - created by archiving change size-aware-tool-timeouts. Update Purpose after archive.
## Requirements
### Requirement: Size-aware timeout derivation

The system SHALL provide a `scaled_timeout(tool, size, regime)` helper that returns an integer subprocess timeout in seconds equal to `max(floor, ceil(margin × predicted_seconds))`, where `predicted_seconds` is the effort model's estimate for `(tool, size, regime)` **adjusted by the current machine's inferred rate for the tool's limiting factor when such a rate is available**, and otherwise the raw reference-machine estimate; `floor` defaults to `TOOL_TIMEOUT_S`, and `margin` defaults to 2. The machine adjustment SHALL only ever be applied to the point estimate fed to the existing `max(floor, margin × ...)` formula; the floor and margin semantics are unchanged.

#### Scenario: Large input scales the cap above the floor
- **WHEN** `scaled_timeout` is called for a modelled tool whose predicted time exceeds `floor / margin`
- **THEN** it returns `ceil(margin × predicted_seconds)`, above the floor

#### Scenario: Small input stays at the floor
- **WHEN** the predicted time is small enough that `margin × predicted_seconds` is below `floor`
- **THEN** `scaled_timeout` returns `floor` unchanged

#### Scenario: Slow machine widens the cap
- **WHEN** the machine's inferred rate for a tool's limiting factor shows it runs slower than the reference machine, and that rate is available before the tool launches
- **THEN** `scaled_timeout` returns a proportionally larger timeout than the raw reference estimate would give, so a legitimately slow tool on a weak or IO-bound machine is not killed prematurely

#### Scenario: No inferred rate yet falls back to the reference estimate
- **WHEN** no tool of the same limiting factor has completed, so no machine rate is inferred
- **THEN** `scaled_timeout` uses the raw reference-machine estimate, preserving today's behaviour

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

