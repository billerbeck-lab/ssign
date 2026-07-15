## MODIFIED Requirements

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
