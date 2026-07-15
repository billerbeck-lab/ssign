## ADDED Requirements

### Requirement: Skip the entire annotation phase with one flag
ssign SHALL provide a `--skip-annotation` flag that skips all optional annotation tools (BLASTp, HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam) at once, equivalent to passing all six individual `--skip-*` flags.

#### Scenario: Skip all annotation tools
- **WHEN** a user runs ssign with `--skip-annotation`
- **THEN** none of the six annotation tools SHALL execute
- **AND** the pipeline SHALL complete through prediction, substrate identification, and reporting, with the annotation columns left empty

#### Scenario: Predictions and substrate calling still run
- **WHEN** `--skip-annotation` is set
- **THEN** the prediction tools (DeepLocPro, DeepSecE, SignalP) and substrate identification SHALL still run (annotation is optional; prediction is core)

#### Scenario: Combines with individual flags
- **WHEN** `--skip-annotation` is combined with an individual annotation flag (e.g. `--no-skip-eggnog` if such an inverse exists)
- **THEN** the individual flag SHALL take precedence for its tool, and the rest SHALL remain skipped
