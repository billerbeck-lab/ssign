# t3ss-detection Specification

## Purpose
Defines ssign's default handling of Type III Secretion Systems: T3SS is detected and substrate-called by default, but DeepSecE is never trusted for T3SS (it cannot distinguish injectisome from flagellar T3SS), so T3SS relies on MacSyFinder/TXSScan detection + DeepLocPro + proximity.

## Requirements

### Requirement: T3SS detected by default

ssign SHALL detect, neighbourhood-extract, and substrate-call Type III Secretion Systems by default. The default `excluded_systems` SHALL be `[Flagellum, Tad]` and SHALL NOT contain `T3SS`.

#### Scenario: Default run includes T3SS

- **WHEN** ssign runs on a genome with a MacSyFinder/TXSScan-validated T3SS and no `--excluded-systems` override
- **THEN** the T3SS system SHALL be validated, its components' neighbourhoods extracted, and its substrates called like any other non-excluded system

#### Scenario: User can still exclude T3SS

- **WHEN** the user passes `--excluded-systems Flagellum,Tad,T3SS`
- **THEN** T3SS SHALL be excluded exactly as before this change (the old default is recoverable)

### Requirement: DeepSecE is never trusted for T3SS

DeepSecE SHALL NOT contribute to a T3SS secretion call anywhere in the pipeline. Any DeepSecE prediction of type `T3SS` SHALL be flagged out of the secretion-evidence count regardless of whether the genome contains a validated T3SS. DeepLocPro (extracellular) is the trusted PLM localization signal for T3SS.

#### Scenario: DSE T3SS call in a genome WITH a validated T3SS

- **WHEN** DeepSecE predicts `T3SS` for a protein in a genome that has a validated T3SS
- **THEN** that DeepSecE prediction SHALL be flagged (`dse_T3SS_flagged = True`) and SHALL NOT count as secretion evidence

#### Scenario: DSE T3SS call in a genome WITHOUT a validated T3SS

- **WHEN** DeepSecE predicts `T3SS` for a protein in a genome with no validated T3SS
- **THEN** that DeepSecE prediction SHALL be flagged and SHALL NOT count as secretion evidence

#### Scenario: Non-T3SS DSE calls are unaffected

- **WHEN** DeepSecE predicts a non-T3SS secretion type (e.g. T1SS, T6SS)
- **THEN** that prediction SHALL be handled by the normal cross-validation logic, unchanged by this requirement

### Requirement: T3SS enrichment is DLP-only

When the enrichment test runs on a genome with a validated T3SS, it SHALL emit a T3SS enrichment row for DeepLocPro (window mode) and SHALL NOT emit a T3SS enrichment row for DeepSecE.

#### Scenario: Enrichment output for a T3SS-bearing genome

- **WHEN** `--enrichment-stats` runs on a genome with a validated T3SS
- **THEN** the enrichment stats SHALL contain a `T3SS` row with `tool = DLP` and `mode = window`
- **AND** SHALL NOT contain any `T3SS` row with `tool = DSE`
