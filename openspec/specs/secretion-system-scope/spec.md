# secretion-system-scope Specification

## Purpose
TBD - created by archiving change figures-v2. Update Purpose after archive.
## Requirements
### Requirement: Non-secretion appendages excluded by default
ssign SHALL NOT substrate-call TXSScan models that are surface/uptake appendages rather than protein secretion systems. The default `excluded_systems` SHALL include `Flagellum, Tad, T4aP, T4bP, MSH, ComM, Archaeal-T4P`. Genuine secretion systems (T1SS, T2SS, T3SS, T4SS/pT4SSi/pT4SSt, T5aSS, T5bSS, T5cSS, T6SS, T9SS) SHALL NOT be excluded by default. The user MAY override the exclusion list via `--excluded-systems`.

#### Scenario: Type IV pilus not substrate-called by default
- **WHEN** ssign runs on a genome with a detected `T4aP` (or `T4bP`/`MSH`/`ComM`/`Archaeal-T4P`) and no `--excluded-systems` override
- **THEN** that system SHALL be excluded and no proteins SHALL be substrate-called as its substrates
- **AND** no figure SHALL show that appendage as a secretion-system type

#### Scenario: Genuine secretion systems remain in scope
- **WHEN** ssign runs with the default exclusion list
- **THEN** T9SS, pT4SSi, and pT4SSt SHALL still be detected and substrate-called

#### Scenario: User can re-include an appendage
- **WHEN** the user passes an `--excluded-systems` list that omits `T4aP`
- **THEN** T4aP SHALL be detected and substrate-called as before this change

