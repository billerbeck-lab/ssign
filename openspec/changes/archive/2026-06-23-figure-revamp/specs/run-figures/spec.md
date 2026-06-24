## ADDED Requirements

### Requirement: Curated per-run figure set
A completed single-genome run SHALL emit a curated set of substrate figures derived from the integrated substrate CSV: substrates per SS type, substrate evidence strength (confidence tier and number of tools agreeing), localization confidence (DLP extracellular probability by SS type), SignalP-positive fraction per SS type, annotation tool coverage, protein length by SS type, and functional categories per SS type. Figure filenames SHALL be zero-padded and numbered in reading order (context, then main result, then drill-down).

#### Scenario: Numbered curated figures emitted
- **WHEN** a run with a non-empty integrated CSV completes
- **THEN** the run's figures directory SHALL contain the numbered curated figures (e.g. `01_*.png` ... `07_*.png`)
- **AND** SHALL contain no `figN_`-prefixed files

#### Scenario: Figure index printed
- **WHEN** the figure step completes
- **THEN** a figure index listing each figure's number and name SHALL be printed to the run log

### Requirement: All figures share one house style
Every auto-generated figure (per-run, pooled cross-genome, and the enrichment null figure) SHALL draw colour, font, and axis style from a single shared style definition. A given SS type SHALL render in the same colour across all figures, and no figure SHALL use an ad-hoc inline palette.

#### Scenario: Consistent SS-type colour
- **WHEN** two different figures in the same run both display an SS type
- **THEN** that SS type SHALL be drawn in the same colour in both, sourced from the shared style module

### Requirement: Physicochemical figure off by default
The physicochemical-properties figure SHALL NOT be emitted by default, because the standard integrated CSV does not carry the physicochemical columns. It MAY be produced only when explicitly enabled AND the required columns are present.

#### Scenario: Default run omits physicochemical figure
- **WHEN** a default run completes
- **THEN** no physicochemical-properties figure SHALL be written

#### Scenario: Enabled but columns absent
- **WHEN** the physicochemical figure is explicitly enabled but its required columns are absent from the CSV
- **THEN** the figure SHALL be skipped with a logged note and no error

### Requirement: Figure generation is non-fatal and column-guarded
Each figure SHALL check for the columns it requires and, if they are absent or empty, SHALL be skipped with a logged note rather than raising. A figure-generation failure SHALL NOT fail the run.

#### Scenario: Missing column skips one figure
- **WHEN** a figure's required column is missing from the integrated CSV
- **THEN** that figure SHALL be skipped and the remaining figures and the run SHALL continue

### Requirement: Pooled cross-genome figures for multi-genome runs
A run over two or more genomes SHALL emit pooled cross-genome figures into the top-level figures output, in addition to the per-genome figures: substrate counts per genome and an SS-type by genome substrate-count view. A single-genome run SHALL NOT emit pooled cross-genome figures.

#### Scenario: Multi-genome pooled figures
- **WHEN** a run over two or more genomes completes
- **THEN** the top-level figures directory SHALL contain the pooled cross-genome figures

#### Scenario: Single genome emits none
- **WHEN** a run over a single genome completes
- **THEN** no pooled cross-genome figures SHALL be emitted

### Requirement: Enrichment null figure renders discrete counts legibly
The circular-shift null figure SHALL render the discrete integer null counts with integer-aligned bins (one bar per achievable count) and SHALL visually de-emphasize panels that are not significant (q >= 0.05), without altering the observed value, null mean, fold, p-value, or q-value displayed.

#### Scenario: Integer-aligned discrete null
- **WHEN** the enrichment figure is rendered
- **THEN** each panel's null SHALL use integer-aligned bins
- **AND** non-significant panels SHALL be visually muted while still showing their statistics
