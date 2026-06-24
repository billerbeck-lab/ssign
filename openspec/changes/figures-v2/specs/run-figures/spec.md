## MODIFIED Requirements

### Requirement: Curated per-run figure set
A completed single-genome run SHALL emit a curated set of substrate figures derived from the integrated substrate CSV. Cargo-secreting types (T1SS, T2SS, T3SS, T4SS, T6SS) and autotransporters (T5aSS, T5cSS) SHALL be handled separately: the substrates-per-type, localization-confidence (DLP extracellular by type), and protein-length figures SHALL cover cargo types only (T5aSS/T5cSS excluded), while T5aSS/T5cSS SHALL be presented by a dedicated autotransporter self-detection figure. SS-type labels in every figure SHALL be the canonical display types (subtypes collapsed: T6SSi/ii/iii -> T6SS, pT4SSt -> T4SS; T5 subtypes kept distinct). SignalP-positive-by-type and annotation-tool-coverage figures SHALL NOT be emitted. Figure filenames SHALL be zero-padded and numbered in reading order.

#### Scenario: Cargo figures exclude autotransporters
- **WHEN** a run with detected T5aSS/T5cSS systems completes
- **THEN** the substrates-per-type, localization, and length figures SHALL NOT include T5aSS or T5cSS
- **AND** the SignalP-by-type and tool-coverage figures SHALL NOT be present

#### Scenario: Canonical SS-type labels
- **WHEN** a figure displays SS types
- **THEN** subtype labels SHALL be collapsed to canonical display types (e.g. `T6SSi` shown as `T6SS`)

### Requirement: Pooled cross-genome figures for multi-genome runs
A run over two or more genomes SHALL emit, into the top-level figures output: (a) the existing cross-genome figures (substrate counts per genome, an SS-type by genome substrate-count view, pooled evidence), AND (b) pooled versions of the per-genome curated set computed over all genomes combined, named `0N_pooled_<name>.png`. A single-genome run SHALL NOT emit pooled figures.

#### Scenario: Pooled per-genome set emitted for a group
- **WHEN** a run over two or more genomes completes
- **THEN** the top-level figures directory SHALL contain `0N_pooled_*` versions of the curated per-genome figures, over all genomes combined
- **AND** SHALL also contain the cross-genome substrates-per-genome and SS-type-by-genome figures

## ADDED Requirements

### Requirement: Autotransporter self-detection figure
A run with detected T5aSS or T5cSS systems SHALL emit a dedicated autotransporter figure showing, per component, its SignalP call, DLP extracellular probability, and DLP outer-membrane probability (the self-detection signal: an autotransporter is the substrate, passing on extracellular OR outer-membrane localization). The integrated CSV SHALL carry `outer_membrane_prob` for this figure.

#### Scenario: Autotransporter figure emitted
- **WHEN** a run with one or more T5aSS/T5cSS components completes
- **THEN** the figures directory SHALL contain the autotransporter self-detection figure showing SignalP, extracellular probability, and outer-membrane probability per component

#### Scenario: Outer-membrane probability available
- **WHEN** the integrated CSV is produced
- **THEN** it SHALL contain an `outer_membrane_prob` column carried from the DLP predictions

### Requirement: Functional-category figures from controlled vocabularies
A run SHALL be able to emit functional-category figures of secreted proteins from controlled vocabularies present in the integrated CSV: COG category (`cog_category`), KEGG (`kegg_ko`), EggNOG description, and the broad-consensus annotation. Each functional figure SHALL be available in two scopes: overall (all secreted proteins) and stacked by SS type. A figure whose source column is absent (e.g. no EggNOG ran) SHALL be skipped with a logged note.

#### Scenario: COG-category figure from existing data
- **WHEN** a run whose integrated CSV carries `cog_category` completes
- **THEN** a COG-category figure of secreted proteins SHALL be emitted, in both overall and per-SS-type scopes

#### Scenario: Functional figure skipped when source absent
- **WHEN** the integrated CSV lacks a functional source column (e.g. `kegg_ko`)
- **THEN** that functional figure SHALL be skipped with a logged note and the run SHALL continue
