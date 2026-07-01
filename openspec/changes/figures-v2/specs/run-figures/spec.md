## MODIFIED Requirements

### Requirement: Curated per-run figure set
A completed run SHALL emit a curated, contiguously-numbered set of figures (`01`-`07`) derived from the integrated substrate CSV: `01` secreted proteins by SS type, `02` autotransporter (T5aSS/T5cSS) self-detection, `03` size & physicochemical properties by SS type (protein length plus each ProtParam property, one violin panel each), and `04`-`07` the functional-category figures. The `01` figure SHALL adapt to scope: a single genome SHALL be shown as one bar per SS type (a count histogram); a group SHALL be shown as one stacked bar per genome (this is also the cross-genome overview). The protein length figure SHALL NOT be a separate figure (it is the first panel of `03`). SS-type labels in every figure SHALL be the canonical display types (subtypes collapsed: T6SSi/ii/iii -> T6SS, pT4SSt -> T4SS; T5 subtypes kept distinct). Secretion-evidence, localization-confidence, SignalP-positive-by-type and annotation-tool-coverage figures SHALL NOT be emitted. Each figure SHALL be column-guarded and skipped with a logged note when its source columns are absent.

#### Scenario: Single-genome vs group counts figure
- **WHEN** a single-genome run completes
- **THEN** figure `01` SHALL show one bar per SS type (a count histogram)
- **AND WHEN** a run over two or more genomes completes, figure `01` SHALL show one stacked-by-SS-type bar per genome

#### Scenario: Length lives with the physicochemical figure
- **WHEN** the run emits figure `03`
- **THEN** it SHALL include a protein-length panel alongside the physicochemical-property panels, and no standalone length figure SHALL be emitted

#### Scenario: Canonical SS-type labels
- **WHEN** a figure displays SS types
- **THEN** subtype labels SHALL be collapsed to canonical display types (e.g. `T6SSi` shown as `T6SS`)

### Requirement: Pooled figures for multi-genome runs
A run over two or more genomes SHALL emit, into the top-level figures output, the curated set computed over all genomes combined, named `0N_pooled_<name>.png`. Figure `01_pooled_*` (one stacked bar per genome) IS the cross-genome overview; no separate `P0N_` cross-genome figures are emitted. A single-genome run SHALL NOT emit pooled figures.

#### Scenario: Pooled set emitted for a group
- **WHEN** a run over two or more genomes completes
- **THEN** the top-level figures directory SHALL contain `0N_pooled_*` versions of the curated figures, over all genomes combined, including `01_pooled_*` as the per-genome stacked overview
- **AND** no `P0N_`-prefixed figures SHALL be emitted

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
A run SHALL emit functional-category figures of secreted proteins, stacked by SS type, from controlled vocabularies present in the integrated CSV: COG category (`cog_category`, `04`), KEGG function (`kegg_ko`, `05`), EggNOG description (`06`), and the curated consensus annotation (`07`). KEGG codes SHALL be shown as human-readable function names via a bundled KO→name table (raw `K#####` IDs only as a fallback). Bounded vocabularies (COG, consensus) SHALL show every present category; high-cardinality vocabularies (KEGG, EggNOG) MAY cap at a top-N with the remainder collapsed to "Other". A figure whose source column is absent (e.g. no EggNOG ran) SHALL be skipped with a logged note. The integrated CSV SHALL also carry readable `cog_category_name` and `kegg_function` columns per protein.

#### Scenario: COG-category figure from existing data
- **WHEN** a run whose integrated CSV carries `cog_category` completes
- **THEN** a COG-category figure of secreted proteins stacked by SS type SHALL be emitted, showing every present category (no "Other" bucket)

#### Scenario: KEGG figure shows function names
- **WHEN** a run whose integrated CSV carries `kegg_ko` completes
- **THEN** the KEGG figure SHALL label categories by human-readable function name (not bare KO IDs), and the integrated CSV SHALL carry a `kegg_function` column

#### Scenario: Functional figure skipped when source absent
- **WHEN** the integrated CSV lacks a functional source column (e.g. `kegg_ko`)
- **THEN** that functional figure SHALL be skipped with a logged note and the run SHALL continue
