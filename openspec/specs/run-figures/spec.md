# run-figures Specification

## Purpose
Defines the figures ssign auto-generates per run, single-genome and pooled multi-genome: their content, the shared house-style conventions (theme, numbering, figure index), column-guarded non-fatal generation, and which are emitted by default.
## Requirements
### Requirement: Curated per-run figure set
A completed run SHALL emit a curated, contiguously-numbered set of figures (`01`-`06`) derived from the integrated substrate CSV: `01` secreted proteins by SS type, `02` size & physicochemical properties by SS type (protein length plus each ProtParam property, one violin panel each), and `03`-`06` the functional-category figures (COG, KEGG, EggNOG, curated consensus). No standalone autotransporter self-detection figure SHALL be emitted (the self signal is carried by the enrichment self-mode result). The `01` figure SHALL adapt to scope: a single genome SHALL be shown as one bar per SS type (a count histogram, without per-bar value labels); a group SHALL be shown as one stacked bar per genome (this is also the cross-genome overview). The protein length figure SHALL NOT be a separate figure (it is the first panel of `02`). SS-type labels in every figure SHALL be the canonical display types (subtypes collapsed: T6SSi/ii/iii -> T6SS, pT4SSt -> T4SS; T5 subtypes kept distinct). Secretion-evidence, localization-confidence, SignalP-positive-by-type and annotation-tool-coverage figures SHALL NOT be emitted. Each figure SHALL be column-guarded and skipped with a logged note when its source columns are absent.

#### Scenario: Single-genome vs group counts figure
- **WHEN** a single-genome run completes
- **THEN** figure `01` SHALL show one bar per SS type (a count histogram) with no per-bar value labels
- **AND WHEN** a run over two or more genomes completes, figure `01` SHALL show one stacked-by-SS-type bar per genome

#### Scenario: No autotransporter self-detection figure
- **WHEN** a run with one or more T5aSS/T5cSS components completes
- **THEN** the figures directory SHALL NOT contain an autotransporter self-detection scatter figure
- **AND** the curated set SHALL be numbered `01`-`06` with no gap

#### Scenario: Canonical SS-type labels
- **WHEN** a figure displays SS types
- **THEN** subtype labels SHALL be collapsed to canonical display types (e.g. `T6SSi` shown as `T6SS`)

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
A run over two or more genomes SHALL emit, into the top-level figures output, the curated set computed over all genomes combined, named `0N_pooled_<name>.png`. Figure `01_pooled_*` (one stacked bar per genome) IS the cross-genome overview; no separate `P0N_` cross-genome figures are emitted. A single-genome run SHALL NOT emit pooled figures.

#### Scenario: Pooled set emitted for a group
- **WHEN** a run over two or more genomes completes
- **THEN** the top-level figures directory SHALL contain `0N_pooled_*` versions of the curated figures, over all genomes combined, including `01_pooled_*` as the per-genome stacked overview
- **AND** no `P0N_`-prefixed figures SHALL be emitted

### Requirement: Enrichment null figure renders discrete counts legibly
The circular-shift null figure SHALL render the discrete integer null counts with integer-aligned bins (one bar per achievable count) and SHALL visually de-emphasize panels that are not significant (q >= 0.05), without altering the observed value, null mean, fold, p-value, or q-value displayed.

#### Scenario: Integer-aligned discrete null
- **WHEN** the enrichment figure is rendered
- **THEN** each panel's null SHALL use integer-aligned bins
- **AND** non-significant panels SHALL be visually muted while still showing their statistics

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

### Requirement: T5a/c self and hitchhiker figure separation
Figures SHALL distinguish the two T5aSS/T5cSS populations: the autotransporter component itself ("self", `substrate_source == "T5SS-self"`) and its secreted-predicted neighbours ("hitchhiker", proximity substrates with a T5a/c nearby type). In the curated count, physicochemical, and functional-category figures the two SHALL be shown as separate categories (e.g. `T5aSS (self)` vs `T5aSS (hitch)`). In the enrichment fold/significance figures (per-tool and combined) each T5a/c type SHALL be drawn as two adjacent x-groups, the self group and the hitchhiker group.

#### Scenario: Annotation figures split self from hitchhiker
- **WHEN** a curated count or functional-category figure includes T5aSS or T5cSS substrates of both kinds
- **THEN** the self components and the hitchhiker neighbours SHALL appear as distinct labelled categories, not merged into one T5a/c bar

#### Scenario: Enrichment figure shows two T5a/c groups
- **WHEN** the enrichment fold/significance figure is rendered for a genome with a T5aSS or T5cSS system
- **THEN** that type SHALL appear as two adjacent x-groups, one for the self result and one for the hitchhiker (window) result

