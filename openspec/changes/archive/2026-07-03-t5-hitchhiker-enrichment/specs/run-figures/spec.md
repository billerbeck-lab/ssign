## MODIFIED Requirements

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

## ADDED Requirements

### Requirement: T5a/c self and hitchhiker figure separation
Figures SHALL distinguish the two T5aSS/T5cSS populations: the autotransporter component itself ("self", `substrate_source == "T5SS-self"`) and its secreted-predicted neighbours ("hitchhiker", proximity substrates with a T5a/c nearby type). In the curated count, physicochemical, and functional-category figures the two SHALL be shown as separate categories (e.g. `T5aSS (self)` vs `T5aSS (hitch)`). In the enrichment fold/significance figures (per-tool and combined) each T5a/c type SHALL be drawn as two adjacent x-groups, the self group and the hitchhiker group.

#### Scenario: Annotation figures split self from hitchhiker
- **WHEN** a curated count or functional-category figure includes T5aSS or T5cSS substrates of both kinds
- **THEN** the self components and the hitchhiker neighbours SHALL appear as distinct labelled categories, not merged into one T5a/c bar

#### Scenario: Enrichment figure shows two T5a/c groups
- **WHEN** the enrichment fold/significance figure is rendered for a genome with a T5aSS or T5cSS system
- **THEN** that type SHALL appear as two adjacent x-groups, one for the self result and one for the hitchhiker (window) result

## REMOVED Requirements

### Requirement: Autotransporter self-detection figure
**Reason**: The autotransporter self-detection signal is now carried statistically by the enrichment `self`-mode result and its fold/significance bars; the standalone DLP extracellular-vs-outer-membrane scatter is redundant and was removed (Teo's call).
**Migration**: Read the T5aSS/T5cSS `mode=self` rows in `*_enrichment_stats.tsv` and the self group in the enrichment fold figure instead of the deleted `02_autotransporter_self_detection.png`.
