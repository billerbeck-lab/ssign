## ADDED Requirements

### Requirement: Circular-shift permutation enrichment test
The enrichment test SHALL quantify, per SS type, whether secreted-predicted proteins are spatially enriched around that type's secretion-system components using a genome-structure-preserving circular-shift null. It SHALL emit an enrichment score (fold = observed count / null-mean count) and a significance score (permutation p-value with Benjamini-Hochberg q-value). The null SHALL be the set of circular rotations of the predictor's genome-ordered positivity vector (all rotations by default, exact; a configured sample of rotations MAY be used for very large genomes). The binomial test SHALL NOT be used as the significance source.

#### Scenario: Per-type fold and significance
- **WHEN** an enrichment-enabled run completes for a genome containing one or more systems of an SS type
- **THEN** the enrichment output SHALL contain, for that SS type and predictor (DLP, DSE), the observed count, the null mean, the fold, the permutation p-value, and the BH q-value

#### Scenario: Permutation null preserves gene order
- **WHEN** the null distribution is generated for a predictor
- **THEN** each replicate SHALL be a circular rotation of that predictor's genome-ordered positivity vector (preserving secreted-gene clustering), and the observed count SHALL correspond to the unrotated (zero-offset) arrangement

### Requirement: Enrichment couples whole-genome predictions and warns on runtime
When the enrichment test is enabled, the run SHALL produce DLP and DSE predictions for every protein in the genome (whole-genome), because the circular-shift null requires each gene's positivity in gene order. The run SHALL surface a note that enabling enrichment increases runtime.

#### Scenario: Enabling enrichment forces whole-genome predictors
- **WHEN** a run is started with the enrichment test enabled
- **THEN** DLP and DSE SHALL run on the whole proteome regardless of the neighborhood-only default
- **AND** the user SHALL be shown a note that enrichment increases runtime

#### Scenario: Enrichment disabled leaves predictor scope unchanged
- **WHEN** a run is started without the enrichment test enabled
- **THEN** DLP and DSE SHALL run on their default scope (neighborhood) and no enrichment runtime note SHALL be shown

### Requirement: Per-type aggregation with autotransporter self-detection
The enrichment test SHALL pool a genome's systems of each SS type. For window types (T1SS, T2SS, T3SS, T4SS, pT4SSt, T5bSS, T6SSi, T6SSii) the observed signal SHALL be secreted-predicted proteins within the configured gene window (default +/-3) of any component of that type. For autotransporter types (T5aSS, T5cSS) the observed signal SHALL be self-detection of the component itself (outer-membrane probability >= the configured threshold OR extracellular probability >= the threshold for DLP; secreted SS-type call for DSE), because the autotransporter is both machinery and substrate. DSE SHALL NOT be tested for T3SS.

#### Scenario: Window type uses neighborhood
- **WHEN** the tested SS type is a window type
- **THEN** the observed count SHALL be secreted-predicted proteins within the configured gene window of that type's components

#### Scenario: Autotransporter type uses self-detection
- **WHEN** the tested SS type is T5aSS or T5cSS
- **THEN** the observed signal SHALL be whether the component protein itself is predicted outer-membrane-or-extracellular (DLP) or secreted-typed (DSE), not a neighborhood window

#### Scenario: T3SS excluded from DSE
- **WHEN** the predictor is DSE
- **THEN** no enrichment row SHALL be produced for T3SS

### Requirement: Per-type null-distribution figure
An enrichment-enabled run SHALL emit one figure summarizing the circular-shift nulls: a grid with one panel per tested SS type, each panel showing the null distribution, the observed value, the null mean, and the fold and p-value, for the DLP and DSE predictors.

#### Scenario: Figure emitted with enrichment
- **WHEN** an enrichment-enabled run completes
- **THEN** a per-SS-type circular-shift null-distribution figure SHALL be written to the run's figures output

## MODIFIED Requirements

### Requirement: Exact background when predictions are whole-genome
Because the enrichment test forces whole-genome DLP/DSE predictions, the genome background positive rate SHALL be computed from ALL non-neighborhood, non-component proteins rather than a random subsample.

#### Scenario: Whole-genome predictor run
- **WHEN** the enrichment test runs
- **THEN** the background rate SHALL be computed over the full set of non-neighborhood, non-component proteins
- **AND** no random subsampling SHALL be applied

## REMOVED Requirements

### Requirement: Enrichment background sample size
**Reason**: The circular-shift test forces whole-genome DLP/DSE predictions, so the background is always computed exactly over the full non-neighborhood pool. The 1000-protein sampled-background path no longer drives per-run enrichment.
**Migration**: No action needed. Enrichment-enabled runs now always use the exact whole-genome background; the null-size knob has no effect on the circular-shift test.
