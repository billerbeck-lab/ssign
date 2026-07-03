# enrichment-stats Specification

## Purpose
Defines ssign's per-run secretion-enrichment statistical test: a per-SS-type circular-shift permutation that quantifies whether secreted-predicted proteins cluster around each secretion system's components, emitting fold enrichment, a permutation p-value, and a BH q-value per type, plus a per-type null-distribution figure. The test uses DLP and DSE only (not PLM-Effector) and forces whole-genome predictions when enabled.
## Requirements
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
The enrichment test SHALL pool a genome's systems of each SS type. For window types (T1SS, T2SS, T3SS, T4SS, pT4SSt, T5bSS, T6SSi, T6SSii) the observed signal SHALL be secreted-predicted proteins within the configured gene window (default +/-3) of any component of that type. For autotransporter types (T5aSS, T5cSS) the test SHALL emit **two** results distinguished by the `mode` column: a `self` result (the component itself is predicted outer-membrane-or-extracellular for DLP, or secreted-typed for DSE, because the autotransporter is both machinery and substrate) AND a `window` "hitchhiker" result (secreted-predicted proteins within the configured gene window of that type's components, testing the hypothesis that neighbours piggyback through the T5 pore). The hitchhiker window SHALL use the same window mask, predictors, and combined rule as other window types. DSE SHALL NOT be tested for T3SS.

#### Scenario: Window type uses neighborhood
- **WHEN** the tested SS type is a window type
- **THEN** the observed count SHALL be secreted-predicted proteins within the configured gene window of that type's components

#### Scenario: Autotransporter type emits self and hitchhiker results
- **WHEN** the tested SS type is T5aSS or T5cSS
- **THEN** the enrichment output SHALL contain a `mode=self` result (component self-detection) AND a `mode=window` hitchhiker result (secreted-predicted proteins within the configured gene window of that type's components), for each applicable predictor

#### Scenario: Hitchhiker combined track uses the window rule
- **WHEN** the combined (union) track is computed for a T5aSS/T5cSS hitchhiker (`mode=window`) result
- **THEN** it SHALL be the DLP-or-DSE union, the same rule as other window types (whereas the `mode=self` combined track remains the DLP-or-SignalP union)

#### Scenario: T3SS excluded from DSE
- **WHEN** the predictor is DSE
- **THEN** no enrichment row SHALL be produced for T3SS

### Requirement: Exact background when predictions are whole-genome
Because the enrichment test forces whole-genome DLP/DSE predictions, the genome background positive rate SHALL be computed from ALL non-neighborhood, non-component proteins rather than a random subsample.

#### Scenario: Whole-genome predictor run
- **WHEN** the enrichment test runs
- **THEN** the background rate SHALL be computed over the full set of non-neighborhood, non-component proteins
- **AND** no random subsampling SHALL be applied

### Requirement: Enrichment predictor set excludes PLM-Effector
The enrichment test SHALL test DLP and DSE only and SHALL NOT include PLM-Effector as an enrichment predictor.

#### Scenario: Enrichment output columns
- **WHEN** the enrichment test runs on any genome
- **THEN** the output SHALL contain rows for tools DLP and DSE only
- **AND** SHALL contain no PLME rows even if a PLM-Effector output file is present

### Requirement: Enrichment fold/significance figures
An enrichment-enabled run SHALL present its per-SS-type circular-shift results as fold-enrichment (observed / null-mean) bar charts with statistical significance. TWO figures SHALL be emitted: (a) a **per-tool** chart with a DeepLocPro and a DeepSecE bar per secretion-system type, and (b) a **combined** chart with a single "DLP or DSE" bar per type. The combined track is a third predictor in the stats: a protein is positive if EITHER valid predictor calls it secreted (DSE dropped where unreliable, e.g. T3SS), run through the same circular-shift test and BH-corrected as its own multiple-testing family. Non-significant bars SHALL be visually de-emphasized and significance SHALL be annotated from the BH q-value. Autotransporter types (T5aSS/T5cSS) SHALL be shown via their self-detection enrichment and labelled as such. The spiky per-type null *histogram grid* SHALL NOT be emitted. The DLP/DSE per-type statistics (observed, null mean, fold, p, q) are unchanged.

A single-genome run SHALL emit both charts from its own stats; a multi-genome run SHALL emit pooled versions computed from the per-genome stats combined.

#### Scenario: Per-genome per-tool and combined bars
- **WHEN** an enrichment-enabled single-genome run completes
- **THEN** the figures output SHALL contain a per-tool fold-enrichment bar chart (DLP and DSE per SS type) AND a separate combined (DLP-or-DSE) fold-enrichment bar chart, both with non-significant types de-emphasized and significance annotated

#### Scenario: Combined track is its own BH family
- **WHEN** the enrichment test runs
- **THEN** the combined (DLP-or-DSE) tests SHALL be BH-corrected separately from the per-tool (DLP/DSE) tests, so the per-tool q-values are unaffected by the combined track

#### Scenario: Pooled fold/significance bars for a genome group
- **WHEN** an enrichment-enabled run over two or more genomes completes
- **THEN** the top-level figures output SHALL contain pooled per-tool AND combined fold-enrichment bar charts combining the per-genome enrichment stats

