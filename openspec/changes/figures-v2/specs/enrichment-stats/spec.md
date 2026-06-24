## MODIFIED Requirements

### Requirement: Per-type null-distribution figure
An enrichment-enabled single-genome run SHALL emit a per-SS-type **fold-enrichment** figure: a bar chart of fold (observed / null-mean) per SS type for the DeepLocPro and DeepSecE predictors, annotated with BH q-value significance, with autotransporters shown by self-detection. The spiky per-type null *histogram* grid SHALL NOT be the per-genome presentation (a single genome's per-type count is a small integer, so its null is uninformative as a histogram). The underlying statistics (observed, null mean, fold, p, q) are unchanged.

#### Scenario: Per-genome fold bars
- **WHEN** an enrichment-enabled single-genome run completes
- **THEN** the figures output SHALL contain a per-SS-type fold-enrichment bar chart (DLP and DSE), with non-significant types visually de-emphasized

## ADDED Requirements

### Requirement: Genome-wide pooled null for genome groups
The enrichment test SHALL additionally compute a genome-wide pooled null per predictor: the count of positives within the configured window of ANY window-type secretion-system component, over all circular rotations (a single pooled mask, not per type). For a multi-genome run, these per-genome genome-wide nulls SHALL be pooled across all genomes and rendered as one smooth DeepLocPro and one DeepSecE null-distribution histogram with the observed value and fold. A single-genome run SHALL NOT require the pooled-group figure.

#### Scenario: Group-level smooth null figure
- **WHEN** an enrichment-enabled run over two or more genomes completes
- **THEN** the top-level figures output SHALL contain a genome-wide pooled circular-shift null figure (one DLP, one DSE histogram) summed across all genomes, with the observed count and fold annotated

#### Scenario: Genome-wide null emitted alongside per-type stats
- **WHEN** the enrichment test runs on a genome
- **THEN** it SHALL emit, in addition to the per-SS-type stats, a genome-wide pooled null (positives near any window-type component) for each predictor
