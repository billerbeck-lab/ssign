## MODIFIED Requirements

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

<!-- Dropped from this change (Teo, 2026-06-24): the planned "Genome-wide pooled
null for genome groups" ADDED requirement (a smooth per-predictor circular-shift
null histogram pooled across genomes). Superseded by the combined fold/significance
bar chart above, which carries the same signal at both single-genome and group
scale without the extra genome-wide-null computation, npz, or histogram. No
genome-wide pooled-null is computed; enrichment_testing keeps its per-(SS-type)
circular-shift stats unchanged. -->
