## ADDED Requirements

### Requirement: SignalP circular-shift enrichment predictor
The enrichment test SHALL evaluate SignalP as a third predictor alongside DeepLocPro and DeepSecE, for every secretion-system type. A protein SHALL be SignalP-positive when SignalP assigns it a signal-peptide class (the exact class set is fixed in design). For window types the SignalP track SHALL use the ±N-gene window mask; for autotransporter types (T5aSS/T5cSS) it SHALL use self-detection. Each (SS type, SignalP) test SHALL emit observed, null mean, fold, permutation p, and BH q, in the same per-tool multiple-testing family as DeepLocPro and DeepSecE.

#### Scenario: SignalP track emitted per type
- **WHEN** an enrichment-enabled run completes
- **THEN** the enrichment stats SHALL contain a SignalP row for each tested secretion-system type, carrying fold, p, and q computed by the same circular-shift test as DeepLocPro/DeepSecE

#### Scenario: Autotransporters enrich for SignalP
- **WHEN** a genome with T5aSS/T5cSS systems is tested
- **THEN** the SignalP track for those types SHALL use self-detection (component positions as the mask), reflecting that autotransporter passengers are Sec-exported and signal-peptide-bearing

#### Scenario: SignalP joins the per-tool BH family
- **WHEN** BH correction is applied
- **THEN** SignalP tests SHALL be corrected together with DeepLocPro and DeepSecE (the per-tool family), separate from the combined DLP-or-DSE family

### Requirement: Whole-genome SignalP forced under enrichment
Enabling `--enrichment-stats` SHALL force SignalP to run on the whole genome (the full proteome), the same way it already forces whole-genome DeepLocPro and DeepSecE. The circular-shift null rotates a whole-genome positivity vector, so every gene needs a SignalP call. SignalP SHALL run locally.

#### Scenario: Enrichment forces whole-genome SignalP
- **WHEN** a run is started with `--enrichment-stats`
- **THEN** SignalP SHALL be run on the whole proteome (not only the ±N neighbourhood), and its whole-genome output SHALL be the SignalP input to the enrichment test

#### Scenario: No webserver dependency
- **WHEN** whole-genome SignalP runs for enrichment
- **THEN** it SHALL run via the local SignalP install (no external web service)

### Requirement: SignalP bar in the per-tool enrichment figure
The per-tool enrichment fold/significance figure SHALL show a SignalP bar per secretion-system type alongside the DeepLocPro and DeepSecE bars, in the shared house style (its own theme colour, non-significant bars de-emphasised, fold and significance annotated). Each bar is that tool's own independent circular-shift score.

#### Scenario: Three bars per type
- **WHEN** the per-tool enrichment figure is rendered for an enrichment-enabled run
- **THEN** each secretion-system type SHALL show DeepLocPro, DeepSecE, and SignalP bars (a type with no valid predictor for one tool shows only the applicable bars)

### Requirement: Combined bar pairs DLP with SignalP for T5SS
The combined (one-bar-per-type) enrichment figure SHALL show, per secretion-system type, a union of the two relevant predictors (a gene counts if either flags it, scored by the same permutation test): DeepLocPro-or-DeepSecE for non-T5 types, and DeepLocPro-or-SignalP for ALL T5SS subtypes (T5aSS, T5bSS, T5cSS). DeepSecE SHALL NOT contribute to the T5SS combined bar (it is unreliable for T5). The combined bar SHALL NOT be an average of the individual folds.

#### Scenario: T5SS combined bar is DLP-or-SignalP
- **WHEN** the combined enrichment figure is rendered for a run with T5SS systems
- **THEN** the combined bar for every T5SS subtype SHALL be scored from the DeepLocPro-or-SignalP union (a locus counts if DeepLocPro OR SignalP is positive), while non-T5 types SHALL show the DeepLocPro-or-DeepSecE union
- **AND** a T5 locus positive by only one of the two (e.g. SignalP-positive but DeepLocPro-negative) SHALL still count toward the combined bar
