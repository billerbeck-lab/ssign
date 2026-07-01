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

### Requirement: Combined bar uses SignalP for autotransporters
The combined (one-bar-per-type) enrichment figure SHALL show, per secretion-system type, the DeepLocPro-or-DeepSecE combined score for window types, and the SignalP score alone for autotransporter types (T5aSS/T5cSS). The combined positivity for autotransporters SHALL NOT be a DLP-or-DSE-or-SignalP union: SignalP-alone has a lower genome background and so yields a stronger, more often significant bar over the few autotransporter loci, whereas the union raises the background and weakens it.

#### Scenario: Autotransporter combined bar is the SignalP score
- **WHEN** the combined enrichment figure is rendered for a run with T5aSS/T5cSS systems
- **THEN** the combined bar for those types SHALL be the SignalP circular-shift score (self-detection), while window types SHALL show the DeepLocPro-or-DeepSecE combined score
- **AND** the combined bar for autotransporters SHALL NOT be a three-way DLP-or-DSE-or-SignalP union
