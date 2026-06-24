## ADDED Requirements

### Requirement: Corrected 66-genome panel

The rerun SHALL be defined over a corrected genome panel of 66 genomes: the 51 input genomes that carry at least one testable proximity-system answer-key entry, plus the 15 T5SS genomes that are absent from the current fleet, and SHALL exclude the 16 input genomes that carry no testable answer-key system. The exact panel MUST be emitted as a committed manifest (genome accession + reason-for-inclusion) so the set is reproducible.

#### Scenario: Panel membership is derived, not hand-listed

- **WHEN** the panel manifest is built
- **THEN** each of the 66 genomes is included because it carries a testable proximity system or a T5SS answer-key effector, and each of the 16 excluded genomes is recorded with the reason "no testable answer-key system"

#### Scenario: Missing T5SS genomes are added

- **WHEN** the panel is compared to the current 67-genome fleet
- **THEN** the 15 T5SS genomes that were never run in the fleet appear in the panel

### Requirement: Extended-tier run configuration

Each panel genome SHALL be run through ssign at extended tier (tier 2) on CX3 via the `cx3-submit` skill with: `--enrichment-stats` (which forces whole-genome DLP+DSE), T3SS detection enabled (the current default), local DeepLocPro/SignalP predictors (CX3 compute nodes cannot reach the DTU webserver), and the extended-tier databases resolved from `$EPHEMERAL`. DeepSecE MUST remain excluded from T3SS calls (critical bug-fix #4).

#### Scenario: Enrichment is computed fleet-wide

- **WHEN** the rerun completes
- **THEN** every genome has a per-genome enrichment-stats table and (for the multi-genome run) a pooled enrichment table, each with per-SS-type fold, permutation p, and BH q

#### Scenario: T3SS is detected without the flagellar blowup

- **WHEN** a genome with a validated T3SS is run
- **THEN** ssign emits a sane number of T3SS substrates (not a flagellar-misclassification blowup) and a T3SS DLP enrichment row, with no T3SS DeepSecE contribution

### Requirement: Reconciliation against the benchmark figures

The rerun outputs SHALL be reconciled against the recall (06), annotation (07/08), and enrichment figures. The four RTX T1SS toxins (HlyA, ApxIA, LtxA, LktA) MUST be confirmed to emit from their real full-assembly runs (retiring the `clean_dataset` staging fix), and the 15 previously-ungraded T5SS effectors MUST be annotation-graded from the rerun.

#### Scenario: Staging fix is retired by real emission

- **WHEN** the RTX T1SS toxins are run on their full assemblies in the rerun
- **THEN** ssign emits them as secreted, and the recall figures no longer depend on the manual `APPLY_T1SS_STAGING_FIX` rescue

#### Scenario: T5SS coverage is completed

- **WHEN** the rerun includes the previously-missing T5SS genomes
- **THEN** the annotation figure grades T5SS effectors from real ssign annotation rather than only the 5 previously available
