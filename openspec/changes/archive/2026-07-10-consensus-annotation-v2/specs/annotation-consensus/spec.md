## ADDED Requirements

### Requirement: Tool-weighted functional consensus
The consensus SHALL assign each substrate one broad functional category by a weighted vote across its per-tool annotation descriptions, where each tool casts one vote weighted by its credibility tier: functional-name tools (BLASTp/Swiss-Prot, EggNOG, Bakta, GBFF) = 3, domain tools (HHpred_Pfam, InterProScan) = 2, structure-only tools (pLM-BLAST/ECOD, HHpred_PDB) = 1. Identical descriptions across tools SHALL be counted once, at the highest contributing weight.

#### Scenario: Functional tool outweighs a structural tool
- **WHEN** one tool's description classifies as Lipase (weight 2) and another classifies as Pore-forming toxin (weight 1)
- **THEN** the broad category is Lipase

#### Scenario: Duplicate descriptions counted once
- **WHEN** two tools emit the identical description string
- **THEN** it contributes a single vote at the higher of the two tool weights

### Requirement: Structural fold names cannot outvote function
A description that is itself a structural fold or superfamily name (matching patterns such as "… superfamily", "…-like domain", "beta-helix/barrel/propeller", "TIM barrel") SHALL have its vote weight capped to 1, regardless of the source tool.

#### Scenario: Fold-superfamily name is down-weighted
- **WHEN** InterProScan reports "Haem peroxidase domain superfamily" and a Tier-1 tool reports a functional category
- **THEN** the fold-derived category does not win over the functional category

### Requirement: Single specificity-ranked category per description
Classifying a description SHALL return exactly one category, chosen by a most-specific-first ordered rule set (first match wins), or none. When no rule matches, the classifier SHALL NOT invent a category from the description text.

#### Scenario: No match yields no minted label
- **WHEN** a description matches no category rule
- **THEN** it contributes no functional vote and never produces a title-case fallback category

### Requirement: Hypothetical and Unclassified are floors
`Hypothetical` and `Unclassified` SHALL never outrank a functional category. A substrate SHALL be reported as `Unclassified` only when no tool yields a functional or apparatus category, and `Hypothetical` only when the sole signal is a hypothetical/uncharacterised annotation.

#### Scenario: A single functional call beats multiple hypothetical calls
- **WHEN** two tools call a protein hypothetical and one tool yields Protease
- **THEN** the broad category is Protease

#### Scenario: Genuinely unannotated protein
- **WHEN** no tool yields a functional, apparatus, or hypothetical category
- **THEN** the broad category is Unclassified

### Requirement: Machinery detected by component identity, competing by weight
Secretion/uptake apparatus and translocator components SHALL be detected by component identity (e.g. VgrG, Hcp, PAAR, Tss*, Gsp*, Vir[BD], Sec/Tat, flagellar gene families with word boundaries; ShlB/FhaC/HecB/TpsB translocators; "haemolysin secretion/activation" context) and mapped to `Apparatus-associated`, which competes in the weighted vote. A protein SHALL NOT be routed to machinery merely because its annotation names a secretion-system type (e.g. "T6SS").

#### Scenario: Effector named after its system is not hijacked
- **WHEN** a T6SS amidase effector is annotated "T6SS amidase effector (Tae)"
- **THEN** the broad category is Peptidoglycan hydrolase, not Apparatus-associated

#### Scenario: Translocator agreed by good tools wins
- **WHEN** EggNOG and Pfam both describe a protein as a "haemolysin secretion protein" (Apparatus) and a structural hit suggests Adhesin
- **THEN** the broad category is Apparatus-associated

### Requirement: Literature-grounded category taxonomy
The category set SHALL cover the effector activity classes of the secretion systems, including Peptidoglycan hydrolase, ADP-ribosyltransferase, Glycosyltransferase, Phosphothreonine lyase, Ubiquitin-pathway, GTPase modulator, Protease/Peptidase, Lipase/Phospholipase, Nuclease, Glycoside hydrolase, Beta-lactamase, Pore-forming toxin, Adhesin, Hemophore/metal-uptake, S-layer, Autotransporter passenger, and Phage/mobile-element. There SHALL be no standalone "RTX toxin" category; RTX proteins SHALL be classified by their actual activity (serralysins → Protease; cytolysins/adenylate-cyclase toxins → Pore-forming toxin).

#### Scenario: RTX metalloprotease classified by activity
- **WHEN** an RTX serralysin metalloprotease is annotated as a peptidase
- **THEN** the broad category is Protease/Peptidase, not a generic RTX bucket

### Requirement: Figure bucket collapse tracks the category set
The figure-side `consensus_bucket` SHALL pass known categories through, route apparatus/flagellar machinery to `Apparatus-associated`, and route anything unrecognised (including blanks) to `Other`, staying consistent with the consensus category set.

#### Scenario: Unknown value routed to Other
- **WHEN** `consensus_bucket` receives a value that is not a known category
- **THEN** it returns `Other`
