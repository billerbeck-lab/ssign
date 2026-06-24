## ADDED Requirements

### Requirement: Full-corpus re-verification coverage

The benchmark SHALL independently re-verify every effector in `data/dataset/positives_all.tsv` (all 337 rows, T1SS through T6SS), not only the proximity-reachable subset covered by the prior audit. Each effector MUST receive a recorded verdict: `confirmed`, `fix` (identifier/citation corrected), or `drop` (not supportable to the same-species experimental-evidence bar).

#### Scenario: Every corpus effector has a verdict

- **WHEN** the re-verification ledger is built
- **THEN** the number of ledger rows equals the number of corpus effectors (337), and every row carries a verdict of `confirmed`, `fix`, or `drop` with no blanks

#### Scenario: Coverage is independent, not inherited

- **WHEN** an effector was already passed by the prior audit (tasks 74-84)
- **THEN** it is still re-checked in this pass and its verdict is recorded from this pass, not copied from the earlier audit

### Requirement: Primary-reference DOI verification

For each effector the benchmark SHALL verify that the cited primary-reference DOI resolves to a real publication AND that the publication actually reports the protein as a secreted substrate of the stated secretion-system type. A DOI that does not resolve, or resolves to a paper that does not support the secretion claim, MUST be recorded as a citation defect and the effector routed to `fix` (replacement reference found) or `drop`.

#### Scenario: Citation supports the secretion claim

- **WHEN** an effector's primary-reference DOI is checked
- **THEN** the ledger records whether the reference resolves and whether it supports the secretion claim, with the supporting quote or page captured in the provenance field

#### Scenario: Unsupported citation is not silently kept

- **WHEN** a DOI fails to resolve or does not support the claim and no replacement reference is found
- **THEN** the effector verdict is `drop`, not `confirmed`

### Requirement: Anti-hallucination procedure

The re-verification SHALL use parallel blind agents plus a manual pass, with a strict output schema that forbids inventing accessions, loci, or DOIs (every claimed identifier must be echoed from a retrieved source). Agent claims that cannot be confirmed against UniProt/NCBI/PubMed MUST be treated as unconfirmed and resolved by the manual pass before entering the ledger.

#### Scenario: Unconfirmable agent claim is quarantined

- **WHEN** a blind agent returns an accession or DOI that cannot be confirmed against a primary source
- **THEN** that claim is flagged unconfirmed and excluded from the applied corrections until manually resolved

### Requirement: Reversible fold into the answer key

Drops and fixes from the ledger SHALL be applied through a committed corrections table consumed by `clean_dataset` (mirroring scripts 55/56), so the change is reversible and every benchmark figure regenerates from the same shared filter. The applied corrections MUST NOT be hard-coded into individual figure scripts.

#### Scenario: Figures reflect the re-verified answer key

- **WHEN** the corrections table is applied and the recall + annotation figures are regenerated
- **THEN** all figures count the same re-verified effector set, and reverting the corrections table restores the prior counts
