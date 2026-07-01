## ADDED Requirements

### Requirement: SignalP is an alternative positive path for T5SS substrate calls
When identifying secreted substrates, a component whose secretion system is a T5SS SHALL be called a substrate when EITHER its DeepLocPro per-component localization rule passes OR SignalP predicts a Sec signal peptide for it. The DeepLocPro per-component rules (`T5SS_COMPONENT_RULES`, including the T5bSS translocator outer-membrane-only rule) SHALL be preserved unchanged; the SignalP path SHALL be evaluated independently and combined by logical OR. This SignalP-as-decider behaviour SHALL apply ONLY to T5SS components; for all other secretion systems SignalP SHALL remain evidence-only.

#### Scenario: SignalP rescues a DeepLocPro-negative T5 component
- **WHEN** a T5SS component fails its DeepLocPro localization rule but SignalP predicts a Sec signal peptide for it
- **THEN** the component SHALL be called a T5SS substrate

#### Scenario: DeepLocPro path still calls a SignalP-negative T5 component
- **WHEN** a T5SS component passes its DeepLocPro per-component rule but SignalP does not predict a signal peptide for it
- **THEN** the component SHALL still be called a T5SS substrate (the DeepLocPro path is unchanged)

#### Scenario: Non-T5 systems keep SignalP evidence-only
- **WHEN** a component belongs to a non-T5SS system and SignalP predicts a signal peptide but the system's own predictors do not call it secreted
- **THEN** SignalP alone SHALL NOT make it a substrate (SignalP stays evidence-only outside T5SS)

### Requirement: T5SS SignalP path accepts only Sec signal peptides
The SignalP positive path for a T5SS substrate call SHALL count only Sec signal peptides (`VALID_SEC_SIGNAL_TYPES`: Sec/SPI and Sec/SPII). Tat (Tat/SPI, Tat/SPII) and PILIN predictions SHALL NOT satisfy the SignalP path, because T5SS passenger export is Sec-dependent.

#### Scenario: Tat signal does not satisfy the SignalP path
- **WHEN** a T5SS component is DeepLocPro-negative and SignalP predicts a Tat (or PILIN) signal, not a Sec signal
- **THEN** the SignalP path SHALL NOT fire, and the component SHALL NOT be called a substrate on the SignalP path alone
