## Why

A T5SS substrate is currently called only when DeepLocPro localizes the component (extracellular/OM per the per-component rule in `T5SS_COMPONENT_RULES`). SignalP is evidence-only: a clear Sec signal peptide cannot, by itself, make a T5 component a substrate. But T5SS passengers are Sec-exported by definition (they carry an N-terminal Sec signal peptide; the Sec-only rule is confirmed by the T5SS biogenesis literature, with no documented Tat-exported T5 passenger). So a genuine T5 substrate whose signal peptide is unambiguous but that DeepLocPro happens to miss is dropped today. The enrichment analysis makes this concrete: SignalP recovers real T5 loci that DeepLocPro scores zero on (e.g. T5cSS in the CX3 5-genome run, where DLP=0 but SignalP marks the component). Adding SignalP as an alternative positive path aligns the substrate caller with T5 biology and with the DLP-or-SignalP combined enrichment track.

## What Changes

- The T5SS per-component substrate rule becomes: a component is a substrate when its system is a T5SS **AND** (**b**) DeepLocPro passes the existing per-component localization rule (`T5SS_COMPONENT_RULES`: T5aSS/T5cSS extracellular-or-OM, T5bSS translocator OM-only) at `conf_threshold` **OR** (**c**) SignalP predicts a Sec signal peptide (`VALID_SEC_SIGNAL_TYPES` = Sec/SPI + Sec/SPII) for the component. Today only (b) exists; this adds (c) as an independent OR path.
- SignalP-positive for the call is **Sec-only** (Tat/PILIN rejected), reusing `is_sec_signal_peptide`. Tat is excluded on biological grounds (Tat exports folded proteins, incompatible with threading the autotransporter β-barrel / TpsB pore).
- The per-component DeepLocPro rules (critical bug-fix #6: T5bSS translocator OM-only, never extracellular) are **preserved** for path (b); SignalP does not relax them, it adds a parallel path.
- Net effect on output: ssign reports **more** T5 substrates (the SignalP-clear, DLP-missed ones); the existing `no_sec_signal` quality flag (non-Sec signal) is unaffected.

## Capabilities

### New Capabilities
- `signalp-t5ss-substrate-call`: SignalP Sec signal peptide as an alternative positive path for calling a T5SS component a secreted substrate, alongside the existing DeepLocPro localization rule.

### Modified Capabilities
<!-- No archived baseline spec captures the T5SS substrate-calling rule as a requirement; it is introduced here as a new capability. The DeepLocPro per-component behavior (bug-fix #6) is preserved, not modified. -->

## Impact

- **Code**: `src/ssign_app/scripts/cross_validate_predictions.py` (the per-component T5 rule check gains a SignalP OR path), possibly `t5ss_handler.py` (keep the Sec-signal quality gate consistent), `ssign_lib/constants.py` (any SignalP threshold constant if a probability floor is chosen).
- **Behavior**: more T5 substrates reported (SignalP-clear components DeepLocPro missed). A regression check on the benchmark genomes should quantify the added calls and confirm they are plausible T5 substrates, not noise.
- **Docs/tests**: CLAUDE.md critical-bug-fix #6 note + Key Parameters; cross_validate unit tests for the new OR path; the T5SS handling tests.
- **Out of scope**: the enrichment track (separate `signalp-enrichment-track` change); DeepLocPro's per-component thresholds; T5bSS TpsA passenger calling via proximity (it is not a TXSScan component; see the design open question).
