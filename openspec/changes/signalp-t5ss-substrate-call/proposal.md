## Why

Today **every** detected T5SS component is reported as a substrate, unconditionally: `t5ss_handler` emits all T5 components and `system_filtering` keeps every `T5SS-self` row. The intended per-component DeepLocPro rule (`T5SS_COMPONENT_RULES`, bug-fix #6) only sets an `is_secreted` column in cross_validate that is never read for inclusion, so it is effectively inert. This over-calls T5: components with no localization evidence and no signal peptide are still reported. This change makes the DeepLocPro rule actually gate T5 substrate calling and pairs it with SignalP: a T5 component is a substrate only if DeepLocPro localizes it OR SignalP predicts a Sec signal peptide. SignalP belongs in the rule because T5SS passengers are Sec-exported by definition (Sec-only, confirmed by the T5SS biogenesis literature, no documented Tat-exported T5 passenger), so a clear Sec signal is direct substrate evidence even when DeepLocPro is under threshold (e.g. T5cSS in the CX3 5-genome run: DLP=0 but SignalP-positive). This aligns the substrate caller with the DLP-or-SignalP combined enrichment track.

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

- **Code**: `src/ssign_app/scripts/system_filtering.py` (the real T5 inclusion gate: a new `_t5_self_has_evidence` predicate on `T5SS-self` rows). NOT cross_validate (its `is_secreted` is inert for inclusion). The per-component DeepLocPro rule (bug-fix #6, incl. T5bSS OM-only) is preserved inside the predicate.
- **Behavior**: **fewer** T5 substrates reported (drops components that are DeepLocPro-negative AND have no Sec signal; previously all were kept). The SignalP OR path ensures Sec-signal-bearing T5 components are not lost to the new DeepLocPro requirement. A benchmark before/after regression quantifies the drop and confirms the dropped components are genuinely evidence-less.
- **Docs/tests**: CLAUDE.md critical-bug-fix #6 note; `system_filtering` unit tests (drop-case, DLP path, SignalP path, T5bSS OM-only regression); regenerate the `t5ass_minimal` golden fixture (its weak T5 substrate is now dropped).
- **Out of scope**: the enrichment track (separate `signalp-enrichment-track` change); DeepLocPro's per-component thresholds; T5bSS TpsA passenger calling via proximity (it is not a TXSScan component; see the design open question).
