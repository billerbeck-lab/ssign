## Context

T5SS substrate calling lives in `cross_validate_predictions.py`: for a component whose `(ss_type, gene_name)` is in `T5SS_COMPONENT_RULES`, the trigger is `max(DeepLocPro probabilities for the rule's columns) >= conf_threshold`. The rules encode bug-fix #6: `T5aSS_PF03797` and `T5cSS_PF03895` pass on `extracellular_prob` OR `outer_membrane_prob`; `T5bSS_translocator` (TpsB pore) is `outer_membrane_prob` only (never extracellular). SignalP today is evidence-only everywhere (`_signalp_supports` feeds `signalp_supports_secretion`, which does not contribute to `is_secreted`); `t5ss_handler` additionally flags a `no_sec_signal` quality flag when a T5 substrate's SignalP class is non-Sec. The T5bSS TpsA passenger is not a TXSScan component and is handled by `proximity_analysis`, not this rule.

The `signalp-enrichment-track` change established the DLP-or-SignalP union as the right pairing for T5 (Sec-dependent, DSE unreliable) and a literature review confirmed T5 export is Sec-only with no credible Tat exception.

## Goals / Non-Goals

**Goals:**
- Call a T5SS component a substrate when EITHER DeepLocPro passes its per-component localization rule OR SignalP predicts a Sec signal peptide for it.
- Reuse the single Sec-signal definition (`is_sec_signal_peptide`, `VALID_SEC_SIGNAL_TYPES`).
- Preserve the per-component DeepLocPro rules (bug-fix #6) unchanged for the DLP path.

**Non-Goals:**
- Making SignalP a global decider for `is_secreted` (this change is scoped to the T5 per-component rule only; SignalP stays evidence-only for all other systems).
- Changing DeepLocPro thresholds or the `no_sec_signal` quality flag.
- Calling T5bSS TpsA passengers (proximity path, not a component rule; see open question).

## Decisions

**1. SignalP-positive for the T5 call = a Sec signal peptide (`is_sec_signal_peptide`, Sec/SPI + Sec/SPII).** Tat (TAT/TATLIPO) and PILIN are rejected. Rationale: the T5SS biogenesis literature is unanimous that passengers cross the inner membrane via Sec (they must stay unfolded to thread the β-barrel/TpsB pore; Tat exports folded proteins). Admitting Tat would add false positives (Tat redox enzymes near a locus) for ~zero true-positive gain. Reuses the shared constant already used by the enrichment track and `t5ss_handler`.

**2. SignalP is an independent OR path; DeepLocPro rules are preserved.** The call becomes `dlp_rule_passes OR signalp_is_sec`. The per-component DeepLocPro columns (bug-fix #6, incl. T5bSS OM-only) are untouched; SignalP does not relax them, it runs in parallel. A component positive by DLP alone, SignalP alone, or both is a substrate.

**3. SignalP decision = the class call, no extra probability floor (default).** SignalP's own class assignment (SP/LIPO) is the decision, consistent with how the enrichment track and `t5ss_handler` already treat it. *Alternative (add a `signalp_probability >= threshold` floor):* deferred; SignalP's class call already encodes its confidence, and adding a second knob risks drift from the enrichment definition. Left as a tunable if false positives appear on the benchmark.

**4. Scope = the T5 TXSScan components in `T5SS_COMPONENT_RULES` (T5aSS/T5cSS self-detection, T5bSS translocator).** These are where "the component is (or gates) the substrate". The T5bSS TpsA passenger is a separate, proximity-detected protein, not a component rule; extending the SignalP-OR to proximity-called TpsA is a possible follow-up but out of scope here (open question below).

## Risks / Trade-offs

- **[More substrates reported]** → the SignalP-OR path adds T5 components DeepLocPro missed. Mitigation: a benchmark regression (before/after counts on the validation genomes) to confirm the added calls are plausible T5 substrates, not noise, before archiving.
- **[SignalP false positives on non-substrates near a T5 locus]** → bounded by scope: the rule only fires for a protein that is itself a T5 TXSScan component (already validated as T5 machinery by MacSyFinder), so a stray Sec signal on an unrelated neighbour cannot trigger it.
- **[Interaction with `no_sec_signal` flag]** → a component called via the DLP path but with a non-Sec SignalP class still gets the `no_sec_signal` quality flag (unchanged); a component called via the SignalP path is by definition Sec-signal-positive, so it is never `no_sec_signal`. Consistent, no special-casing needed.

## Open Questions

- **T5bSS TpsA passenger.** Should the SignalP-OR also apply to proximity-called TpsA passengers (Sec-dependent, signal-peptide-bearing), or stay limited to the TXSScan components? Leaning: separate follow-up, since TpsA calling lives in `proximity_analysis`, not the component rule.
