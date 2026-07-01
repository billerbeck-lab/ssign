## 1. T5SS substrate rule

- [ ] 1.1 In `cross_validate_predictions.py`, locate the per-component T5 check that applies `T5SS_COMPONENT_RULES` (`max(DeepLocPro probs for the rule columns) >= conf_threshold`). Refactor it to compute `dlp_pass` (the existing rule, unchanged) and `signalp_pass = is_sec_signal_peptide(signalp_prediction)` for the same component, and call it a substrate when `dlp_pass or signalp_pass`.
- [ ] 1.2 Ensure the SignalP path only affects T5SS components (guard on the component being in `T5SS_COMPONENT_RULES` / ss_type being a T5 subtype); SignalP stays evidence-only for every other system.
- [ ] 1.3 Confirm the `no_sec_signal` quality flag in `t5ss_handler` stays consistent: a SignalP-path call is Sec-positive so never `no_sec_signal`; a DLP-path call with a non-Sec class keeps the flag.

## 2. Tests

- [ ] 2.1 `test_cross_validate_predictions.py`: a T5aSS/T5cSS component that is DeepLocPro-negative but SignalP Sec-positive is called a substrate (SignalP rescue); a DeepLocPro-positive SignalP-negative component is still called (DLP path); a Tat-signal DeepLocPro-negative component is NOT called (Sec-only).
- [ ] 2.2 A non-T5 component with a Sec signal but no other positive predictor is NOT called (SignalP stays evidence-only outside T5).
- [ ] 2.3 T5bSS translocator OM-only DeepLocPro rule still holds on the DLP path (bug-fix #6 regression).
- [ ] 2.4 Full unit suite green; `openspec validate signalp-t5ss-substrate-call --strict`.

## 3. Docs

- [ ] 3.1 CLAUDE.md: update critical bug-fix #6 note and the T5SS Key Parameters to describe the DLP-or-SignalP(Sec) T5 substrate rule; note Tat is excluded (Sec-only, literature-backed).
- [ ] 3.2 `docs/reference/output_files.md` / relevant how-to: mention that T5 substrates can be called via a Sec signal peptide alone.

## 4. Validation

- [ ] 4.1 Benchmark regression: run the validation genomes before/after; report the number of added T5 substrate calls per genome and spot-check that they are plausible T5 substrates (Sec-signal-bearing T5 components DeepLocPro missed), not noise. Record in NOTES.md.
- [ ] 4.2 CX3 single + multi-genome run confirming the added T5 calls appear in `combined_results.csv` and the counts are sane.
