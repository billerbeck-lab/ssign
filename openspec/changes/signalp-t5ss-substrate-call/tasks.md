## 1. T5SS substrate gate

> Discovery during implementation: the T5 DLP rule in `cross_validate._dlp_flag`
> only sets an `is_secreted` column that is never read for inclusion, and
> `t5ss_handler` emits every detected T5 component, which `system_filtering`
> kept unconditionally (the `or substrate_source == "T5SS-self"` bypass). So the
> real (and only) inclusion gate is in `system_filtering`, and adding the
> DLP-or-SignalP rule there is a *new gate* that drops evidence-less T5
> components (not a rescue of dropped ones).

- [x] 1.1 `system_filtering.py`: add `_t5_self_has_evidence(row, conf)` = DeepLocPro localizes it (T5aSS/T5cSS extracellular-or-OM; T5bSS OM-only, bug-fix #6) OR `is_sec_signal_peptide(signalp_prediction)`. Gate T5SS-self rows on it in the filter loop; drop + log the evidence-less ones.
- [x] 1.2 Scope: only `substrate_source == "T5SS-self"` rows are gated; proximity/non-T5 rows are untouched, so SignalP stays evidence-only outside T5.
- [x] 1.3 `no_sec_signal` quality flag in `t5ss_handler` unchanged (a SignalP-path call is Sec-positive so never `no_sec_signal`; a DLP-path call with a non-Sec class keeps the flag).

## 2. Tests

- [x] 2.1 `test_system_filtering.py`: T5 kept via Sec signal (DLP-negative); T5 kept via DeepLocPro (no signal); T5 with neither is dropped.
- [x] 2.2 Non-T5 rows unaffected (the gate only touches T5SS-self rows; existing proximity tests still green).
- [x] 2.3 T5bSS translocator OM-only rule preserved on the DLP path (bug-fix #6 regression: extracellular-high/OM-low/no-signal T5bSS is dropped).
- [x] 2.4 Full unit suite green (1413). `openspec validate signalp-t5ss-substrate-call --strict`.
- [ ] 2.5 Regenerate the `t5ass_minimal` golden fixture: its T5 substrate (`dlp_extracellular_prob=0.60` < 0.8, `signalp=OTHER`) is now correctly dropped by the gate, so the frozen results change. Regenerate per `tests/fixtures/golden/REGENERATE.md` on a full install (the e2e test is integration-scoped and skips without tools). Alternatively bump the fixture input so the T5 component has evidence, to keep the T5 path exercised.

## 3. Docs

- [x] 3.1 CLAUDE.md critical bug-fix #6 rewritten: the T5 substrate gate (DLP-or-SignalP-Sec), enforced in system_filtering, Sec-only, and the note that it was previously inert.
- [x] 3.2 `docs/reference/output_files.md`: note that a T5 component needs DeepLocPro localization or a Sec signal peptide to be reported. (Already present at the "T5SS substrates" callout, lines 55-59; verified accurate.)

## 4. Validation

- [ ] 4.1 Benchmark regression: run the validation genomes before/after; report how many T5 components get dropped per genome and spot-check they are genuinely evidence-less (DLP < threshold AND no Sec signal), not real substrates. Record in NOTES.md.
- [ ] 4.2 CX3 single + multi-genome run: confirm the T5 evidence-gate log line fires and `combined_results.csv` T5 counts are sane.
