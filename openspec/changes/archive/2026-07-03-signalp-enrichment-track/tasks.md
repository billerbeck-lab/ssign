## 1. Trigger whole-genome SignalP under enrichment

- [x] 1.1 `runner.py` `Config.__post_init__`: add `self.sp_whole_genome = True` to the existing block that forces `dlp_whole_genome`/`dse_whole_genome` when `enrichment_stats` is on. Update the log line to mention SignalP.
- [x] 1.2 Confirm `_step_signalp` already resolves the full proteome via `_resolve_step_input_fasta(self.config.sp_whole_genome)` and writes a whole-genome `signalp` output the runner can locate (`self.files["signalp"]`).

## 2. SignalP track in the enrichment test

- [x] 2.1 `enrichment_testing.py`: add a `--signalp` arg (whole-genome SignalP TSV); load it keyed by locus_tag.
- [x] 2.2 Add `is_signalp_positive(row)` = `signalp_prediction` maps to a class in `VALID_SEC_SIGNAL_TYPES` (Sec/SPI / Sec/SPII). Reuse the same parse `cross_validate_predictions` uses; keep it in one place.
- [x] 2.3 `positivity_vectors`: build a `signalp` vector over the gene order (missing locus -> negative, like DLP/DSE).
- [x] 2.4 `run_enrichment`: emit a SignalP row for every type (window types -> window mask, autotransporters -> self mask); no skip rule. Add `"SignalP"` to `ENRICH_TOOLS` in `constants.py` so it joins the per-tool BH family via `bh_fdr_by_family`.
- [x] 2.5 Verify `write_nulls` dumps the SignalP null arrays and `pool_enrichment_stats` pools the SignalP `(type, SignalP)` cells unchanged. (Covered by the 2-genome integration test: `T2SS__SignalP` in the npz and pooled `SignalP` observed = 4.)
- [x] 2.6 COMBINED track: make its positivity type-dependent — `signalp` vector when `dt in ENRICH_AUTOTRANSPORTER_TYPES`, else `max(dlp_track, dse_track)` (existing DSE-dropped-for-T3SS rule). The COMBINED row for T5a/c is therefore the SignalP self-detection score; window types are unchanged.

## 3. Figure

- [x] 3.1 `plot_style.py`: add a `THEME["SignalP"]` colour distinct from DLP/DSE/COMBINED.
- [x] 3.2 `run_enrichment_figure.py`: per-tool figure shows a third SignalP bar per type (offset/width logic already generalizes over `len(tools)`); legend gains a SignalP patch. The combined figure renders the per-type COMBINED rows as-is (window types = DLP-or-DSE; T5a/c = SignalP, per task 2.6); the autotransporter (`mode == "self"`) bar is drawn in the SignalP colour with a legend/title note so it reads true.

## 4. Tests + docs

- [x] 4.1 `test_enrichment_testing.py`: SignalP rows emitted for all types incl. T3SS; SignalP-positive uses the Sec-class rule; a Sec-using planted type shows fold > 1; SignalP joins the per-tool BH family. (Also pins the autotransporter combined-uses-SignalP rule via a T5cSS that is SignalP-positive but DLP-self-negative.)
- [x] 4.2 `test_run_enrichment_figure.py`: per-tool figure renders with three tools (DLP/DSE/SignalP) including the inf-fold and self-mode rows.
- [x] 4.3 Integration: a 2-genome enrichment run forces whole-genome SignalP and emits SignalP stats + pooled SignalP rows.
- [x] 4.4 Docs: CLAUDE.md enrichment paragraph (three predictor tracks; SignalP forced whole-genome local), CLI `--enrichment-stats` help. Also fixed the stale CLAUDE.md "remote/DTU default" lines for DLP/SignalP (now: local by default, remote is fallback). README enrichment description is generic ("per-tool and combined") and stays accurate.
- [x] 4.5 Full unit suite green (1410 passed); run `openspec validate signalp-enrichment-track --strict`.

## 5. Validation

- [x] 5.1 Validate on a CX3 single-genome + multi-genome enrichment run: SignalP fold high for Sec-using types (T2SS/T5SS), near/below 1 for Sec-bypassing types (T3SS/T6SS); per-tool figure shows three bars; pooled SignalP renders.
