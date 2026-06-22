## 1. Constants + core circular-shift module

- [x] 1.1 Add to `ssign_lib/constants.py`: `ENRICH_PERM_REPS`/exact-rotation cap, `ENRICH_WINDOW=3`, `ENRICH_CONF=0.8`, `ENRICH_WINDOW_TYPES` (T1SS,T2SS,T3SS,T4SS,pT4SSt,T5bSS,T6SSi,T6SSii), `ENRICH_AUTOTRANSPORTER_TYPES` (T5aSS,T5cSS), `ENRICH_DSE_NO_WINDOW` (T3SS).
- [x] 1.2 Port `rotation_counts(pos_vec, win_mask)` (FFT cross-correlation, `np.rint().astype(int)`, `c[0]`==observed) from `fleet_67/04_circular_shift_enrichment.py` into `enrichment_testing.py`.
- [x] 1.3 Add `ordered_positivity_vectors(gene_order, dlp, dse, conf)`: build per-genome circular gene-order, return DLP/DSE positivity vectors over ALL proteins (window-type positivity) and the locus→ordinal map; multi-contig ordered contig-then-start.
- [x] 1.4 Add self-positivity for autotransporters: DLP = (outer_membrane_prob ≥ conf OR dlp_extracellular_prob ≥ conf), DSE = `is_dse_positive`.

## 2. Per-type test + output

- [x] 2.1 Replace `score_scope`/binomial `main` with a per-type circular-shift test: for each SS type present, build `win_mask` (±W neighborhood minus components for window types; component positions for autotransporters) and run `rotation_counts` for DLP and DSE (skip DSE for T3SS).
- [x] 2.2 Compute observed=`c[0]`, null=`c[1:]` (subsample if above the exact cap), `null_mean`, `fold=observed/null_mean`, `p_perm=(#{null≥observed}+1)/(len+1)`; BH-correct across the non-skipped rows with existing `bh_fdr`.
- [x] 2.3 New OUT_FIELDS: `sample_id, ss_type, tool, mode, observed, n_mask, null_mean, fold, p_perm, qvalue, significant, n_rotations`; write stats TSV.
- [x] 2.4 Dump per-(type,tool) null arrays + observed/fold/p to `*_enrichment_nulls.npz` for the figure.
- [x] 2.5 Update the CLI: drop `--null-ids`; require whole-genome `--dlp`/`--dse`; keep `--ss-components`, `--gene-order`, `--window`, `--conf-threshold`, `--sample`, `--out`; add `--nulls-out`.
- [x] 2.6 Remove dead binomial code (`binom_pvalue`, sampled-background logic) and update the module docstring.

## 3. Figure generator

- [x] 3.1 Add `scripts/run_enrichment_figure.py`: read `*_enrichment_nulls.npz`, render a grid (one panel per tested SS type × predictor) with null histogram, observed line, null mean, fold + p annotation; figure5 styling; obey `~/.claude/skills/publication-plots` defaults.
- [x] 3.2 Save as a numbered figure into the run's figures dir; print it in the run's figure index.

## 4. Runner wiring

- [x] 4.1 In config normalization, when `enrichment_stats` is set, force `dlp_whole_genome=dse_whole_genome=True` and emit a one-line runtime-increase note (cite ~13 min/genome).
- [x] 4.2 Update `_step_enrichment` to call the new `enrichment_testing.py` CLI (whole-genome DLP/DSE inputs, `--nulls-out`); remove the `_step_sample_null_proteins` dependency from the enrichment flow.
- [x] 4.3 Add the enrichment figure step (calls `run_enrichment_figure.py`) gated on `enrichment_stats`, placed in the figures phase.
- [x] 4.4 Verify resume/checkpoint + core-vs-optional failure handling still hold for the changed steps.

## 5. Tests

- [x] 5.1 Unit: `rotation_counts` equals brute-force rotation counts on random vectors/masks (ported assertion).
- [x] 5.2 Unit: window-type observed equals direct ±W neighborhood count; `c[0]` == observed.
- [x] 5.3 Unit: autotransporter self-detection observed equals count of self-positive components; DSE-T3SS row absent.
- [x] 5.4 Unit: per-type aggregation pools multiple systems of a type; BH skips no real rows and flags significance correctly.
- [x] 5.5 Integration: small genome end-to-end with `--enrichment-stats` produces the stats TSV + `.npz` + figure; assert columns and that whole-genome DLP/DSE were forced.

## 6. Docs

- [x] 6.1 Update CLAUDE.md Key Parameters (enrichment now circular-shift; `--enrichment-stats` forces whole-genome DLP/DSE + runtime note; n_null_proteins retired for this path).
- [x] 6.2 Update README enrichment section + `docs/design_decisions.md` (why circular-shift over binomial; autotransporter self-detection).
- [x] 6.3 NOTES.md: mark the #69 circular-shift decision shipped; carry the PAO1-duplicate dedup reminder for any validation re-run.

## 7. Combined cross-genome (pooled) figure

- [x] 7.1 `pool_enrichment_stats` takes `nulls_output` and dumps the pooled Monte-Carlo null arrays (same `enrich_null_key` scheme) alongside the pooled stats TSV.
- [x] 7.2 Multi-genome call site (Home.py) renders `pooled_enrichment_null_distributions.png` from the pooled TSV + nulls via `run_enrichment_figure.py`.
- [x] 7.3 Tests: unit (pooled nulls npz written with correct keys/shape) + integration (pooled figure renders end-to-end).
- [x] 7.4 Docs: output_files.md notes the pooled TSV/npz/figure.
