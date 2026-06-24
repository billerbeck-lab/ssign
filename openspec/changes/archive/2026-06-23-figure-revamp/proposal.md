## Why

ssign's auto-generated run figures don't meet the v1.0.0 publication bar. The 7 regular per-run figures use ad-hoc palettes (`Set2`/`Set3`/`steelblue`), unnumbered `figN_` filenames, no shared theme, and no figure index, none of the publication-plots house rules. One figure (physicochemical violins) reads `gravy`/`mw_da`/`isoelectric_point`/`instability_index`, columns the integrated CSV does not carry unless ProtParam ran and merged, so it usually ships empty. Multi-genome runs emit no pooled cross-genome figures at all (only the pooled enrichment grid). And the set ignores the richest signal in the output: `confidence_tier`, `n_tools_agreeing`, `dlp_extracellular_prob`, and `signalp_prediction` (which ties directly to the Family A/B secretion biology).

## What Changes

- Add a shared plot-style module (`THEME` + `rcParams`, palette-by-data-type rules, numbered filenames, figure index, stale-file cleanup) per the publication-plots house rules, imported by both figure scripts (centralizes the `THEME` that currently lives inline only in `run_enrichment_figure.py`).
- **Curate the regular per-run figure set** (content + style):
  - Drop the usually-empty physicochemical figure from the default.
  - Merge the two near-duplicate functional figures (category bar + functional stacked-bar) into one.
  - Add three biology-grounded figures: substrate evidence strength (`confidence_tier` / `n_tools_agreeing`), localization confidence (`dlp_extracellular_prob` distribution), and SignalP-positive fraction per SS type (the Family A vs B signal: T2/T5 expected +, T3/T6 expected -).
  - Restyle the survivors (SS-type distribution, tool coverage, protein length) onto the shared theme.
- Zero-padded numbered filenames (`01_…`, `02_…`), printed figure index, stale unnumbered `.png` cleanup before each run.
- **New: pooled cross-genome regular figures for multi-genome runs** (substrates-per-genome, SS-type x genome substrate-count heatmap, pooled confidence/evidence), alongside the existing pooled enrichment figure.
- Enrichment figure polish: import the shared theme, render the discrete small-count nulls legibly (integer-aligned bins instead of the current spiky 40-bin histogram), and de-emphasize non-significant panels. Presentational only, the stats shown (observed, null mean, fold, p, q) are unchanged.
- Update the `PipelineConfig` figure toggles and CLI `--no-*` flags to match the curated set.
- **BREAKING (output only)**: figure filenames change (`figN_` -> `0N_`) and the physicochemical figure is no longer emitted by default. No API, data-schema, or dependency change.

## Capabilities

### New Capabilities
- `run-figures`: the set of figures ssign auto-generates per run, single-genome and pooled multi-genome, their content, house-style conventions (theme, numbering, index), and which are emitted by default.

### Modified Capabilities
<!-- none: the enrichment-stats "Per-type null-distribution figure" requirement is unchanged at spec level; the enrichment-figure restyle is presentational. -->

## Impact

- Code: new `src/ssign_app/scripts/ssign_lib/plot_style.py` (shared `THEME`/`rcParams`/helpers); `generate_figures.py` (curate + restyle); `run_enrichment_figure.py` (use shared theme, discrete-null rendering); `core/runner.py` `_step_figures` (toggle/arg wiring); `core/multi_runner.py` (new pooled regular-figure pass); `PipelineConfig` `fig_*` fields + `cli.py` flags.
- Docs: CLAUDE.md figure notes; any README/tutorial figure references.
- Tests: figure scripts run headless (`Agg`) and assert expected output filenames per the curated set; pooled pass asserts cross-genome figures on a 2-genome fixture.
- Dependencies: none added (matplotlib/seaborn already base-tier).
