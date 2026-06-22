## 1. Shared style module

- [x] 1.1 Add `src/ssign_app/scripts/ssign_lib/plot_style.py`: `THEME` (semantic keys), `apply_house_style()` (sets `rcParams` once), palette selectors (`ss_type_palette(types)` returning a stable type->colour dict, `SEQUENTIAL_CMAP`, `DIVERGING_CMAP`), and figure helpers (`numbered_path`/`pooled_path`, `clear_figure_set(outdir)` removing only owned `fig*_`/`0N_`/`P0N_` PNGs so the enrichment figure survives, `print_figure_index(entries)`). Follow `~/.claude/skills/publication-plots/` house rules.
- [x] 1.2 Unit-test `plot_style`: `ss_type_palette` stable/deterministic + variant inheritance; `numbered_path`/`pooled_path` zero-pad; `clear_figure_set` removes owned figs, keeps the enrichment PNG. (7 tests pass)

## 2. Enrichment figure onto the shared style

- [x] 2.1 `run_enrichment_figure.py`: delete the inline `THEME`/`rcParams`, import from `plot_style`; keep tool colours via the shared theme.
- [x] 2.2 Replace the 40-bin histogram with integer-aligned bins (`np.arange(nmin, nmax+2) - 0.5`) so each achievable count is one bar; keep observed/null-mean lines and the fold/p/q annotations.
- [x] 2.3 Visually mute panels with `q >= 0.05` (muted hist + observed line + box) while still showing their statistics; split legend (top-left) from fold box (top-right) to stop the pre-existing collision.
- [x] 2.4 Re-rendered against the CX3 PAO1 reference stats (nulls `.npz` not in the tar, so synthesized Poisson nulls from the reported `null_mean`, seed=0); viewed the PNG: discrete bars + n.s. muting + no legend/box overlap confirmed.

## 3. Curate + restyle the regular per-run figures

- [x] 3.1 `generate_figures.py`: import `plot_style`, call `apply_house_style()`, route all colours through `THEME`/`ss_type_palette`; `clear_figure_set` at start and `print_figure_index` at end.
- [x] 3.2 Restyle and renumber the kept figures: `01_substrates_per_type` (was fig1), `05_tool_coverage` (was fig2), `06_protein_length` (was fig3); merge old fig5+fig7 into `07_functional_categories`.
- [x] 3.3 Add `02_secretion_evidence` (from the `tool` column, the secretion predictors that flagged each protein, NOT `confidence_tier` which is annotation-tool agreement), `03_localization_confidence` (`dlp_extracellular_prob` by SS type, with the call-threshold line), `04_signalp_by_type` (`signalp_prediction` positive-fraction per SS type, captioned as the predictor call with Family-B negative expected); each column-guarded.
- [x] 3.4 Physicochemical is now opt-in (`--physicochemical`), column-guarded, off the default emit path.

## 4. Pooled cross-genome figures (multi-genome)

- [x] 4.1 `generate_figures.py`: added `--mode {per_genome,pooled}` (default `per_genome`). Pooled emits `P01_substrates_per_genome`, `P02_sstype_by_genome` (viridis heatmap), `P03_evidence_basis`; no-ops for <2 genomes (verified). Per-genome + pooled render clean on synthetic fixtures.
- [x] 4.2 `core/multi_runner.py`: after the per-genome pass, when >=2 integrated CSVs exist, call `generate_figures.py --mode pooled` with every genome's CSV, writing into the top-level `figures/` dir (mirrors `pool_and_plot_enrichment`).

## 5. Config, CLI, and runner wiring

- [x] 5.1 `PipelineConfig.fig_*` now the curated 7: `fig_ss_comp/fig_evidence/fig_localization/fig_signalp/fig_tool_heatmap/fig_length/fig_func_summary`; dropped `fig_category`/`fig_substrate_count`; physicochemical opt-in only.
- [x] 5.2 `generate_figures.py` `--no-*` flags + `runner._step_figures` toggle map + `cli.py` flags/dict + `Home.py` checkboxes/kwargs all updated in lockstep. Grep confirms no `fig_category`/`fig_substrate_count` stragglers in `src/`; imports + `--help` parse.
- [x] 5.3 Confirmed `generate_report.py` embeds no `figN_` filenames (no change needed).

## 6. Docs

- [x] 6.1 CLAUDE.md: added a figures Key-Params bullet (curated `01`-`07`, shared `plot_style`, pooled `P01`-`P03`, physicochemical opt-in, fig 02 reads `tool` not `confidence_tier`).
- [x] 6.2 README: replaced "five summary figures" with the curated set + pooled `P01`-`P03`.

## 7. Tests and validation

- [x] 7.1 `tests/unit/test_generate_figures.py`: numbered set emitted + non-empty + no `figN_`; missing column skips only its figure; toggle skips a named figure. Headless `Agg`.
- [x] 7.2 Pooled-mode test: 2-genome fixture emits `P01`-`P03`; 1-genome input emits none.
- [x] 7.3 Updated the golden e2e `_EXPECTED_FIGURES` to the seven new numbered names (all seven render on the golden integrated CSV).
- [x] 7.4 Full unit suite green: 1376 passed (was 1362; +7 plot_style, +6 generate_figures, +1 net).
- [x] 7.5 Validated the figure step on the golden integrated CSV (all seven emit) and realistic synthetic single/2-genome CSVs (seed-pinned); viewed every per-genome (01-07), pooled (P01-P03), and the restyled enrichment PNG against the house-rule checklist. A full MacSyFinder->prediction e2e is the opt-in golden integration test (filenames now updated).
