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
- [ ] 4.2 `core/multi_runner.py`: after the per-genome pass, when >=2 genomes, call `generate_figures.py --mode pooled` with every genome's integrated CSV, writing into the top-level figures dir (mirror `pool_and_plot_enrichment`).

## 5. Config, CLI, and runner wiring

- [ ] 5.1 Update `PipelineConfig.fig_*` fields to the curated set (add `fig_evidence`, `fig_localization`, `fig_signalp`; drop fields with no remaining figure; default physicochemical off).
- [ ] 5.2 Update `generate_figures.py` `--no-*` flags and `runner._step_figures` arg-wiring in lockstep; grep the old flag names to confirm no stragglers (cross-file-drift guard).
- [ ] 5.3 Confirm `generate_report.py` does not embed figure paths by the old `figN_` names; if it does, update the references.

## 6. Docs

- [ ] 6.1 CLAUDE.md: update the figures description (curated set, numbered filenames, pooled multi-genome figures, physicochemical off by default).
- [ ] 6.2 README / any tutorial: update figure references and example filenames.

## 7. Tests and validation

- [ ] 7.1 Unit tests for `generate_figures.py`: on a small fixture integrated CSV, assert the expected numbered filenames exist, are non-empty, and that a missing column skips only its figure (run headless via `Agg`).
- [ ] 7.2 Unit test for pooled mode: on a 2-genome fixture, assert `P0*` files appear; on a 1-genome input, assert none.
- [ ] 7.3 Update any golden/fixture or test asserting old `figN_` filenames.
- [ ] 7.4 Run the full unit suite (`pytest tests/unit/ -q`); all green.
- [ ] 7.5 Validation: run a small single-genome fixture and a 2-genome run end-to-end; view every emitted PNG against the house-rule checklist (no overlap/clipping, theme colours, populated panels).
