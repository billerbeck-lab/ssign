## 1. Shared style module

- [ ] 1.1 Add `src/ssign_app/scripts/ssign_lib/plot_style.py`: `THEME` (semantic keys), `apply_house_style()` (sets `rcParams` once), palette selectors (`ss_type_palette(types)` returning a stable type->colour dict, `SEQUENTIAL_CMAP`, `DIVERGING_CMAP`), and figure helpers (`numbered_path(outdir, n, name)`, `clear_unnumbered(outdir)`, `print_figure_index(entries)`). Follow `~/.claude/skills/publication-plots/` house rules.
- [ ] 1.2 Unit-test `plot_style`: `ss_type_palette` is stable/deterministic for a given type list; `numbered_path` zero-pads; `clear_unnumbered` removes only non-numbered `*.png`.

## 2. Enrichment figure onto the shared style

- [ ] 2.1 `run_enrichment_figure.py`: delete the inline `THEME`/`rcParams`, import from `plot_style`; keep tool colours via the shared theme.
- [ ] 2.2 Replace the 40-bin histogram with integer-aligned bins (`np.arange(nmin, nmax+2) - 0.5`) so each achievable count is one bar; keep observed/null-mean lines and the fold/p/q annotations.
- [ ] 2.3 Visually mute panels with `q >= 0.05` (muted title + observed line) while still showing their statistics.
- [ ] 2.4 Re-render against the CX3 reference stats (`validation_sweeps/cx3_enrichment/...enrichment_stats.tsv` + `..._nulls.npz`); view the PNG and confirm discrete bars + n.s. muting read correctly.

## 3. Curate + restyle the regular per-run figures

- [ ] 3.1 `generate_figures.py`: import `plot_style`, call `apply_house_style()`, route all colours through `THEME`/`ss_type_palette`; add `clear_unnumbered` at start and `print_figure_index` at end.
- [ ] 3.2 Restyle and renumber the kept figures: `01_substrates_per_ss_type` (was fig1), `05_tool_coverage` (was fig2), `06_protein_length_by_type` (was fig3); merge old fig5+fig7 into `07_functional_categories_by_type`.
- [ ] 3.3 Add `02_evidence_strength` (`confidence_tier` x `n_tools_agreeing`), `03_localization_confidence` (`dlp_extracellular_prob` by SS type), `04_signalp_positive_by_type` (`signalp_prediction` fraction per SS type); each column-guarded and skipped-with-note if absent. Title `04` to state it is the predictor's SignalP call (Family-B SignalP-negative is expected).
- [ ] 3.4 Move physicochemical to an opt-in, column-guarded figure that is off by default; remove it from the default emit path.

## 4. Pooled cross-genome figures (multi-genome)

- [ ] 4.1 `generate_figures.py`: add `--mode {per_genome,pooled}` (default `per_genome`). In `pooled` mode emit `P01_substrates_per_genome`, `P02_sstype_by_genome_heatmap` (sequential cmap), `P03_pooled_evidence` into the given `--outdir`; no-op when fewer than 2 genomes are present in the combined CSVs.
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
