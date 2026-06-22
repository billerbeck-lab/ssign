## Context

Two scripts produce ssign's run figures. `generate_figures.py` emits 7 "regular" per-run figures (`fig1`–`fig7`) from the integrated substrate CSV; it predates the publication-plots house rules and uses inline `Set2`/`Set3`/`steelblue` palettes, `figN_` filenames, no shared theme, no figure index. `run_enrichment_figure.py` emits the circular-shift null grid; it already defines an inline `THEME` + `rcParams` and is mostly house-styled. The runner wires them as `_step_figures` (gated on an integrated CSV) and `_step_enrichment_figure` (gated on `--enrichment-stats`); `multi_runner` runs both per-genome and additionally pools only the enrichment grid via `pool_and_plot_enrichment` (no pooled regular figures).

The integrated CSV schema (golden fixture `t5ass_minimal_integrated.csv`) carries: `nearby_ss_types`, `tool`, `dlp_extracellular_prob`, `predicted_localization`, `dse_ss_type`, `signalp_prediction`/`signalp_probability`, `aa_length`, `product`, `broad_annotation`, `confidence_tier`, `n_tools_agreeing`, `concordance_ratio`, `annotation_tools`, `sample_id`. It does **not** carry `gravy`/`mw_da`/`isoelectric_point`/`instability_index` (those exist only if ProtParam ran and merged), which is why `fig4_physicochemical` usually renders empty.

Constraint: figures run headless (`matplotlib.use("Agg")`); no new dependencies; figure failures must stay non-core (a plot error must never fail a genome, as `_step_enrichment_figure` already guarantees).

## Goals / Non-Goals

**Goals:**
- One shared style module so every figure (regular + enrichment, single + pooled) draws from the same `THEME`, `rcParams`, palette rules, numbering, and figure index.
- A curated regular-figure set that uses the output's real signal (confidence, tool agreement, localization prob, SignalP-by-type) and drops what ships empty or duplicated.
- Pooled cross-genome regular figures for multi-genome runs, not just pooled enrichment.
- Enrichment grid renders discrete small-count nulls legibly and de-emphasizes non-significant panels, without changing the statistics shown.

**Non-Goals:**
- Changing the enrichment statistics, the substrate-calling logic, or any data schema.
- The HTML report (`generate_report.py`) layout, beyond filenames it references.
- New figure formats (SVG/PDF) or interactivity. PNG at configured DPI stays.
- T3SS-on-by-default (separate change `t3ss-on-by-default-dlp-only`).

## Decisions

**1. Single shared style module `ssign_lib/plot_style.py`.** Exposes `THEME` (semantic keys), an `apply_house_style()` that sets `rcParams` once, palette selectors (`ss_type_palette(types)`, `sequential_cmap`, `diverging_cmap`), and figure-management helpers (`numbered_path(outdir, n, name)`, `clear_unnumbered(outdir)`, `print_figure_index(entries)`). Both `generate_figures.py` and `run_enrichment_figure.py` import it; the inline `THEME` in the enrichment script is deleted. *Alternative — keep per-script themes:* rejected; that is exactly the cross-file-drift this codebase has been bitten by (the two scripts already disagree on palette today).

**2. Curated regular set (per-run), 7 → 7 but re-aimed.** Order = context → main result → drill-down (house rule):
| # | Figure | Key columns | Status |
|---|---|---|---|
| 01 | Substrates per SS type | `nearby_ss_types` | keep, restyle (was fig1) |
| 02 | Evidence strength (confidence tier x tool agreement) | `confidence_tier`, `n_tools_agreeing` | **new** |
| 03 | Localization confidence (DLP extracellular prob, by SS type) | `dlp_extracellular_prob`, `nearby_ss_types` | **new** |
| 04 | SignalP-positive fraction per SS type | `signalp_prediction`, `nearby_ss_types` | **new** (Family A/B biology) |
| 05 | Annotation tool coverage | `*_top1`/`signalp_prediction`/… | keep, restyle (was fig2) |
| 06 | Protein length by SS type | `aa_length`, `nearby_ss_types` | keep, restyle (was fig3) |
| 07 | Functional categories per SS type | `broad_annotation`, `nearby_ss_types` | keep, **merge** old fig5+fig7 |

Dropped from default: `fig4_physicochemical` (columns absent in the standard CSV). It may remain as an opt-in figure guarded on the columns being present, but is off by default. *Alternative — keep all 7 and only restyle:* rejected per the user's "curate" decision; an empty figure and two near-duplicate functional figures are not publication-grade.

**3. SS-type colour is one house qualitative palette, fixed across all figures.** SS types are categorical (up to ~10), so a qualitative palette is allowed by the house rules, but it must be the single house palette applied identically everywhere (so T6SS is the same colour in fig 01, 03, 04, 06). `ss_type_palette()` returns a stable type→colour dict. Sequential data (counts in the pooled heatmap) uses `viridis`; the enrichment fold is not a matrix so no diverging map is needed.

**4. Pooled regular figures via a `--mode pooled` path in `generate_figures.py`.** `generate_figures.py` already accepts `--master-csvs` (nargs +). Add `--mode {per_genome,pooled}` (default `per_genome`). In `pooled` mode it emits cross-genome figures (`P01_substrates_per_genome`, `P02_sstype_by_genome_heatmap`, `P03_pooled_evidence`) into the top-level figures dir. `multi_runner` calls it once with every genome's integrated CSV after the per-genome pass, mirroring how `pool_and_plot_enrichment` already pools the enrichment grid. *Alternative — a separate `pool_figures.py`:* rejected; reuses the same style module and data-loading, and keeping one script avoids a second place to drift. Pooled pass no-ops for <2 genomes.

**5. Enrichment grid: integer-aligned discrete nulls + n.s. de-emphasis.** Null counts are small integers, so the current 40-bin histogram is spiky and reads as noise. Replace with integer-aligned bins (`bins = np.arange(nmin, nmax+2) - 0.5`) so each achievable count is one bar. Panels with `q >= 0.05` get a muted title/observed-line treatment so the eye lands on the significant rows. The observed/null-mean lines, fold, p, and q annotations are unchanged. This is presentational; the enrichment-stats spec's figure requirement still holds, so no spec modification.

**6. Toggle surface follows the curated set.** `PipelineConfig.fig_*` fields and `generate_figures.py` `--no-*` flags are renamed/added to match (e.g. `fig_evidence`, `fig_localization`, `fig_signalp`; drop `fig_*` that no longer map). Pre-release, so renaming is acceptable; `runner._step_figures` arg-wiring updated in lockstep (cross-file-drift guard: grep the flag names).

## Risks / Trade-offs

- **[Figure references a column an older/odd CSV lacks]** → every figure guards its required columns and skips with a logged note (existing pattern in `generate_figures.py`); never raises.
- **[Breaking filename change `figN_`→`0N_`]** → no external users yet; update CLAUDE.md/README and any test/golden asserting the old names. The HTML report does not embed figure filenames by path (verify in task 1); if it does, update together.
- **[Pooled pass on a single genome]** → gated to >1 genome, consistent with `pool_and_plot_enrichment`.
- **[SignalP-by-type figure misread as ground truth]** → caption/title states it is the predictor's call, and Family-B types (T3/T6) being SignalP-negative is expected, not a failure (per `docs/secretion_substrate_localization.md`).
- **[Plot regressions invisible in CI]** → tests assert filenames + non-empty PNGs on fixtures; the author views every PNG before declaring done (house-rule checklist).

## Migration Plan

1. Add `plot_style.py`; restyle `run_enrichment_figure.py` onto it (smallest, already themed) and verify the enrichment PNGs against the current CX3 outputs.
2. Curate + restyle `generate_figures.py` (per-genome mode); numbered names, index, stale cleanup.
3. Add pooled mode + wire `multi_runner`; update `PipelineConfig`/CLI toggles + `_step_figures`.
4. Update docs (CLAUDE.md, README) and tests/golden filenames.
5. Validate: run a small single-genome fixture and a 2-genome run end-to-end; view every PNG.
6. Rollback: revert the change; no data or schema migration is involved.

## Open Questions

- None blocking. Whether to keep `fig4_physicochemical` as an opt-in (column-guarded) figure or delete it outright is a minor call resolved during implementation; default-off either way.
