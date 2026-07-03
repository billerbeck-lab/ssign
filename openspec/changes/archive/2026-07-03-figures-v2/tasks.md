## 1. Data plumbing + substrate scope (foundational)

- [x] 1.1 Carry `outer_membrane_prob` into the integrated CSV. Done at the substrate producers (`proximity_analysis.py`, `t5ss_handler.py`) rather than `integrate_annotations` directly: producers write it onto each substrate row, `system_filtering` unions row keys, integrate uses `substrates_filtered` as its base frame (columns preserved), and runner `priority_cols` orders it. Verified end-to-end by code inspection + new test. `dlp_extracellular_prob`/`cog_category`/`kegg_ko` confirmed present.
- [x] 1.2 Default `excluded_systems` = `[Flagellum, Tad, T4aP, T4bP, MSH, ComM, Archaeal-T4P]` (Teo approved 2026-06-30). Wired via the `DEFAULT_EXCLUDED_SYSTEMS` constant across cli/Home(x2)/runner/system_filtering/validate_macsyfinder_systems (no more hardcoded `["Flagellum","Tad"]`). Verified labels match the emitted `model_fqn` `ss_type` strings (`T4aP` etc.); validate regex extended (specific-before-generic). Home options list gained the appendages so the default is a valid multiselect subset.
- [x] 1.3 Tests: integrated CSV has `outer_membrane_prob` (`test_integrate_annotations.py::test_dlp_localization_probs_preserved`); default set matches the constant (`test_runner.py`); T4aP/T4bP/MSH/ComM/Archaeal-T4P substrates dropped by default (`test_system_filtering.py::test_appendage_substrate_dropped_under_default_exclusions`).

## 2. Per-run figure restructure

- [x] 2.1 Remove T5aSS/T5cSS from `01_substrates_per_type`, `03_localization_confidence`, `04_protein_length` (cargo-only). (Label collapse + dropping SignalP/tool-coverage already shipped 5152c0a.)
- [x] 2.2 New autotransporter self-detection figure: per T5aSS/T5cSS component, SignalP call + DLP `extracellular_prob` + DLP `outer_membrane_prob`. Column-guarded.
- [x] 2.3 Wire the new figure into `PER_GENOME_FIGS` + config toggle + runner/cli/Home; renumber the per-run set coherently; remove the now-dead `fig_signalp`/`fig_tool_heatmap` toggles.

## 3. Functional-category figures (candidates; Teo prunes)

- [x] 3.1 COG-category figure: map the single-letter `cog_category` codes to the standard COG **category names** (e.g. M = Cell wall/membrane/envelope biogenesis), bars by name. Overall + per-SS-type scopes. (`functional_vocab.COG_CATEGORY_NAMES`; multi-letter codes like `MU` split per letter.)
- [x] 3.2 COG-detailed figure: folded into the EggNOG-description figure (no finer per-COG column exists in the integrated CSV; design allows the fallback). Logged note.
- [x] 3.3 KEGG figure (`kegg_ko` -> KO) in both scopes. Pathway names unavailable offline; plots raw KO IDs (logged note).
- [x] 3.4 EggNOG-annotation figure (`eggnog_description`) in both scopes.
- [x] 3.5 Broad-consensus figure (`broad_annotation`), curated via `consensus_bucket` (reuses `annotation_consensus.CATEGORY_PATTERNS`); machinery (Secretion system / Flagellar) routed to an "Apparatus-associated" bucket, fallback-noise to "Other". Both scopes.
- [x] 3.6 `log()` the full candidate set (label/file/desc + folding notes) so Teo can review and name the keep-set.

## 4. Enrichment redesign

- [~] 4.1 DROPPED (Teo redirect 2026-06-24): no genome-wide pooled null computed. The combined fold/significance bar chart supersedes the smooth genome-wide-null histogram; `enrichment_testing` keeps its per-(SS-type) circular-shift stats unchanged.
- [x] 4.2 `run_enrichment_figure.py`: per-genome figure becomes a per-type fold-enrichment bar chart (fold + BH-q stars, DLP+DSE, autotransporters self), replacing the histogram grid. (Redirected per Teo: ONE combined DLP+DSE fold/significance bar chart per SS type; replaces both the grid and the planned genome-wide pooled null figure. See task #6 spec update.)
- [x] 4.3 `pool_and_plot_enrichment` pools the per-genome enrichment stats into one combined `pooled_enrichment_fold.png` (same combined fold/significance bar chart at group scale). Genome-wide-null histogram dropped with 4.1.

## 5. Genome-group pooled per-genome figures

- [x] 5.1 `0N_pooled_<name>.png` for the genome group: implemented by extending `generate_figures.py` pooled mode to also run the curated per-genome set over the combined frame (renamed `0N_pooled_*`) alongside P01-P03. One `clear_figure_set` / one subprocess, so `multi_runner` needs no change (it already calls pooled mode).
- [x] 5.2 Tests: `test_pooled_mode_two_genomes` asserts the `0N_pooled_*` set + P01-P03; `test_pooled_mode_single_genome_noop` asserts a 1-genome run emits none.

## 6. Docs, tests, validation

- [x] 6.1 CLAUDE.md (excluded_systems default, 01-07 figure set, combined enrichment, KEGG bundle, readable CSV columns) + README updated for the round-2/3 final state.
- [x] 6.2 Tests updated for the new per-run set + autotransporter + functional figures + pooled-per-genome (`test_generate_figures.py`, new `test_functional_vocab.py`, `test_integrate_annotations.py`). Full unit suite green (1392).
- [x] 6.3 Validate on a small local + CX3 multi-genome run. Teo reviews the functional candidates (06-13) and names the keep-set; rejected ones deleted in follow-up. (Genome-wide-null check dropped with 4.1.)

## 7. Round 2/3 redesign (Teo feedback, 2026-06-30)

The figure set was iterated twice after first review. Final state (supersedes the
01-13 candidate set in sections 2-3 above):

- [x] 7.1 Curated per-genome set renumbered to a clean **01-07**: 01 secreted-by-SS-type
  (single genome = one bar per type histogram; group = one stacked bar per genome,
  merging old 01 + P01), 02 autotransporter self-detection, 03 size & physicochemical
  (length folded in as the first panel, + GRAVY/MW/pI/instability/aromaticity/charge),
  04-07 functional by SS type (COG / KEGG / EggNOG / consensus). Dropped: the
  secretion-evidence, localization, P02 heatmap, P03 evidence, and all "overall"
  functional figures. No more `P0N_` figures (01 is the cross-genome view).
- [x] 7.2 Functional "Other" reduced: COG + consensus (bounded vocabularies) show every
  category; KEGG + EggNOG cap at top-15. KEGG labels are real function names via a
  bundled `data/kegg_ko_names.tsv.gz` (KEGG list/ko), loaded lazily in `functional_vocab`.
- [x] 7.3 Integrated CSV carries readable `cog_category_name` + `kegg_function` per protein
  (`integrate_annotations._add_functional_names`).
- [x] 7.4 Enrichment: KEPT the per-tool DLP+DSE figure and ADDED a separate **combined
  DLP-or-DSE** figure (`enrichment_testing` emits a `COMBINED` track via the same
  circular-shift test, BH'd as its own family; `run_enrichment_figure --combined`;
  runner + pooling emit `*_enrichment_fold_combined.png`).
- [x] 7.5 DEFERRED (Teo): revisit the consensus-function keyword grouping
  (`annotation_consensus.CATEGORY_PATTERNS`) vocabulary/granularity. Recorded in NOTES.
