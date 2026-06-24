## 1. Data plumbing + substrate scope (foundational)

- [ ] 1.1 `integrate_annotations.py`: carry `outer_membrane_prob` from the predictions file into the integrated CSV (DLP already computes it; cross_validate already emits it). Confirm `dlp_extracellular_prob`/`cog_category`/`kegg_ko` already present.
- [ ] 1.2 Default `excluded_systems` -> exclude everything that is not a secretion system: `[Flagellum, Tad, T4aP, T4bP, MSH, ComM, Archaeal-T4P]` in all entry points (cli, Home x2, runner, system_filtering, validate_macsyfinder_systems, `DEFAULT_EXCLUDED_SYSTEMS`). Verify the labels match the emitted `ss_type` strings. Grep-guard for stragglers.
- [ ] 1.3 Tests: integrated CSV has `outer_membrane_prob`; default `excluded_systems` matches the new set across config + constant; a T4aP-bearing genome yields no T4aP substrates by default.

## 2. Per-run figure restructure

- [ ] 2.1 Remove T5aSS/T5cSS from `01_substrates_per_type`, `03_localization_confidence`, `04_protein_length` (cargo-only). (Label collapse + dropping SignalP/tool-coverage already shipped 5152c0a.)
- [ ] 2.2 New autotransporter self-detection figure: per T5aSS/T5cSS component, SignalP call + DLP `extracellular_prob` + DLP `outer_membrane_prob`. Column-guarded.
- [ ] 2.3 Wire the new figure into `PER_GENOME_FIGS` + config toggle + runner/cli/Home; renumber the per-run set coherently; remove the now-dead `fig_signalp`/`fig_tool_heatmap` toggles.

## 3. Functional-category figures (candidates; Teo prunes)

- [ ] 3.1 COG-category figure: map the single-letter `cog_category` codes to the standard COG **category names** (e.g. M = Cell wall/membrane/envelope biogenesis), bars by name. Overall + per-SS-type scopes.
- [ ] 3.2 COG-detailed figure: finer per-COG-function descriptions where available, else fall back to category names. Both scopes.
- [ ] 3.3 KEGG figure (`kegg_ko` -> KO / pathway) in both scopes.
- [ ] 3.4 EggNOG-annotation figure (`eggnog_description`/`preferred_name`) in both scopes.
- [ ] 3.5 Broad-consensus figure (`broad_consensus_annotation`), with a keyword-matching revamp toward a curated vocabulary (Transporter/Hydrolase/Protease/...); machinery-associated secreted proteins (consensus = secretion-system/translocator) routed to a labeled "apparatus-associated" bucket. Both scopes.
- [ ] 3.6 `log()` the full candidate set so Teo can review and name the keep-set (rejected ones deleted in follow-up).

## 4. Enrichment redesign

- [ ] 4.1 `enrichment_testing.py`: additionally compute a genome-wide pooled null per predictor (positives near ANY window-type component, single mask, all rotations); emit to the stats/npz output.
- [ ] 4.2 `run_enrichment_figure.py`: per-genome figure becomes a per-type fold-enrichment bar chart (fold + BH-q stars, DLP+DSE, autotransporters self), replacing the histogram grid.
- [ ] 4.3 `multi_runner.py` / `pool_and_plot_enrichment`: pool the per-genome genome-wide nulls across genomes into one smooth DLP + one DSE histogram (the `figure5` reproduction) at the group level.

## 5. Genome-group pooled per-genome figures

- [ ] 5.1 `multi_runner.py`: after the per-genome pass, run the per-genome figure functions once over all genomes' integrated CSVs, writing `0N_pooled_<name>.png` into the top-level figures dir, alongside the existing P01-P03.
- [ ] 5.2 Test: a 2-genome fixture emits `0N_pooled_*` for each curated per-genome figure; a 1-genome run emits none.

## 6. Docs, tests, validation

- [ ] 6.1 CLAUDE.md (excluded_systems default, figure set, autotransporter/functional figures), README.
- [ ] 6.2 Update `generate_figures`/golden tests for the new per-run set + autotransporter + functional figures; full unit suite green.
- [ ] 6.3 Validate on a small local + CX3 multi-genome run; verify the genome-wide pooled null is smooth and per-genome enrichment matches the single-genome null_mean (the earlier multi-genome background fix). Teo reviews the functional candidates and names the keep-set.
