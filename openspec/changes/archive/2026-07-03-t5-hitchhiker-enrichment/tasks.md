## 1. Enrichment: T5a/c dual self+hitchhiker

- [x] 1.1 In `enrichment_testing.py:run_enrichment`, for T5aSS/T5cSS emit BOTH a `mode=self` row-set (unchanged: `dlp_self` + DSE + SignalP, COMBINED=DLP-or-SignalP) AND a `mode=window` hitchhiker row-set (window_mask + `dlp` + DSE + SignalP, COMBINED=DLP-or-DSE). Key the emitted rows by `(ss_type, tool, mode)` so nulls don't collide.
- [x] 1.2 Confirm `enrich_null_key` (and the pooled `pool_enrichment_stats` accumulation) disambiguates by mode so a T5a/c self and window null are pooled separately, not merged. Adjust the key if it only uses `(ss_type, tool)`.
- [x] 1.3 Keep T5bSS window-only, T3SS DSE-dropped, and all other types unchanged. Verify PLM-Effector stays excluded.
- [x] 1.4 Unit tests (`test_enrichment_testing.py`): a T5aSS fixture emits a `mode=self` and a `mode=window` row per predictor; hitchhiker COMBINED == DLP-or-DSE union; self COMBINED == DLP-or-SignalP union; window observed uses the +/-3 mask.

## 2. Enrichment figure: two adjacent T5a/c groups

- [x] 2.1 In `run_enrichment_figure.py`, group the x-axis by `(display_type, mode)` so each T5a/c type renders two adjacent groups labelled `T5aSS (self)` / `T5aSS (hitch)` (and T5cSS); non-T5a/c types stay single-group. Apply to both the per-tool and combined charts, pooled and per-genome.
- [x] 2.2 Style: self group keeps the SignalP-colour DLP-or-SignalP COMBINED; hitch group uses the DLP-or-DSE COMBINED styling. Preserve the faded-not-grey n.s. treatment and the "significance" titles/axes.
- [x] 2.3 Unit test (`test_run_enrichment_figure.py`): a stats TSV with T5aSS self+window rows renders without error and produces two T5aSS x-groups.

## 3. Annotation figures: split self vs hitchhiker

- [x] 3.1 Add a helper mapping a T5a/c substrate row to a split display label (`T5aSS (self)` when `substrate_source=="T5SS-self"`, else `T5aSS (hitch)`); other types pass through `display_type`. Make palette lookups strip the ` (self)`/` (hitch)` suffix so both share the T5a/c colour.
- [x] 3.2 Apply the split in `fig01_secreted_by_genome` and the physicochemical + functional figures so T5a/c self and hitchhiker appear as distinct categories.
- [x] 3.3 Unit test: an integrated CSV with T5aSS self and proximity rows yields split categories in the count and functional figures.

## 4. Delete fig-02 + renumber curated set

- [x] 4.1 Remove `fig02_autotransporter`, its `PER_GENOME_FIGS` entry, the `--no-autotransporter` arg, and (if now unused) the `CONF_THRESHOLD`/`AUTOTRANSPORTER_TYPES`/`ENRICH_AUTOTRANSPORTER_TYPES` imports and `_signalp_positive`/`_SIGNALP_NEGATIVE` helpers.
- [x] 4.2 Renumber physicochemical/COG/KEGG/EggNOG/consensus `03`-`07` -> `02`-`06`: `numbered_path` indices, `_functional_sources` n-values (4-7 -> 3-6), the module docstring figure index (drop `02 autotransporter`, renumber, `01`-`07` -> `01`-`06`), and the figure-index print.
- [x] 4.3 Remove the `fig_autotransporter` toggle across `cli.py` (arg + mapping), `Home.py` (checkbox + kwarg), `core/runner.py` (Config field + `--no-autotransporter` mapping).
- [x] 4.4 Update `test_generate_figures.py`: drop `02_autotransporter_self_detection.png` expectations, assert the renumbered `02_physicochemical`..`06_consensus`, and the self/hitch split.
- [x] 4.5 Grep the repo for stale references (`fig_autotransporter`, `02_autotransporter`, `01`-`07`, `--fig-autotransporter`) and fix drift.

## 5. Docs + concept

- [x] 5.1 Add a short "hitchhikers" definition (what: secreted-predicted neighbours of T5a/c; why: pore-piggyback hypothesis; how tested: proximity-window enrichment) to `docs/explanation/design_decisions.md` and reference it from `docs/reference/output_files.md`.
- [x] 5.2 Update `output_files.md`: enrichment table now has self+window rows for T5a/c; curated set is `01`-`06`, no fig-02.
- [x] 5.3 Update `CLAUDE.md`: enrichment paragraph (T5a/c dual self+hitchhiker, window COMBINED=DLP-or-DSE); figures paragraph (drop `02`, renumber `01`-`06`, note the self/hitch split); drop `--fig-autotransporter` from `docs/reference/cli.md`.

## 6. Validate + archive ordering

- [x] 6.1 Run `pytest tests/unit/ -q` (enrichment, enrichment-figure, generate-figures suites green) and re-run the simplify pass on the diff.
- [x] 6.2 Regenerate figures from an existing 5-genome integrated CSV and eyeball: two T5a/c enrichment groups, split annotation categories, no `02_autotransporter`, contiguous `01`-`06`.
- [x] 6.3 Note archival order: this change's `run-figures` delta supersedes figures-v2's autotransporter-figure requirement, so archive `figures-v2` (and `signalp-enrichment-track`) before archiving this change.
