## Why

A T5aSS/T5cSS autotransporter is a single self-contained polypeptide (its passenger threads through its own outer-membrane beta-barrel), so the component IS the substrate: this is the SELF signal we already test. Separately, secreted-predicted proteins cluster significantly in the +/-3 gene window around T5a/c components ("hitchhikers", hypothesized to piggyback through the T5 pore). Those hitchhikers already exist as proximity substrates but have no enrichment statistic, because the enrichment test hard-branches T5a/c to self-mode only. This change gives T5a/c both tests and splits the two populations in the figures.

## What Changes

- T5aSS and T5cSS each emit **two** enrichment results, distinguished by the existing `mode` column: `self` (unchanged) and a new `window` "hitchhiker" test (+/-3 window, DLP+DSE+SignalP, COMBINED = DLP-or-DSE).
- Enrichment figures render T5a/c as **two adjacent x-groups** per type (`T5aSS (self)` / `T5aSS (hitch)`), in both the per-tool and combined charts.
- Annotation figures split T5a/c substrates into **self vs hitchhiker** categories (`substrate_source`), in the count figure and the physicochemical + functional figures.
- **BREAKING (figure set):** delete figure 02 (autotransporter self-detection scatter) for single and grouped runs; renumber the curated set `03`-`07` -> `02`-`06`; remove the `fig_autotransporter` toggle (CLI/GUI/runner).
- Document the "hitchhiker" concept (definition + pore-piggyback hypothesis), which currently appears nowhere in the repo.

## Capabilities

### New Capabilities
<!-- none: this modifies two existing capabilities -->

### Modified Capabilities
- `enrichment-stats`: for autotransporter types (T5aSS/T5cSS) the test SHALL emit both a self-detection result and a proximity-window (hitchhiker) result, instead of self-detection only. The hitchhiker window uses the same window mask, predictors, and COMBINED = DLP-or-DSE rule as other window types.
- `run-figures`: the curated set SHALL be `01`-`06` (the autotransporter self-detection figure is removed); T5a/c substrates SHALL be shown split by self vs hitchhiker in the count/physicochemical/functional figures; the enrichment figures SHALL show T5a/c as two adjacent self/hitchhiker groups.

## Impact

- Code: `scripts/enrichment_testing.py` (dual self+window emission for T5a/c), `scripts/run_enrichment_figure.py` (two-group rendering), `scripts/generate_figures.py` (self/hitch split, fig-02 removal, renumber), `scripts/ssign_lib/constants.py` (any new type-class set), `cli.py` / `Home.py` / `core/runner.py` (remove `fig_autotransporter` toggle).
- Specs superseded: `figures-v2` (its autotransporter-figure requirement + `01`-`07` set + T5a/c-split design) and `signalp-enrichment-track` (its T5a/c self-only enrichment).
- Docs: `docs/reference/output_files.md`, `docs/explanation/design_decisions.md` (hitchhiker definition), `CLAUDE.md` (enrichment + figures paragraphs), `docs/reference/cli.md` (drop `--fig-autotransporter`).
- Tests: `tests/unit/test_enrichment_testing.py`, `tests/unit/test_run_enrichment_figure.py`, `tests/unit/test_generate_figures.py`.
- No change to substrate-calling: proximity already emits the hitchhiker rows; what counts as a substrate is unchanged. T5bSS, the T5 self evidence gate, and PLM-Effector exclusion are untouched.
