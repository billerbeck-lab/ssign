## Why

ssign's opt-in per-run enrichment test is a binomial that assumes each neighborhood protein is independently secreted at the genome background rate. Secreted genes cluster (operons, genomic islands), so that independence assumption makes the null too narrow and over-states significance. The fleet validation (task #70, `validation_sweeps/benchmark/analysis/fleet_67/04_circular_shift_enrichment.py`) showed the binomial gives the SAME fold as a genome-structure-preserving circular-shift null, but a p-value that is anti-conservative (null SD 14.6 observed vs 10.7 under independence). For a publication-grade per-run significance claim, ssign should report the circular-shift result, not the binomial.

## What Changes

- Replace the binomial enrichment test with a **circular-shift permutation** test as the per-run enrichment analysis. Each run emits, per SS type: an **enrichment score** (fold = observed / null-mean) and a **significance score** (permutation p + BH q over 10,000 shifts). **BREAKING** (output schema): the per-run enrichment table changes columns/semantics.
- Keep the test **opt-in** via the existing `--enrichment-stats` flag. When on, **auto-force whole-genome DLP/DSE** (the circular-shift null needs every gene's positivity in gene order) and surface a CLI note that this increases runtime (~13 min/genome from runtime calibration).
- **Per-type design**: pool a genome's systems of each type. Window types (T1SS, T2SS, T3SS, T4SS, pT4SSt, T5bSS, T6SSi, T6SSii) use ±3-gene windows; autotransporter types (T5aSS, T5cSS) use **self-detection** (the component itself OM≥0.8 or extracellular≥0.8 for DLP, secreted-type for DSE) because the protein is both machinery and substrate. T3SS excluded from DSE; PLM-Effector excluded entirely (unchanged).
- **Auto-generate one figure per enrichment-enabled run**: a grid of per-SS-type circular-shift null histograms (each panel: null distribution, observed line, null mean, fold + p annotation) for DLP and DSE, styled like the prior `figure5_null_distributions.png`.
- Retire the binomial significance path in `enrichment_testing.py`.

## Capabilities

### New Capabilities
<!-- none; this modifies the existing enrichment-stats capability -->

### Modified Capabilities
- `enrichment-stats`: test method changes from binomial-with-sampled-background to circular-shift permutation (fold + permutation p + BH q); `--enrichment-stats` now couples whole-genome DLP/DSE; per-type aggregation with autotransporter self-detection is specified; a per-type null-distribution figure is emitted; the 1000-protein background-sample requirement is superseded by the whole-genome coupling.

## Impact

- **Code**: `src/ssign_app/scripts/enrichment_testing.py` (swap binomial for circular-shift; keep `is_dlp_positive`/`is_dse_positive`/`bh_fdr`), `src/ssign_app/core/runner.py` (`_step_enrichment`, whole-genome coupling, figure step), `src/ssign_app/scripts/ssign_lib/constants.py` (REPS, WINDOW, CONF, window/autotransporter type sets), a new figure generator, CLI help/note.
- **Output**: per-run enrichment TSV schema changes; a new figure artifact per enrichment-enabled run.
- **Runtime**: enrichment-enabled runs get ~13 min/genome slower (whole-genome DLP/DSE), as noted to the user. Default runs (flag off) unchanged.
- **Tests**: new unit tests (rotation_counts vs brute force, autotransporter self-detection, per-type aggregation, BH-skip handling) + a small-genome integration test.
- **Docs**: CLAUDE.md key parameters, README, `docs/design_decisions.md`.
- Porting source: `validation_sweeps/benchmark/analysis/fleet_67/04_circular_shift_enrichment.py` (validated logic).
