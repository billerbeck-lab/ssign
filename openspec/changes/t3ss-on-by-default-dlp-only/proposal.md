## Why

T3SS is currently excluded from ssign by default (`excluded_systems = [Flagellum, Tad, T3SS]`). The exclusion was a blunt fix for a real problem: DeepSecE over-predicts T3SS (1,808 calls across 74 genomes, mostly flagellar misclassification, vs 0 MacSyFinder-validated T3SS in that set). But MacSyFinder + TXSScan detect *bona fide* injectisome T3SS reliably; it is only DeepSecE that cannot tell injectisome from flagellar. Excluding the whole system throws away real T3SS detection to avoid one untrustworthy predictor. The targeted fix is to turn T3SS detection on and stop trusting DeepSecE for T3SS specifically (task #68).

## What Changes

- **BREAKING (default behavior):** drop `T3SS` from the default `excluded_systems`, so the default becomes `[Flagellum, Tad]`. T3SS systems are now detected, neighbourhood-extracted, and substrate-called by default.
- DeepSecE **never** contributes to a T3SS call. In `cross_validate_predictions.py` the existing conditional guard (`dse_T3SS_flagged` only when the genome lacks a validated T3SS) becomes **unconditional**: any DeepSecE `T3SS` prediction is flagged out of the secretion-evidence count regardless of whether the genome has a validated T3SS. This is what prevents the flagellar false positives from returning once T3SS detection is on.
- DLP (DeepLocPro extracellular) remains the trusted PLM localization signal for T3SS; SignalP is unaffected (T3SS effectors lack Sec signal peptides, so SignalP-negative is the correct expectation); PLM-Effector is unaffected (off by default).
- The per-run enrichment test surfaces T3SS as **DLP-window-only** (already pre-wired via `ENRICH_WINDOW_TYPES ∋ T3SS` and `ENRICH_DSE_NO_WINDOW = {T3SS}`); it was simply never reached because T3SS was excluded upstream. No T3SS DSE enrichment row is emitted.

## Capabilities

### New Capabilities
- `t3ss-detection`: T3SS is detected and substrate-called by default, with DeepSecE excluded as a T3SS evidence source across cross-validation and enrichment (DLP-only for T3SS).

### Modified Capabilities
<!-- none: enrichment-stats already specifies DLP/DSE handling and is mid-revision under enrichment-circular-shift-per-run; T3SS-specific behaviour is owned by the new t3ss-detection capability. -->

## Impact

- Code: `src/ssign_app/core/runner.py` (default `excluded_systems`), `src/ssign_app/cli.py` (flag default + help), `src/ssign_app/Home.py` (two GUI defaults), `src/ssign_app/scripts/system_filtering.py` (CLI default), `src/ssign_app/scripts/cross_validate_predictions.py` (`_dse_flag` unconditional T3SS flag; `_genome_has_t3ss` may become unused for flagging).
- Docs: `CLAUDE.md` (excluded_systems default + Key Parameters; reframe critical-bug-fix #4 from "T3SS excluded" to "T3SS on, DSE excluded for T3SS"), `README`, `docs/design_decisions.md`.
- Tests: cross-validate T3SS-DSE-always-flagged (genomes with and without a validated T3SS); new default `excluded_systems`; enrichment emits a T3SS DLP row and no T3SS DSE row for a T3SS-bearing genome.
- Validation: re-run Salmonella LT2 (SPI-1/SPI-2) on CX3; confirm a sane T3SS substrate count (not the old 1,808-style flagellar blowup) and a T3SS DLP enrichment row.
- No new dependencies. Users who relied on the old default can restore it with `--excluded-systems Flagellum,Tad,T3SS`.
