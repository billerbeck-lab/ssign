## Why

The CX3 1/4/20-genome validation runs (2026-06-23) surfaced the figure set's real-world gaps, reviewed against the lab's existing Xanthobacter figures. T3SS detection and the `figure-revamp` styling work, but: (1) the per-SS-type enrichment histograms are spiky and uninformative at single-genome scale (a count of 0–4 over ~5000 shifts piles on the left), where the readable presentation is a per-type fold-bar chart plus a *group-pooled* genome-wide null; (2) T5aSS/T5cSS are autotransporters where the component *is* the substrate, so showing them in the generic cargo figures (substrate count, localization, length) conflates two different quantities; (3) the functional-category figure uses raw `broad_annotation` (too granular, and "secretion system" appears because some secreted proteins genuinely are apparatus), where the target is curated controlled vocabularies (COG, KEGG); (4) multi-genome runs emit only the cross-genome P01–P03 figures, not pooled versions of the per-genome 01–05; and (5) Type IV pili and other non-secretion appendages (T4aP/T4bP/MSH/ComM/Archaeal-T4P) are being detected and called as substrate-bearing "systems," which is biologically wrong (they are not protein secretion systems).

## What Changes

**Data plumbing (foundational):**
- Surface `outer_membrane_prob` in the integrated CSV (DLP already computes it; it's in the predictions file but dropped at integrate). `dlp_extracellular_prob`, `cog_category`, `kegg_ko` already reach the integrated CSV.

**Substrate scope (BREAKING default):**
- Exclude non-secretion appendages by default: add `T4aP, T4bP, MSH, ComM, Archaeal-T4P` to `excluded_systems` (joining `Flagellum, Tad`). These are pili / competence / motility machinery, not protein secretion systems. `T9SS`, `pT4SSi`, `pT4SSt` stay (they are secretion). Opt back in via `--excluded-systems`.

**Per-run figures (restructure):**
- (Already shipped this session) drop SignalP-by-type + tool-coverage; collapse SS variant labels via `display_type`.
- Remove T5aSS/T5cSS from the generic cargo figures (substrates-per-type, localization, length), leaving those cargo-only.
- New autotransporter self-detection figure for T5aSS/T5cSS: per component, SignalP call + DLP extracellular prob + DLP outer-membrane prob.
- **Functional figures** (candidates; user prunes after review), each in BOTH scopes (overall + stacked-by-SS-type): COG category, COG detailed, KEGG, EggNOG annotation, broad-consensus annotation (with possible keyword-matching revamp).

**Enrichment figures (redesign):**
- Per genome: replace the spiky per-type null-histogram grid with a per-type **fold-enrichment bar chart** (fold + BH-q stars), like the lab's `figure4`.
- Genome group (multi-genome): add a **genome-wide pooled null** figure (all window-type components pooled across all genomes, one DLP + one DSE smooth histogram), like the lab's `figure5`. Requires `enrichment_testing.py` to emit a genome-wide pooled null.

**Genome-group pooled figures:**
- Multi-genome runs also emit pooled versions of the per-genome set as `01_pooled_*` … `05_pooled_*` (all genomes combined), alongside the existing P01–P03.

## Capabilities

### New Capabilities
- `secretion-system-scope`: which TXSScan systems ssign treats as substrate-bearing secretion systems by default (excludes non-secretion appendages).

### Modified Capabilities
- `run-figures`: the per-run + pooled figure set (T5a/T5c split, autotransporter figure, functional figures, genome-group pooled 01–05, OM-prob/COG/KEGG data use).
- `enrichment-stats`: per-genome figure becomes per-type fold-bars; a genome-wide pooled null is emitted at group level.

## Impact

- Code: `integrate_annotations.py` (carry `outer_membrane_prob`); `excluded_systems` defaults (cli/Home/runner/system_filtering/constants); `generate_figures.py` (T5a/c split, autotransporter fig, functional figs, pooled 01–05); `run_enrichment_figure.py` (fold-bars + group null); `enrichment_testing.py` (genome-wide pooled null); `multi_runner.py` (pooled 01–05 wiring); config toggles.
- Docs: CLAUDE.md (excluded_systems, figure set), README.
- Tests: figure emission per the new set; appendage exclusion default; OM-prob in integrated CSV; genome-wide pooled null.
- No new dependencies (COG/KEGG already from EggNOG; matplotlib/seaborn base-tier).
