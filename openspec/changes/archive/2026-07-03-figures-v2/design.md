## Context

Builds on `figure-revamp` (shared `plot_style`, curated 01–07 set) and the CX3 validation. Verified facts driving the design:
- DLP computes `extracellular_prob`, `outer_membrane_prob`, `periplasmic_prob`, `cytoplasmic_prob`. `cross_validate_predictions` carries `outer_membrane_prob` in its predictions output (FIELDNAMES line 196, populated line 277), but `integrate_annotations` does NOT merge it into the integrated CSV the figures read. `dlp_extracellular_prob` is carried.
- `cog_category` and `kegg_ko` (from `run_eggnog`) DO reach the integrated CSV when EggNOG runs.
- TXSScan models include non-secretion appendages (`T4aP, T4bP, MSH, ComM, Archaeal-T4P`) that currently pass `excluded_systems = [Flagellum, Tad]` and get substrate-called.
- The lab's `figure5_null_distributions` is **group-pooled** (observed ~94, null ~26): big counts → smooth. A single genome's per-type count is a small integer → spiky. This is scope, not a bug.
- `enrichment_testing` currently computes per-(SS-type) nulls only; the genome-wide pooled null (all window components in one mask) is a separate computation (the fleet script `04_circular_shift_enrichment.py` does it).

## Goals / Non-Goals

**Goals:**
- Figures that distinguish cargo-secreting systems (T1/T2/T3/T4/T6) from autotransporters (T5a/c, self-secreting), and use controlled functional vocabularies.
- Enrichment figures that are readable at single-genome scale (fold bars) and reproduce the smooth group-pooled null at multi-genome scale.
- Pooled per-genome figures for a genome group.
- Stop calling non-secretion appendages as substrate-bearing systems.

**Non-Goals:**
- Changing the enrichment statistics themselves (fold/p/q math unchanged; only an additional genome-wide pooled null is computed).
- Changing substrate-calling logic beyond the appendage exclusion default.
- Re-running annotation tools (COG/KEGG already produced by EggNOG).

## Decisions

**1. Non-secretion appendage exclusion (default).** `excluded_systems` default becomes `[Flagellum, Tad, T4aP, T4bP, MSH, ComM, Archaeal-T4P]`. Rationale: T4 pili, MSH pilus, ComM competence machinery, and archaeal pili are surface/uptake appendages, not protein secretion systems; calling neighbourhood proteins their "substrates" is biologically meaningless. T9SS / pT4SSi / pT4SSt are genuine secretion systems and stay. *Alternative — figure-only filter:* rejected; the proteins shouldn't be substrate-called at all, not just hidden. CONFIRM the exact set with Teo (he named T4aP).

**2. Surface `outer_membrane_prob`.** Add it to the columns `integrate_annotations` carries from the predictions file. No new computation. Enables the autotransporter figure to plot extracellular vs OM.

**3. T5a/T5c split.** Autotransporters are removed from `01_substrates_per_type`, `03_localization_confidence`, `04_protein_length` (cargo-only), and get a dedicated `autotransporter self-detection` figure: per T5aSS/T5cSS component, SignalP positive/negative, DLP `extracellular_prob`, DLP `outer_membrane_prob` (the `T5SS_COMPONENT_RULES` view: pass on extracellular OR OM). T5b stays a window type (TpsA cargo is separate). *Alternative — keep in both:* rejected per Teo (remove from the cargo figures).

**4. Functional figures: five sources, two scopes, prune later.** Build candidates and let Teo keep a subset:
- `cog` (COG category: the 1-letter code → 25 standard groups, full names on axis)
- `cog_detailed` (per-COG-category with the specific COG descriptions, or finer grouping — TBD in impl)
- `kegg` (KEGG KO / pathway from `kegg_ko`)
- `eggnog` (`eggnog_description` / `preferred_name`)
- `consensus` (`broad_consensus_annotation` / `broad_annotation`, possibly with revamped keyword matching)
Each in BOTH scopes: overall (one bar chart over all secreted proteins, Xanthobacter `figure9` style) and stacked-by-SS-type. Machinery handling: secreted proteins whose consensus annotation IS a secretion-system/translocator are bucketed/flagged so the chart reflects cargo (decision point in impl: filter out vs labeled "apparatus-associated" bucket). *Note:* this is intentionally over-produced; the keep-set is decided after Teo reviews real output.

**5. Enrichment redesign, split by scope.** *(Redirected 2026-06-24, Teo: collapse to ONE combined DLP+DSE fold/significance bar chart per SS type for both single-genome and pooled scale; the genome-wide pooled-null computation + smooth histogram below were DROPPED. See specs/enrichment-stats.)* A single genome's per-type null is inherently a low-count spike, so:
- Per genome: `run_enrichment_figure` emits a per-type **fold bar chart** (fold = observed / null-mean, BH-q stars, DLP+DSE, autotransporters as self), not the histogram grid. The per-type stats already exist.
- Genome group: `enrichment_testing` additionally computes a **genome-wide pooled null** (positivity near ANY window-type component, one mask, per tool), and `multi_runner` pools these across genomes into one smooth DLP + one DSE histogram (the `figure5` reproduction). *Alternative — keep per-type histograms:* rejected; uninformative at single-genome scale.

**6. Genome-group pooled 01–05.** The per-genome figure functions already accept a combined frame. `multi_runner` calls them once on all genomes' integrated CSVs, writing `01_pooled_*`…`05_pooled_*` into the top-level figures dir, ALONGSIDE the existing P01–P03 (which are genuinely group-specific: substrates-per-genome, SS-type×genome heatmap, pooled evidence). Naming per Teo: `0N_pooled_<name>.png`.

## Risks / Trade-offs

- **[Appendage exclusion changes substrate counts for existing genomes]** → BREAKING default, documented escape hatch (`--excluded-systems`); no external users. Regression test pins the new default set.
- **[Over-many functional figures clutter the output]** → intentional + temporary; `log()` the candidate set; Teo prunes and we delete the rejected ones in a follow-up.
- **[Genome-wide pooled null still smallish for a single genome]** → it's only emitted/meaningful at group level (≥2 genomes); single-genome runs get the fold-bar chart instead.
- **[COG/KEGG absent when EggNOG didn't run]** → column-guarded like every figure; skipped-with-note (base tier has no EggNOG).
- **[Machinery-vs-cargo classification is fuzzy]** → translocators are legitimately secreted; default to a labeled "apparatus-associated" bucket rather than silent drop, so nothing disappears unexplained.

## Migration Plan

1. Data plumbing (`outer_membrane_prob`) + appendage exclusion default + tests.
2. Per-run figure restructure (T5a/c split + autotransporter fig) + functional candidates.
3. Enrichment redesign (fold bars + genome-wide pooled null) + multi_runner pooling.
4. Genome-group pooled 01–05.
5. Docs + full unit suite; validate on a small local/CX3 multi-genome run; Teo reviews, prunes the functional set.
6. Rollback: revert; no data migration.

## Open Questions (resolved 2026-06-24)

- **Appendage set:** RESOLVED, exclude everything that is not a secretion system: `Flagellum, Tad, T4aP, T4bP, MSH, ComM, Archaeal-T4P`. Keep T1–T6SS, T9SS, pT4SSi, pT4SSt. (Verify the exact `ss_type` label strings the pipeline emits during apply, e.g. the heatmap shows `T4aP`, while `validate_macsyfinder_systems` parsing may normalize, match the emitted labels.)
- **COG figures:** RESOLVED, the COG-category figure uses the human-readable COG **category names** (single-letter → name via the standard COG table), not the letters. `cog_detailed` = the finer per-COG-function descriptions where available; if too sparse, fall back to category names (implementer's call).
- **Machinery-associated secreted proteins:** RESOLVED, route to a labeled "apparatus-associated" bucket in the functional figures (not filtered out).
