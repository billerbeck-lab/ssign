## Context

`enrichment_testing.py:run_enrichment` branches each SS type into one of two modes: window types get a +/-W mask over neighbours; autotransporter types (T5aSS/T5cSS) get a self-mask over the component positions (`is_auto -> mode="self"`). T5a/c therefore produce **only** a self result. But `proximity_analysis.py` builds a +/-W window around every component including T5a/c, so secreted-predicted neighbours of T5a/c are already emitted as `substrate_source=="proximity"` substrates ("hitchhikers"). The self and hitchhiker populations are biologically distinct (the autotransporter itself vs its neighbours) and currently only the first is tested and the second is silently merged into T5a/c in the figures.

The positivity machinery already computes every vector this needs: `positivity_vectors` returns `dlp` (extracellular), `dlp_self` (OM-or-ext), `dse`, and `signalp`. The window path (mask + DLP/DSE/SignalP + COMBINED=DLP-or-DSE) is exactly what T5bSS already runs.

## Goals / Non-Goals

**Goals:**
- T5aSS and T5cSS each emit both a `self` and a `window` (hitchhiker) enrichment result, reusing existing masks/vectors.
- Figures separate the two populations: enrichment charts show two adjacent groups per T5a/c type; annotation figures split T5a/c substrates by self vs hitchhiker.
- Retire the now-redundant autotransporter self-detection scatter (fig-02) and renumber the curated set to `01`-`06`.
- Name and document the "hitchhiker" concept.

**Non-Goals:**
- No change to substrate-calling (proximity already emits hitchhiker rows).
- No change to T5bSS, T3SS, the T5 self evidence gate, or PLM-Effector exclusion.
- Not persisting gate-dropped components (out of scope; hitchhikers are a separate population from dropped self-components).

## Decisions

**1. T5a/c become dual-mode; keep `ss_type`, distinguish by `mode`.** `run_enrichment` emits, for T5aSS/T5cSS, both a self row-set and a window row-set (same `ss_type`, `mode` in {self, window}). *Alternative — a distinct `T5aSS-hitch` pseudo-type:* rejected; it fragments the palette/display order and the `mode` column already carries the distinction. Downstream (figures, pooling via `enrich_null_key`) keys on `(ss_type, tool, mode)` so the two nulls don't collide.

**2. Hitchhiker window follows the non-T5 window convention.** Predictors DLP(extracellular)+DSE+SignalP; COMBINED = DLP-or-DSE. Rationale: a hitchhiker is an ordinary neighbour substrate, not the Sec-dependent autotransporter, so it uses the same rule as T1/T2/T6 windows. The `self` result keeps DLP-self(OM-or-ext)+SignalP and COMBINED = DLP-or-SignalP (unchanged from `signalp-t5ss-substrate-call`). So per T5a/c type: self-COMBINED is Sec-flavoured, hitch-COMBINED is DSE-flavoured.

**3. BH families unchanged in definition, larger in membership.** The per-tool family and the COMBINED family gain the new T5a/c window rows. This is correct (more genuine tests) but shifts q-values slightly across the family; noted as a trade-off, not a bug.

**4. Enrichment figure: group by `(display_type, mode)` for T5a/c.** Two adjacent x-groups labelled `T5aSS (self)` / `T5aSS (hitch)`. Self group keeps the SignalP-colour COMBINED styling; hitch group uses the DLP-or-DSE COMBINED styling. Non-T5a/c types stay single-group. Figure width already scales with group count.

**5. Annotation split by `substrate_source`.** A small helper maps a T5a/c substrate row to a split display label (`T5aSS (self)` when `substrate_source=="T5SS-self"`, else `T5aSS (hitch)`); window/other types pass through `display_type` unchanged. Palette lookups strip the ` (self)`/` (hitch)` suffix so both share the T5a/c colour. Applied in fig 01 and figs 03-06 (physicochemical + functional).

**6. fig-02 deletion + renumber.** Remove `fig02_autotransporter`, its `PER_GENOME_FIGS` entry, the `--no-autotransporter` arg, and the `fig_autotransporter` toggle across CLI/GUI/runner. Renumber physicochemical/COG/KEGG/EggNOG/consensus `03`-`07` -> `02`-`06` (numbered_path indices + `_functional_sources` n-values + docstring index + figure-index print). Drop imports/helpers that only fig-02 used (`CONF_THRESHOLD`, `AUTOTRANSPORTER_TYPES`, `_signalp_positive`, `_SIGNALP_NEGATIVE`) after confirming no other user.

## Risks / Trade-offs

- **BH q-values shift for the whole family** (more T5a/c rows) → expected and correct; document in the enrichment note.
- **Figure crowding** from up to 4 extra T5a/c groups → mitigated by the existing width-scaling on group count; verify legibility on the 5-genome pooled figure.
- **Palette KeyError on split labels** → mitigated by suffix-stripping in the palette lookup; unit-test a T5a/c split render.
- **Superseding an active spec (`figures-v2`, 22/24)** → update the `run-figures` main spec via this change's delta and record the supersede; figures-v2 archival reconciles against the updated main spec.
- **`_signalp_positive` (denylist) removal** also retires the figure-02-only looser SignalP rule that a prior simplify flagged as inconsistent with `is_sec_signal_peptide`; no remaining user, so the inconsistency is resolved by deletion.

## Open Questions

None outstanding; all four scoping questions resolved with Teo (hitch predictors = DLP/DSE/SignalP; two adjacent figure groups; annotation split; standalone change that documents the concept).
