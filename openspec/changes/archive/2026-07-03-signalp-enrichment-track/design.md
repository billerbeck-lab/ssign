## Context

The circular-shift enrichment test (`enrichment_testing.py`) builds, per genome, whole-genome gene-ordered positivity vectors for DeepLocPro and DeepSecE, then for each SS type counts positives within a mask (±N-gene window for window types; component positions for autotransporter self-detection) and compares the observed count to the exact set of circular rotations. `--enrichment-stats` already forces whole-genome DLP/DSE (`Config.__post_init__`); the per-tool figure shows a DLP and a DSE bar per type, BH-corrected as one family. A separate `COMBINED` (DLP-or-DSE) track exists, BH'd as its own family.

SignalP already runs in the pipeline (`run_signalp.py`, output `signalp_prediction` ∈ {`OTHER`, `SP(Sec/SPI)`, `LIPO(Sec/SPII)`, `TAT(Tat/SPI)`, `TATLIPO(Tat/SPII)`, `PILIN(Sec/SPIII)`}), with a whole-genome mode (`Config.sp_whole_genome` + `_step_signalp` resolving the full proteome) that enrichment does not yet trigger. `enrichment_testing.py` takes `--dlp`/`--dse` but no `--signalp`. SignalP runs **local** (no webserver).

## Goals / Non-Goals

**Goals:**
- A SignalP fold/p/q per SS type, computed by the same circular-shift test, for all types.
- SignalP forced whole-genome under `--enrichment-stats`, reusing existing machinery; local only.
- A third SignalP bar in the per-tool enrichment figure.
- Reuse the pipeline's existing definition of a secretion-supporting signal peptide.

**Non-Goals:**
- Changing the DLP/DSE statistics or the figures-v2 figure set.
- A new SignalP install/runner (the whole-genome step already exists).
- Webserver execution (removed in a prior version).

## Decisions

**1. SignalP-positive = a Sec signal peptide (`VALID_SEC_SIGNAL_TYPES`).** A protein is positive when `signalp_prediction` is Sec/SPI or Sec/SPII (the `("SP", "LIPO")` set in `constants.py` that `cross_validate_predictions` already uses for `signalp_supports_secretion`). Tat/PILIN are excluded: Tat exports folded proteins by a different route, PILIN marks pili (a default-excluded appendage). Reusing the existing constant keeps one definition of "secretion signal" across the pipeline. *Alternative — any signal-peptide class (what the autotransporter figure's marker uses):* rejected for the statistic; "any class" is fine for a visual marker but the enrichment should test the Sec-secretion signal specifically.

**2. SignalP is valid for every type (no skip rule).** Unlike DSE (skipped for T3SS), SignalP is emitted for all types. T3SS/T6SS effectors bypass Sec, so their SignalP fold should sit near or below 1 — that low value IS the informative result (a clean negative contrast), not a reason to skip.

**3. SignalP joins the per-tool BH family.** Added to `ENRICH_TOOLS = ("DLP", "DSE", "SignalP")`, so it is BH-corrected together with DLP/DSE and shown as a third per-type bar. The existing per-tool-vs-combined family split (`bh_fdr_by_family`) keeps working unchanged.

**4. COMBINED bar = DLP-or-DSE for non-T5 types, DLP-or-SignalP for ALL T5SS (T5a/b/c).** Two figures: a **per-tool** figure (a separate DLP, DSE, and SignalP bar per type, each its own circular-shift score) and a **combined** figure (one bar per type). The combined bar is a per-type *union* of the two relevant predictors: a gene counts if either flags it, then the same circular-shift test is run fresh on that combined positivity (it is NOT an average of the two folds, and NOT "take the taller bar"). For non-T5 types the combined stays `np.maximum(dlp, dse)` (with DSE dropped for T3SS → DLP only). For **all** T5SS subtypes (T5aSS, T5bSS, T5cSS via `ENRICH_T5SS_TYPES`) the combined is `np.maximum(dlp_track, signalp)`; DSE is unreliable for T5 and dropped.

Rationale: T5 substrates are Sec-dependent (signal-peptide-bearing; the Sec-only rule is confirmed by the T5SS biogenesis literature, with no documented Tat-exported T5 passenger), so the two informative predictors for a T5 locus are localization (DLP extracellular/OM) and the Sec signal peptide (SignalP), not DSE. Pairing DLP with SignalP mirrors the substrate-caller intent. The union does raise the genome background vs a single track (lowering the fold a little), exactly as DLP-or-DSE does for the other systems, which is the accepted behaviour of the combined bar. *Implementation:* the combined positivity vector is type-dependent, `np.maximum(dlp_track, signalp)` when `dt in ENRICH_T5SS_TYPES`, `dlp_track` alone for T3SS (DSE dropped), else `np.maximum(dlp_track, dse)`. (Superseded the earlier "SignalP-alone for T5a/c" proposal after CX3 validation: DLP self-detection is strongly significant for T5aSS, so unioning DLP with SignalP is more faithful than SignalP alone.)

**5. Plumbing reuses what exists.** `Config.__post_init__` adds `sp_whole_genome = True` to the existing enrichment force block. The runner passes the whole-genome SignalP TSV to `enrichment_testing.py` via a new `--signalp` arg. `positivity_vectors` gains a `signalp` vector; `run_enrichment` emits the SignalP track with the same masks; pooling (`pool_enrichment_stats`) already iterates all `(type, tool)` keys, so SignalP pools for free.

## Risks / Trade-offs

- **[Extra whole-genome SignalP pass per genome]** → adds runtime under `--enrichment-stats`; accepted by Teo, local execution, no rate limits.
- **[SignalP baseline rate is high]** → ~10–20% of a proteome is Sec-positive, so the null mean is non-trivial; folds for Sec-using systems are still strongly > 1 (near-100% vs background), and the test stays exact. Not a correctness issue, just calibrates expectations.
- **[`ENRICH_TOOLS` now has three entries]** → the figure widens to 3 bars/type; the existing offset/width logic already generalizes over `len(tools)`. Theme needs a SignalP colour.
- **[Whole-genome SignalP must cover every gene for the rotation null]** → if SignalP output is missing a locus, treat it as negative (matches DLP/DSE handling of missing rows); a fully-missing SignalP TSV yields an all-negative vector (fold 0), surfaced rather than silently wrong.
