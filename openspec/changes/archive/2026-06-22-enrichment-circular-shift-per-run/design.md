## Context

The per-run enrichment test (`enrichment_testing.py`, opt-in via `--enrichment-stats`) is a one-sided binomial: per system it counts DLP/DSE-positive proteins in the ±3 neighborhood and compares to a background rate estimated from a sampled null (`sample_null_proteins.py`, default 1000). Task #70's fleet validation (`validation_sweeps/benchmark/analysis/fleet_67/04_circular_shift_enrichment.py`) showed the binomial fold is correct but its p-value is anti-conservative: it assumes neighborhood proteins are independently secreted, while secreted genes cluster (operons, islands), so the real null is wider (SD 14.6 vs 10.7). The circular-shift permutation preserves that clustering and gives an honest p-value. This change makes the circular-shift the per-run test.

The validated logic already exists in script 04 (FFT `rotation_counts`, unit-tested against brute force). This design ports it into the pipeline.

## Goals / Non-Goals

**Goals:**
- Per enrichment-enabled run, emit per SS type: fold (enrichment) + permutation p + BH q, from a circular-shift null.
- Auto-generate one per-SS-type null-distribution figure (DLP + DSE), figure5-style.
- Keep it opt-in; couple whole-genome DLP/DSE behind the same flag with a runtime note.
- Handle autotransporters (T5aSS/T5cSS) by self-detection, not neighborhood windows.

**Non-Goals:**
- Not flipping enrichment on by default (stays opt-in; the runtime cost is real).
- No PLM-Effector enrichment (unchanged exclusion).
- No per-system-instance rows in the primary output (per-type pooling only; a single genome usually has ≤1 system per type anyway).
- No pooling across genomes (that was the fleet-validation mode; a per-run test is one genome).

## Decisions

**1. Exact all-rotations FFT null, not Monte-Carlo sampling.** For a single genome the exact permutation null is the full set of circular rotations of the predictor's positivity vector. `rotation_counts(pos_vec, win_mask)` via FFT cross-correlation returns the count at every offset in O(n log n); `c[0]` is the observed, `c[1:]` is the exact null. So we report an exact permutation p = (#{c[1:] ≥ observed} + 1)/(n) and fold = observed / mean(c[1:]). The 10,000-shift Monte-Carlo (what the prior figure used) is unnecessary; we keep a configurable cap (`PERM_MAX_EXACT`, e.g. 200k) above which we subsample offsets for very large genomes. *Alternative considered:* 10,000 random shifts like the reference script, rejected as strictly noisier than exact for the same cost.

**2. One FFT path for both window types and autotransporters.** The window mask differs by type class:
- Window types (T1SS, T2SS, T3SS, T4SS, pT4SSt, T5bSS, T6SSi, T6SSii): `win_mask` marks positions within ±W of any component of that type (minus the component positions, matching `get_neighborhood_proteins`); `pos_vec` is `is_dlp_positive` (extracellular ≥ conf) / `is_dse_positive`.
- Autotransporter types (T5aSS, T5cSS): `win_mask` marks the component positions *themselves* (width 1, no window); `pos_vec` for DLP is self-positivity (outer-membrane ≥ conf OR extracellular ≥ conf), for DSE is the component's secreted type-call. Then "observed = how many of this type's components are self-secreted" and the rotation null answers "how often would that many randomly-placed marks hit self-positive genes" — i.e. detection rate vs the genome's OM/extracellular background.

This unifies both as `rotation_counts(pos_vec, win_mask)`; only the two vectors differ. *Alternative considered:* a separate hypergeometric for autotransporters, rejected because the rotation null gives the same answer with one code path and the same figure.

**3. The rotation null replaces the sampled-background step.** `mean(c[1:])` equals rate × (mask size) exactly, which is the background expectation the binomial estimated by sampling. So `_step_sample_null_proteins` / `sample_null_proteins.py` are not needed on the circular-shift path and are dropped from the enrichment flow. Fewer moving parts, exact background.

**4. Whole-genome coupling lives in config normalization.** When `enrichment_stats` is set, the runner sets `dlp_whole_genome = dse_whole_genome = True` (already-supported flags) and logs a one-line runtime note. The circular-shift needs every gene's positivity in gene order; a neighborhood-only run physically cannot build `pos_vec`. Coupling it to the flag prevents a silently-wrong null.

**5. Keep `enrichment_testing.py` as the stats core; add a figure wrapper.** `enrichment_testing.py` is refactored to the circular-shift method (keep `is_dlp_positive`, `is_dse_positive`, `bh_fdr`, `load_systems`, `broad_type`; drop `binom_pvalue`, `score_scope`, `--null-ids`). It writes the stats TSV plus a compact `*_enrichment_nulls.npz` (per type×tool: null array or histogram bins + observed + fold + p). A new `run_enrichment_figure.py` renders the per-type grid from that dump, wired as a figure step in the runner (matching the `run_*.py` + figure-step convention). *Alternative:* render inline in `enrichment_testing.py`, rejected to keep stats and plotting separable (and testable without matplotlib).

**6. Output schema (BREAKING, opt-in pre-1.0).** New per-type TSV columns: `sample_id, ss_type, tool, mode (window|self), observed, n_mask, null_mean, fold, p_perm, qvalue, significant, n_rotations`. Drops binomial-only fields (`p_bg`, `M`/`k` renamed to `n_mask`/`observed`, `n_null`).

**7. Multi-contig.** Order genes contig-then-start, concatenate into one circular sequence (matching scripts 02/04). Rotations can cross contig boundaries; acceptable for genomes dominated by a few large replicons, documented.

## Risks / Trade-offs

- [Single-genome power for rare types] A type with one small system has a coarse null; fold is reported but p can be weak. → Report honestly; the figure shows the null spread so weak signals are visually obvious. BH q across the few per-type tests is expected to be low-power, not a bug.
- [Autotransporter counts are small] Few T5a/T5c components → wide null, rarely significant per single genome. → Self-detection is reported as a detection rate with its null; the fleet aggregate (benchmark) is where significance accrues. Documented.
- [Multi-contig wrap] Rotations cross replicon boundaries. → Order contig-then-position; large contigs dominate; documented as an approximation (same one the validated script used).
- [FFT float rounding] → `np.rint(...).astype(int)`; unit test asserts FFT counts equal brute-force counts (ported from script 04).
- [Runtime] Whole-genome DLP/DSE adds ~13 min/genome. → Opt-in only; explicit CLI note; default runs untouched.
- [BREAKING output schema] → Opt-in feature, pre-1.0, no external consumers; README/docs updated.

## Migration Plan

1. Land the refactor behind the existing `--enrichment-stats` flag (default off): no change to default runs.
2. Couple whole-genome predictors + note in config normalization.
3. Replace binomial internals; drop null-sampling step from the enrichment flow.
4. Add figure wrapper + runner figure step.
5. Update tests + docs. Rollback = revert the change; the flag-off path is unchanged throughout.

## Open Questions

- Histogram detail in the `.npz`: full null array (simplest, a few MB) vs pre-binned counts (smaller). Lean full-array unless size is a problem. (Resolve in implementation.)
- Whether to also emit per-system-instance rows as a secondary table for debugging. Default no; add only if the per-type view proves too coarse in practice.
