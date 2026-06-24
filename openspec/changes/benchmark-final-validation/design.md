## Context

The benchmark (change `effector-recovery-benchmark`) is feature-complete: recall figure 06 reports 45/60 reachable system instances, annotation figures 07/08 grade emitted secreted proteins, and the enrichment machinery (`enrichment-circular-shift-per-run`) is in product. Two gaps block a publication-grade claim:

1. **Answer-key confidence.** The anti-hallucination audit (tasks 74-84, then a re-audit 81-84) only covered the proximity-reachable proteins that reach the recall figure (~245). The full corpus is 337 effectors. Unreachable and non-emitted effectors still define the benchmark denominators (the reachable/unreachable split and the per-type totals), so a fabricated "unreachable" effector silently distorts the result just as much as a fabricated reachable one.
2. **Measurement provenance.** Every number was read off one extended-tier fleet (`/tmp/ssign_fleet_67`) with three defects: only 5 of 20 T5SS genomes were annotated, the four RTX T1SS toxins were rescued by a manual `clean_dataset` staging flip rather than a real run, and `--enrichment-stats` was never run fleet-wide.

Constraints: CX3 is user-driven (password + 2FA, no agent SSH); genome inputs (`inputs_gb/`, 537 MB) are gitignored; the concurrent running-job cap is ~2; extended-tier databases on `$EPHEMERAL` are purged after ~30 days idle.

## Goals / Non-Goals

**Goals:**
- Independent verdict (`confirmed`/`fix`/`drop`) for all 337 corpus effectors and every primary-reference DOI, folded reversibly into `clean_dataset`.
- One clean extended-tier CX3 run over a corrected 66-genome panel with enrichment, T3SS-on, and full T5SS coverage.
- Retire the manual T1SS staging fix by confirming the RTX toxins emit from real full-assembly runs.
- Close benchmark tasks 70 (enrichment working) and 71 (annotation-tool benchmark).

**Non-Goals:**
- No ssign product/runtime code changes (analysis + validation only).
- Not the full database tier (that is task 72); this matches the fleet's extended tier for comparability.
- Not a re-derivation of the machinery answer key (locus coordinates stay; only the effector corpus is re-verified).

## Decisions

**D1. Re-verify the full 337-row corpus, not just the figure proteins.** The recall denominators depend on every effector's reachability classification, so coverage must be 100% or the "could-never-reach" bars remain unverified. Alternative (re-check only the ~245 figure proteins again) was rejected: it re-audits what is already audited and leaves the denominator unproven.

**D2. Genome panel = 51 + 15 - 16 = 66, emitted as a derived manifest.** Include the 51 input genomes with a testable proximity system and the 15 T5SS genomes missing from the fleet; drop the 16 inputs with no testable answer-key system. The 16 only ever served as enrichment background, and the circular-shift null is **per-genome** (each predictor's own gene-ordered positivity vector, rotated), so removing them does not weaken any enrichment test. Rerunning all 67 was rejected as wasted GPU on genomes with no recall/annotation signal.

**D3. Match the fleet's extended tier; drive via `cx3-submit`.** Keeps the rerun directly comparable to the current figures. Use `submit_batched_overnight.sh --small --enrichment-stats` so the user pastes one short line per batch (raw `qsub` gets mangled on paste). `--small` (32c/64gb) places faster than the 64c default.

**D4. Batch the 66 genomes; per-genome enrichment is the primary deliverable.** At ~30-48 min/genome and a ~2-job cap, 66 genomes is multi-day. Submit in overnight batches. The benchmark consumes per-genome enrichment + per-genome recall/annotation; pooled enrichment is per-batch and treated as a bonus, not a headline number.

**D5. Retire the staging fix with real emission, keep the flag as fallback.** Once the RTX toxins emit from full-assembly runs, recall no longer needs `APPLY_T1SS_STAGING_FIX`; leave the flag (default False) so the pre-rerun figures remain reproducible.

## Risks / Trade-offs

- **The 15 missing T5SS genomes may not be staged in `inputs_gb/`.** → Verify first (early task); if absent, fetch via `09_fetch_refseq.py` before submitting. This is the top open question.
- **`$EPHEMERAL` database purge kills the run at step 2.** → Run `ssign doctor --tier extended` + `fetch_databases.sh --tier extended` on the login node before submitting.
- **Re-verification drops more effectors → recall/annotation counts shift.** → Reversible corrections table; report the before/after delta explicitly rather than silently moving the numbers.
- **Agent content-filter trips on toxin/effector literature (seen before).** → Resolve those rows deterministically via UniProt/NCBI/PubMed in the manual pass.
- **Multi-day CX3 wallclock.** → Batched overnight submission; `qdel` stale jobs before each resubmit so they do not hold the ~2-job cap.

## Open Questions

- Are all 15 missing T5SS genomes present in `inputs_gb/` (just never run), or do some need fetching? Resolve before building the panel manifest.
- Do we want pooled enrichment across all 66 genomes (requires one big multi-genome run or a post-hoc pooled recompute), or are per-genome enrichment tables sufficient for the headline? Default: per-genome sufficient.
