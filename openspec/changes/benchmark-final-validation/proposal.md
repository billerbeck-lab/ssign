## Why

The effector-recovery benchmark now reports 45/60 reachable system instances and graded annotation correctness, but two confidence gaps remain before it can back a publication claim. First, the anti-hallucination audit (tasks 74-84) only re-verified the proximity-reachable subset (~245 proteins); the other corpus rows and many primary-reference DOIs were never independently checked, so the answer key is not provably 100% real. Second, the numbers were measured against a single extended-tier fleet (`/tmp/ssign_fleet_67`) that has known defects: T5SS was annotated for only 5 of 20 genomes, the four RTX T1SS toxins were rescued by a manual staging fix rather than a real run, and enrichment was not computed fleet-wide. We need one final, fully-verified corpus and one clean tier-2 rerun so every headline figure is reproducible from a real ssign run on a defensible answer key.

## What Changes

- **Re-verify 100% of the effector corpus.** Independently check all 337 `positives_all.tsv` effectors (T1:25, T2:54, T3:136, T4:48, T5:19, T6:55) and every primary-reference DOI, using the proven parallel blind-agent + manual method with the strict anti-hallucination output schema. Produce a verification ledger (per-effector verdict + provenance) and fold any new drops/fixes into `clean_dataset` reversibly, the same way scripts 55/56 did.
- **Rerun the benchmark genome panel at extended tier (tier 2) on CX3.** Replace the defective fleet with a corrected 66-genome panel: the 51 genomes carrying a testable proximity system plus the 15 T5SS genomes missing from the current fleet, dropping the 16 input genomes that carry no testable answer-key system. Run with `--enrichment-stats` (forces whole-genome DLP+DSE), T3SS detection on (now default), and proper T5SS handling, via the `cx3-submit` skill.
- **Reconcile the rerun against the recall + annotation figures.** Confirm the four RTX T1SS toxins now emit from real full-assembly runs (retiring the staging fix), grade the 15 currently-ungraded T5SS effectors, and verify the per-SS-type circular-shift enrichment emits sane fold/p/q (closing task 70).
- This is **analysis + validation only**: no ssign product/runtime code changes. The benchmark scripts, `clean_dataset`, and the summary figures are the only things edited.

## Capabilities

### New Capabilities
- `effector-corpus-reverification`: the 100%-coverage re-verification of every corpus effector and its primary-reference DOI (blind-agent + manual procedure, output schema, verification ledger, and the reversible fold into `clean_dataset`).
- `tier2-benchmark-rerun`: the corrected 66-genome panel definition, the extended-tier CX3 run procedure (enrichment + T3SS-on + T5SS), and the reconciliation of its outputs against the recall, annotation, and enrichment figures.

### Modified Capabilities
<!-- None. This change adds analysis/validation capabilities only and touches no existing ssign product spec. -->

## Impact

- **Re-verification** reads `data/dataset/positives_all.tsv` and the corpus citation columns; writes a new ledger under `validation_sweeps/benchmark/data/` and (if drops/fixes result) a new `removed_*`/corrections table consumed by `clean_dataset`. May shift recall/annotation counts; all figures regenerate from the shared `clean_dataset` filter.
- **CX3 rerun** consumes the gitignored `validation_sweeps/benchmark/inputs_gb/` genomes (66 of the 67), needs the extended-tier databases on `$EPHEMERAL` (re-fetch if purged), and produces per-genome + pooled enrichment outputs retrieved to the laptop. GPU compute: ~30-48 min/genome x 66 (the `--enrichment-stats` whole-genome DLP/DSE is the cost driver).
- **External dependencies**: PubMed/web + UniProt/NCBI for the verification agents; no change to ssign's runtime code or dependencies.
- **Closes** open benchmark tasks 70 (enrichment working?) and the annotation-tool benchmark (71); supersedes the manual T1SS staging fix in `clean_dataset`.
