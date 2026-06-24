## 1. Phase A: full-corpus re-verification

- [ ] 1.1 Build the re-verification input: dump all 337 `positives_all.tsv` effectors with their accession, locus, ss_type, gene, organism, family, citation_quote, and primary-reference DOI into one ledger-input TSV.
- [ ] 1.2 Define the strict anti-hallucination output schema (echoed-identifier-only; verdict ∈ confirmed/fix/drop; provenance quote required) and the per-effector agent prompt, reusing the tasks 74-84 method.
- [ ] 1.3 Run the parallel blind-agent re-verification over all 337 effectors (batched), capturing per-agent verdicts + provenance; resolve content-filter-tripped toxin rows deterministically via UniProt/NCBI/PubMed.
- [ ] 1.4 Verify every primary-reference DOI resolves AND supports the secretion claim; flag citation defects, find replacement references where possible.
- [ ] 1.5 Manual reconciliation pass: adjudicate agent disagreements and unconfirmed claims into a single verdict per effector; write the final verification ledger (337 rows, no blanks).
- [ ] 1.6 Emit the corrections table (drops + fixes) in the scripts 55/56 format and wire it into `clean_dataset` reversibly (new `removed_*`/corrections file, no figure-script hard-coding).
- [ ] 1.7 Regenerate recall (01-06) + annotation (07/08) figures from the updated `clean_dataset`; record the before/after count delta in NOTES.

## 2. Phase B: build the corrected genome panel

- [x] 2.1 Confirm which T5SS genomes are already staged in `inputs_gb/`; fetch any absent ones (top open question RESOLVED, `scripts/57_build_rerun_panel.py`): 13 of 20 T5SS genomes already in 10 staged units; 7 are cache-not-staged (stage from `refseq_cache`, **0 need a network fetch**). Accession drift (RefSeq NC_ vs INSDC) resolved via the Datasets API + eutils `assemblyacc`.
- [x] 2.2 Derive the panel manifest (`data/phase2/rerun_panel_manifest.tsv`): the real numbers are **59 genomes** (51 proximity + 8 net T5SS = 52 staged-to-run + 7 stage-from-cache), with **15** staged units dropped, NOT the 66/16 estimated in the proposal. `AE002098` (N. meningitidis MC58, listed as `NC_003112.2` in the T5SS table) was rescued from the drop list because it bears a T5SS effector. Cross-checked: no dropped unit secretly carries a T5SS effector.
- [ ] 2.2b Stage the 7 cache-not-staged T5SS genomes from `refseq_cache` into `inputs_gb/` `.gbff` units (extend script 22's staging) so the panel is submission-ready on CX3.
- [ ] 2.3 Dry-run `submit_batched_overnight.sh --dry-run --small --enrichment-stats` over a sample of the panel to confirm `gpu_type=RTX6000`, colon-joined genomes, and the enrichment flag.

## 3. Phase B: run on CX3 (user-driven)

- [ ] 3.1 On CX3: `git pull` the branch, run `ssign doctor --tier extended`, re-fetch any purged `$EPHEMERAL` databases.
- [ ] 3.2 Submit the 66-genome panel in overnight batches via `cx3-submit`; monitor each batch to `Pipeline complete` / per-genome `N/N steps succeeded`.
- [ ] 3.3 Retrieve per-genome + pooled enrichment TSVs, figures, `combined_results.csv`, and run logs to the laptop.

## 4. Phase B: reconcile + close

- [ ] 4.1 Confirm the four RTX T1SS toxins (HlyA, ApxIA, LtxA, LktA) emit as secreted from their real full-assembly runs; if so, mark the `clean_dataset` staging fix retired (flag stays as fallback).
- [ ] 4.2 Grade annotation correctness for the 15 previously-ungraded T5SS effectors from the rerun; extend `annotation_correctness.py` and regenerate 07/08.
- [ ] 4.3 Verify the per-SS-type circular-shift enrichment emits sane fold/p/q across the panel (including a T3SS DLP row, no T3SS DSE row); close task 70.
- [ ] 4.4 Reconcile rerun recall against figure 06 (does the real run match the strict-emission counts?); document any discrepancies in NOTES.
- [ ] 4.5 Update the benchmark docs/methods with the re-verified answer key and the tier-2 rerun provenance; mark tasks 70 + 71 done.
