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
- [x] 2.2 Derive the enriched panel manifest (`data/phase2/rerun_panel_manifest.tsv`): **59 run-units = 57 distinct assemblies** (50 staged + 7 stage-from-cache), **15** dropped, NOT the 66/16 estimated. Per-genome columns added: NCBI-authoritative organism, `genome_group` (biosample + accession-alias union), systems, n_secreted_proteins, n_found_by_ssign, and the secreted-protein list (`*` = emitted by ssign). Findings: `AE002098` (N. meningitidis MC58 = `NC_003112.2`) rescued from drop; 2 same-assembly merges (PAO1 `AE004091`+`NC_002516.2`; A. fabrum C58 `NC_003063.2`+`NC_003065.3`); corpus organism mislabels corrected (`NC_006351`=B. pseudomallei not B. thailandensis, `NC_013716`=Citrobacter not EPEC, `NC_007508`=Xanthomonas not P. aeruginosa).
- [x] 2.2b Stage rerun inputs (`scripts/58_stage_rerun_inputs.py` -> `data/phase2/rerun_inputs.txt`, 57 inputs): 49 reused (already on CX3) + 8 NEW (7 cached T5SS staged, incl. P. mirabilis 2-replicon merge; + the C58 2-replicon merge). PAO1 collapsed to one copy (`NC_002516.2`; the two units are the same chromosome). All 8 new `.gbff` validated parseable (merged ones have 2 records).
- [x] 2.3 Submit driver written (`scripts/cx3/SUBMIT_rerun.sh`). DECISION (Teo): annotation ON, enrichment ON, NO `--use-input-annotations` (Bakta re-annotates everything, clean slate, all 57 identical). Confirmed `runner.py:1320`: GenBank input without `--use-input-annotations` re-annotates with Bakta by default, so the generic `submit_batched_overnight.sh` wrapper (extended tier, annotation on, no skip flags) is the right tool. One job per genome (per-genome ceiling ~4.4h from calibration: pLM-BLAST ~4h + EggNOG ~3h dominate), trickled 2 at a time, idempotent.
- [ ] 2.4 (Teo, on CX3) rsync `inputs_gb/` to the repo path + `git pull` + `ssign doctor --tier extended` (re-fetch purged DBs) + run `SUBMIT_rerun.sh` under tmux.

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
