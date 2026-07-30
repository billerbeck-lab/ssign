# Output files

What ssign writes to your `--outdir` after a successful run.

## Single-genome layout

Everything below is relative to your `--outdir`:

| Path | What it is |
|---|---|
| `<sample-id>_results.csv` | Main results, three chunks (see below) |
| `<sample-id>_results_raw.csv` | All annotations, no filtering, no column pruning |
| `<sample-id>_summary.txt` | Combined report text + enrichment summary |
| `figures/<sample-id>/01_secreted_by_genome.png` | Secreted proteins by SS type |
| `figures/<sample-id>/02_physicochemical.png` | Size & physicochemical properties |
| `figures/<sample-id>/03_cog_category_by_sstype.png` | COG functional category |
| `figures/<sample-id>/04_kegg_function_by_sstype.png` | KEGG function |
| `figures/<sample-id>/05_eggnog_description_by_sstype.png` | EggNOG description |
| `figures/<sample-id>/06_consensus_function_by_sstype.png` | Curated consensus function |
| `.ssign/<sample-id>_progress.json` | Resume manifest (used by `--resume`) |

Per-step intermediate files (proteins.faa, gene_info.tsv, individual tool
outputs, etc.) are written to a temporary work directory during the run and
removed on success. On a failure they are kept under `/tmp/ssign_*` and the
log line points at the location.

## Multi-genome layout

When you pass several genomes to one `ssign run`, each genome is written to its
own subdirectory under `--outdir` (the single-genome layout above, one folder
per genome). At the `--outdir` root, alongside those per-genome folders:

| Path | What it is |
|---|---|
| `<genome>/` | Per-genome files, one folder per genome (the single-genome layout above) |
| `combined_results.csv` | Secreted proteins pooled across all genomes, with a leading `source_genome` column |
| `combined_summary.txt` | Aggregated report: pooled counts, per-type totals, per-tool coverage |
| `cross_genome_orthologs.csv` | Every pooled substrate with its `sample_id` and cross-genome `ortholog_group` (BLAST+ required) |
| `cross_genome_ortholog_groups.csv` | Per group: `n_members`, `n_genomes`, `genomes`, `members`, `mean_pident` |
| `figures/0N_pooled_*.png` | Curated set over all genomes (01–06 numbering) |
| `figures/07_cross_genome_orthologs.png` | Ortholog conservation, 2-panel (BLAST+ required) |
| `figures/pooled_enrichment_fold[_combined].png` | Fold-enrichment charts (with `--enrichment-stats`) |

There is no combined raw CSV; each genome keeps its own
`<genome>/<genome>_results_raw.csv`.

The cross-genome ortholog step pools every genome's substrates, runs one
all-vs-all BLASTp, and clusters the result. It only runs for two or more genomes
and only when NCBI BLAST+ (`blastp` + `makeblastdb`) is on PATH; otherwise it is
skipped. Each genome's integrated CSV also gains an `xg_ortholog_group` column
linking its substrates to the cross-genome groups.

## `<sample-id>_results.csv` (main results)

The file opens with a short `#`-prefixed overview block naming the three
sections; skip it when parsing. Below it are up to three
chunks separated by blank lines, each with a `#`-prefixed header. Empty
chunks are omitted (e.g. genomes with no "other" systems will not have a
chunk 3):

1. `# Secreted Proteins`, one row per predicted substrate.
2. `# Secretion Systems (with secreted proteins)`, one row per system or
   component, for systems whose neighbourhoods contained at least one
   substrate.
3. `# Secretion Systems (other)`, systems detected without high-confidence
   substrates.

> **T5SS substrates:** a detected T5SS component is reported as a substrate only
> if it has evidence: DeepLocPro localizes it (extracellular/OM; the T5bSS
> translocator is outer-membrane only) OR SignalP predicts a Sec signal peptide.
> A T5 component with neither is not reported (it is still counted in the domain
> audit / secretion-system detection, just not called a substrate).

### Chunk 1 column reference (Secreted Proteins)

Columns appear in this order when present; missing columns indicate the
producing step was skipped or had no output.

| Group | Columns |
|---|---|
| Identity | `locus_tag`, `sample_id` |
| Annotation consensus | `broad_consensus_annotation`, `broad_annotation`, `detailed_annotation`, `detailed_consensus_annotation`, `evidence_keywords`, `n_tools_agreeing`, `n_tools_with_hits`, `concordance_ratio`, `confidence_tier` |
| Physicochemical | `aa_length`, `gravy`, `mw_da`, `isoelectric_point`, `charge_ph7`, `instability_index`, `aromaticity` |
| Secretion-system context | `nearby_ss_types`, `secretion_evidence`, `is_secreted`, `n_prediction_tools_agreeing` (0-2: how many of DeepLocPro + DeepSecE flagged secretion; SignalP is evidence-only and does not count) |
| DeepLocPro | `predicted_localization`, `dlp_extracellular_prob`, `dlp_max_localization`, `dlp_max_probability`, `periplasmic_prob`, `outer_membrane_prob`, `cytoplasmic_prob` |
| DeepSecE | `dse_ss_type`, `dse_max_prob` |
| SignalP | `signalp_prediction`, `signalp_probability`, `signalp_cs_position` |
| Original GenBank | `gbff_annotation` |
| BLASTp | `blastp_hit_accession`, `blastp_hit_description`, `blastp_pident`, `blastp_qcov`, `blastp_evalue` |
| HHpred Pfam | `pfam_top1_id`, `pfam_top1_description`, `pfam_top1_probability`, `pfam_top1_evalue`, `pfam_top1_score` |
| HHpred PDB | `pdb_top1_id`, `pdb_top1_description`, `pdb_top1_probability`, `pdb_top1_evalue`, `pdb_top1_score` |
| InterProScan | `interpro_domains`, `interpro_go_terms`, `interpro_pfam_ids`, `interpro_descriptions` |
| Ortholog groups | `ortholog_group`, `og_n_members`, `og_mean_pident` |
| Tool inventory | `annotation_tools` |
| T5aSS annotation source | `t5_annotation_source` (`passenger` or `full`; empty for non-T5aSS rows). See [`design_decisions.md` § 4.3](../explanation/design_decisions.md). |
| T5aSS whole-protein second pass (opt-in via `--t5ass-annotate-whole`) | `t5ass_whole_eggnog_*`, `t5ass_whole_blastp_*`, `t5ass_whole_ecod_top1_*`, `t5ass_whole_pfam_top1_*`, `t5ass_whole_pdb_top1_*`, `t5ass_whole_<protparam>_*`. Each mirrors the corresponding default-pass column on T5aSS substrates only. |
| Sequence | `sequence` (always last when present) |

Any tool-specific column not listed above (e.g. EggNOG, pLM-BLAST fields)
appears alphabetically after the last priority group, before `sequence`.

### Chunk 2 + 3 column reference (Secretion Systems)

ssign writes its own condensed schema (not a MacSyFinder table passthrough). Each
row is either a whole system or one of its components, tagged by `record_type`:

| `record_type` | Columns |
|---|---|
| `system` | `record_type`, `sample_id`, `sys_id`, `ss_type`, `wholeness`, `n_components`, `excluded` |
| `component` | `record_type`, `sample_id`, `sys_id`, `ss_type`, `locus_tag`, `gene_name`, `gene_status`, `wholeness`, `excluded` |

Excluded systems (default: Flagellum, Tad, and the type-IV pili / uptake
appendages T4aP, T4bP, MSH, ComM, Archaeal-T4P) and their components do not
appear in either chunk. T3SS is not excluded by default, so it does appear.

## `<sample-id>_results_raw.csv` (full annotations)

Every column ssign computed for every protein that reached the integration
step, with no filtering. Used for further downstream analysis or for
debugging why a particular protein was or was not called as a substrate.

Same columns as Chunk 1 of `_results.csv` plus any tool-specific columns
that did not make the priority list (the `_results.csv` priority order is
purely cosmetic; nothing is dropped from raw).

## `<sample-id>_summary.txt`

Plain text concatenation of:

1. The HTML report's text version (substrate counts, per-SS breakdowns, tool
   contribution summary).
2. The enrichment-analysis table (only with `--enrichment-stats`): one row per
   (`sample_id`, `ss_type`, `tool`, `mode`) with `observed`, `n_mask`,
   `null_mean`, `fold`, `p_perm`, `qvalue`, `significant`, `n_null`. The
   predictor column is `tool` (DeepLocPro / DeepSecE / SignalP / `COMBINED`),
   and `sample_id` is the leading column. `mode` is `window` for ordinary
   types (secreted-predicted proteins clustering near the components) or `self`
   for autotransporter self-detection. **T5aSS/T5cSS emit two results**: a `self`
   row-set (the autotransporter component detecting itself) AND a `window`
   "hitchhiker" row-set (secreted-predicted neighbours that may piggyback through
   the T5 pore; see `docs/explanation/design_decisions.md` § 5.2). Predictors are
   DeepLocPro, DeepSecE, and SignalP, plus a `COMBINED` row that pools the two
   relevant predictors per type: DLP-or-SignalP for the Sec-dependent T5 results
   (autotransporter self + T5bSS), DLP-or-DSE for every other window (including
   the T5a/c hitchhiker window), DLP-only for T3SS.
   `n_null` is the null sample size: the exact per-genome rotations (n genes give
   n-1 offsets) for a single genome, or 10000 Monte-Carlo draws when pooled. The
   exact pooled null (every joint combination of per-genome offsets) is too large
   to enumerate, so random joint rotations are drawn instead; each genome still
   contributes its own exact rotation set, only the join is sampled.

> **Reading significance:** the test's power scales with how many genomes and
> loci contribute. A **single-genome** run often shows `significant = False`
> even for real, correctly-detected systems, simply because a few loci against a
> whole-genome null can't reach q < 0.05. **A non-significant bar does not mean
> the system is absent or wrong**: it just means there wasn't enough statistical
> power. Pool several genomes (multi-genome run) for a powered test.

With `--enrichment-stats` the run also writes
`<sample-id>_enrichment_stats.tsv` (the table above),
`<sample-id>_enrichment_nulls.npz` (per-type nulls, used to pool the fold/p
across genomes), and two bar charts under `figures/<sample-id>/`:
`<sample-id>_enrichment_fold.png` (one bar per predictor per SS type) and
`<sample-id>_enrichment_fold_combined.png` (one combined bar per type). Bars are
fold enrichment (observed / expected) annotated with BH q-value significance
stars; non-significant bars keep their colour but are faded.

For a multi-genome run, the genomes' results are additionally pooled into
`pooled_enrichment_stats.tsv` and the same charts at
`figures/pooled_enrichment_fold[_combined].png`, computed over all genomes.

## `figures/<sample-id>/*.png`

Summary figures rendered at `--dpi` (default 300). Toggle individual figures
with the `--fig-*` flags listed in [`reference/cli.md`](cli.md). These are
summary-quality.

## `.ssign/<sample-id>_progress.json`

Resume manifest. Records every successful step plus the temp work-dir path,
so that `--resume` can skip already-completed steps after a partial failure.
Not meant to be read by the user.
