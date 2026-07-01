# Output files

What ssign writes to your `--outdir` after a successful run.

## Single-genome layout

```
results/
├── <sample-id>_results.csv          Main results, three chunks (see below)
├── <sample-id>_results_raw.csv      All annotations, no filtering, no column pruning
├── <sample-id>_summary.txt          Combined report text + enrichment summary
├── figures/
│   └── <sample-id>/
│       ├── fig_category.png
│       ├── fig_ss_comp.png
│       ├── fig_tool_heatmap.png
│       ├── fig_substrate_count.png
│       └── fig_func_summary.png
└── .ssign/
    └── <sample-id>_progress.json    Resume manifest (used by --resume)
```

Per-step intermediate files (proteins.faa, gene_info.tsv, individual tool
outputs, etc.) are written to a temporary work directory during the run and
removed on success. On a failure they are kept under `/tmp/ssign_*` and the
log line points at the location.

## Multi-genome layout (GUI batch run)

When the GUI processes several genomes in one session, the per-genome files
above are produced for each, and three combined files are added at the
output-directory root:

```
results/
├── <each genome's files as above>
├── ssign_results.csv                Combined results across all genomes
├── ssign_results_raw.csv            Combined raw annotations
└── ssign_summary.txt                Combined summary
```

## `<sample-id>_results.csv` (main results)

Up to three chunks separated by blank lines, each with a `#`-prefixed
header. Empty chunks are omitted (e.g. genomes with no "other" systems
will not have a chunk 3):

1. `# Secreted Proteins` — one row per predicted substrate.
2. `# Secretion Systems (with secreted proteins)` — one row per system or
   component, for systems whose neighbourhoods contained at least one
   substrate.
3. `# Secretion Systems (other)` — systems detected without high-confidence
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
| Secretion-system context | `nearby_ss_types`, `secretion_evidence`, `is_secreted` |
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
| T5aSS annotation source | `t5_annotation_source` (`passenger` or `full`; empty for non-T5aSS rows). See `design_decisions.md` § 4.3. |
| T5aSS whole-protein second pass (opt-in via `--t5ass-annotate-whole`) | `t5ass_whole_eggnog_*`, `t5ass_whole_blastp_*`, `t5ass_whole_ecod_top1_*`, `t5ass_whole_pfam_top1_*`, `t5ass_whole_pdb_top1_*`, `t5ass_whole_<protparam>_*`. Each mirrors the corresponding default-pass column on T5aSS substrates only. |
| Sequence | `sequence` (always last when present) |

Any tool-specific column not listed above (e.g. EggNOG, pLM-BLAST, PLM-Effector
fields) appears alphabetically after the last priority group, before
`sequence`.

### Chunk 2 + 3 column reference (Secretion Systems)

Columns mirror MacSyFinder's output table with two added at the front:

| Column | Source |
|---|---|
| `record_type` | `system` or `component` (added by ssign for chunked-CSV navigation). |
| `sample_id` | Genome ID. |
| `ss_type`, `wholeness_score`, `model_fqn`, `replicon`, `genes` etc. | MacSyFinder columns; see [MacSyFinder docs](https://macsyfinder.readthedocs.io/) for the full list. |

Excluded systems (default: Flagellum, Tad, T3SS) and their components do not
appear in either chunk.

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
   (SS type, predictor) with `mode` (window/self), `observed`, `n_mask`,
   `null_mean`, `fold`, `p_perm`, `qvalue`, `significant`, `n_null`. Predictors
   are DeepLocPro, DeepSecE, and SignalP, plus a `COMBINED` row that pools the
   two relevant predictors per type (DLP-or-DSE, or DLP-or-SignalP for T5SS).
   `n_null` is the null sample size: the exact per-genome rotations (~one per
   gene) for a single genome, or 10000 Monte-Carlo draws when pooled.

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
summary-quality; paper-grade figures are regenerated separately from scripts
in the top-level `figures/` directory.

## `.ssign/<sample-id>_progress.json`

Resume manifest. Records every successful step plus the temp work-dir path,
so that `--resume` can skip already-completed steps after a partial failure.
Not meant to be read by the user.
