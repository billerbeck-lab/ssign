# Output files

What ssign writes to your `--outdir` after a successful run.

## Single-genome layout

Everything below is relative to your `--outdir`:

| Path | What it is |
|---|---|
| `<sample-id>_results.csv` | Main results |
| `<sample-id>_results_raw.csv` | Raw output from tools |
| `<sample-id>_summary.txt` | Text summary of run |
| `figures/<sample-id>/01_secreted_by_genome.png` | Secreted proteins by SS type |
| `figures/<sample-id>/02_physicochemical.png` | Size & physicochemical properties |
| `figures/<sample-id>/03_cog_category_by_sstype.png` | COG functional category |
| `figures/<sample-id>/04_kegg_function_by_sstype.png` | KEGG function |
| `figures/<sample-id>/05_eggnog_description_by_sstype.png` | EggNOG description |
| `figures/<sample-id>/06_consensus_function_by_sstype.png` | Consensus function |
| `.ssign/<sample-id>_progress.json` | Resume manifest (used by `--resume`) |

## Multi-genome layout

When you pass several genomes to one `ssign run`, each genome is written to its
own subdirectory under `--outdir` (the single-genome layout above, one folder
per genome) alongside these per-genome folders:

| Path | What it is |
|---|---|
| `<genome>/` | Per-genome files (the single-genome layout above) |
| `combined_results.csv` | Secreted proteins pooled across all genomes, with a leading `source_genome` column |
| `combined_summary.txt` | Aggregated report, pooled counts, per-type totals, per-tool coverage |
| `cross_genome_orthologs.csv` | Every pooled substrate with its `sample_id` and cross-genome `ortholog_group` (BLAST+ required) |
| `cross_genome_ortholog_groups.csv` | Per group: `n_members`, `n_genomes`, `genomes`, `members`, `mean_pident` |
| `figures/0N_pooled_*.png` | Same figures as single-genome runs |
| `figures/07_cross_genome_orthologs.png` | Ortholog conservation, 2-panel (BLAST+ required) |
| `figures/pooled_enrichment_fold[_combined].png` | Fold-enrichment charts (with `--enrichment-stats`) |

There is no combined raw CSV; each genome keeps its own
`<genome>/<genome>_results_raw.csv`.

The cross-genome ortholog step pools every genome's substrates and runs one
all-vs-all BLASTp. Each genome's integrated CSV also gains an `xg_ortholog_group` column
linking its substrates to the cross-genome groups.

## `<sample-id>_results.csv` (main results)

The file opens with a short `#`-prefixed overview block naming the three
sections

1. `# Secreted Proteins`, one row per predicted secreted protein.
2. `# Secretion Systems (with secreted proteins)`, one row per system or
   component, for systems whose neighbourhoods contained at least one
   substrate.
3. `# Secretion Systems (other)`, systems detected without high-confidence
   substrates.

> **T5SS substrates:** a detected T5SS component is reported as a substrate only
> if it has evidence: DeepLocPro localizes it (extracellular/OM) OR SignalP predicts a Sec signal peptide.
> A T5 component with neither is not reported (it is still counted in the domain
> audit / secretion-system detection, just not called a substrate).

### Chunk 1 column reference (Secreted Proteins)

Columns appear in this order.

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

| `record_type` | Columns |
|---|---|
| `system` | `record_type`, `sample_id`, `sys_id`, `ss_type`, `wholeness`, `n_components`, `excluded` |
| `component` | `record_type`, `sample_id`, `sys_id`, `ss_type`, `locus_tag`, `gene_name`, `gene_status`, `wholeness`, `excluded` |
