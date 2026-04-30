# Output Files

## Directory Structure

```
results/
├── per_genome/
│   └── <sample_id>/
│       ├── proteins.faa            # Extracted protein sequences
│       ├── gene_info.tsv           # Gene metadata (position, product, etc.)
│       ├── gene_order.tsv          # Sorted gene order per contig
│       ├── macsyfinder/            # Raw MacSyFinder output
│       ├── ss_components.tsv       # Validated SS components
│       ├── valid_systems.tsv       # Complete systems (wholeness >= threshold)
│       ├── deeplocpro.tsv          # Subcellular localization predictions
│       ├── deepsece.tsv            # SS type predictions
│       ├── signalp.tsv             # Signal peptide predictions (if run)
│       ├── predictions_validated.tsv  # Cross-validated predictions
│       ├── proximity_results.tsv   # Proteins near SS components
│       ├── t5ss_substrates.tsv     # T5SS self-substrates
│       ├── substrates_filtered.tsv # Final filtered substrate list
│       ├── substrates_scored.tsv   # Substrates with confidence scores
│       └── annotations/
│           ├── blastp.csv          # BLASTp hits (if run)
│           ├── hhsuite.csv         # HH-suite hits (if run)
│           ├── interproscan.csv    # InterProScan domains (if run)
│           ├── protparam.csv       # Physicochemical properties (if run)
│           └── integrated.csv      # All annotations merged
├── master_substrates.csv           # Combined substrates across all genomes
├── master_all_proteins.csv         # All proteins with all annotations
├── enrichment_results.csv          # Statistical enrichment per SS type
├── enrichment_summary.txt          # Human-readable enrichment summary
├── figures/
│   ├── fig1_ss_type_distribution.png
│   ├── fig2_tool_coverage.png
│   ├── fig3_protein_lengths.png
│   └── fig4_physicochemical.png
├── ssign_report.html               # Interactive HTML report
└── ssign_report.txt                # Plain text report
```

## Key Output Files

### master_substrates.csv

The primary output — one row per predicted substrate across all genomes.

| Column                   | Description                                           |
| ------------------------ | ----------------------------------------------------- |
| `locus_tag`              | Protein identifier                                    |
| `sample_id`              | Source genome                                         |
| `tool`                   | How identified (DLP, DSE, T5SS-self, DLP+DSE)         |
| `nearby_ss_types`        | Secretion system(s) this substrate is associated with |
| `dlp_extracellular_prob` | DeepLocPro extracellular probability                  |
| `dse_ss_type`            | DeepSecE predicted SS type                            |
| `dse_max_prob`           | DeepSecE confidence                                   |
| `signalp_prediction`     | Signal peptide type (if SignalP run)                  |
| `blastp_hit_description` | Best BLASTp hit (if run)                              |
| `interpro_domains`       | InterPro domain IDs (if run)                          |
| `pfam_top1_*`            | HHpred Pfam results (if run)                          |
| `pdb_top1_*`             | HHpred PDB results (if run)                           |
| `gravy`, `mw_da`, etc.   | Physicochemical properties (if run)                   |

### substrates_filtered.tsv (per genome)

Per-genome substrate list after all filtering:

- System completeness >= threshold
- Component localization >= fraction threshold
- Excluded systems removed
- DeepSecE cross-genome leakage corrected
