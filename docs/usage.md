# Usage

## Easy Mode (GUI)

```bash
pip install ssign
ssign                     # opens Streamlit GUI in browser
```

The GUI lets you upload genomes, configure thresholds, select annotation tools,
and run the pipeline with a progress tracker. Results are saved to your chosen
output directory.

## Power Mode (Nextflow CLI)

```bash
# Single genome (GenBank format — recommended)
nextflow run ssign/main.nf \
    --input genome.gbff \
    -profile local

# Multiple genomes via samplesheet
nextflow run ssign/main.nf \
    --input samplesheet.csv \
    -profile local

# Resume a failed run
nextflow run ssign/main.nf \
    --input genome.gbff \
    -profile local \
    -resume
```

## Input Formats

### Single Genome File

| Format        | Extension               | Description                                                                |
| ------------- | ----------------------- | -------------------------------------------------------------------------- |
| GenBank       | `.gbff`, `.gbk`, `.gb`  | Annotated genome (recommended — contains protein sequences and gene order) |
| GFF3 + FASTA  | `.gff3` + `.fasta`      | Provide both files via samplesheet                                         |
| FASTA contigs | `.fasta`, `.fna`, `.fa` | Raw contigs — Prodigal will predict ORFs                                   |

### Samplesheet (CSV)

For multiple genomes, create a CSV with columns:

```csv
sample,input_1,input_2
Genome_A,/path/to/genomeA.gbff,
Genome_B,/path/to/genomeB.gff3,/path/to/genomeB.fasta
Genome_C,/path/to/contigs.fasta,
```

- `sample`: Unique sample name (used in output file names)
- `input_1`: Path to primary input file (required)
- `input_2`: Path to secondary file (only for GFF3+FASTA pairs)

## Execution Profiles

| Profile     | Container   | Scheduler | Use case                             |
| ----------- | ----------- | --------- | ------------------------------------ |
| `local`     | Docker      | Local     | Desktop/laptop                       |
| `hpc_pbs`   | Singularity | PBS       | HPC with PBS (e.g. Imperial College) |
| `hpc_slurm` | Singularity | SLURM     | HPC with SLURM                       |
| `cloud_aws` | Docker      | AWS Batch | Cloud execution                      |
| `test`      | Docker      | Local     | CI testing with small genomes        |

## Annotation Tool Modes

Each annotation tool supports **local** and/or **remote** mode:

| Tool         | Local                                               | Remote                                 | Recommended              | Notes                                         |
| ------------ | --------------------------------------------------- | -------------------------------------- | ------------------------ | --------------------------------------------- |
| BLASTp       | `--blastp_db /path` (required)                      | **removed in v1.0.0**                  | Local only               | NR: ~390GB, Swiss-Prot: ~1.5GB                |
| HH-suite     | `--hhsuite_mode local` + DB paths                   | `--hhsuite_mode remote` (default)      | Remote for <100 proteins | Pfam-A: ~3GB, PDB70: ~20GB, UniClust30: ~25GB |
| InterProScan | `--interproscan_mode local --interproscan_db /path` | `--interproscan_mode remote` (default) | Remote for <200 proteins | Local install: ~80GB                          |
| pLM-BLAST    | `--plmblast_ecod_db /path`                          | N/A (MPI ECOD70 is broken)             | Local only               | ECOD70: ~5-10GB                               |

### Rate Limits (Remote Mode)

| Tool                   | Rate Limit              | Practical Limit               |
| ---------------------- | ----------------------- | ----------------------------- |
| HH-suite (MPI Toolkit) | 200 jobs/hr, 2000/day   | ~100 proteins/run comfortable |
| InterProScan (EBI)     | 30 req/sec, 25k seq/day | ~200 proteins/run comfortable |

## Skipping Tools

Any tool can be skipped:

```bash
nextflow run ssign/main.nf \
    --input genome.gbff \
    --deeplocpro_path /path/to/deeplocpro \
    --skip_signalp \
    --skip_hhsuite \
    --skip_plmblast \
    --skip_structure \
    -profile local
```

## System Filtering

By default, Flagellum, Tad, and T3SS are excluded from substrate identification:

```bash
# Default (recommended)
--excluded_systems 'Flagellum,Tad,T3SS'

# Include T3SS
--excluded_systems 'Flagellum,Tad'

# Include everything
--excluded_systems ''
```

**Why T3SS is excluded by default:** MacSyFinder found zero T3SS across 74 Xanthomonas
genomes, yet DeepSecE predicted 1,808 T3SS substrates — mostly flagellar protein
misclassifications. T3SS filtering can be disabled if working with organisms known
to have T3SS.

## GUI (Easy Mode)

The GUI is the primary interface for Easy Mode. It launches automatically via
the `ssign` command. You can also start it directly:

```bash
ssign                        # recommended
ssign --port 8502            # custom port
ssign --no-browser           # don't auto-open browser
```
