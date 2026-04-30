# Usage

## GUI

```bash
pip install ssign
ssign                     # opens Streamlit GUI in browser
```

The GUI lets you upload genomes, configure thresholds, select annotation tools,
and run the pipeline with a progress tracker. Results are saved to your chosen
output directory.

## CLI

```bash
ssign run input.gbff --outdir results
ssign run input.gbff --outdir results --skip-hhsuite --skip-plmblast
ssign run input.gbff --outdir results --resume   # resume after a failed step
```

For HPC and batch use, the CLI is the recommended interface. The GUI is for
single-genome interactive runs.

## Input Formats

| Format        | Extension               | Description                                                                |
| ------------- | ----------------------- | -------------------------------------------------------------------------- |
| GenBank       | `.gbff`, `.gbk`, `.gb`  | Annotated genome (recommended — contains protein sequences and gene order) |
| GFF3 + FASTA  | `.gff3` + `.fasta`      | Pass both with `--input-gff` and `--input-fasta`                           |
| FASTA contigs | `.fasta`, `.fna`, `.fa` | Raw contigs — Bakta (or Prodigal as fallback) predicts ORFs                |

## Annotation Tool Modes

Each annotation tool runs locally against a database the user installs once via
`scripts/fetch_databases.sh` (Phase 4b). Remote API modes were removed in v1.0.0
to make runs reproducible offline.

| Tool         | Database                                                   | Approx. size                |
| ------------ | ---------------------------------------------------------- | --------------------------- |
| BLASTp       | `--blastp-db /path` (NR or Swiss-Prot)                     | NR ~390 GB, Swiss-Prot ~2 GB |
| HH-suite     | `--hhsuite-pfam`, `--hhsuite-pdb70`, `--hhsuite-uniref30`  | Pfam ~3 GB, PDB70 ~22 GB, UniRef30 ~25 GB |
| InterProScan | `--interproscan-path /path/to/install`                     | ~24 GB unpacked             |
| EggNOG       | `--eggnog-db /path`                                        | ~50 GB                      |
| pLM-BLAST    | `--plmblast-ecod-db /path`                                 | ECOD70 ~10 GB               |

See `docs/install.md` for the install-tier table (base / extended / full).

## Skipping Tools

Any annotation tool can be skipped on the command line:

```bash
ssign run input.gbff --outdir results \
    --skip-signalp \
    --skip-hhsuite \
    --skip-plmblast
```

## System Filtering

By default, Flagellum, Tad, and T3SS are excluded from substrate identification:

```bash
ssign run input.gbff --outdir results \
    --excluded-systems 'Flagellum,Tad,T3SS'   # default
```

**Why T3SS is excluded by default:** MacSyFinder found zero T3SS across 74 Xanthomonas
genomes, yet DeepSecE predicted 1,808 T3SS substrates — mostly flagellar protein
misclassifications. Override with a different list if you're working with organisms
known to have T3SS.

## GUI options

```bash
ssign                        # default port 8501, opens browser
ssign --port 8502            # custom port
ssign --no-browser           # don't auto-open browser
```
