# Configuration

All parameters can be passed on the CLI (`ssign run --foo bar`) or in the
Streamlit GUI sidebar.

## Input/Output

| Parameter     | Default     | Description                          |
| ------------- | ----------- | ------------------------------------ |
| `--input`     | (required)  | Path to genome file                  |
| `--outdir`    | `./results` | Output directory                     |
| `--sample-id` | (auto)      | Sample name used in output filenames |

## Phase 1: Input Processing

| Parameter                  | Default | Description                                                     |
| -------------------------- | ------- | --------------------------------------------------------------- |
| `--use-input-annotations`  | `false` | Skip Bakta re-annotation of GenBank input (preserve incoming annotations) |
| `--bakta-db`               | `null`  | Path to Bakta database (~2 GB light, ~30 GB full)               |

## Phase 2: Secretion System Detection

| Parameter               | Default              | Description                             |
| ----------------------- | -------------------- | --------------------------------------- |
| `--wholeness-threshold` | `0.8`                | MacSyFinder system completeness minimum |
| `--excluded-systems`    | `Flagellum,Tad,T3SS` | Comma-separated SS types to exclude     |

## Phase 3: Secreted Protein Prediction

| Parameter           | Default  | Description                                                |
| ------------------- | -------- | ---------------------------------------------------------- |
| `--conf-threshold`  | `0.8`    | DeepLocPro extracellular probability minimum               |
| `--deeplocpro-path` | `null`   | Path to local DeepLocPro install (omit to use remote API)  |
| `--signalp-path`    | `null`   | Path to local SignalP 6.0 install (omit to use remote API) |
| `--skip-signalp`    | `false`  | Skip SignalP                                               |
| `--skip-deepsece`   | `true`   | Skip DeepSecE (off by default — needs ~7 GB ESM model)     |

## Phase 4: Substrate Identification

| Parameter                     | Default | Description                                   |
| ----------------------------- | ------- | --------------------------------------------- |
| `--proximity-window`          | `3`     | +/- N genes per SS component                  |
| `--required-fraction-correct` | `0.8`   | Fraction of SS components correctly localized |

## Phase 5: Annotation Tools

| Parameter                | Default | Description                                      |
| ------------------------ | ------- | ------------------------------------------------ |
| `--skip-blastp`          | `false` | Skip BLASTp                                      |
| `--blastp-db`            | `null`  | Path to local BLAST database (required to run)   |
| `--blastp-evalue`        | `1e-5`  | E-value threshold                                |
| `--blastp-min-pident`    | `80`    | Minimum percent identity                         |
| `--blastp-min-qcov`      | `80`    | Minimum query coverage                           |
| `--blastp-exclude-taxid` | `null`  | Taxonomy ID to exclude (query organism)          |
| `--skip-hhsuite`         | `true`  | Skip HH-suite (off by default — needs large DBs) |
| `--hhsuite-pfam`         | `null`  | Pfam-A database path                             |
| `--hhsuite-pdb70`        | `null`  | PDB70 database path                              |
| `--hhsuite-uniref30`     | `null`  | UniRef30 database path                           |
| `--skip-eggnog`          | `false` | Skip EggNOG-mapper                               |
| `--eggnog-db`            | `null`  | EggNOG database path                             |
| `--skip-interproscan`    | `false` | Skip InterProScan                                |
| `--interproscan-path`    | `null`  | InterProScan install directory                   |
| `--skip-plmblast`        | `true`  | Skip pLM-BLAST (off by default)                  |
| `--plmblast-ecod-db`     | `null`  | ECOD70 database path                             |
| `--skip-plm-effector`    | `false` | Skip PLM-Effector                                |
| `--skip-protparam`       | `false` | Skip ProtParam physicochemical property compute  |

## Phase 6: Reporting

| Parameter | Default | Description       |
| --------- | ------- | ----------------- |
| `--dpi`   | `300`   | Figure resolution |

## Resources

| Parameter   | Default        | Description                                               |
| ----------- | -------------- | --------------------------------------------------------- |
| `--threads` | (CPU count)    | Threads used by parallel-friendly steps (BLAST, HH-suite) |

## HPC

The CLI runs as a single process; HPC use is straightforward via job submission:

### PBS (e.g. Imperial CX3)

```bash
qsub -l select=1:ncpus=16:mem=64gb:walltime=24:00:00 -- bash -c '
  module load anaconda3
  source activate ssign
  ssign run input.gbff --outdir results --threads 16
'
```

### SLURM

```bash
sbatch --cpus-per-task=16 --mem=64G --time=24:00:00 --wrap='
  ssign run input.gbff --outdir results --threads 16
'
```

For a Singularity-image workflow on HPC (recommended for reproducibility), see
`docs/how-to/run_on_hpc.md`.
