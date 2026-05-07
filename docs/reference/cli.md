# Configuration

All parameters can be passed on the CLI (`ssign run INPUT --flag value …`)
or in the Streamlit GUI sidebar. Boolean flags use Python 3.9+
`argparse.BooleanOptionalAction`, so each `--<flag>` has a `--no-<flag>`
inverse (e.g. `--skip-blastp` / `--no-skip-blastp`).

## Input/Output

| Parameter        | Default     | Description                                                      |
| ---------------- | ----------- | ---------------------------------------------------------------- |
| `INPUT_PATH`     | (required)  | Positional — path to the genome file (GenBank or FASTA contigs). |
| `--outdir`       | `./results` | Output directory                                                 |
| `--sample-id`    | (auto)      | Sample name used in output filenames; defaults to input stem.    |
| `--resume`       | `false`     | Skip steps that already succeeded in a previous run.             |

## Phase 1: Input Processing

| Parameter                  | Default | Description                                                     |
| -------------------------- | ------- | --------------------------------------------------------------- |
| `--use-input-annotations`  | `false` | Skip Bakta re-annotation of GenBank input (preserve incoming annotations) |
| `--bakta-db`               | `null`  | Path to Bakta database (~2 GB light, ~30 GB full)               |

## Phase 2: Secretion System Detection

| Parameter               | Default              | Description                             |
| ----------------------- | -------------------- | --------------------------------------- |
| `--wholeness-threshold` | `0.8`                  | MacSyFinder system completeness minimum                |
| `--excluded-systems`    | `Flagellum Tad T3SS`   | Space-separated SS types to exclude (`nargs="+"`)      |

## Phase 3: Secreted Protein Prediction

| Parameter           | Default  | Description                                                |
| ------------------- | -------- | ---------------------------------------------------------- |
| `--conf-threshold`  | `0.8`    | DeepLocPro extracellular probability minimum               |
| `--deeplocpro-path` | `null`   | Path to local DeepLocPro install (omit to use remote API)  |
| `--signalp-path`    | `null`   | Path to local SignalP 6.0 install (omit to use remote API) |
| `--skip-signalp`    | `false`  | Skip SignalP                                               |
| `--skip-deepsece`   | `false`  | Skip DeepSecE (`--skip-deepsece` to opt out; `--no-skip-deepsece` to opt back in) |

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
| `--hhsuite-pfam-db`      | `null`  | Pfam-A database path (or `$SSIGN_HHSUITE_PFAM`)    |
| `--hhsuite-pdb70-db`     | `null`  | PDB70 database path (or `$SSIGN_HHSUITE_PDB70`)    |
| `--hhsuite-uniclust-db`  | `null`  | UniClust database path (or `$SSIGN_HHSUITE_UNICLUST`) |
| `--skip-eggnog`          | `true`  | Skip EggNOG-mapper (off by default)              |
| `--eggnog-db`            | `null`  | EggNOG database path                             |
| `--skip-interproscan`    | `false` | Skip InterProScan                                |
| `--interproscan-db`      | `null`  | InterProScan install directory                   |
| `--skip-plmblast`        | `true`  | Skip pLM-BLAST (off by default)                  |
| `--plmblast-db`          | `null`  | ECOD70 database path                             |
| `--skip-plm-effector`    | `true`  | Skip PLM-Effector (off by default — GPU-heavy)   |
| `--skip-protparam`       | `false` | Skip ProtParam physicochemical property compute  |

## Phase 6: Reporting

| Parameter | Default | Description       |
| --------- | ------- | ----------------- |
| `--dpi`   | `300`   | Figure resolution |

## Resources

| Parameter           | Default     | Description                                                            |
| ------------------- | ----------- | ---------------------------------------------------------------------- |
| `--cpu-per-genome`  | (CPU count) | CPUs available to per-genome subtools (e.g. macsyfinder -w, BLAST -num_threads) |
| `--bakta-threads`   | `4`         | Threads passed specifically to Bakta.                                  |

## HPC

The CLI runs as a single process; HPC use is straightforward via job submission:

### PBS (e.g. Imperial CX3)

```bash
qsub -l select=1:ncpus=16:mem=64gb:walltime=24:00:00 -- bash -c '
  module load anaconda3
  source activate ssign
  ssign run input.gbff --outdir results --cpu-per-genome 16
'
```

### SLURM

```bash
sbatch --cpus-per-task=16 --mem=64G --time=24:00:00 --wrap='
  ssign run input.gbff --outdir results --cpu-per-genome 16
'
```

For a Singularity-image workflow on HPC (recommended for reproducibility),
see the `containers/` README and the v1.0.0 release notes — a dedicated
HPC how-to lands as part of the Phase 5 docs overhaul.
