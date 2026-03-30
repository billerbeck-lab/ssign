# Configuration

## All Parameters

### Input/Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Path to genome file or samplesheet CSV |
| `--outdir` | `./results` | Output directory |

### Phase 1: Input Processing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_bakta` | `false` | Use Bakta instead of Prodigal for unannotated inputs |
| `--bakta_db` | `null` | Path to Bakta database (~30GB) |

### Phase 2: Secretion System Detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wholeness_threshold` | `0.8` | MacSyFinder system completeness minimum |
| `--excluded_systems` | `Flagellum,Tad,T3SS` | Comma-separated SS types to exclude |

### Phase 3: Secreted Protein Prediction

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--conf_threshold` | `0.8` | DeepLocPro extracellular probability minimum |
| `--deeplocpro_path` | `null` | Path to local DeepLocPro install (not needed for remote mode) |
| `--signalp_path` | `null` | Path to local SignalP 6.0 install (not needed for remote mode) |
| `--skip_signalp` | `false` | Skip SignalP |

### Phase 4: Substrate Identification

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--proximity_window` | `3` | +/- N genes per SS component |
| `--required_fraction_correct` | `0.8` | Fraction of SS components correctly localized |

### Phase 5: Annotation Tools

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_blastp` | `false` | Skip BLASTp |
| `--blastp_mode` | `remote` | `local` or `remote` |
| `--blastp_db` | `null` | Path to local BLAST database |
| `--blastp_evalue` | `1e-5` | E-value threshold |
| `--blastp_min_pident` | `80` | Minimum percent identity |
| `--blastp_min_qcov` | `80` | Minimum query coverage |
| `--blastp_exclude_taxid` | `null` | Taxonomy ID to exclude (query organism) |
| `--skip_hhsuite` | `true` | Skip HH-suite (off by default — needs large DBs) |
| `--hhsuite_mode` | `remote` | `local` or `remote` |
| `--hhsuite_pfam_db` | `null` | Pfam-A database path |
| `--hhsuite_pdb70_db` | `null` | PDB70 database path |
| `--hhsuite_uniclust_db` | `null` | UniClust30 database path |
| `--skip_interproscan` | `false` | Skip InterProScan |
| `--interproscan_mode` | `remote` | `local` or `remote` |
| `--interproscan_db` | `null` | InterProScan install path |
| `--skip_foldseek` | `false` | Skip Foldseek |
| `--foldseek_db` | `null` | Foldseek database path |
| `--foldseek_evalue` | `0.001` | Foldseek E-value threshold |
| `--foldseek_tmscore_threshold` | `0.8` | Minimum TM-score (query-normalized) |
| `--skip_plmblast` | `true` | Skip pLM-BLAST (off by default) |
| `--plmblast_ecod_db` | `null` | ECOD70 database path |
| `--skip_protparam` | `false` | Skip ProtParam |
| `--skip_structure` | `false` | Skip structure prediction |
| `--plddt_threshold` | `70` | Minimum mean pLDDT for structure acceptance |

### Phase 6: Reporting

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dpi` | `300` | Figure resolution |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `240.h` | Maximum wall time per process |

## HPC Configuration

### PBS (e.g. Imperial College)

```bash
nextflow run ssign/main.nf \
    --input samplesheet.csv \
    -profile hpc_pbs
```

The PBS profile uses Singularity containers and submits jobs via `qsub`.
Customize in `conf/hpc_pbs.config`:

```groovy
process {
    executor = 'pbs'
    queue = 'your_queue_name'
    clusterOptions = '-l walltime=48:00:00'
}
singularity {
    cacheDir = '/path/to/singularity/cache'
}
```

### SLURM

```bash
nextflow run ssign/main.nf \
    --input samplesheet.csv \
    -profile hpc_slurm
```

Customize in `conf/hpc_slurm.config`.

## Custom Container Paths

Override container images in your Nextflow config:

```groovy
process {
    withLabel: 'process_annotation' {
        container = '/path/to/custom/ssign-annotation.sif'
    }
}
```
