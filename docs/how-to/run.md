# Running ssign

This page covers the basic invocation, how to
pick flags for common situations, and how to run on an HPC cluster.

For install steps (both the container and native routes) see
[`install.md`](install.md). For the full flag list see
[`reference/cli.md`](../reference/cli.md); for environment variables see
[`reference/env_vars.md`](../reference/env_vars.md).

## Basic invocation

```bash
ssign run input.gbff --outdir results
```

That runs the whole pipeline on one genome and writes to `results/`:

The sample id defaults to the input filename's stem; override it with
`--sample-id`. See
[`reference/output_files.md`](../reference/output_files.md) for the full output
layout.

Two more subcommands help around a run:

```bash
ssign doctor                             # verify the install before running
ssign fetch-databases --tier extended    # download reference databases
```

### Resuming a failed run

```bash
ssign run input.gbff --outdir results --resume
```

`--resume` skips steps that already have a successful entry in the progress
manifest at `<outdir>/.ssign/<sample-id>_progress.json`. If you changed flags
between runs, ssign treats the changed config as a new run and starts over.

## Single vs multiple genomes

Pass one path for a single-genome run, or several paths to run them as one
batched job:

```bash
# single
ssign run genome.gbff --outdir results

# multiple (batched)
ssign run a.gbff b.gbff c.gbff --outdir results
```

In a batched run, per-genome outputs land in `<outdir>/<sample_id>/` and a
top-level `combined_results.csv` aggregates every genome's substrates with a
`source_genome` column. Batching saves on compute costs and increases statistical power.

## Input formats

`.gbff`, `.gbk`, `.gb` (GenBank), `.fasta`, `.fna`, `.fa` (nucleotide FASTA),
`.faa` (protein FASTA), `.gff` (GFF3).

GFF3 carries no sequence, so ssign needs a companion nucleotide FASTA with the same
stem next to it (e.g. `genome.gff` + `genome.fna`). The container launcher `ssign-run`
finds and binds that companion automatically.

## Choosing flags for your situation

`--use-input-annotations`

This flag tells ssign to trust the input GenBank's `/product` qualifiers and
skip Bakta re-annotation. Don't pass it unless you have a specific reason to.

Skipping annotation tools

Any annotation tool can be skipped on the command line:

```bash
ssign run input.gbff --outdir results \
    --skip-signalp \
    --skip-hhsuite \
    --skip-plmblast
```

The pipeline still runs everything else and records "skipped" for the matching
columns. `--skip-annotation` drops the whole functional-annotation phase (BLASTp,
HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam) at once.

Filtering or including secretion-system types

By default Flagellum, Tad, and the type-IV pili / uptake appendages (T4aP, T4bP,
MSH, ComM, Archaeal-T4P) are excluded from substrate identification. Override with
a different space-separated list:

```bash
# the default set, written out explicitly:
ssign run input.gbff --outdir results \
    --excluded-systems Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P
```

Tuning the proximity window

`--proximity-window N` (default `3`) controls how many genes on either side of
each detected SS component are treated as candidate substrates. ssign's approach
detects substrates that are **co-located** with their secretion system in the
chromosome; widening the window captures more candidates near each component but
does not help with substrates encoded far from the machinery. It also increases false positive rates.

```bash
ssign run input.gbff --outdir results --proximity-window 5
```

Per-type enrichment statistics (opt-in)

`--enrichment-stats` emits a per-SS-type circular-shift enrichment test (fold +
permutation p + BH q for DeepLocPro, DeepSecE, and SignalP) plus enrichment
figures. It forces whole-genome DeepLocPro + DeepSecE + SignalP, which adds compute requirements. Statistical power scales with genome
count: single-genome runs rarely reach q<0.05 even for real systems, so pool
across genomes. Method details are in
[`explanation/design_decisions.md`](../explanation/design_decisions.md).

## Running on HPC

For running ssign on more than a handful of genomes with an HPC account. The
**recommended** cluster path is the self-contained Apptainer/Singularity image
[`install.md`](install.md).

### PBS template

```bash
#!/bin/bash
#PBS -N ssign-ecoli
#PBS -l select=1:ncpus=16:mem=64gb:ngpus=1
#PBS -l walltime=06:00:00
#PBS -j oe

module load apptainer
cd $PBS_O_WORKDIR          # the ssign clone, so scripts/ssign-run resolves

scripts/ssign-run /path/to/genome.gbff $EPHEMERAL/ssign-out/${PBS_JOBID%%.*} \
    --tier extended \
    --sif $EPHEMERAL/ssign.sif \
    --db-root $EPHEMERAL/ssign-databases \
    --stage-image \
    --scratch $TMPDIR \
    --max-ram 60 \
    -- --cpu-per-genome 16
```

Submit with `qsub ssign-job.sh`.

Three of those flags matter more on a cluster than on a laptop:

- `--max-ram` is **required** on PBS. The scheduler hides the job's allocation from
  the container, so without it ssign sizes tool memory to the whole node and can be
  killed for overrunning. Pass a little under the `mem=` request. SLURM is read
  automatically from `SLURM_MEM_PER_NODE`.
- `--stage-image` copies the `.sif` to node-local disk first, turning a slow
  network-filesystem startup into one sequential copy.
- `--scratch` points Bakta's working files at real disk. The container's default
  `/tmp` is a 64 MiB tmpfs, which Bakta overflows.

Drop `ngpus=1` and add `--cpu` to run without a GPU; the predictors work but are much
slower. Anything after `--` is passed straight through to `ssign run`.

For the native/pip route instead, activate your environment and call `ssign run`
directly, exporting the database paths the fetcher printed. See
[`install-secondary.md`](install-secondary.md).

### Make `--outdir` unique per job

The moment you submit multiple jobs that could land concurrently (array jobs,
overnight bursts, parameter sweeps), each needs its own output directory. If two
jobs share the same `--outdir`, both ssign processes read and write the same
intermediate files; the pipeline emits no error but the final results table
silently mixes data from both runs. Use the scheduler's job id in the path.

```bash
# PBS: strip everything after the first dot in case PBS reports JOBID.server-name
--outdir $EPHEMERAL/ssign-out/${PBS_JOBID%%.*}
```

## Troubleshooting

If your error is not here,
[open a GitHub issue](https://github.com/billerbeck-lab/ssign/issues)

### Out of memory

DeepSecE wraps a large ESM-1b model (about 7.3 GB resident). On a small node or
laptop this shows up as `MemoryError` during checkpoint loading. Free up memory
or use a machine with at least 16 GB RAM.

`RuntimeError: CUDA out of memory` during pLM-BLAST means the GPU lacks VRAM for
the ProtT5 encoder.

### EggNOG is slow on a node with < 64 GB RAM

No action needed. EggNOG's in-RAM database load (`--dbmem`) auto-disables when the
job's RAM share is under 50 GB, so it mmaps the database from local SSD instead,
which is fine for the small substrate set ssign annotates. Force it either way with
`-- --eggnog-dbmem` / `-- --no-eggnog-dbmem`.

### Bakta "No space left on device" (container scratch too small)

Inside the Apptainer/Singularity image, `/tmp` is often a small tmpfs that Bakta's
working files overflow. Point scratch at real disk:

```bash
ssign run input.gbff --outdir results --scratch-dir /path/to/big/scratch
```

`--scratch-dir` keeps `$TMPDIR` when it has adequate free space and otherwise
falls back to a directory under `--outdir`. On a cluster you can instead bind a
real-disk path over `/tmp` when launching the container.

### SignalP or DeepLocPro not found (DTU local install)

If a local install is not on `PATH` (and you have not passed `--signalp-mode
remote` / `--deeplocpro-mode remote`), ssign stops with install instructions. Point it at the install directory with
`--signalp-path` / `--deeplocpro-path`, or set `$SSIGN_SIGNALP_PATH` /
`$SSIGN_DEEPLOCPRO_PATH`.

SignalP 6.0 has to live in its own environment, so this error usually means the
install step was skipped or the environment is not where ssign looks. Setup steps
and why the environment must be separate are in
[`install-secondary.md`](install-secondary.md#signalp-60).

### `hhblits binary not found`

You enabled HH-suite (`--no-skip-hhsuite`) but it is not installed:

```bash
conda install -c bioconda hhsuite
```

If `which hhblits` finds the binary but the error persists, it is not on `PATH`
inside your job script. Add the conda env activation to the script.
