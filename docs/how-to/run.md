# Running ssign

ssign is a command-line tool. This page covers the basic invocation, how to
pick flags for common situations, how to run on an HPC cluster, and how to fix
the problems people actually hit.

For install steps (both the container and native routes) see
[`install.md`](install.md). For the full flag list see
[`reference/cli.md`](../reference/cli.md); for environment variables see
[`reference/env_vars.md`](../reference/env_vars.md).

## Basic invocation

```bash
ssign run input.gbff --outdir results
```

That runs the whole pipeline on one genome and writes to `results/`:

- `<sample-id>_results.csv` secreted proteins + secretion systems (three sections)
- `<sample-id>_summary.txt` readable summary (and enrichment stats if enabled)
- `figures/<sample-id>/` figures `01`-`06`

The sample id defaults to the input filename's stem; override it with
`--sample-id` (single-genome runs only). See
[`reference/output_files.md`](../reference/output_files.md) for the full output
layout.

Two more subcommands help around a run:

```bash
ssign doctor                             # verify the install before running
ssign fetch-databases --tier extended    # download reference databases
```

Running `ssign` with no subcommand just prints a usage hint. `pip install ssign`
lands at the v1.0.0 release; until then install from source (git clone +
`pip install -e .`), see [`install.md`](install.md).

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
`source_genome` column. Batching pays heavy startup costs (InterProScan JVM,
EggNOG DB load, pLM-BLAST embeddings) once per batch instead of once per genome.
The per-genome sample id is derived from each filename's stem, so input
basenames must be distinct; `--sample-id` is rejected for multi-genome runs.
Turn the combined table off with `--no-combined-summary`.

For hundreds of genomes, prefer an HPC job array (one job per genome) over a
single long batched job; see the HPC section below.

## Input formats

| Format        | Extensions              | Notes |
|---------------|-------------------------|-------|
| GenBank       | `.gbff`, `.gbk`, `.gb`  | Recommended. Carries protein sequences and gene order; ssign re-annotates with Bakta by default unless `--use-input-annotations` is set. |
| FASTA contigs | `.fasta`, `.fna`, `.fa` | Bakta (or pyrodigal as fallback) predicts ORFs. |
| GFF3          | `.gff`                  | Works when a same-stem nucleotide FASTA sits alongside it (e.g. `x.gff` + `x.fna`). |

## Choosing flags for your situation

The defaults are tuned for the common case (Bakta re-annotation, local DTU
tools, proximity window 3). The recipes below cover the non-obvious choices.

### `--use-input-annotations`: leave it off for most inputs

This flag tells ssign to trust the input GenBank's `/product` qualifiers and
skip Bakta re-annotation. Don't pass it unless you have a specific reason to:
most public GenBank records (including NCBI's K-12 reference) lack the
Pfam-domain annotations MacSyFinder needs to detect secretion systems, so
substrate recall drops sharply.

On the K-12 reference genome:

| Configuration                 | Substrates detected |
|-------------------------------|---------------------|
| Default (Bakta re-annotation) | 17 |
| `--use-input-annotations`     | 8  |

The same drop applies to GFF3-derived inputs (NCBI's GFF3 also lacks Pfam IDs in
the form ssign needs). Only use `--use-input-annotations` when you have manually
curated a GenBank record with Pfam-tagged CDS qualifiers and want them preserved.

### Skipping annotation tools

Any annotation tool can be skipped on the command line:

```bash
ssign run input.gbff --outdir results \
    --skip-signalp \
    --skip-hhsuite \
    --skip-plmblast
```

The pipeline still runs everything else and records "skipped" for the matching
columns. `--skip-annotation` drops the whole functional-annotation phase (BLASTp,
HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam) at once; predictions and
substrate calls still run. A per-tool `--no-skip-<tool>` (e.g. `--no-skip-eggnog`)
re-enables one tool on top of that master switch.

### Filtering or including secretion-system types

By default Flagellum, Tad, and the type-IV pili / uptake appendages (T4aP, T4bP,
MSH, ComM, Archaeal-T4P) are excluded from substrate identification. T3SS is
**not** excluded (it is detected and substrate-called by default). Override with
a different space-separated list:

```bash
# the default set, written out explicitly:
ssign run input.gbff --outdir results \
    --excluded-systems Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P
```

To also exclude T3SS (restore the pre-2026 behaviour), append it to the list.
T3SS stays on by default, but DeepSecE is never trusted for it (it produces a
high false-positive rate, mostly flagellar misclassifications); T3SS calls rely
on TXSScan detection + DeepLocPro + proximity instead. Full rationale, including
the 74-genome benchmark, is in
[`explanation/design_decisions.md`](../explanation/design_decisions.md).

### Tuning the proximity window

`--proximity-window N` (default `3`) controls how many genes on either side of
each detected SS component are treated as candidate substrates. ssign's approach
detects substrates that are **co-located** with their secretion system in the
chromosome; widening the window recovers more candidates near each component but
does not help with substrates encoded far from the machinery.

```bash
ssign run input.gbff --outdir results --proximity-window 5
```

On the K-12 reference: window 3 gives 17 substrates, window 5 gives 20. Compute
cost grows roughly linearly with the window (the neighborhood feeds DeepLocPro,
DeepSecE, and SignalP), and the false-positive count grows with it too:
substrates further from a component are more often co-located by chance than by
secretion biology.

### SignalP and DeepLocPro: local vs webserver

ssign is offline-first: the default path is a local install of SignalP and
DeepLocPro. If the binaries are on `PATH`, no flags are needed. If they live
elsewhere, point ssign at them:

```bash
ssign run input.gbff --outdir results \
    --signalp-path /path/to/signalp6/bin \
    --deeplocpro-path /path/to/deeplocpro
```

If no local install is found, ssign stops with install instructions; it does
**not** silently upload to DTU. Without a DTU licence you can opt into the
webserver fallback (no licence needed, internet required):

```bash
ssign run input.gbff --outdir results \
    --signalp-mode remote \
    --deeplocpro-mode remote
```

The webserver is a convenience for users without a licence and for first-time
trials. It depends on DTU continuing to host the service, which we cannot
guarantee long-term; for paper-grade or long-running cohorts, install SignalP
and DeepLocPro locally. See [`install.md`](install.md) for the local-install
steps.

### Per-type enrichment statistics (opt-in)

`--enrichment-stats` emits a per-SS-type circular-shift enrichment test (fold +
permutation p + BH q for DeepLocPro, DeepSecE, and SignalP) plus enrichment
figures. It forces whole-genome DeepLocPro + DeepSecE + SignalP (the rotation
null needs every gene's positivity in gene order), which adds roughly 13 min per
genome, so size HPC walltime accordingly. Statistical power scales with genome
count: single-genome runs rarely reach q<0.05 even for real systems, so pool
across genomes. Method details are in
[`explanation/design_decisions.md`](../explanation/design_decisions.md).

## Running on HPC

For running ssign on more than a handful of genomes with an HPC account. The
**recommended** cluster path is the self-contained Apptainer/Singularity image
([`install.md`](install.md), see its container "On HPC" section). This
section documents the **native pip** alternative: pip-install in a venv, fetch
databases to scratch, submit a SLURM or PBS job that calls `ssign run`. Use it
when the container does not fit your cluster.

### What is different on a cluster

- You have no `sudo`. System binaries (BLAST+, HH-suite) come from the cluster's
  `module` system or from conda in your home directory.
- The login node is for setup, not for running pipelines. Real work runs as a
  job on a compute node, submitted via SLURM (`sbatch`) or PBS (`qsub`).
- Filesystems are split. Code in `$HOME` (small quota), databases in `$SCRATCH`
  or `$WORK` (large, fast). A 50 GB EggNOG download in `$HOME` will likely hit
  your quota.
- GPUs are requested explicitly in the job script. Forget that flag and DeepSecE
  and pLM-BLAST both fall back to CPU (10-100x slower).
- Jobs have walltime limits. A bigger cohort or the full tier may need to be
  split across multiple jobs.

Do not run from a JupyterHub or interactive session. They look like Linux shells
but are typically cgroup-throttled to about one CPU's worth of compute. `nproc`
reports the true quota; if it returns 1 or 2 you are not on a compute node. On a
1-CPU session DeepSecE on a single genome takes about 90 min instead of seconds
on GPU or about 10 min on a 16-core CPU node. Always submit a proper job for the
actual run.

### 1. Set up the environment (login node)

```bash
# Load Python from the cluster's module system. The exact name varies;
# common patterns: python/3.11, anaconda3, miniforge.
module avail python                     # see what's available
module load python/3.11

python -m venv ~/.ssign-env
source ~/.ssign-env/bin/activate
which python                            # verify the venv is active

pip install ssign[extended]             # or: pip install ssign  (base only)
ssign --version                         # verify install
```

If the cluster firewalls outbound HTTPS, install from a pre-staged wheel or via a
configured pip proxy. Most institutional clusters allow PyPI; ask your admin if
`pip install` hangs.

### 2. Stage databases on scratch

```bash
# Pick a directory on scratch (NOT home). The exact env var depends on your
# cluster: $SCRATCH on some, $TMPDIR or /scratch/$USER on others. Check the docs.
export SSIGN_DBS=$SCRATCH/ssign-databases
mkdir -p $SSIGN_DBS && cd $SSIGN_DBS

# Fetch the tier matching your pip install.
ssign fetch-databases --tier extended --target $SSIGN_DBS
```

Disk sizes per tier are in [`install.md`](install.md) (about 2.5 / 100 / 500 GB
for base / extended / full). Rough fetch time on a fast cluster link: base 15-30 min,
extended 1-3 h, full 6-12 h. The fetcher prints a list of database-path exports
at the end; copy them into your shell rc file or job script so the paths are set
when ssign runs. See [`reference/env_vars.md`](../reference/env_vars.md) for the
variables ssign reads.

### 3. Submit a job

#### SLURM template

```bash
#!/bin/bash
#SBATCH --job-name=ssign-ecoli
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.log
# For pLM-BLAST add: #SBATCH --gres=gpu:1

module load python/3.11
source ~/.ssign-env/bin/activate

# Database paths (printed by the fetcher in step 2)
export SSIGN_HHSUITE_PFAM=$SCRATCH/ssign-databases/hhsuite/pfam
export SSIGN_HHSUITE_PDB70=$SCRATCH/ssign-databases/hhsuite/pdb70
export SSIGN_HHSUITE_UNICLUST=$SCRATCH/ssign-databases/hhsuite/uniref30

ssign run /path/to/genome.gbff \
    --outdir $SCRATCH/ssign-out/$SLURM_JOB_ID \
    --bakta-db $SCRATCH/ssign-databases/bakta/db-light \
    --eggnog-db $SCRATCH/ssign-databases/eggnog \
    --interproscan-db $SCRATCH/ssign-databases/interproscan/interproscan-5.77-108.0 \
    --plmblast-db $SCRATCH/ssign-databases/plm_blast/ECOD30 \
    --cpu-per-genome 16
```

Submit with `sbatch ssign-job.sh`. Monitor with `squeue -u $USER` and tail the
log with `tail -f ssign-ecoli-*.log`.

#### PBS template

```bash
#!/bin/bash
#PBS -N ssign-ecoli
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=06:00:00
#PBS -j oe
# For pLM-BLAST add: #PBS -l select=1:ncpus=16:mem=64gb:ngpus=1

module load anaconda3
source ~/.ssign-env/bin/activate

# Database paths (printed by the fetcher in step 2)
export SSIGN_HHSUITE_PFAM=$EPHEMERAL/ssign-databases/hhsuite/pfam
# ... other database exports ...

ssign run /path/to/genome.gbff \
    --outdir $EPHEMERAL/ssign-out/${PBS_JOBID%%.*} \
    --bakta-db $EPHEMERAL/ssign-databases/bakta/db-light \
    --cpu-per-genome 16
```

Submit with `qsub ssign-job.sh`.

#### Make `--outdir` unique per job

The moment you submit multiple jobs that could land concurrently (array jobs,
overnight bursts, parameter sweeps), each needs its own output directory. If two
jobs share the same `--outdir`, both ssign processes read and write the same
intermediate files; the pipeline emits no error but the final results table
silently mixes data from both runs. Use the scheduler's job id in the path, as
both templates above do:

```bash
# SLURM
--outdir $SCRATCH/ssign-out/$SLURM_JOB_ID

# PBS: strip everything after the first dot in case PBS reports JOBID.server-name
--outdir $EPHEMERAL/ssign-out/${PBS_JOBID%%.*}
```

For job arrays, also include the array index (`$SLURM_ARRAY_TASK_ID` /
`$PBS_ARRAY_INDEX`) so per-task outputs do not overlap.

### GPU access

DeepSecE and pLM-BLAST both benefit from a GPU:

- **DeepSecE:** seconds on GPU, about 10 min on a 16-core CPU node, about 90 min
  on a 1-CPU session. Auto-detected; no flag.
- **pLM-BLAST:** the ProtT5 embedding is about 100x slower on CPU than GPU.

Request a GPU in the job script (`#SBATCH --gres=gpu:1` for SLURM, add `:ngpus=1`
to the PBS `select` line). After the job starts, verify it is visible before
ssign launches:

```bash
nvidia-smi
python -c "import torch; print(torch.cuda.is_available())"
```

If `nvidia-smi` is missing or reports no GPU, ssign falls back to CPU silently.
If `nvidia-smi` works but `torch.cuda.is_available()` prints `False`, PyTorch was
installed CPU-only; reinstall the CUDA build for your driver version.

### Installing DTU tools on HPC

SignalP 6.0 and DeepLocPro both come from DTU but through different channels:
SignalP via DTU HealthTech's portal, DeepLocPro via the maintainer's GitHub
([Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro)). Two
HPC-specific gotchas:

1. SignalP's portal URL is an Apache directory listing, not a file. Append the
   filename `signalp-6.0i.fast.tar.gz` to the URL when you `wget`, otherwise you
   get a 1-2 KB HTML index page instead of the ~1.5 GB tarball. DeepLocPro avoids
   this because it is a `git clone`.
2. Do not `mamba init` your shell on HPC. Cluster shell profiles are often shared
   or sandboxed, and `init` writes to your rc file. Use the conda env's absolute
   paths instead (`PYBIN=~/.conda/envs/<env>/bin`, then `$PYBIN/pip ...`).

Step-by-step instructions are in [`install.md`](install.md). Both tools go into
isolated conda envs; ssign points at them via `--signalp-path` /
`--deeplocpro-path` (the directory holding the console script):

```bash
ssign run input.gbff --outdir results \
    --signalp-mode local --signalp-path ~/.conda/envs/signalp6/bin \
    --deeplocpro-mode local --deeplocpro-path ~/.conda/envs/deeplocpro/bin
```

### Walltime considerations

Approximate per-genome wall time at extended tier on a 16-core compute node:

| Step | Time |
|---|---|
| Bakta re-annotation | 10-20 min |
| MacSyFinder | 5-10 min |
| DeepLocPro (local) | 2-5 min, GPU-accelerated |
| SignalP 6 (local) | 2-5 min, CPU-only |
| DeepLocPro + SignalP (DTU webserver, opt-in) | 5-15 min, network-bound |
| DeepSecE (GPU / 16-core CPU / 1-core throttled) | seconds / ~10 min / ~90 min |
| HH-suite (Pfam + PDB70) | 10-30 min |
| InterProScan | 5-15 min |
| BLASTp Swiss-Prot | a few min (full-tier default; NR opt-in is far slower) |
| EggNOG-mapper | 5-15 min |
| pLM-BLAST (GPU) | 5-15 min |
| Reporting | 1-2 min |

A typical extended-tier run on one genome fits in 4-6 h. The full tier
(`--tier full`) adds BLASTp-vs-Swiss-Prot and HH-suite-vs-UniRef30 on the
substrate set. HH-suite runs one MSA per substrate, so on a large pooled cohort
(hundreds to a thousand-plus substrates) it is the tail; Swiss-Prot BLASTp
finishes in minutes. Size `--walltime` to 24 h+ for a full-tier panel, and
validate on a 2-4 genome smoke test first to measure the real per-substrate
HH-suite time before committing a long job.

Full tier also needs its databases present: Swiss-Prot at `<db_root>/blast_swissprot/`
(sets `SSIGN_BLAST_SWISSPROT`) and UniRef30 at `<db_root>/hhsuite/uniref30/`
(sets `SSIGN_HHSUITE_UNICLUST`). With both present, ssign resolves BLASTp's `-db`
and HH-suite's UniRef30 automatically. Swiss-Prot is the default because full NR
(~800 GB) is impractically slow on a real substrate set; to use NR, fetch it and
pass `--blastp-db <nr-dir>/nr`.

### Cohort runs as a job array

```bash
#!/bin/bash
#SBATCH --job-name=ssign-cohort
#SBATCH --array=1-100
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00

module load python/3.11
source ~/.ssign-env/bin/activate

GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $HOME/genomes.txt)

ssign run $GENOME \
    --outdir $SCRATCH/ssign-out/$(basename $GENOME .gbff) \
    --bakta-db $SCRATCH/ssign-databases/bakta/db-light \
    --cpu-per-genome 16
```

`genomes.txt` is a flat file with one input path per line. SLURM schedules them
in parallel up to your queue's array limit.

## Troubleshooting

Only the problems that actually come up in practice. If your error is not here,
[open a GitHub issue](https://github.com/billerbeck-lab/ssign/issues) with the
full log, the exact command, and your environment (`pip freeze`,
`python --version`, OS).

### Out of memory

DeepSecE wraps a large ESM-1b model (about 7.3 GB resident). On a small node or
laptop this shows up as `MemoryError` during checkpoint loading. Free up memory
or use a machine with at least 12 GB RAM; on a 16 GB laptop, close other heavy
applications first. On HPC, size `--mem` to at least 64 GB for an extended-tier
run and set `SSIGN_MAX_RAM_GB` if the scheduler does not report the real cap.

`RuntimeError: CUDA out of memory` during pLM-BLAST means the GPU lacks VRAM for
the ProtT5 encoder (wants about 16 GB). Options: run pLM-BLAST on CPU (slow but
works), skip it with `--skip-plmblast` (or `--skip-annotation` for all annotation
tools), or use a GPU node with more VRAM.

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
remote` / `--deeplocpro-mode remote`), ssign stops with install instructions
rather than uploading to DTU. Point it at the install directory with
`--signalp-path` / `--deeplocpro-path`, or set `$SSIGN_SIGNALP_PATH` /
`$SSIGN_DEEPLOCPRO_PATH`.

SignalP 6.0 pins PyTorch < 2.0 while every other ssign tool needs PyTorch 2.x, so
installing it into the main environment breaks DeepSecE. Install SignalP into its
own conda or venv environment, then point ssign at it via `--signalp-path`.
Detailed steps are in [`install.md`](install.md).

### `hhblits binary not found`

You enabled HH-suite (`--no-skip-hhsuite`) but it is not installed:

```bash
conda install -c bioconda hhsuite
```

If `which hhblits` finds the binary but the error persists, it is not on `PATH`
inside your job script. Add the conda env activation to the script.
