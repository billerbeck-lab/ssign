# Running ssign on an HPC cluster

For users with an HPC account who want to run ssign on more than a handful
of genomes. Verified path: pip install in a venv, fetch databases to
scratch, submit a SLURM or PBS job that calls `ssign run`.

A Singularity container path is planned with v1.0.0's Docker image and
will be added to this page when the image lands; for now use the pip
path.

## Overview: what is different on a cluster

Compared to a laptop install:

- You have no `sudo`. System binaries (BLAST+, HH-suite) come from the
  cluster's `module` system or from conda in your home directory.
- The login node is for setup, not for running pipelines. Real work runs
  as a job on a compute node, submitted via SLURM (`sbatch`) or PBS
  (`qsub`).
- Filesystems are split. Code in `$HOME` (small quota), databases in
  `$SCRATCH` or `$WORK` (large, fast). A 50 GB EggNOG download in
  `$HOME` will likely hit your quota.
- GPUs are requested explicitly in the job script. Forget that flag and
  DeepSecE, PLM-Effector, and pLM-BLAST all fall back to CPU (10-100x
  slower).
- Jobs have walltime limits. A bigger cohort or extended tier may need
  to be split across multiple jobs.

> **Don't run from a JupyterHub or interactive session.** They look like
> Linux shells but are typically cgroup-throttled to ~1 CPU's worth of
> compute (even when `Cpus_allowed_list` says you can touch 4 cores).
> `nproc` reports the true quota — if it returns 1 or 2 you are not on a
> compute node. On a 1-CPU session DeepSecE on a single genome takes
> ~90 min instead of seconds on GPU or ~10 min on a 16-core CPU node.
> Always submit a proper SLURM / PBS job for the actual run.

## 1. Set up the environment (login node)

```bash
# Load Python from the cluster's module system. The exact name varies;
# common patterns: python/3.11, anaconda3, miniforge.
module avail python                     # see what's available
module load python/3.11

# Create a venv in your home directory. Activate it.
python -m venv ~/.ssign-env
source ~/.ssign-env/bin/activate
which python                            # verify the venv is active

# Install ssign at the tier you want.
pip install ssign[extended]             # or: pip install ssign  (base only)
ssign --version                         # verify install
```

If the cluster firewalls outbound HTTPS, install from a pre-staged wheel
or via a configured pip proxy. Most institutional clusters allow PyPI;
ask your admin if `pip install` hangs.

## 2. Stage databases on scratch

```bash
# Pick a directory on scratch (NOT home). The exact env var depends on
# your cluster: $SCRATCH on Imperial CX3, $TMPDIR on some clusters,
# /scratch/$USER on others. Check the cluster docs.
export SSIGN_DBS=$SCRATCH/ssign-databases
mkdir -p $SSIGN_DBS && cd $SSIGN_DBS

# Run the fetcher at the tier matching your pip install.
bash $(python -c "import ssign_app, os; print(os.path.dirname(ssign_app.__file__))")/../../scripts/fetch_databases.sh \
    --tier extended \
    --target $SSIGN_DBS
```

Disk usage by tier (approximate):

| Tier | Disk | Time to fetch |
|---|---|---|
| `base` | 17 GB | 15-30 min |
| `extended` | 130 GB | 1-3 h |
| `full` | 630 GB | 6-12 h |

The fetcher prints a list of `export SSIGN_*` lines at the end; copy them
into your shell rc file or your job script so the database paths are set
when ssign runs.

## 3. Submit a job

### SLURM template

```bash
#!/bin/bash
#SBATCH --job-name=ssign-ecoli
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.log
# For PLM-Effector or pLM-BLAST add: #SBATCH --gres=gpu:1

module load python/3.11
source ~/.ssign-env/bin/activate

# Database paths (printed by the fetcher in step 2)
export SSIGN_HHSUITE_PFAM=$SCRATCH/ssign-databases/pfam
export SSIGN_HHSUITE_PDB70=$SCRATCH/ssign-databases/pdb70_from_mmcif_2026-02-20
export SSIGN_HHSUITE_UNICLUST=$SCRATCH/ssign-databases/UniRef30_2023_02

ssign run /path/to/genome.gbff \
    --outdir $SCRATCH/ssign-out \
    --bakta-db $SCRATCH/ssign-databases/bakta_light \
    --eggnog-db $SCRATCH/ssign-databases/eggnog \
    --interproscan-db $SCRATCH/ssign-databases/interproscan-5.77-108.0 \
    --plmblast-db $SCRATCH/ssign-databases/ECOD70 \
    --cpu-per-genome 16
```

Submit with `sbatch ssign-job.sh`. Monitor with `squeue -u $USER` and tail
the log file with `tail -f ssign-ecoli-*.log`.

### PBS template

```bash
#!/bin/bash
#PBS -N ssign-ecoli
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=06:00:00
#PBS -j oe
# For PLM-Effector or pLM-BLAST add: #PBS -l select=1:ncpus=16:mem=64gb:ngpus=1

module load anaconda3
source ~/.ssign-env/bin/activate

# Database paths (printed by the fetcher in step 2)
export SSIGN_HHSUITE_PFAM=$EPHEMERAL/ssign-databases/pfam
# ... other SSIGN_* exports ...

ssign run /path/to/genome.gbff \
    --outdir $EPHEMERAL/ssign-out \
    --bakta-db $EPHEMERAL/ssign-databases/bakta_light \
    --cpu-per-genome 16
```

Submit with `qsub ssign-job.sh`.

## 4. GPU access

DeepSecE, PLM-Effector, and pLM-BLAST all benefit from a GPU:

- **DeepSecE:** seconds on GPU, ~10 min on a 16-core CPU node, ~90 min on
  a 1-CPU session. Auto-detected; no flag.
- **PLM-Effector:** seconds on GPU, ~1-2 min per protein on CPU. Default
  is `cuda`; skipped automatically if no GPU is visible.
- **pLM-BLAST:** ProtT5 embedding is ~100x slower on CPU than GPU.

Request a GPU in your job script:

```bash
# SLURM:
#SBATCH --gres=gpu:1
# PBS:
#PBS -l select=1:ncpus=16:mem=64gb:ngpus=1
```

After the job starts, verify the GPU is visible inside the job before
ssign launches:

```bash
nvidia-smi    # should print one or more GPUs
python -c "import torch; print(torch.cuda.is_available())"
```

If `nvidia-smi` is missing or reports no GPU, ssign falls back to CPU
silently for DeepSecE and skips PLM-Effector entirely.

## 5. Walltime considerations

Approximate per-genome wall time at extended tier on a 16-core compute
node:

| Step | Time |
|---|---|
| Bakta re-annotation | 10-20 min |
| MacSyFinder | 5-10 min |
| DeepLocPro + SignalP (local) | 2-5 min, GPU-accelerated |
| DeepLocPro + SignalP (DTU webserver fallback) | 5-15 min, network-bound |
| DeepSecE (GPU / 16-core CPU / 1-core throttled) | seconds / ~10 min / ~90 min |
| HH-suite (Pfam + PDB70) | 10-30 min |
| InterProScan | 5-15 min |
| BLASTp NR | 30-90 min depending on n_proteins and NR vintage |
| EggNOG-mapper | 5-15 min |
| pLM-BLAST (GPU) | 5-15 min |
| Reporting | 1-2 min |

A typical extended-tier run on a single bacterial genome fits in 4-6 h.
Cohorts spanning hundreds of genomes need a job array (one job per
genome) rather than a single long job.

## 6. Cohort runs as a job array

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
    --bakta-db $SCRATCH/ssign-databases/bakta_light \
    --cpu-per-genome 16
```

`genomes.txt` is a flat file with one input path per line. SLURM will
schedule them in parallel up to your queue's array limit.

## 7. When the Docker image lands

v1.0.0 will publish a SHA-pinned Docker image to the GitHub Container
Registry that bundles every Python dependency. On a cluster with
Singularity (the standard HPC container runtime):

```bash
# Pre-v1.0.0: this section is a stub. The path lights up when the image
# is built and tagged on Docker Hub / GHCR.
singularity pull docker://billerbeck-lab/ssign:1.0.0
singularity run --bind $SCRATCH/ssign-databases:/data ssign_1.0.0.sif \
    run /path/to/genome.gbff --outdir /scratch/out
```

The Singularity path will be filled in once the image is published.
