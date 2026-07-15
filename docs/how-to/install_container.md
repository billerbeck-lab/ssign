# Installing ssign (container, the recommended way)

The Singularity/Apptainer image is the recommended way to run ssign on Linux and
HPC. It is **self-contained and reproducible**: the entire free toolchain and all
model weights are baked in and version-locked (digest-pinned base, `uv.lock`,
pinned tool commits, baked ESM/ProtT5 weights), so a run reproduces over time. The
only things not included are the DTU-licensed **SignalP 6** and the **reference
databases**.

This guide covers the **extended** tier (secretion-system detection + secreted
protein identification + functional/structural annotation). Base and full differ
only in which databases you fetch, see [Tiers](#tiers).

The image is for **non-commercial research use** (it bundles DeepLocPro,
CC BY-NC-SA, and EggNOG-mapper, AGPL).

> Building or publishing the image (maintainers) is documented separately in
> [`containers/README.md`](../../containers/README.md).

## Hardware requirements (extended)

| Resource | Recommended | Notes |
| --- | --- | --- |
| GPU | CUDA NVIDIA, ~16-24 GB VRAM (RTX 6000 / A40 / etc.) | DeepLocPro, DeepSecE, and pLM-BLAST embeddings use it. CPU-only works but the predictors are much slower. |
| RAM | 32 GB works; **>= 64 GB** recommended | 64 GB only buys the EggNOG in-RAM database speedup (`--dbmem`), which auto-disables below a 50 GB job share, so smaller nodes run fine (just slower EggNOG). See [Troubleshooting](#troubleshooting). |
| CPU | 8+ cores (the more the better) | Bakta, HMMER, and the parallel annotation group all scale with cores. |
| Disk | image ~20 GB + databases ~100 GB + **fast node-local scratch ~100 GB** | The scratch holds the staged image, the cached EggNOG DB, and InterProScan temp; put it on local SSD, not a network mount. |

## Prerequisites

- **`apptainer`** (or `singularity`). On HPC it is almost always preinstalled as a
  module: `module avail apptainer` (or `singularity`), then `module load apptainer`,
  then `apptainer --version` to confirm. If your cluster genuinely lacks it, email HPC
  support to add it (it is a standard rootless cluster tool; users normally should not
  install a system container runtime themselves). User-space fallback on clusters that
  allow user namespaces: `conda install -c conda-forge apptainer`.
- A licensed **SignalP 6** tarball. Register (free for academics) at
  <https://services.healthtech.dtu.dk/services/SignalP-6.0/> and download once.

## Install (extended)

```bash
# 0a. Load apptainer (on HPC it is a module) and confirm it runs
module load apptainer && apptainer --version

# 0b. Get the launcher scripts (they wrap apptainer, so pip does not put them on
#     PATH; a clone gives you scripts/ and the docs):
git clone https://github.com/billerbeck-lab/ssign && cd ssign

# 0c. Get the image (.sif) and point ssign-run at it. At release:
#       apptainer pull oras://ghcr.io/billerbeck-lab/ssign:1.0.0
#     until then, build it once (see containers/README.md).
export SSIGN_SIF=$HOME/ssign.sif

# 1. Reference databases, fetched FROM the image (no host tools needed; ~100 GB, once)
apptainer run --writable-tmpfs --containall -B $HOME/ssign-databases:$HOME/ssign-databases \
  "$SSIGN_SIF" fetch-databases --tier extended --target $HOME/ssign-databases

# 2. SignalP 6, the only licence-gated tool (register + download the tarball first, then:)
bash scripts/ssign-setup-dtu ~/signalp-6.0i.fast.tar.gz --signalp-only

# 3. Run (one line per genome)
bash scripts/ssign-run $HOME/mygenome.gbff $HOME/ssign_out --tier extended \
  --sif "$SSIGN_SIF" --db-root $HOME/ssign-databases --max-ram 60
```

Everything else (Bakta + its toolchain, EggNOG-mapper, DeepLocPro + ESM2, DeepSecE
+ ESM-1b, pLM-BLAST + ProtT5, InterProScan's Java, BLAST+ 2.17) is inside the
image, so a first run touches the network only for the one database fetch.

- `--signalp-env` is auto-detected from your conda envs after step 2; the flag only
  overrides it. No SignalP licence? Drop step 2 and add `-- --signalp-mode remote`
  (uses the DTU webserver, which uploads your sequences).
- `--max-ram <GB>` = your job's RAM. Required on schedulers that hide the allocation
  from the container (PBS), else ssign sizes tool memory to the whole node. On SLURM
  it is read from `SLURM_MEM_PER_NODE` automatically.
- `--stage-image` copies the `.sif` to fast node-local disk before the run, turning a
  slow network-filesystem startup into one sequential copy. Recommended on HPC.

## What's baked vs what you provide

| Component | In the image | You provide |
| --- | :---: | --- |
| CUDA runtime, Python 3.10, `ssign[extended]` deps (pinned `uv.lock`) | yes | |
| MacSyFinder + TXSScan 1.1.4 profiles, `hmmsearch` (pyhmmer) | yes | |
| Bakta 1.12.0 + toolchain (BLAST+ 2.17, real HMMER, DIAMOND, tRNAscan-SE, aragorn, INFERNAL, PilerCR, AMRFinderPlus) | yes | |
| EggNOG-mapper (`emapper.py`), Java 11, rsync | yes | |
| DeepLocPro + ESM2 backbone | yes | |
| DeepSecE + its ESM-1b checkpoint | yes | |
| pLM-BLAST code + ProtT5 encoder | yes | |
| **SignalP 6** (DTU academic licence) | no | conda env via `ssign-setup-dtu`, or `--signalp-mode remote` |
| **Reference databases** (Bakta DB, EggNOG DB, InterProScan install+DBs, ECOD; full adds HH-suite + Swiss-Prot) | no | `fetch-databases --tier <t>` |

InterProScan is the one in-between: its ~24 GB software+member-DB directory is part
of the database set (pulled by `fetch-databases`) and that same directory also
supplies `interproscan.sh`. The Java to run it is baked.

## How a run works

- **"Mounting" (or "binding") a folder** just means "let the box see this folder on
  the cluster". `ssign-run` mounts your genome, the database folder, and the SignalP
  folder automatically.
- **A conda environment** is a self-contained folder holding one tool plus its
  dependencies (kept separate so tools don't clash). The only one you install is
  SignalP 6 (via `ssign-setup-dtu`), and `ssign-run` finds it for you. Everything
  else is baked into the image.
- **The GPU** is turned on automatically (`--nv`) when one is present. Nothing to do.
- **Where does the genome go?** Anywhere you can read it on the cluster (home,
  scratch, wherever). You just give `ssign-run` the path to the file. There is no
  special location, and you do not copy it into the container.

## On HPC

Build the image off-cluster (an HPC compute node usually can't, see
[`containers/README.md`](../../containers/README.md)), copy the `.sif` up, then
fetch the DBs + install SignalP once and run per job:

```bash
scp ssign.sif you@cluster:$EPHEMERAL/ssign.sif
```

On PBS/CX3 the ready wrapper `scripts/cx3/run_container_extended.pbs` stages the
image to node-local disk, pins the RAM budget to the job's `mem`, adds `--nv`, and
mounts only SignalP + the DB root:

```bash
qsub -v GENOME=$HOME/g.gbff,SIF=$EPHEMERAL/ssign.sif scripts/cx3/run_container_extended.pbs
```

See [`run_on_hpc.md`](run_on_hpc.md) for scheduler templates, walltime, and job arrays.

## Verify

```bash
apptainer run --writable-tmpfs --containall "$SSIGN_SIF" doctor --tier extended
```

## Tiers

Same image every tier; the tier is the database set you fetch plus `--tier`.

| Tier | `fetch-databases` pulls | Approx DB size |
| --- | --- | --- |
| `base` | NCBI taxdump + Bakta light | ~4 GB |
| `extended` | + EggNOG + InterProScan + ECOD30 (pLM-BLAST) | ~100 GB |
| `full` | + Bakta full + HH-suite (Pfam/PDB70/UniRef30) + BLASTp Swiss-Prot | ~500 GB |

**macOS / base tier without a container:** every base-tier tool is pip-installable, so
macOS base needs no image, `python3 -m pip install -c containers/requirements-base.lock.txt ssign`.
SignalP/DeepLocPro still come from DTU (install locally, or `--*-mode remote`).
Extended/full stay on Linux/HPC via the `.sif`.

## Troubleshooting

- **< 64 GB RAM:** no action needed. EggNOG's in-RAM database load (`--dbmem`)
  auto-disables when the job's RAM share is under 50 GB, so it mmaps the DB from
  local SSD instead (fine for the small substrate set ssign annotates). Force it
  either way with `-- --eggnog-dbmem` / `-- --no-eggnog-dbmem` if you need to.
- **No GHCR access / air-gapped:** build the image on an internet-connected machine
  and copy the `.sif`, or get it from a colleague. Nothing at run time needs the
  network except the one-time database fetch.
- **`image not found: ./ssign.sif`:** pass `--sif /path/to/your.sif` or
  `export SSIGN_SIF=...`; `ssign-run` looks in the current directory by default.

## Without `ssign-run` (raw apptainer)

`ssign-run` only assembles the command below; run it by hand if you can't use the
wrapper. `scripts/cx3/run_container_extended.pbs` is the fully-worked reference
(all binds + env vars); the essentials are:

```bash
apptainer run --nv --writable-tmpfs --containall \
  -B $HOME/ssign-databases:$HOME/ssign-databases \
  -B $HOME/.conda/envs/signalp6:$HOME/.conda/envs/signalp6:ro \
  -B genome.gbff:/work/in.gbff:ro -B out:/work/out \
  --env SSIGN_MAX_RAM_GB=60 \
  --env EGGNOG_DATA_DIR=$HOME/ssign-databases/eggnog \
  "$SSIGN_SIF" run /work/in.gbff --outdir /work/out --tier extended \
    --signalp-mode local --signalp-path $HOME/.conda/envs/signalp6/bin
```

`ssign-run` additionally auto-detects the ECOD / InterProScan sub-trees under
`--db-root` and sets `SSIGN_ECOD_DB` / `SSIGN_INTERPROSCAN_PATH`; set them yourself
here (see [`../reference/env_vars.md`](../reference/env_vars.md)) if you go the raw
route with those tools on.

## See also

- [`containers/README.md`](../../containers/README.md), building + publishing the
  image (maintainers).
- [`install.md`](install.md), the native/pip install of the optional tools.
- [`run_on_hpc.md`](run_on_hpc.md), cluster specifics (scheduler templates, walltime).
