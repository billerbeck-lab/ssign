# Installing ssign

There are two ways to install ssign:

- **Container (recommended).** A single Singularity/Apptainer image that bakes
  the entire free toolchain and every model weight, version-locked. It is
  self-contained and reproducible: a run behaves the same today and in three
  years. This is the durable choice for paper-grade and HPC work. Start here.
- **Native / from source (secondary).** Lighter-weight (clone the repo and
  `pip install -e .`, plus a few optional tools), but you assemble and maintain
  the external toolchain yourself, and those tool versions can drift over time. See
  [Secondary install options](install-secondary.md).

---

# 1. Container install (recommended)

The image is **self-contained and reproducible**: the entire free toolchain and
all model weights are baked in and version-locked (digest-pinned CUDA base,
`uv.lock`, pinned tool commits, baked ESM/ProtT5 weights), so a run reproduces
over time.

Two things are **host-provided** rather than baked into the image:

- **SignalP 6** (DTU academic licence, cannot be redistributed). Install it
  locally with `scripts/ssign-setup-dtu`, or run with `--signalp-mode remote`.
- **InterProScan** (the scan engine plus its member databases). The image ships
  only a Java runtime to run a host-mounted InterProScan install; the engine and
  its ~24 GB of member databases come with the reference-database set.

Reference databases are likewise always fetched to the host (via
`fetch-databases`), never baked. Everything else (Bakta and its toolchain,
EggNOG-mapper, BLAST+, MacSyFinder, HH-suite, DeepLocPro + its ESM2 weights, the
DeepSecE checkpoint, pLM-BLAST + ProtT5) lives inside the image.

The image is for **non-commercial research use**: it bundles DeepLocPro
(CC BY-NC-SA 4.0) and EggNOG-mapper (AGPL). ssign's own code is GPL-3.0-or-later.

This guide covers the **extended** tier (secretion-system detection + secreted
protein identification + functional/structural annotation). Base and full differ
only in which databases you fetch, see [Tiers](#tiers).

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

- **`apptainer`** (or `singularity`). On HPC it is almost always preinstalled as
  a module: `module avail apptainer` (or `singularity`), then
  `module load apptainer`, then `apptainer --version` to confirm. If your cluster
  genuinely lacks it, email HPC support to add it (it is a standard rootless
  cluster tool; users normally should not install a system container runtime
  themselves). User-space fallback on clusters that allow user namespaces:
  `conda install -c conda-forge apptainer`.
- A licensed **SignalP 6** tarball (only if you want signal-peptide calls).
  Register (free for academics) at
  <https://services.healthtech.dtu.dk/services/SignalP-6.0/> and download once.

## Install (extended)


1. Confirm apptainer runs:
```
module load apptainer 2>/dev/null; apptainer --version
```
2. Get the launcher scripts:
```
git clone https://github.com/billerbeck-lab/ssign && cd ssign
```
3. Download the image:
```
wget -O $HOME/ssign.sif https://zenodo.org/records/21441318/files/ssign_1.0.0.sif
export SSIGN_SIF=$HOME/ssign.sif
```
4. Reference databases, can decide tier now (this is for extended):
```
mkdir -p $HOME/ssign-databases
apptainer run --writable-tmpfs --containall -B $HOME/ssign-databases:$HOME/ssign-databases \
  "$SSIGN_SIF" fetch-databases --tier extended --target $HOME/ssign-databases
```
5. SignalP 6; register + download first from https://services.healthtech.dtu.dk/services/SignalP-6.0/ then:
```
scripts/ssign-setup-dtu ~/signalp-6.0i.fast.tar.gz --signalp-only
```
6. Run
```
scripts/ssign-run $HOME/mygenome.gbff $HOME/ssign_out --tier extended \
  --sif "$SSIGN_SIF" --db-root $HOME/ssign-databases --max-ram 60
```

Everything else (Bakta + its toolchain, EggNOG-mapper, DeepLocPro + ESM2,
DeepSecE + ESM-1b, pLM-BLAST + ProtT5, InterProScan's Java, BLAST+ 2.17) is
inside the image, so a first run touches the network only for the one database
fetch.

- `--signalp-env` is auto-detected from your conda envs after step 2; the flag
  only overrides it. No SignalP licence? Drop step 2 and add
  `-- --signalp-mode remote` (uses the DTU webserver, which uploads your
  sequences), or omit SignalP entirely to run without signal-peptide calls.
- `--max-ram <GB>` = your job's RAM. Required on schedulers that hide the
  allocation from the container (PBS), else ssign sizes tool memory to the whole
  node. On SLURM it is read from `SLURM_MEM_PER_NODE` automatically.
- `--stage-image` copies the `.sif` to fast node-local disk before the run,
  turning a slow network-filesystem startup into one sequential copy.
  Recommended on HPC.

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
| HH-suite (`hhblits`/`hhsearch`, full tier) | yes | |
| **SignalP 6** (DTU academic licence) | no | conda env via `ssign-setup-dtu`, or `--signalp-mode remote` |
| **InterProScan** engine + member databases | no (only its Java is baked) | fetched by `fetch-databases` |
| **Reference databases** (Bakta DB, EggNOG DB, ECOD; full adds HH-suite DBs + Swiss-Prot) | no | `fetch-databases --tier <t>` |

InterProScan is the in-between one: its ~24 GB software + member-DB directory is
part of the database set (pulled by `fetch-databases`) and that same directory
supplies `interproscan.sh`. Only the Java to run it is baked.

## How a run works

- **"Mounting" (or "binding") a folder** just means "let the container see this
  folder on the host". `ssign-run` mounts your genome, the database folder, and
  the SignalP folder automatically.
- **A conda environment** is a self-contained folder holding one tool plus its
  dependencies (kept separate so tools don't clash). The only one you install is
  SignalP 6 (via `ssign-setup-dtu`), and `ssign-run` finds it for you.
  Everything else is baked into the image.
- **The GPU** is turned on automatically (`--nv`) when one is present. Nothing to
  do.
- **Where does the genome go?** Anywhere you can read it on the host (home,
  scratch, wherever). You give `ssign-run` the path to the file; there is no
  special location, and you do not copy it into the container.

## On HPC

Get the `.sif` onto the cluster: `apptainer pull` it on the login node, or
download it on your laptop and `scp` it up. Then fetch the DBs + install SignalP
once and run per job:

```bash
scp ssign.sif you@cluster:/path/to/scratch/ssign.sif   # your cluster's scratch/work dir
```

On PBS/SLURM, submit one job per genome that calls `ssign-run` with
`--stage-image` (copies the `.sif` to node-local disk first) and
`--max-ram <job GB>`. See [`run.md`](run.md) for ready SLURM and PBS job
templates, GPU requests, walltime, and per-job output dirs.

## Verify

```bash
apptainer run --writable-tmpfs --containall "$SSIGN_SIF" doctor --tier extended
```

`ssign doctor` checks every dependency and reports what is missing with the fix
command; it exits non-zero on failure so you can chain
`... doctor && ... run ...` in a job script.

## Tiers

Same image every tier; the tier is the database set you fetch plus `--tier`.

| Tier | `fetch-databases` pulls | Approx DB size |
| --- | --- | --- |
| `base` | Bakta light | ~4 GB |
| `extended` | + EggNOG + InterProScan + ECOD30 (pLM-BLAST) | ~100 GB |
| `full` | + Bakta full + HH-suite (Pfam/PDB70/UniRef30) + BLASTp Swiss-Prot | ~500 GB |

> **macOS or base tier?** You do not need the container. Every base-tier tool is
> pip-installable, so macOS and base-tier runs use the native path (see
> [Secondary install options](install-secondary.md)): clone the repo and
> `pip install -e . -c containers/requirements-base.lock.txt` (the lock pins the
> base deps). SignalP 6 is DTU-licensed (install locally, or
> `--signalp-mode remote`); DeepLocPro is open source (install from its GitHub
> repo, or `--deeplocpro-mode remote`). Extended/full stay on Linux/HPC via the
> `.sif`.

## Troubleshooting

- **< 64 GB RAM:** no action needed. EggNOG's in-RAM database load (`--dbmem`)
  auto-disables when the job's RAM share is under 50 GB, so it mmaps the DB from
  local SSD instead (fine for the small substrate set ssign annotates). Force it
  either way with `-- --eggnog-dbmem` / `-- --no-eggnog-dbmem` if you need to.
- **No GHCR access / air-gapped:** download the image on an internet-connected
  machine and copy the `.sif`, or get it from a colleague. Nothing at run time
  needs the network except the one-time database fetch.
- **`image not found: ./ssign.sif`:** pass `--sif /path/to/your.sif` or
  `export SSIGN_SIF=...`; `ssign-run` looks in the current directory by default.

## Without `ssign-run` (raw apptainer)

`ssign-run` only assembles the command below; run it by hand if you can't use the
wrapper. The essentials are:

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
`--db-root` and sets `SSIGN_ECOD_DB` / `SSIGN_INTERPROSCAN_PATH`; set them
yourself here (see [`../reference/env_vars.md`](../reference/env_vars.md)) if you
go the raw route with those tools on.

## See also

- [`containers/README.md`](../../containers/README.md), building + publishing the
  image (maintainers).
- [`run.md`](run.md), running the pipeline: configuration flags and cluster
  specifics (scheduler templates, walltime).
- [`install-secondary.md`](install-secondary.md), native / pip install without
  the container (per-tool setup).
