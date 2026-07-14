# Installing ssign (container, the recommended way)

The Singularity/Apptainer image is the recommended way to run ssign on Linux and
HPC. It is **self-contained and reproducible**: the entire free toolchain and all
model weights are baked in and version-locked (digest-pinned base, `uv.lock`,
pinned tool commits, baked ESM/ProtT5 weights), so a run is reproducable over time. 
The things that are not included are: the DTU-licensed **SignalP 6**, and the **reference
databases**.

This guide covers the **extended** tier (secretion-system detection + secreted
protein identification + functional/structural annotation). Base and full differ
only in which databases you fetch, see [Tiers](#tiers).

The image is for **non-commercial research use** (it bundles DeepLocPro,
CC BY-NC-SA, and EggNOG-mapper, AGPL).

## Hardware requirements (extended)

| Resource | Recommended | Notes |
| --- | --- | --- |
| GPU | CUDA NVIDIA, ~16-24 GB VRAM (RTX 6000 / A40 / etc.) | DeepLocPro, DeepSecE, and pLM-BLAST embeddings use it. CPU-only works but the predictors are much slower. |
| RAM | **>= 64 GB** | EggNOG loads its database into RAM and the annotation tools run in parallel. 32 GB works only with `--no-eggnog-dbmem` (see [Troubleshooting](#troubleshooting)). |
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

# 0b. Get the launcher scripts (ssign-run, ssign-setup-dtu are bash wrappers, so pip
#     does NOT install them). Until this is merged to main, clone the branch:
git clone -b enrichment-circular-shift-per-run https://github.com/billerbeck-lab/ssign && cd ssign

# 0c. Get the image (.sif) and point ssign-run at it. At release:
#       apptainer pull oras://ghcr.io/billerbeck-lab/ssign:1.0.0
#     Until then, download the .sif a colleague shares (see "Getting the image" below)
#     into a folder with space:
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
  from the container (PBS), else ssign sizes tool memory to the whole node.

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

So after the install, running one genome is a single line:

```bash
scripts/ssign-run path/to/genome.gbff path/to/output_dir --tier extended \
  --sif "$SSIGN_SIF" --db-root /data/ssign-databases --max-ram 60
```

`genome.gbff` is your input (GenBank, GFF3, or FASTA), `output_dir` is where results
land.

## Verify

```bash
apptainer run --writable-tmpfs --containall "$SSIGN_SIF" doctor --tier extended
```

## Getting the image

The `.sif` is a single file, roughly 13-20 GB. Ways to get it onto a cluster,
simplest first:

- **A shared link** (OneDrive, Google Drive, Dropbox, or a large-file service such
  as WeTransfer Pro / Smash / Filemail). Upload the `.sif` once, share the link; the
  recipient downloads it onto their cluster, browser then `scp`, or directly on the
  cluster with `wget`/`gdown`/`rclone` if it has outbound internet. (Institutional
  OneDrive usually has room for 20 GB; a free Google Drive is only 15 GB total.)
- **Globus** for HPC-to-HPC: the standard for large research-data transfers, and most
  university clusters have a Globus endpoint. Reliable and resumable.
- **Build it** from the repo (`containers/build_sif.sh`): no transfer, fully
  reproducible; needs ~1 hour and internet on the build machine.
- **Zenodo** (after release): a public download link / DOI anyone can `wget`.

## Tiers

Same image every tier; the tier is the database set you fetch plus `--tier`.

| Tier | `fetch-databases` pulls | Approx size |
| --- | --- | --- |
| `base` | NCBI taxdump + Bakta light | ~22 GB |
| `extended` | + EggNOG + InterProScan + ECOD30 (pLM-BLAST) | ~100 GB |
| `full` | + Bakta full + HH-suite (Pfam/PDB70/UniRef30) + BLASTp Swiss-Prot | ~500 GB |

## Troubleshooting

- **< 64 GB RAM:** EggNOG's in-RAM database load can OOM. Add `-- --no-eggnog-dbmem`
  to the `ssign-run` command; it mmaps the DB from local SSD instead (fine for the
  small substrate set ssign annotates).
- **No GHCR access / air-gapped:** build the image on an internet-connected machine
  and copy the `.sif`, or get it from a colleague. Nothing at run time needs the
  network except the one-time database fetch.
- **`image not found: ./ssign.sif`:** pass `--sif /path/to/your.sif` or
  `export SSIGN_SIF=...`; `ssign-run` looks in the current directory by default.

## See also

- [`containers/README.md`](../../containers/README.md), per-tier mounts, the
  baked-vs-provided matrix, and the raw `apptainer` invocation.
- [`containers/PUBLISH.md`](../../containers/PUBLISH.md), how the release image is
  published (GHCR + Zenodo).
