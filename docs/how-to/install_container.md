# Installing ssign (container, the recommended way)

The Singularity/Apptainer image is the recommended way to run ssign on Linux and
HPC. It is **self-contained and reproducible**: the entire free toolchain and all
model weights are baked in and version-locked (digest-pinned base, `uv.lock`,
pinned tool commits, baked ESM/ProtT5 weights), so a run produces the same result
today and in five years and does not break when an upstream tool changes. You
provide only two things: the DTU-licensed **SignalP 6**, and the **reference
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

- `apptainer` (or `singularity`). On HPC it is usually a module: `module load apptainer`.
- A licensed **SignalP 6** tarball. Register (free for academics) at
  <https://services.healthtech.dtu.dk/services/SignalP-6.0/> and download once.

## Install (extended)

```bash
# 0. Get the launcher scripts + the image. The scripts (ssign-run, ssign-setup-dtu)
#    are bash wrappers, so pip does NOT install them; git clone gives you scripts/.
git clone https://github.com/billerbeck-lab/ssign && cd ssign
apptainer pull oras://ghcr.io/billerbeck-lab/ssign:1.0.0     # at release
#    (until release, build it once: containers/build_sif.sh, see below)
export SSIGN_SIF=$PWD/ssign_1.0.0.sif                        # so ssign-run finds it

# 1. Reference databases, fetched FROM the image (no host tools needed; ~100 GB)
apptainer run --writable-tmpfs --containall -B /data/ssign-databases:/data/ssign-databases \
  "$SSIGN_SIF" fetch-databases --tier extended --target /data/ssign-databases

# 2. SignalP 6, the only licence-gated tool
scripts/ssign-setup-dtu ~/Downloads/signalp-6.0i.fast.tar.gz --signalp-only

# 3. Run
scripts/ssign-run genome.gbff out --tier extended \
  --db-root /data/ssign-databases \
  --signalp-env ~/.conda/envs/signalp6 \
  --max-ram 60
```

Everything else (Bakta + its toolchain, EggNOG-mapper, DeepLocPro + ESM2, DeepSecE
+ ESM-1b, pLM-BLAST + ProtT5, InterProScan's Java, BLAST+ 2.17) is inside the
image, so a first run touches the network only for the one database fetch.

- `--signalp-env` is auto-detected from your conda envs after step 2; the flag only
  overrides it. No SignalP licence? Drop step 2 and add `-- --signalp-mode remote`
  (uses the DTU webserver, which uploads your sequences).
- `--max-ram <GB>` = your job's RAM. Required on schedulers that hide the allocation
  from the container (PBS), else ssign sizes tool memory to the whole node.

## How a run works (plain English)

The `.sif` image is a sealed box that already contains every program ssign needs.
A container is walled off from the rest of the machine by default (that is what
makes it reproducible), so it cannot see your files unless you explicitly hand them
over. `ssign-run` does all of that handing-over for you; you never type the plumbing.

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

So after the one-time install, running one genome is a single line:

```bash
scripts/ssign-run path/to/genome.gbff path/to/output_dir --tier extended \
  --sif "$SSIGN_SIF" --db-root /data/ssign-databases --max-ram 60
```

`genome.gbff` is your input (GenBank, GFF3, or FASTA), `output_dir` is where results
land. That is the whole run; ssign-run handles the mounts, the GPU, and finding SignalP.

## Verify

```bash
apptainer run --writable-tmpfs --containall "$SSIGN_SIF" doctor --tier extended
```

## On HPC

Build the image off-cluster (compute nodes usually can't reach the package hosts),
copy the `.sif` up, then fetch DBs once and run per job.

```bash
scp ssign.sif you@cluster:$EPHEMERAL/ssign.sif
```

- Add `--stage-image` to `ssign-run` so it copies the `.sif` to node-local SSD
  before running (a `.sif` on a network filesystem otherwise adds a long
  import-startup).
- **PBS/CX3:** use the ready wrapper `scripts/cx3/run_container_extended.pbs`, it
  stages the image, pins the RAM budget to the job's `mem`, and mounts only SignalP
  + the DB root:
  ```bash
  qsub -v GENOME=$HOME/g.gbff,SIF=$EPHEMERAL/ssign.sif scripts/cx3/run_container_extended.pbs
  ```

## Building the image yourself

Until the release image is published, build once on a machine with open internet
(laptop, workstation, CI); this is normal, HPC compute nodes usually cannot build
it. Needs `apptainer` and a real-disk build tmp (not a small tmpfs). ~1 hour; the
image is ~20 GB.

```bash
containers/build_sif.sh          # build + offline smoke test -> ssign.sif
```

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
