# containers/

Reproducible install artifacts for **ssign**. The Singularity/Apptainer image is
the default way to run ssign on Linux/HPC: one archivable `.sif`, every tier, no
root. It is **self-contained**. The entire free toolchain and all model weights
are baked in; you provide only two things: the DTU-licensed **SignalP 6**, and
the **reference databases**.

- **`ssign.def`** the image definition (digest-pinned CUDA base + frozen `uv.lock`).
- **`build_sif.sh`** build + offline smoke-test helper.
- **`scripts/ssign-run`** the one-command launcher (recommended). It builds the
  `apptainer run` line for you: binds the DBs + SignalP, sets the tier and RAM
  budget, stages the image, adds `--nv`.
- **`scripts/ssign-setup-dtu`** installs SignalP 6 into a conda env.

The image bundles DeepLocPro (CC BY-NC-SA) and EggNOG-mapper (AGPL), so it is for
**non-commercial** research use.

macOS: the **base** tier installs via pinned pip (no container); extended/full
are Linux/HPC only. Windows: use WSL2 (the Linux path).

## Quickstart (extended / tier 2)

Three commands once you have the image. `ssign-run` and `ssign-setup-dtu` live in
`scripts/` in the repo, and inside the image at `/opt/ssign/scripts/`.

```bash
# 0. Get the launcher scripts + the image. The scripts wrap apptainer, so pip does
#    not install them; git clone gives you scripts/ (and the docs):
git clone https://github.com/billerbeck-lab/ssign && cd ssign
#    Image: apptainer pull oras://ghcr.io/billerbeck-lab/ssign:1.0.0 (at release),
#    or build once (see "Build" below) -> ssign.sif

# 1. Reference databases, fetched FROM the image (no host tools needed; ~100 GB extended)
apptainer run --writable-tmpfs --containall -B /data/ssign-databases:/data/ssign-databases \
  ssign.sif fetch-databases --tier extended --target /data/ssign-databases

# 2. SignalP 6, the ONLY licence-gated tool. Register + download the tarball once
#    (https://services.healthtech.dtu.dk/services/SignalP-6.0/), then:
scripts/ssign-setup-dtu ~/Downloads/signalp-6.0i.fast.tar.gz --signalp-only

# 3. Run
scripts/ssign-run genome.gbff out --tier extended \
  --sif ssign.sif \
  --db-root /data/ssign-databases \
  --signalp-env ~/.conda/envs/signalp6 \
  --max-ram 60
```

- `--sif` is your image path (default `./ssign.sif`); or `export SSIGN_SIF=/path/to/ssign.sif`
  once instead of passing the flag each run.
- `--signalp-env` is auto-detected from your conda envs after step 2, so the flag
  is only needed to override it.
- `--max-ram GB`: pass your job's RAM on schedulers that hide the allocation from
  the container (e.g. PBS). Without it ssign reads the whole physical node and can
  over-size tool heaps / EggNOG `--dbmem` and OOM. On SLURM it is read from
  `SLURM_MEM_PER_NODE` automatically.
- Add `--stage-image` on HPC to copy the `.sif` to fast node-local disk first
  (turns a slow network-filesystem startup into one sequential copy).
- No SignalP licence? Skip step 2 and add `-- --signalp-mode remote` (uses the
  DTU webserver, which uploads your sequences), or run without signal-peptide calls.

Everything else (Bakta + its toolchain, EggNOG-mapper, DeepLocPro + ESM2, DeepSecE
+ its checkpoint, pLM-BLAST + ProtT5, BLAST+ 2.17, Java 11) is inside the image,
so a first run touches the network only for the database fetch.

## What's baked vs what you provide

| Component | In the image | You provide |
| --- | :---: | --- |
| CUDA 12.4 runtime, Python 3.10, `ssign[extended]` deps (pinned `uv.lock`) | yes | |
| MacSyFinder + TXSScan 1.1.4 profiles, `hmmsearch` (pyhmmer) | yes | |
| Bakta 1.12.0 + toolchain (BLAST+ 2.17, real HMMER, DIAMOND, tRNAscan-SE, aragorn, INFERNAL, PilerCR, AMRFinderPlus) | yes | |
| EggNOG-mapper (`emapper.py`), Java 11, rsync | yes | |
| DeepLocPro + ESM2 backbone | yes | |
| DeepSecE + its ESM-1b checkpoint | yes | |
| pLM-BLAST code + ProtT5 encoder | yes | |
| **SignalP 6** (DTU academic licence) | no | conda env via `ssign-setup-dtu`, or `--signalp-mode remote` |
| **Reference databases** (Bakta DB, EggNOG DB, InterProScan install+DBs, ECOD; full adds HH-suite + Swiss-Prot) | no | `fetch-databases --tier <t>` |

InterProScan is the one in-between: its ~24 GB software+member-DB directory is
part of the database set (pulled by `fetch-databases`), and that same directory
also supplies `interproscan.sh`. The Java to run it is baked.

## Tiers

Same image for every tier; the tier is the database set you fetch, plus `--tier`.

| Tier | `fetch-databases` pulls | Approx size |
| --- | --- | --- |
| `base` | NCBI taxdump + Bakta light | ~22 GB |
| `extended` | + EggNOG + InterProScan + ECOD30 (pLM-BLAST) | ~100 GB |
| `full` | + Bakta full + HH-suite (Pfam/PDB70/UniRef30) + BLASTp Swiss-Prot | ~500 GB |

## Build the image (maintainers)

Build on a machine with **open internet** (laptop, workstation, or CI). An HPC
compute node usually **can't** build it: clusters allowlist package hosts and
block the general Ubuntu servers the base image's `apt` needs. Build once, copy
the `.sif` up.

```bash
# needs apptainer + a real-disk build tmp (NOT a small tmpfs). ~15 GB image, ~30-60 min.
containers/build_sif.sh          # build + offline smoke test
# ...or by hand:
export APPTAINER_TMPDIR=$HOME/.ssign_build/tmp APPTAINER_CACHEDIR=$HOME/.ssign_build/cache TMPDIR=$APPTAINER_TMPDIR
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
apptainer build --fakeroot ssign.sif containers/ssign.def
```

- The build downloads ~7.5 GB of model weights (ESM2, the DeepSecE checkpoint from
  a slow academic mirror, ProtT5). Keep the connection up: a container build has no
  resume, and a dropped download restarts that file from 0 and can abort the build.
  The DeepSecE mirror is the slow part; everyone who pulls the finished image skips
  it entirely.
- `--fakeroot` builds unprivileged (standard `/etc/subuid` mapping; no system root).
- The base is digest-pinned; re-pin per the `ssign.def` comment only if a base CVE
  forces it.

## On HPC

```bash
# 1. build off-cluster (above), then copy the image up:
scp ssign.sif you@cluster:$EPHEMERAL/ssign.sif
# 2. fetch the DBs once (from the image) + install SignalP once, then run per job.
```

On PBS/CX3, use the ready wrapper `scripts/cx3/run_container_extended.pbs`: it
stages the image to node-local disk, pins the RAM budget to the job's `mem`,
`--nv`, and mounts only SignalP + the DB root (everything else is baked).

```bash
qsub -v GENOME=$HOME/g.gbff,SIF=$EPHEMERAL/ssign.sif scripts/cx3/run_container_extended.pbs
```

## DTU webserver fallback (no SignalP licence)

```bash
scripts/ssign-run genome.gbff out --tier extended --db-root /data/ssign-databases \
  -- --signalp-mode remote
```

Needs outbound internet (so not air-gapped nodes) and uploads your sequences to DTU.

## macOS (base tier)

Every base-tier tool is pip-installable, so macOS base needs no container:

```bash
python3 -m pip install -c containers/requirements-base.lock.txt ssign
```

SignalP/DeepLocPro still come from DTU (install locally, or `--*-mode remote`).
Extended/full stay on Linux/HPC via the `.sif`.

## Without `ssign-run` (raw apptainer)

`ssign-run` only assembles the command below; run it by hand if you can't use the
wrapper. `scripts/cx3/run_container_extended.pbs` is the fully-worked reference
(all binds + env vars); the essentials are:

```bash
apptainer run --nv --writable-tmpfs --containall \
  -B /data/ssign-databases:/data/ssign-databases \
  -B $HOME/.conda/envs/signalp6:$HOME/.conda/envs/signalp6:ro \
  -B genome.gbff:/work/in.gbff:ro -B out:/work/out \
  --env SSIGN_MAX_RAM_GB=60 \
  --env EGGNOG_DATA_DIR=/data/ssign-databases/eggnog \
  ssign.sif run /work/in.gbff --outdir /work/out --tier extended \
    --signalp-mode local --signalp-path $HOME/.conda/envs/signalp6/bin
```

`ssign-run` additionally auto-detects the ECOD / InterProScan sub-trees under
`--db-root` and sets `SSIGN_ECOD_DB` / `SSIGN_INTERPROSCAN_PATH`; set them yourself
here (see `docs/reference/env_vars.md`) if you go the raw route with those tools on.
