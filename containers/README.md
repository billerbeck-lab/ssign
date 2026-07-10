# containers/

Reproducible install artifacts for **ssign**.

- **`ssign.def`** — a Singularity/Apptainer definition. This is the **default**
  reproducible artifact: Linux + HPC, all tiers, one archivable `.sif`, runs
  with no root. Built from a digest-pinned CUDA base and the frozen `uv.lock`.
- **`build_sif.sh`** — build + offline smoke-test helper.
- **macOS** users install the **base** tier with a pinned pip install (every
  base tool is pip-installable; see below). Extended/full are Linux/HPC only
  (InterProScan, HH-suite, BLAST+ have no macOS build) and run via the `.sif`.
- **Windows** is not supported natively (several tools are Linux/macOS only);
  use WSL2, which is just the Linux path.

Databases and model weights are **not** baked in (licence + up to ~500 GB);
they are fetched on the host and bind-mounted at run time.

## Build the image

Build on a machine with **open internet** (a laptop, workstation, or CI). It is
normal for an **HPC compute node to be unable to build it**: clusters typically
allow only an allowlist of package hosts (PyPI, conda, Docker Hub, NVIDIA) and
block the general Ubuntu package servers the base image's `apt` needs. Build
once, then copy the `.sif` to the cluster (see "On HPC" below).

```bash
# from the repo root; needs apptainer + a real-disk build tmp (NOT a small tmpfs)
export APPTAINER_TMPDIR=$HOME/.ssign_build/tmp APPTAINER_CACHEDIR=$HOME/.ssign_build/cache TMPDIR=$APPTAINER_TMPDIR
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
apptainer build --fakeroot ssign.sif containers/ssign.def
# ...or the helper (also runs the offline smoke test):
containers/build_sif.sh /path/to/dir/containing/deeplocpro
```

Notes:
- The `.sif` is ~5 GB. The build pulls the CUDA base + the torch stack, so it
  overflows a RAM-backed `/tmp` (`[Errno 122] disk quota exceeded`) unless
  `APPTAINER_TMPDIR` points at real disk (as above).
- `--fakeroot` builds unprivileged (needs a `/etc/subuid` mapping, standard on
  most Linux). Nothing here needs system root.
- The base is pinned by digest; re-pin per the comment in `ssign.def` only if a
  base-image CVE forces it.

## Run the image

ssign needs, per run: (1) an install **tier** (`--tier base|extended|full`);
(2) the DTU predictors (DeepLocPro / SignalP), bind-mounted from a host install;
(3) the databases for the tier, bind-mounted; (4) the ESM model weights.

Three rules the container makes explicit (all validated on Imperial CX3):

1. **Always pass `--tier`.** The image has no tier marker, so it defaults to
   `extended` and will expect host-provided tools (EggNOG, etc.). Pick the tier
   matching the databases you mount.
2. **Mount DTU tools at their real host path.** DeepLocPro/SignalP are usually a
   conda env; bind-mount the *whole env at its own absolute path* so its bundled
   Python and libraries resolve, then point `--deeplocpro-path <env>/bin`.
3. **Provide the ESM weights.** DeepLocPro/DeepSecE load a ~2.5 GB ESM2 model. If
   you have run ssign natively once, mount your `~/.cache/torch` into the
   container; on an air-gapped node it cannot download it fresh.

Worked example (base tier, offline; the exact command validated on CX3):

```bash
DLPENV=$HOME/.conda/envs/deeplocpro          # host DeepLocPro conda env
OUT=$(mktemp -d)
apptainer run --nv --writable-tmpfs --containall \
  -B "$DLPENV":"$DLPENV":ro \
  -B $HOME/.cache/torch:/opt/ssign_home/.cache/torch \
  -B $PWD/genome.gbff:/work/in.gbff:ro \
  -B "$OUT":/work/out \
  ssign.sif run /work/in.gbff --outdir /work/out --tier base \
    --deeplocpro-mode local --deeplocpro-path "$DLPENV/bin"
```

- `--nv` passes the host GPU (drop it for CPU-only nodes; DeepLocPro then runs on
  CPU, fine for small inputs).
- `--writable-tmpfs` gives the read-only image a scratch overlay (caches, `$HOME`).
- `--containall` isolates it from your host environment; explicit `-B` mounts
  still apply. Apptainer passes host `SSIGN_*` env vars through, so you can wire
  databases with `export SSIGN_*_DB=...` instead of flags if you prefer.

### Per-tier database mounts

Same image for every tier; the tier is what you mount + `--tier`.

| Tier       | Mount into the container                                          | Extra host tools           |
| ---------- | ---------------------------------------------------------------- | -------------------------- |
| `base`     | (none — base needs no reference DB)                              | DeepLocPro, SignalP        |
| `extended` | EggNOG DB, InterProScan install, ECOD (pLM-BLAST DB), Bakta-light | + EggNOG-mapper, InterProScan install |
| `full`     | extended + HH-suite (Pfam/PDB70/UniRef30), Swiss-Prot, Bakta-full | + HH-suite install         |

Bind each host DB dir and point ssign at it via its `--*-db` flag or `SSIGN_*`
env var (see `docs/reference/env_vars.md` for the full list). What the image now
bakes for extended/full: **Java 11** (InterProScan's launcher needs it),
**pLM-BLAST code** (pinned clone; only its ECOD DB is mounted), **`ncbi-blast+`**.
Still host-provided: the **EggNOG-mapper** and **InterProScan** *installs* (mount
their conda env / install dir the same way as the DTU tools) and all reference
databases. EggNOG-mapper is host-provided because it hard-pins `biopython==1.76`
against ssign's `>=1.78`; mount it or `--skip-eggnog`.

Worked example (extended tier; DTU env, ESM cache, and the extended DBs mounted):

```bash
DLPENV=$HOME/.conda/envs/deeplocpro
EGGENV=$HOME/.conda/envs/eggnog             # host EggNOG-mapper env (emapper.py)
DBROOT=$EPHEMERAL/ssign-databases           # host DB tree
OUT=$(mktemp -d)
apptainer run --nv --writable-tmpfs --containall \
  -B "$DLPENV":"$DLPENV":ro \
  -B "$EGGENV":"$EGGENV":ro \
  -B $HOME/.cache/torch:/opt/ssign_home/.cache/torch \
  -B "$DBROOT":"$DBROOT" \
  -B $PWD/genome.gbff:/work/in.gbff:ro -B "$OUT":/work/out \
  --env EGGNOG_DATA_DIR="$DBROOT/eggnog" \
  --env SSIGN_ECOD_DB="$(ls -d $DBROOT/plm_blast/ECOD* | head -1)" \
  --env SSIGN_INTERPROSCAN_PATH="$(ls -d $DBROOT/interproscan/interproscan-* | head -1)" \
  ssign.sif run /work/in.gbff --outdir /work/out --tier extended \
    --deeplocpro-mode local --deeplocpro-path "$DLPENV/bin" \
    --eggnog-mode local --interproscan-mode local
```

Notes on this recipe:
- `$DBROOT` is mounted **read-write** (no `:ro`): InterProScan's launcher writes
  scratch files under its own install dir, so a read-only mount makes it fail
  (native runs it from writable `$EPHEMERAL`, so this matches native).
- pLM-BLAST needs no mount (baked); only its ECOD DB comes via `SSIGN_ECOD_DB`.
  The image sets `NUMBA_CACHE_DIR` to a writable tmpdir so pLM-BLAST's numba
  code imports in the read-only image (no action needed).
- Java is baked, so InterProScan's `interproscan.sh` finds `java` in-container.
- EggNOG-mapper: mount the host env at its real path (above) and prepend its
  `bin` to the container PATH so `emapper.py` resolves:
  `export APPTAINERENV_PREPEND_PATH="$EGGENV/bin"` before the run. `--skip-eggnog`
  if you have no EggNOG install.

### DTU webserver fallback

No local DeepLocPro/SignalP and no licence? Use the DTU webserver (needs
outbound internet, so not air-gapped HPC nodes):

```bash
apptainer run --writable-tmpfs --containall -B $PWD:/work \
  ssign.sif run /work/genome.gbff --outdir /work/out --tier base \
    --deeplocpro-mode remote --signalp-mode remote
```

## On HPC (build elsewhere, run here)

HPCs run Apptainer but usually can't build the image (allowlisted internet). So:

```bash
# 1. build on a machine with open internet (above), then copy the .sif up:
scp ssign.sif you@cluster:$EPHEMERAL/ssign.sif
# 2. on the cluster, inside a job, run as in the worked example, pointing at
#    $EPHEMERAL/ssign.sif and your host DeepLocPro env + ~/.cache/torch.
```

The archival release image (Zenodo) will bake the ESM weights in so it is fully
self-contained on air-gapped nodes.

## macOS (base tier)

Every base-tier tool is pip-installable, so macOS base needs no container:

```bash
# pinned to the same locked versions as the container (base tier = core deps)
python3 -m pip install -c containers/requirements-base.lock.txt ssign
```

DeepLocPro/SignalP still come from DTU (install locally, or `--*-mode remote`).
Extended/full stay on Linux/HPC via the `.sif`.

## What's bundled vs provided

| Component                                   | In the image | Provided at run time                    |
| ------------------------------------------- | ------------ | --------------------------------------- |
| CUDA 12.4 runtime, Python 3.10, ssign[extended] deps (pinned via `uv.lock`) | Yes | — |
| MacSyFinder + TXSScan 1.1.4 profiles, `hmmsearch` (pyhmmer), `ncbi-blast+` | Yes | — |
| Java 11 runtime (InterProScan launcher), pLM-BLAST code (pinned `4dddea3`) | Yes | — |
| ESM2 model weights (~2.5 GB)                | No (mount now; baked in the archival image) | host `~/.cache/torch` |
| DeepLocPro, SignalP (DTU licence)           | No           | bind-mount host env, or `--*-mode remote` |
| EggNOG-mapper + InterProScan installs, ECOD/EggNOG/HH-suite DBs | No | host install/DBs, bind-mounted (extended/full) |
| Reference databases                         | No           | `scripts/fetch_databases.sh` + bind-mount |

EggNOG-mapper is separate because it hard-pins `biopython==1.76` while ssign +
Bakta need `>=1.78`; install it on the host (`conda install -c bioconda
eggnog-mapper`) and mount it, or leave it off (`--skip-eggnog`).
