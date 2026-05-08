# containers/

The `Dockerfile` here builds the official **ssign** bundle image: a
SHA-pinned CUDA + Python stack with the ssign Python package installed,
ready to run the pipeline against user-supplied genomes. Model weights
and reference databases are **not** baked in; they are fetched on the
host and bind-mounted at run time. This keeps the image at a few GB
instead of 600+ GB and avoids licence-redistribution friction.

## Prerequisites

- Linux host with NVIDIA GPU (≥ 16 GB VRAM recommended for the prediction tools)
- NVIDIA driver ≥ 550 (matches CUDA 12.4 runtime in the image)
- Docker with the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html), or Singularity ≥ 3.8 with `--nv`

## Fetch the assets first

```bash
# Reference databases (pick a tier: ~3 GB / ~150 GB / ~630 GB)
bash scripts/fetch_databases.sh --tier extended

# Model weights (~18 GB)
bash scripts/fetch_weights.sh
```

Both default to `~/.ssign/{databases,models}`. Override with `--target` if
you need a different host path.

## Build

```bash
docker build -f containers/Dockerfile -t ssign:1.0.0 .
```

The Dockerfile pins its CUDA base by SHA digest. The placeholder digest
in the file is replaced at release time; see the FRAGILE comment at the
top of `Dockerfile`.

## Run (Docker)

```bash
docker run --gpus all --rm \
    -v $HOME/.ssign/databases:/home/ssign/.ssign/databases:ro \
    -v $HOME/.ssign/models:/home/ssign/.ssign/models:ro \
    -v $PWD:/work \
    ssign:1.0.0 run /work/genome.gbff --outdir /work/output
```

Notes:

- `--gpus all` is required for the GPU-accelerated steps (DLP / DSE /
  SignalP / PLM-Effector / pLM-BLAST embedding).
- The two `:ro` mounts let the container read your fetched assets without
  modifying them.
- `/work` is the convention for input + output; mount whichever host
  directory contains your genomes.

## Run (Singularity, HPC)

Most academic HPC environments use Singularity instead of Docker. Pull
and run the same image:

```bash
singularity pull docker://ghcr.io/billerbeck-lab/ssign:1.0.0
singularity run --nv \
    --bind $HOME/.ssign/databases:/home/ssign/.ssign/databases:ro \
    --bind $HOME/.ssign/models:/home/ssign/.ssign/models:ro \
    --bind $PWD:/work \
    ssign_1.0.0.sif run /work/genome.gbff --outdir /work/output
```

Singularity uses `--nv` instead of `--gpus all`. Bind-mounts default to
read-write; the `:ro` suffix above mirrors the safety posture of the
Docker example for the asset directories.

## uid / gid alignment

The image creates a `ssign` user with uid 1000. If your host user has a
different uid, output written under `/work` will be owned by uid 1000
inside the container, which may show as a different owner on the host.
Two ways to fix this:

```bash
# Option 1: run as your host uid
docker run --user $(id -u):$(id -g) --gpus all ... ssign:1.0.0 ...

# Option 2: chown after the run
sudo chown -R $(id -u):$(id -g) ./output
```

Singularity automatically maps the host user into the container, so this
caveat doesn't apply there.

## What's bundled and what isn't

| Component                                    | Bundled? | Fetched by                        |
| -------------------------------------------- | -------- | --------------------------------- |
| CUDA 12.4 runtime + cuDNN                    | Yes      | base image                        |
| Python 3.10 + ssign[extended] deps           | Yes      | `pip install` during build (see note below)        |
| MacSyFinder + TXSScan profiles               | Yes      | pip install (macsyfinder package) |
| `ncbi-blast+` (`blastp` binary)              | Yes      | apt-get during build              |
| Model weights (DeepSecE / ProtT5 / ESM / …)  | No       | `scripts/fetch_weights.sh`        |
| Reference DBs (Bakta / EggNOG / HH-suite / …) | No       | `scripts/fetch_databases.sh`      |
| DeepLocPro + SignalP binaries                | No       | User-installed locally; DTU webserver as fallback * |
| EggNOG database                              | No       | `download_eggnog_data.py` *       |

The image pip-installs `ssign[extended]`. The `full` tier ships the same
Python dependency set; only the database tier picked by
`fetch_databases.sh` differs.

\* The DTU and EggNOG entries default to the conservative path because
their licences don't currently permit redistribution inside a public
image. If those terms change, the Dockerfile gets a one-line edit per
plan addendum E.6 and the entries flip to `Yes`.

## Troubleshooting

- **`Could not select device driver "" with capabilities: [[gpu]]`:**
  the NVIDIA Container Toolkit isn't installed or the daemon needs a
  restart. See the toolkit install guide linked above.
- **`CUDA error: no kernel image is available for execution on the device`:**
  your driver is older than CUDA 12.4 expects. Either upgrade the driver
  or rebuild against an older CUDA base (edit the `CUDA_BASE` ARG in
  `Dockerfile`).
- **`ssign run` exits with `database not found`:** the bind mount paths
  don't match the image's expected layout. Confirm the host directories
  exist (`ls $HOME/.ssign/databases`) and that they were populated by
  `fetch_databases.sh`.
- **First run downloads 7 GB of ESM weights:** expected if
  `fetch_weights.sh` was skipped. Either let it finish (it caches into
  `~/.ssign/models` for next time) or run the fetch script first.
