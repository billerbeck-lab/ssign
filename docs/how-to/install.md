# Installing ssign

There are two ways to install ssign:

- **Container (recommended).** A single Singularity/Apptainer image that bakes
  the entire free toolchain and every model weight, version-locked. It is
  self-contained and reproducible: a run behaves the same today and in three
  years. This is the durable choice for paper-grade and HPC work. Start here.
- **Native / pip (secondary).** Lighter to get going (`pip install ssign` plus a
  few optional tools), but you assemble and maintain the external toolchain
  yourself, and those tool versions can drift out of compatibility over time.

Both paths run the same command-line pipeline. ssign is CLI-only:

```bash
ssign run <genome.gbff> --outdir <out>   # run the pipeline
ssign doctor                             # verify the install
ssign fetch-databases --tier <tier>      # download reference databases
```

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

```bash
# 0a. Load apptainer (on HPC it is a module) and confirm it runs
module load apptainer && apptainer --version

# 0b. Get the launcher scripts (ssign-run / ssign-setup-dtu wrap apptainer, so
#     pip does not put them on PATH; a clone gives you scripts/ and the docs):
git clone https://github.com/billerbeck-lab/ssign && cd ssign

# 0c. Get the image (.sif) and point ssign-run at it. At release:
#       apptainer pull oras://ghcr.io/billerbeck-lab/ssign:1.0.0
#     until then, build it once: containers/build_sif.sh (see containers/README.md).
export SSIGN_SIF=$HOME/ssign.sif

# 1. Reference databases, fetched FROM the image (no host tools needed; ~100 GB, once)
apptainer run --writable-tmpfs --containall -B $HOME/ssign-databases:$HOME/ssign-databases \
  "$SSIGN_SIF" fetch-databases --tier extended --target $HOME/ssign-databases

# 2. SignalP 6, the only licence-gated tool (register + download the tarball first, then:)
scripts/ssign-setup-dtu ~/signalp-6.0i.fast.tar.gz --signalp-only

# 3. Run (one line per genome)
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

Build the image off-cluster (an HPC compute node usually can't, see
[`containers/README.md`](../../containers/README.md)), copy the `.sif` up, then
fetch the DBs + install SignalP once and run per job:

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
> pip-installable, so macOS and base-tier runs use the native path below:
> `python3 -m pip install -c containers/requirements-base.lock.txt ssign`
> (available from PyPI at the v1.0.0 release; until then install from source,
> git clone + `pip install -e .`). SignalP 6 is DTU-licensed (install locally, or
> `--signalp-mode remote`); DeepLocPro is open source (install from its GitHub
> repo, or `--deeplocpro-mode remote`). Extended/full stay on Linux/HPC via the
> `.sif`.

## Troubleshooting

- **< 64 GB RAM:** no action needed. EggNOG's in-RAM database load (`--dbmem`)
  auto-disables when the job's RAM share is under 50 GB, so it mmaps the DB from
  local SSD instead (fine for the small substrate set ssign annotates). Force it
  either way with `-- --eggnog-dbmem` / `-- --no-eggnog-dbmem` if you need to.
- **No GHCR access / air-gapped:** build the image on an internet-connected
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

---

# 2. Native / pip install (secondary)

This path is lighter-weight: `pip install ssign` plus a handful of optional
tools. The trade-off is that you install and maintain the external toolchain
yourself (Bakta, EggNOG-mapper, BLAST+, HH-suite, InterProScan, SignalP 6,
DeepLocPro), and those tool versions can drift out of compatibility over time. If
you want reproducible, paper-grade runs that behave identically years from now,
use the [container](#1-container-install-recommended) instead.

`pip install ssign` gets you the base pipeline: secretion-system detection,
secreted-protein prediction (DeepLocPro, DeepSecE, and SignalP, see the
DeepLocPro and SignalP sections below for installing the DTU tools locally, or
for the opt-in webserver fallback when you don't have a DTU licence), proximity
analysis, and reporting. The tools below extend what the pipeline reports and how
it runs.

The simplest path is one of the named tiers:

```bash
pip install ssign                 # base
pip install ssign[extended]       # base + Bakta + pLM-BLAST + extended-tier pins
pip install ssign[full]           # extended deps + full database tier
```

> Available from PyPI at the v1.0.0 release; until then install from source
> (git clone + `pip install -e .`).

After pip install, fetch the matching database bundle:

```bash
bash scripts/fetch_databases.sh --tier base       # or: extended / full
```

### Database tier sizes

Real on-disk numbers, measured 2026-06-03 against the fetched databases.
Cumulative, each row includes everything from the rows above it:

| Tier | Cumulative size | Adds |
|---|---|---|
| base | ~4 GB | Bakta light DB (3.9 GB) |
| extended | ~100 GB | + EggNOG (47 GB), InterProScan (24 GB), pLM-BLAST ECOD30 (11 GB) |
| full | ~500 GB | + Bakta full DB (84 GB, replaces light), HH-suite Pfam + PDB70 + UniRef30 (340 GB total), BLASTp-vs-Swiss-Prot (0.3 GB) |

These are databases only. Model weights are separate and shared across all tiers,
~18 GB total (base auto-downloads the DeepSecE checkpoint and ESM backbones;
ProtT5 is extended-tier only; DeepLocPro and SignalP weights ship with their own
installs). The "Acquired by" column says which weights auto-download on first run
and which come from a user-provided local install; the container bakes all of them:

| Model | Size | Used by | Acquired by |
|---|---|---|---|
| DeepSecE checkpoint (fine-tuned ESM-1b) | ~2.5 GB | DeepSecE | auto-download on first run |
| ESM-1b | ~7.3 GB | DeepSecE backbone | auto-download on first run |
| ESM-2 | ~2.5 GB | DeepLocPro backbone | auto-download on first run |
| ProtT5 (half) | ~2.4 GB | pLM-BLAST query embedding (extended tier) | auto-download on first run |
| DeepLocPro weights | ~2 GB | DeepLocPro (CC BY-NC-SA 4.0, cloned from public GitHub at a pinned commit; no acquisition step) | ships with the local DeepLocPro install |
| SignalP 6.0 weights | ~1.5 GB | SignalP (DTU academic licence, user-acquired) | ships with the SignalP tarball |

DeepLocPro/DeepSecE pull their ESM backbones via `esm.pretrained` (torch-hub
cache), pLM-BLAST pulls ProtT5 via `huggingface_hub`, and the DeepSecE checkpoint
downloads via `run_deepsece._ensure_checkpoint()`, so those need no separate
fetch step. The container bakes all of them (container runs are fully offline).
The two weights you provide yourself are SignalP 6 (DTU tarball) and DeepLocPro
(its GitHub clone); their weights ship with those local installs (see below).

HH-suite (extracted) is the long pole for full-tier disk. BLASTp defaults to
Swiss-Prot (tiny, curated); full NR (~800 GB) is opt-in only (too slow to blast a
real substrate set), fetch it via `fetch_blast_nr` and pass
`--blastp-db <nr-dir>/nr`.

If a tier doesn't fit your storage budget, pick individual tools below.

## Verify the install

After installing the tools and fetching databases, run:

```bash
ssign doctor --tier extended      # or: base / full to match what you installed
```

`ssign doctor` checks every dependency ssign needs and reports what's missing
with the exact fix command for each: Python packages, external binaries on PATH
(Bakta, EggNOG-mapper, HH-suite, BLAST+, InterProScan), on-disk databases (read
from `~/.ssign/db_root` written by `fetch_databases.sh`, or `SSIGN_*` env vars if
you set them), and model weights. Exits non-zero on failure so you can chain
`ssign doctor && ssign run …` in scripts.

If `ssign doctor` is green, the pipeline can run.

For environment variables (mirror URLs, database paths, dev-only flags), see
[`reference/env_vars.md`](../reference/env_vars.md).

---

## DeepSecE (base tier, GPU recommended)

DeepSecE predicts secretion-system type per protein and runs as a second opinion
alongside DeepLocPro. It ships in the base tier: `pip install ssign` already
includes it (PyTorch plus the ESM protein language model add ~7.3 GB to the
install), so there is no separate install step. The legacy `ssign[deepsece]`
extra is a no-op alias kept only so old invocations keep working; it installs
nothing.

CUDA GPU strongly recommended; CPU runs are slow.

---

## Bakta (pip extra + database)

Bakta provides annotation-grade gene calling and functional descriptions. ssign
re-annotates inputs with Bakta by default; if Bakta is not installed, the
pipeline falls back to a pyrodigal-only call without functional annotations.

**Important:** `pip install ssign[bakta]` only installs the Bakta Python wrapper.
Bakta also depends on several binary tools (AMRFinderPlus, DIAMOND, HMMER,
tRNAscan-SE, aragorn) that aren't pip-installable. Without them, **even
`bakta_db download` fails** because Bakta's startup runs the same dependency
check as `bakta` itself.

Building that binary set piecemeal always misses one (`ERROR: tRNAscan-SE not
found`, then `ERROR: PilerCR not found`, ...). Install all of Bakta from conda in
one shot, which pulls every binary (tRNAscan-SE, aragorn, PilerCR, AMRFinderPlus,
DIAMOND, HMMER, BLAST+) at compatible versions:

```bash
mamba create -n bakta-deps -c bioconda bakta -y     # or -n base to skip a new env
export PATH=~/.conda/envs/bakta-deps/bin:$PATH       # only for a dedicated env
```

Then download the light database (~3.9 GB extracted):

```bash
bakta_db download --output ~/bakta_db --type light
```

Pass `--bakta-db ~/bakta_db`, or set `SSIGN_BAKTA_DB=~/bakta_db` (read by
`scripts/fetch_databases.sh`). The full database (~84 GB extracted) is the
`--type full` variant; `fetch_databases.sh --tier full` pulls it.

---

## EggNOG-mapper (separate install + database)

EggNOG-mapper provides ortholog-based functional annotation (COG, KEGG, GO, PFAM)
for substrate proteins. It is invoked as a subprocess (`emapper.py`), not
imported by ssign, and is **not** included in the `[extended]` / `[full]` pip
extras: upstream eggnog-mapper hard-pins `biopython==1.76` while ssign and Bakta
need `biopython>=1.78`, which makes the two unsatisfiable in a single pip
resolution. Install it separately.

The conda path is recommended; bioconda has shipped eggnog-mapper against modern
biopython for years without breakage:

```bash
conda install -c bioconda eggnog-mapper
```

If you don't use conda, `--no-deps` skips the upstream pin and lets ssign's
biopython (`>=1.80`) satisfy eggnog-mapper at runtime. The only API
incompatibility (`Bio.Alphabet`, removed in biopython 1.78) is already guarded
with `try/except` in eggnog-mapper itself, so this works in practice; bioconda
relies on the same:

```bash
pip install --no-deps eggnog-mapper
```

Verify the install:

```bash
which emapper.py
emapper.py --version
```

> **HPC users (Imperial CX3 and similar):** ssign's bundled `bakta-deps` conda
> env does **not** ship eggnog-mapper. Either install it into your active ssign
> venv (the `pip install --no-deps` line above), or install it into `bakta-deps`
> explicitly with `conda install -n bakta-deps -c bioconda eggnog-mapper` so the
> runner's auto-discovery picks it up. ssign now hard-fails at pre-flight if
> `emapper.py` is missing and `--skip-eggnog` is not set, so a missing install no
> longer wastes a 1-hour PBS allocation.

Fetch the database (~47 GB extracted):

```bash
scripts/fetch_databases.sh --tier extended --target ~/.ssign/databases
```

This wgets the three required files (`eggnog.db`, `eggnog.taxa.db`,
`eggnog_proteins.dmnd`) directly from `eggnog5.embl.de`. We don't use
`download_eggnog_data.py` because eggnog-mapper 2.1.13 (the latest on bioconda as
of 2026-05) still hardcodes the retired `eggnogdb.embl.de` hostname and produces
0-byte files with exit-code 0. 2.1.14 fixed it upstream but never reached PyPI,
so the fetch script bypasses that breakage. Tell ssign where the database lives:

```bash
ssign run input.gbff --outdir results --eggnog-db ~/.ssign/databases/eggnog
```

EggNOG annotation is on by default at the extended (default) and full tiers; it
is off only at `--tier base`. Pass `--skip-eggnog` to disable it.

> **HPC / shared scratch users:** `--eggnog-dbmem` (loads `eggnog.db` into RAM,
> ~44 GB resident) is now **auto** by default: on only when the job's RAM share
> is >= 50 GB, else the on-disk SQLite is memory-mapped. The runner stages that
> DB to local scratch first, so the old NFS-mmap hang no longer applies. Force
> either way with `--eggnog-dbmem` / `--no-eggnog-dbmem`.

---

## BLAST+ (system binary)

`blastp` is needed for cross-genome ortholog grouping. Download is about 200 MB.

```bash
sudo apt install ncbi-blast+        # Debian / Ubuntu
brew install blast                  # macOS
conda install -c bioconda blast     # cross-platform
```

ssign auto-detects `blastp` on the next launch.

---

## InterProScan

InterProScan (EBI) scans proteins against a panel of member databases (Pfam,
TIGRFAM, HAMAP, SMART, PIRSF, SUPERFAMILY, Gene3D, ProSite, CDD) to annotate
domains, family memberships, and GO terms. ssign uses it to add domain-level
annotation to substrate proteins. Java required.

The bundle is large (~24 GB extracted) but installs as a single tarball:

```bash
# 1. Java 11+ on PATH (Ubuntu / Debian):
sudo apt install openjdk-17-jre-headless

# 2. Fetch + extract via the helper script (pinned version, currently
#    5.77-108.0; bump scripts/fetch_databases.sh `IPS_VERSION` when EBI
#    rotates it off, old version dirs get 404'd after a few releases):
scripts/fetch_databases.sh --tier extended --target ~/.ssign/databases
```

Point ssign at the install directory (the one containing `interproscan.sh`):

```bash
export SSIGN_INTERPROSCAN_PATH=~/.ssign/databases/interproscan/interproscan-5.77-108.0
```

ssign runs InterProScan with the bacterial-relevant member DBs by default
(PANTHER, the slowest member and eukaryote-leaning, is excluded). Per-protein
scan time is typically 5-30 s; a whole-genome run on ~5,000 proteins is 30-90
minutes. The first run also queries EBI's precalculated-match lookup service for
a 5-10x speedup on known sequences; add `-dp` behaviour via your own wrapper
script if you need air-gapped operation.

---

## HH-suite (system binary + databases)

`hhsearch` and `hhblits` (Steinegger / Söding labs) detect remote evolutionary
relationships by comparing HMM-vs-HMM profiles. ssign uses HH-suite to annotate
substrate proteins with Pfam domain families and PDB structural homologs.

Install the binaries:

```bash
conda install -c bioconda hhsuite      # recommended
# or build from source: https://github.com/soedinglab/hh-suite
```

Download the three databases. Two canonical mirrors host them; **prefer
Tübingen** (fresher, recommended by Söding lab issue #382) and fall back to GWDG
only if Tübingen is unreachable.

```bash
# Pfam (domain families). Tübingen has v38 (2024+); GWDG has v35 (2021).
mkdir -p $HHSUITE_DBS && cd $HHSUITE_DBS
wget http://ftp.tuebingen.mpg.de/pub/ebio/protevo/toolkit/databases/hhsuite_dbs/PfamA_v38_2.tar.gz
tar -xzf PfamA_v38_2.tar.gz && rm PfamA_v38_2.tar.gz

# PDB70 (structural homology). Tübingen has 2026-02 build; GWDG has 2022-03.
wget http://ftp.tuebingen.mpg.de/pub/ebio/protevo/toolkit/databases/hhsuite_dbs/pdb70_from_mmcif_2026-02-20.tar.gz
tar -xzf pdb70_from_mmcif_2026-02-20.tar.gz && rm pdb70_from_mmcif_2026-02-20.tar.gz

# UniRef30 (clustered UniProt for hhblits MSA generation). Only at GWDG.
wget https://wwwuser.gwdg.de/~compbiol/uniclust/2023_02/UniRef30_2023_02_hhsuite.tar.gz
tar -xzf UniRef30_2023_02_hhsuite.tar.gz && rm UniRef30_2023_02_hhsuite.tar.gz
```

Total disk after extraction: ~340 GB (Pfam 8 GB + PDB70 70 GB + UniRef30 261 GB).
Tell ssign where to find them:

```bash
export SSIGN_HHSUITE_PFAM=$HHSUITE_DBS/pfam
export SSIGN_HHSUITE_PDB70=$HHSUITE_DBS/pdb70_from_mmcif_2026-02-20
export SSIGN_HHSUITE_UNICLUST=$HHSUITE_DBS/UniRef30_2023_02
```

These three are read at run time as fallbacks for the matching CLI flags
(`--hhsuite-pfam-db`, `--hhsuite-pdb70-db`, `--hhsuite-uniclust-db`).

### Mirror caveats

The Söding lab is in maintenance mode; per
[hh-suite issue #382](https://github.com/soedinglab/hh-suite/issues/382), no
funding for support, but the files at GWDG are stable. Tübingen has fresher
builds. If you see a 404 on the GWDG host, probe with `curl -IL` (capital L
follows the 302 redirect to the sibling domain). If both mirrors are down for an
extended period, pre-stage the DBs on your HPC scratch directory.

---

## pLM-BLAST (clone + database)

pLM-BLAST does embedding-based remote homology against ECOD30 (default;
ECOD50/70/90 also work). Not on PyPI; clone the upstream repo and point ssign at
`plmblast.py`:

```bash
git clone https://github.com/labstructbioinf/pLM-BLAST.git ~/pLM-BLAST
export SSIGN_PLMBLAST_SCRIPT=~/pLM-BLAST/scripts/plmblast.py
```

Verify the install:

```bash
test -f "$SSIGN_PLMBLAST_SCRIPT" && echo OK || echo NOT FOUND
test -f "$(dirname "$(dirname "$SSIGN_PLMBLAST_SCRIPT")")/embeddings.py" && echo OK || echo NOT FOUND
```

> **HPC / persistent shells:** add the `export SSIGN_PLMBLAST_SCRIPT=...` line to
> `~/.bashrc` (or the equivalent for your shell) AND to any PBS / SLURM batch
> script that runs ssign. The variable does not survive across new login shells
> unless persisted. ssign now hard-fails at pre-flight if both the env var and
> the PATH entry are missing, so a missing install no longer wastes hours of job
> time.

Pre-built ECOD30 database (~10 GB compressed, ~11 GB extracted):

```bash
mkdir -p ~/pLM-BLAST/db && cd ~/pLM-BLAST/db
wget http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs/ecod30db_20240417.tar.gz
tar -xzf ecod30db_20240417.tar.gz && rm ecod30db_20240417.tar.gz
export SSIGN_ECOD_DB=~/pLM-BLAST/db/ECOD30
```

Other cluster levels (ECOD50/70/90) live at the same FTP path with names
`ecod50db_20240417.tar.gz` etc., if you need more redundancy. ECOD30 still has
>=1 representative per F-group so annotation labels are identical across cluster
levels.

GPU strongly recommended: ProtT5 embedding is ~5-10 sec per 500-aa protein on CPU
vs ~0.1 sec on a modern GPU. A whole-genome run on CPU (~5,000 proteins) is 10+
hours just for embedding.

---

## SignalP 6.0

ssign is offline-first: the canonical path is a local SignalP install. **DTU
confirmed on 2026-05-07 that SignalP 6.0 cannot be redistributed**, so each user
acquires it from the DTU portal directly (free academic licence).

If you do not have a DTU licence (or just want to try ssign without one), you can
opt into the DTU webserver fallback instead with `--signalp-mode remote`. The
webserver requires no licence on your part and works on any machine with internet
access. Treat this as a convenience for first-time users and small pilots: the
route depends on DTU continuing to host the service, which is outside our control
and can rate-limit, change, or disappear over time. For anything you intend to
publish or repeat, install locally.

### When the webserver fallback is fine

- Single-genome analyses or few-dozen-genome pilots
- Quick first runs to evaluate ssign before installing DTU tools
- Machines without a DTU licence

### When to install locally (the canonical path)

- Cohorts of >100 genomes where webserver throughput becomes the bottleneck
- Air-gapped environments with no outbound HTTPS to DTU
- Reproducible / paper-ready runs that need every tool fully offline
- Any long-running project where you don't want a third-party webserver on the
  critical path

### Install (Linux / macOS)

The quickest path is `scripts/ssign-setup-dtu <tarball>`, which does every step
below for you and drops SignalP in `~/.conda/envs/signalp6` (a location ssign
auto-detects). The manual recipe follows for reference.

SignalP 6.0 pins **Python <= 3.10** and **PyTorch < 2.0**, while ssign itself
runs on Python 3.11+ with PyTorch 2.x. Installing SignalP into your ssign env
will downgrade PyTorch and break DeepSecE / DeepLocPro / pLM-BLAST. **Install
SignalP into its own env** and point ssign at the binary.

```bash
# 1. Register and request a download at
#    https://services.healthtech.dtu.dk/services/SignalP-6.0/
#    (academic email required). DTU emails / displays a one-time URL.
#    The URL points at an Apache *directory listing*, NOT a file. The
#    directory contains:
#        signalp-6.0_license.txt
#        signalp-6.0i.fast.tar.gz   <-- the ~1.5 GB tarball to wget

# 2. Append the filename to the directory URL and wget on the install
#    machine directly. If you get back a tiny (1-2 KB) HTML file instead
#    of the tarball, your URL is missing the trailing filename.
mkdir -p ~/build && cd ~/build
wget -O signalp6.tar.gz "https://services.healthtech.dtu.dk/download/<your-token>/signalp-6.0i.fast.tar.gz"
ls -lh signalp6.tar.gz        # ~1.5 GB
file signalp6.tar.gz          # "gzip compressed data"

# 3. Create a dedicated Python 3.10 env. Any conda-family tool works:
#    mamba, micromamba, conda. On HPC you typically `module load` one.
mamba create -n signalp6 -c conda-forge python=3.10 pip "numpy<2" -y

# 4. Use the env's binaries by absolute path, avoids needing `mamba init`
#    (which would permanently modify your shell rc). Works identically
#    on a laptop and inside an HPC JupyterHub / batch job.
#
#    The CPU torch wheel is deliberate. We tested swapping in
#    `torch==1.13.1+cu117` on a CUDA 13 / A40 host: SignalP did NOT
#    actually move inference to the GPU because the device is baked
#    into the JIT-compiled `signalp-6-package` model at install time
#    (no runtime `.to(cuda)`). Result: ~10% slower than the CPU wheel
#    from CUDA-lib startup overhead, with zero speedup. Stick with cpu.
PYBIN=~/.conda/envs/signalp6/bin     # adjust if your conda envs live elsewhere
$PYBIN/pip install "torch<2.0" --index-url https://download.pytorch.org/whl/cpu
$PYBIN/python -c "import torch, numpy; print(torch.__version__, numpy.__version__)"
# expected: 1.13.x+cpu  1.26.x

# 5. Extract + install the SignalP package.
tar xzf signalp6.tar.gz       # creates signalp6_fast/ with signalp-6-package inside
cd signalp6_fast
$PYBIN/pip install ./signalp-6-package/

# 6. Copy the model weights into the installed package. `pip install`
#    only installs ~10 MB of Python code; the ~1.4 GB of weights ship
#    separately in the tarball at signalp-6-package/models/.
SIGNALP_DIR=$($PYBIN/python -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
cp -r signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"

# 7. Verify
$PYBIN/signalp6 --version
```

### Wire ssign to the local install

If `signalp6` is on PATH (true after `mamba activate signalp6`), ssign
auto-detects local mode, no flags needed:

```bash
ssign run input.gbff --outdir results
```

If the binary lives outside PATH, point ssign at the install dir via the env var
(preferred, set once per shell) or the CLI flag:

```bash
export SSIGN_SIGNALP_PATH=~/.conda/envs/signalp6/bin    # or:
ssign run input.gbff --outdir results --signalp-path ~/.conda/envs/signalp6/bin
```

ssign invokes SignalP with `--organism other --mode fast --format txt`
(Gram-negatives use the `other` group in v6; v5's `gram-` was removed). Pass
`--signalp-mode remote` to force the DTU webserver even when a local install is
available; `--signalp-mode local` forces local and errors out if no binary is
found.

---

## DeepLocPro

ssign is offline-first, so the canonical path is a local DeepLocPro install
(~5 GB of model files, GPU recommended for cohort speed). Unlike SignalP,
DeepLocPro is not DTU-licensed: it is open source (CC BY-NC-SA 4.0,
non-commercial), cloned directly from the maintainer's public GitHub repository
[Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro) at a pinned
commit. There is no acquisition or licence step; DTU's portal URL is only the web
prediction service, not the source of the local install.

`scripts/ssign-setup-dtu --deeplocpro-only` runs the recipe below for you and
installs into `~/.conda/envs/deeplocpro` (a location ssign auto-detects).

If you don't want to install locally, opt into the DTU webserver fallback with
`--deeplocpro-mode remote` (no install needed on your part, internet required).
Same caveat as SignalP: this is a convenience path that depends on DTU keeping
the service alive, so use it for trial runs and install locally for production /
paper work.

DeepLocPro's only hard pin is Python >= 3.6, much more permissive than SignalP
6.0's torch<2 constraint. It could technically live in the ssign venv, but we
keep it in its own conda-family env for the same reasons we isolate SignalP:
insulate ssign from any transformers / torch constraint DLP might add in a future
release, and keep the "DTU tool != ssign env" convention consistent.

```bash
# 1. Clone the upstream repo
cd ~/build
git clone https://github.com/Jaimomar99/deeplocpro
cd deeplocpro

# 2. Create a dedicated conda-family env. Any version >=3.6 works;
#    3.11 picks up the latest stable wheels for torch / transformers.
mamba create -n deeplocpro -c conda-forge python=3.11 pip -y
PYBIN=~/.conda/envs/deeplocpro/bin

# 3. Install into the env via absolute-path pip (avoids needing
#    `mamba init`, same pattern as the SignalP install above).
$PYBIN/pip install .

# 4. Verify (DeepLocPro has no --version; --help is the canonical check)
$PYBIN/deeplocpro --help | head
ls $PYBIN/deeplocpro          # note this path for the --deeplocpro-path flag below
```

If the GitHub install grows extra steps over time (model-weights download,
license click-through, etc.), follow whatever the README at
[Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro) says. Our
recipe above is the minimum that worked at the v1.0 release.

### Wire ssign to the local install

If `deeplocpro` is on PATH, ssign auto-detects local mode:

```bash
ssign run input.gbff --outdir results
```

Otherwise point at the install dir via env var or flag:

```bash
export SSIGN_DEEPLOCPRO_PATH=~/.conda/envs/deeplocpro/bin    # or:
ssign run input.gbff --outdir results --deeplocpro-path ~/.conda/envs/deeplocpro/bin
```

`--deeplocpro-mode local` / `remote` force the choice when you don't want
auto-detection.

---

## Verifying the install

After adding any tool, ssign's pre-flight check reports what it found:

```bash
ssign run --help                 # confirms ssign itself is on PATH
ssign run input.gbff --outdir /tmp/test --resume   # pre-flight prints all detected tools
```

The pre-flight log lists every external tool plus its detected version (or a
warning if not found). A missing tool is non-fatal: the corresponding step is
skipped at run time and the rest of the pipeline continues.

## Canonical command, extended tier (every annotation tool on)

After `fetch_databases.sh --tier extended` has run and every tool from the
sections above is installed, the canonical command is just:

```bash
ssign run input.gbff --outdir results
```

That's it. `fetch_databases.sh` records the tier at `~/.ssign/tier`; ssign reads
it and enables exactly the tools the extended bundle ships (EggNOG, InterProScan,
pLM-BLAST) while leaving BLASTp and HH-suite off; BLASTp-vs-Swiss-Prot and
HH-suite UniRef30 are full-tier only. The database paths come from
`~/.ssign/db_root` (also written by `fetch_databases.sh`), no per-DB env var
exports needed for the common case.

If your databases live somewhere other than what's recorded in
`~/.ssign/db_root`, set the env vars below to point at them:

```bash
DBROOT=/path/to/your/databases
export BAKTA_DB=$DBROOT/bakta/db-light
export SSIGN_INTERPROSCAN_PATH=$DBROOT/interproscan/interproscan-5.77-108.0
export EGGNOG_DATA_DIR=$DBROOT/eggnog
export SSIGN_ECOD_DB=$DBROOT/plm_blast/ECOD30
# Full-tier-only DBs, uncomment if you've fetched them:
# export SSIGN_HHSUITE_PFAM=$DBROOT/hhsuite/pfam
# export SSIGN_HHSUITE_PDB70=$DBROOT/hhsuite/pdb70
# export SSIGN_HHSUITE_UNICLUST=$DBROOT/hhsuite/uniref30
```

To deviate from the tier default for one tool (e.g. you have NR locally even on
an extended-tier install, or you want to skip EggNOG on this specific run), pass
the per-tool override:

```bash
ssign run input.gbff --outdir results --no-skip-blastp        # opt-in to BLASTp
ssign run input.gbff --outdir results --skip-eggnog           # opt-out of EggNOG
ssign run input.gbff --outdir results --tier base             # force a tier
```

For HPC users who installed each tool into a separate conda env, prepend them to
PATH **before** activating the ssign venv:

```bash
export PATH=~/.conda/envs/bakta/bin:~/.conda/envs/hhsuite/bin:~/.conda/envs/eggnog/bin:$SSIGN_INTERPROSCAN_PATH:$PATH
source ~/.ssign-env/bin/activate     # activate LAST so the venv's python wins
```

Order matters: if you prepend the conda envs *after* activating the venv, one of
their `python` binaries (which lacks ssign's deps, including torch) shadows the
venv's interpreter and `ssign` fails with `ModuleNotFoundError`.
