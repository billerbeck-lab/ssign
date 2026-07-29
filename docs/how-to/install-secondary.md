# Secondary install options (native / pip)

Install ssign from the cloned repo with `pip`, plus
a handful of optional tools. The trade-off is that you install and maintain the
external toolchain yourself (Bakta, EggNOG-mapper, BLAST+, HH-suite, InterProScan,
SignalP 6, DeepLocPro), and those tool versions can drift out of compatibility
over time. For reproducible runs use the recommended [container](install.md#1-container-install-recommended) instead.

Clone the repo and install the tier you want:

```bash
git clone https://github.com/billerbeck-lab/ssign && cd ssign
pip install -e .                 # base
pip install -e '.[extended]'     # base + Bakta + pLM-BLAST + extended-tier pins
pip install -e '.[full]'         # extended deps + full database tier
```

After pip install, fetch the matching database bundle:

```bash
bash scripts/fetch_databases.sh --tier base       # or: extended / full
```

### Database tier sizes

| Tier | Cumulative size | Adds |
|---|---|---|
| base | ~4 GB | Bakta light DB (3.9 GB) |
| extended | ~100 GB | + EggNOG (47 GB), InterProScan (24 GB), pLM-BLAST ECOD30 (11 GB) |
| full | ~500 GB | + Bakta full DB (84 GB, replaces light), HH-suite Pfam + PDB70 + UniRef30 (340 GB total), BLASTp-vs-Swiss-Prot (0.3 GB) |

These are databases only. Model weights are separate and shared across all tiers,
~18 GB total.

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
downloads via `run_deepsece._ensure_checkpoint()`. The two weights you provide yourself are SignalP 6 (DTU tarball) and DeepLocPro
(its GitHub clone).

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
you set them), and model weights.

For environment variables (mirror URLs, database paths, dev-only flags), see
[`reference/env_vars.md`](../reference/env_vars.md).

---

## DeepSecE

DeepSecE predicts secretion-system type per protein and runs
alongside DeepLocPro.

CUDA GPU strongly recommended; CPU runs are slow.

---

## Bakta (pip extra + database)

Bakta provides annotation-grade gene calling and functional descriptions. ssign
re-annotates inputs with Bakta by default; if Bakta is not installed, the
pipeline falls back to a pyrodigal-only call without functional annotations.

**Important:** installing the `[bakta]` extra only installs the Bakta Python wrapper.
Bakta also depends on several binary tools (AMRFinderPlus, DIAMOND, HMMER,
tRNAscan-SE, aragorn) that aren't pip-installable. Without them, **even
`bakta_db download` fails** because Bakta's startup runs the same dependency
check as `bakta` itself.

Install all of Bakta from conda (tRNAscan-SE, aragorn, PilerCR, AMRFinderPlus,
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
for substrate proteins.

```bash
conda install -c bioconda eggnog-mapper
```

If you don't use conda, `--no-deps` skips the upstream pin and lets ssign's
biopython (`>=1.80`) satisfy eggnog-mapper at runtime.

```bash
pip install --no-deps eggnog-mapper
```

Verify the install:

```bash
which emapper.py
emapper.py --version
```

Fetch the database (~47 GB extracted):

```bash
scripts/fetch_databases.sh --tier extended --target ~/.ssign/databases
```

Tell ssign where the database lives:

```bash
ssign run input.gbff --outdir results --eggnog-db ~/.ssign/databases/eggnog
```

EggNOG annotation is on by default at the extended (default) and full tiers; it
is off only at `--tier base`. Pass `--skip-eggnog` to disable it.

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
Tübingen** and fall back to GWDG if Tübingen is unreachable.

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

Pre-built ECOD30 database (~10 GB compressed, ~11 GB extracted):

```bash
mkdir -p ~/pLM-BLAST/db && cd ~/pLM-BLAST/db
wget http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs/ecod30db_20240417.tar.gz
tar -xzf ecod30db_20240417.tar.gz && rm ecod30db_20240417.tar.gz
export SSIGN_ECOD_DB=~/pLM-BLAST/db/ECOD30
```

GPU strongly recommended

---

## SignalP 6.0

**DTU confirmed on 2026-05-07 that SignalP 6.0 cannot be redistributed**, so each user
acquires it from the DTU portal directly (free academic licence).

If you do not have a DTU licence (or just want to try ssign without one), you can
opt into the DTU webserver fallback instead with `--signalp-mode remote`.

### Install (Linux / macOS)

The quickest path is `scripts/ssign-setup-dtu <tarball>`, which does every step
below for you and drops SignalP in `~/.conda/envs/signalp6` (a location ssign
auto-detects). The manual recipe follows for reference.

SignalP 6.0 pins **Python <= 3.10** and **PyTorch < 2.0**, while ssign itself
runs on Python 3.10+ with PyTorch 2.x. Installing SignalP into your ssign env
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
#This path is lighter-weight: 
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

---

## DeepLocPro

Unlike SignalP, DeepLocPro is not DTU-licensed: it is open source (CC BY-NC-SA 4.0,
non-commercial), cloned directly from the maintainer's public GitHub repository
[Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro) at a pinned
commit.

`scripts/ssign-setup-dtu --deeplocpro-only` runs the recipe below for you and
installs into `~/.conda/envs/deeplocpro` (a location ssign auto-detects).

If you don't want to install locally, opt into the DTU webserver fallback with
`--deeplocpro-mode remote`.

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

## Extended tier (every tool on)

```bash
ssign run input.gbff --outdir results
```

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

To deviate from the tier default for one tool pass the per-tool override:

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
