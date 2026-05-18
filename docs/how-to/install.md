# Installing optional tools

`pip install ssign` gets you the base pipeline: secretion-system detection,
secreted-protein prediction (DeepLocPro, SignalP, and PLM-Effector — see
the DeepLocPro and SignalP sections below for installing the DTU tools
locally, or for the opt-in webserver fallback when you don't have a DTU
licence), proximity analysis, and reporting. The tools below extend what
the pipeline reports and how it runs.

The simplest path is one of the named tiers:

```bash
pip install ssign                 # base
pip install ssign[extended]       # base + DeepSecE + Bakta + pLM-BLAST + extended-tier pins
pip install ssign[full]           # extended deps + full database tier
```

After pip install, fetch the matching database bundle:

```bash
bash scripts/fetch_databases.sh --tier base       # or: extended / full
```

If a tier doesn't fit your storage budget, pick individual tools below.

For environment variables (mirror URLs, database paths, dev-only flags), see
[`reference/env_vars.md`](../reference/env_vars.md).

---

## DeepSecE (pip extra)

DeepSecE predicts secretion-system type per protein. Runs as a second opinion
alongside DeepLocPro. About 7.3 GB once installed (PyTorch plus the ESM
protein language model).

```bash
source ~/.ssign-env/bin/activate
pip install ssign[deepsece]
```

ssign auto-detects the install on the next launch. CUDA GPU strongly
recommended; CPU runs are slow.

---

## Bakta (pip extra + database)

Bakta provides annotation-grade gene calling and functional descriptions.
ssign re-annotates inputs with Bakta by default; if Bakta is not installed,
the pipeline falls back to a pyrodigal-only call without functional
annotations. Install + download the light database (~2 GB):

```bash
source ~/.ssign-env/bin/activate
pip install ssign[bakta]
bakta_db download --output ~/bakta_db --type light
```

Pass `--bakta-db ~/bakta_db` on the command line, or set
`SSIGN_BAKTA_DB=~/bakta_db` (read by `scripts/fetch_databases.sh`).

The full Bakta database (~30 GB) is the `--type full` variant. The full
tier in `fetch_databases.sh --tier full` pulls it.

---

## EggNOG-mapper (separate install + database)

EggNOG-mapper provides ortholog-based functional annotation (COG, KEGG, GO,
PFAM) for substrate proteins. It is invoked as a subprocess (`emapper.py`),
not imported by ssign, and is **not** included in the `[extended]` /
`[full]` pip extras: upstream eggnog-mapper hard-pins `biopython==1.76`
while ssign and Bakta need `biopython>=1.78`, which makes the two
unsatisfiable in a single pip resolution. Install it separately.

The conda path is recommended; bioconda has shipped eggnog-mapper against
modern biopython for years without breakage:

```bash
conda install -c bioconda eggnog-mapper
```

If you don't use conda, `--no-deps` skips the upstream pin and lets ssign's
biopython (`>=1.80`) satisfy eggnog-mapper at runtime. The only API
incompatibility (`Bio.Alphabet`, removed in biopython 1.78) is already
guarded with `try/except` in eggnog-mapper itself, so this works in
practice; bioconda relies on the same:

```bash
pip install --no-deps eggnog-mapper
```

After installing eggnog-mapper, fetch the database (~50 GB):

```bash
download_eggnog_data.py -y --data_dir ~/.ssign/databases/eggnog
```

`scripts/fetch_databases.sh --tier {extended,full}` runs that step for you
once `download_eggnog_data.py` is on PATH. Tell ssign where the database
lives:

```bash
ssign run input.gbff --outdir results --eggnog-db ~/.ssign/databases/eggnog
```

EggNOG annotation is off by default (`--skip-eggnog` defaults to `true`).
Pass `--no-skip-eggnog` to enable it.

---

## BLAST+ (system binary)

`blastp` is needed for cross-genome ortholog grouping. Download is about
200 MB.

```bash
sudo apt install ncbi-blast+        # Debian / Ubuntu
brew install blast                  # macOS
conda install -c bioconda blast     # cross-platform
```

ssign auto-detects `blastp` on the next launch.

---

## InterProScan

InterProScan (EBI) scans proteins against a panel of member databases
(Pfam, TIGRFAM, HAMAP, SMART, PIRSF, SUPERFAMILY, Gene3D, ProSite, CDD)
to annotate domains, family memberships, and GO terms. ssign uses it to
add domain-level annotation to substrate proteins. Java required.

The bundle is large (~24 GB extracted) but installs as a single tarball:

```bash
# 1. Java 11+ on PATH (Ubuntu / Debian):
sudo apt install openjdk-17-jre-headless

# 2. Download the latest InterProScan release (replace 5.74-105.0 with
#    the current version listed at https://www.ebi.ac.uk/interpro/download/):
mkdir -p ~/interproscan && cd ~/interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.74-105.0/interproscan-5.74-105.0-64-bit.tar.gz
tar -xzf interproscan-5.74-105.0-64-bit.tar.gz
cd interproscan-5.74-105.0

# 3. Initialise HMMs + indexes (writes ~1 GB of derived files):
python3 setup.py interproscan.properties
```

Point ssign at the install directory (the one containing
`interproscan.sh`):

```bash
export SSIGN_INTERPROSCAN_PATH=~/interproscan/interproscan-5.74-105.0
# or pass --interproscan-db ~/interproscan/interproscan-5.74-105.0
```

ssign runs InterProScan with the bacterial-relevant member DBs by default
(PANTHER, the slowest member and eukaryote-leaning, is excluded). Per-
protein scan time is typically 5-30 s; a whole-genome run on ~5,000
proteins is 30-90 minutes. The first run also queries EBI's precalculated-
match lookup service for a 5-10x speedup on known sequences; add `-dp`
behaviour via your own wrapper script if you need air-gapped operation.

---

## HH-suite (system binary + databases)

`hhsearch` and `hhblits` (Steinegger / Söding labs) detect remote evolutionary
relationships by comparing HMM-vs-HMM profiles. ssign uses HH-suite to
annotate substrate proteins with Pfam domain families and PDB structural
homologs.

Install the binaries:

```bash
conda install -c bioconda hhsuite      # recommended
# or build from source: https://github.com/soedinglab/hh-suite
```

Download the three databases. Two canonical mirrors host them; **prefer
Tübingen** (fresher, recommended by Söding lab issue #382) and fall back to
GWDG only if Tübingen is unreachable.

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

Total disk after extraction: ~55 GB. Tell ssign where to find them:

```bash
export SSIGN_HHSUITE_PFAM=$HHSUITE_DBS/pfam
export SSIGN_HHSUITE_PDB70=$HHSUITE_DBS/pdb70_from_mmcif_2026-02-20
export SSIGN_HHSUITE_UNICLUST=$HHSUITE_DBS/UniRef30_2023_02
```

These three are read at run time as fallbacks for the matching CLI flags
(`--hhsuite-pfam-db`, `--hhsuite-pdb70-db`, `--hhsuite-uniclust-db`).

### Mirror caveats

The Söding lab is in maintenance mode; per [hh-suite issue #382](https://github.com/soedinglab/hh-suite/issues/382),
no funding for support, but the files at GWDG are stable. Tübingen has fresher
builds. If you see a 404 on the GWDG host, probe with `curl -IL` (capital L
follows the 302 redirect to the sibling domain). If both mirrors are down
for an extended period, pre-stage the DBs on your HPC scratch directory.

---

## pLM-BLAST (clone + database)

pLM-BLAST does embedding-based remote homology against ECOD70. Not on PyPI;
clone the upstream repo and point ssign at `plmblast.py`:

```bash
git clone https://github.com/labstructbioinf/pLM-BLAST.git ~/pLM-BLAST
export SSIGN_PLMBLAST_SCRIPT=~/pLM-BLAST/scripts/plmblast.py
```

Pre-built ECOD70 database (~21 GB compressed, ~24 GB extracted):

```bash
mkdir -p ~/pLM-BLAST/db && cd ~/pLM-BLAST/db
wget http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs/ecod70db_20240417.tar.gz
tar -xzf ecod70db_20240417.tar.gz && rm ecod70db_20240417.tar.gz
export SSIGN_ECOD70_DB=~/pLM-BLAST/db/ECOD70
```

GPU strongly recommended: ProtT5 embedding is ~5-10 sec per 500-aa protein on
CPU vs ~0.1 sec on a modern GPU. A whole-genome run on CPU (~5,000 proteins)
is 10+ hours just for embedding.

---

## PLM-Effector (vendored; weights download only)

PLM-Effector source code is shipped with ssign under `src/ssign_app/scripts/plm_effector/`.
Only the trained weights and pretrained protein language models need
downloading.

```bash
# Trained models (~1.7 GB)
wget https://www.mgc.ac.cn/PLM-Effector/download/sourcecode.zip
unzip sourcecode.zip
mv sourcecode/trained_models /path/to/plm_effector_weights/

# Pretrained PLMs from HuggingFace (~17 GB total)
hf download Rostlab/prot_t5_xl_uniref50 \
    --include "config.json" "spiece.model" "tokenizer_config.json" \
    "special_tokens_map.json" "pytorch_model.bin" \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/prot_t5_xl_uniref50
hf download facebook/esm1b_t33_650M_UR50S \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/esm1b_t33_650M_UR50S
hf download facebook/esm2_t33_650M_UR50D \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/esm2_t33_650M_UR50D
hf download Rostlab/prot_bert \
    --local-dir /path/to/plm_effector_weights/transformers_pretrained/prot_bert

export SSIGN_PLM_EFFECTOR_WEIGHTS=/path/to/plm_effector_weights
```

CUDA GPU required for practical runtime; CPU is ~100x slower. ssign's
PLM-Effector test hard-skips on systems without a GPU.

---

## SignalP 6.0

ssign is offline-first: the canonical path is a local SignalP install. **DTU
confirmed on 2026-05-07 that SignalP 6.0 cannot be redistributed**, so each
user acquires it from the DTU portal directly (free academic licence).

If you do not have a DTU licence (or just want to try ssign without one),
you can opt into the DTU webserver fallback instead with
`--signalp-mode remote`. The webserver requires no licence on your part
and works on any machine with internet access. Treat this as a
convenience for first-time users and small pilots: the route depends on
DTU continuing to host the service, which is outside our control and
can rate-limit, change, or disappear over time. For anything you intend
to publish or repeat, install locally.

### When the webserver fallback is fine

- Single-genome analyses or few-dozen-genome pilots
- Quick first runs to evaluate ssign before installing DTU tools
- Machines without a DTU licence

### When to install locally (the canonical path)

- Cohorts of >100 genomes where webserver throughput becomes the bottleneck
- Air-gapped environments with no outbound HTTPS to DTU
- Reproducible / paper-ready runs that need every tool fully offline
- Any long-running project where you don't want a third-party webserver
  on the critical path

### Install (Linux / macOS)

SignalP 6.0 pins **Python ≤ 3.10** and **PyTorch < 2.0**, while ssign itself
runs on Python 3.11+ with PyTorch 2.x. Installing SignalP into your ssign env
will downgrade PyTorch and break DeepSecE / DeepLocPro / PLM-Effector /
pLM-BLAST. **Install SignalP into its own env** and point ssign at the
binary.

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

# 4. Use the env's binaries by absolute path — avoids needing `mamba init`
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

```bash
ssign run input.gbff --outdir results \
    --signalp-mode local \
    --signalp-path ~/.conda/envs/signalp6/bin
```

Pass the directory containing the `signalp6` console script (typically
`<env>/bin`) via `--signalp-path`, or set it in the GUI's local-install
field. ssign invokes SignalP with `--organism other --mode fast --format txt`
(Gram-negatives use the `other` group in v6; v5's `gram-` was removed).

---

## DeepLocPro

ssign is offline-first, so the canonical path is a local DeepLocPro
install (~5 GB of model files, GPU recommended for cohort speed).
Unlike SignalP, DeepLocPro is not distributed through DTU's download
portal — that URL is only the web prediction service. The local install
is from the maintainer's GitHub repository
[Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro).

If you don't want to install locally, opt into the DTU webserver
fallback with `--deeplocpro-mode remote` (no install needed on your
part, internet required). Same caveat as SignalP: this is a convenience
path that depends on DTU keeping the service alive, so use it for trial
runs and install locally for production / paper work.

DeepLocPro's only hard pin is Python ≥ 3.6 — much more permissive than
SignalP 6.0's torch<2 constraint. It could technically live in the ssign
venv, but we keep it in its own conda-family env for the same reasons
we isolate SignalP: insulate ssign from any transformers / torch
constraint DLP might add in a future release, and keep the "DTU tool ≠
ssign env" convention consistent.

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

If the GitHub install grows extra steps over time (model-weights
download, license click-through, etc.), follow whatever the README at
[Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro) says
— our recipe above is the minimum that worked at the v1.0 release.

### Wire ssign to the local install

```bash
ssign run input.gbff --outdir results \
    --deeplocpro-mode local \
    --deeplocpro-path ~/.conda/envs/deeplocpro/bin
```

`--deeplocpro-path` takes the directory containing the `deeplocpro`
console script.

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
