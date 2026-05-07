# Installing optional tools

`pip install ssign` gets you the base pipeline: secretion-system detection,
secreted-protein prediction (DeepLocPro and SignalP via DTU's webserver,
PLM-Effector locally), proximity analysis, and reporting. The tools below
extend what the pipeline reports and how it runs.

The simplest path is one of the named tiers:

```bash
pip install ssign                 # base
pip install ssign[extended]       # base + DeepSecE + Bakta + EggNOG + pLM-BLAST + extended-tier pins
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

## SignalP 6.0 local install (optional)

ssign uses DTU's free SignalP webserver by default. The webserver requires no
licence on your part; it is the lowest-friction path and works out of the
box on any machine with internet access. Most users do not need anything
beyond the default.

A local install is faster (no upload latency, parallelisable) and works
offline, but requires a DTU academic licence. **DTU confirmed on 2026-05-07
that SignalP 6.0 cannot be redistributed.** Each user has to acquire it from
the DTU portal directly.

### When to install locally

- Cohorts of >100 genomes where webserver throughput becomes the bottleneck
- Air-gapped environments with no outbound HTTPS to DTU
- Reproducibility-paranoid runs that want every tool fully offline

For a single-genome analysis or a few-dozen-genome pilot, the default
webserver mode is fine.

### Install (Linux / macOS)

SignalP 6.0 pins **Python ≤ 3.10** and **PyTorch < 2.0**, while ssign itself
runs on Python 3.11+ with PyTorch 2.x. Installing SignalP into your ssign env
will downgrade PyTorch and break DeepSecE / DeepLocPro / PLM-Effector /
pLM-BLAST. **Install SignalP into its own env** and point ssign at the
binary.

```bash
# 1. Register and download from DTU (academic email required):
#    https://services.healthtech.dtu.dk/services/SignalP-6.0/
#    Pick the "fast" model variant (signalp-6.0i.fast.tar.gz, ~1.5 GB).

# 2. Create a dedicated Python 3.10 env (micromamba example):
micromamba create -n signalp6 -c conda-forge python=3.10 pip "numpy<2" -y
micromamba activate signalp6

# 3. Pre-install CPU-only torch<2 BEFORE the package, so torch 2.x isn't
#    pulled in transitively (also keeps the wheel small at ~150 MB):
pip install "torch<2.0" --index-url https://download.pytorch.org/whl/cpu

# 4. Extract and install:
tar xzf ~/Downloads/signalp-6.0i.fast.tar.gz -C ~/build
cd ~/build/signalp6_fast
pip install ./signalp-6-package/

# 5. Copy the model weights into the installed package directory:
SIGNALP_DIR=$(python -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
cp -r signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"

# 6. Verify:
signalp6 --version
```

### Wire ssign to the local install

```bash
ssign run input.gbff --outdir results \
    --signalp-mode local \
    --signalp-path ~/micromamba/envs/signalp6/bin
```

Pass the directory containing the `signalp6` console script (typically
`<env>/bin`) via `--signalp-path`, or set it in the GUI's local-install
field. ssign invokes SignalP with `--organism other --mode fast --format txt`
(Gram-negatives use the `other` group in v6; v5's `gram-` was removed).

---

## DeepLocPro local install (optional)

ssign uses DTU's free DeepLocPro webserver by default. Same shape as SignalP:
the webserver requires no licence on your part and is the lowest-friction
path.

A local install is faster but requires a DTU academic licence (~5 GB of model
files; GPU recommended). DTU's licence terms for redistribution are pending
clarification with Ole, the DeepLocPro maintainer (status as of 2026-05-08).
For now, treat as user-acquires-it.

```bash
# 1. Register at https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/
# 2. Download and install per DTU's instructions
# 3. Wire ssign to the install:
ssign run input.gbff --outdir results \
    --deeplocpro-mode local \
    --deeplocpro-path /path/to/deeplocpro
```

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
