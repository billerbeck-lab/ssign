# Troubleshooting

Common errors and what to do about them. If your error is not listed,
[open a GitHub issue](https://github.com/billerbeck-lab/ssign/issues)
with the full log and the command you ran.

## Install and dependency errors

### `command not found: ssign`

Your venv is not active. Run:

```bash
source ~/.ssign-env/bin/activate
ssign --version
```

If the venv exists but `ssign --version` still fails, the install
finished in a different env. Reinstall in the current one:
`pip install ssign`.

### `ModuleNotFoundError: No module named 'deepsece'` (or `bakta`, `plmblast`, etc.)

You are running a step whose pip extra is not installed. Either install
the matching extra or skip the step:

```bash
pip install ssign[deepsece]            # add DeepSecE
# or
ssign run input.gbff --skip-deepsece   # skip it
```

Same pattern for `bakta`, `plmblast`, `plm-effector`.

### `emapper.py: command not found` (EggNOG-mapper)

EggNOG-mapper is **not** in the ssign pip extras: its hard-pinned
`biopython==1.76` clashes with `bakta>=1.78`, so the two cannot
co-resolve. Install separately:

```bash
conda install -c bioconda eggnog-mapper      # recommended
# or
pip install --no-deps eggnog-mapper          # no conda available
```

Then run `bash scripts/fetch_databases.sh --tier extended` to download
the EggNOG database, or skip the step with `--skip-eggnog`. See
[`how-to/install.md`](install.md#eggnog-mapper-separate-install--database)
for the full recipe.

### `transformers >=5.0` or `numpy >=2.0` import errors during pLM-BLAST or PLM-Effector

Upstream removed APIs that pLM-BLAST and PLM-Effector still depend on.
Pin to compatible versions:

```bash
pip install "transformers>=4.38,<5.0" "numpy>=1.26,<2.0"
```

These bounds are already in `pip install ssign[extended]`. If they got
overridden by a later install, reinstall `ssign[extended]`.

### SignalP 6.0 install downgrades PyTorch and breaks DeepSecE

SignalP pins PyTorch < 2.0 while every other ssign tool needs PyTorch
2.x. **Install SignalP into its own conda or venv environment**, then
point ssign at the binary via `--signalp-path`. Detailed steps in
[`how-to/install.md`](install.md#signalp-60-local-install-optional).

## Database and path errors

### `RuntimeError: hhblits not found on PATH`

You enabled HH-suite (`--no-skip-hhsuite`) but HH-suite is not installed.

```bash
conda install -c bioconda hhsuite
```

If `which hhblits` finds the binary but the error persists, the binary
may not be on PATH inside your job script. Add the conda env activation
to the script.

### `FileNotFoundError: <db-path>` for HH-suite, EggNOG, BLAST, etc.

The CLI flag (or its `SSIGN_*` env var fallback) points at a path that
does not exist. Check:

```bash
ls $SSIGN_HHSUITE_PFAM      # should list a directory of pfam_a3m_*.* files
ls $SSIGN_EGGNOG_DB         # should list eggnog.db, eggnog_proteins.dmnd
```

If the path is correct but the listing is empty, the database download
was incomplete. Re-run `bash scripts/fetch_databases.sh --tier <tier>`.

### Bakta complains "tRNAscan-SE not found"

`tRNAscan-SE` is a Bakta dependency that bioconda packages alongside
Bakta. If you `pip install bakta` without going through conda, you may
miss it. Reinstall via conda:

```bash
conda install -c bioconda bakta
```

Or skip Bakta and trust the input GenBank's annotations:

```bash
ssign run input.gbff --use-input-annotations
```

## DTU webserver errors (DeepLocPro, SignalP)

### `DTU server returned HTTP 503` or `HTTP 504`

The DTU webserver is overloaded or briefly down. Wait a few minutes and
retry; ssign retries with exponential backoff already.

If it keeps failing, check
[https://services.healthtech.dtu.dk](https://services.healthtech.dtu.dk)
for status. As a fallback, run with `--skip-signalp` and `--skip-deeplocpro`,
or install both locally (each requires a free DTU academic licence;
see [install.md](install.md)).

### `DTU job <id> timed out`

The job exceeded the polling window (default ~1 hour for large
submissions). For genomes with thousands of proteins, the DTU
webserver can queue behind other users for >30 minutes per batch. Two
options:

- Re-run with `--resume` so finished steps do not repeat.
- Install DeepLocPro and SignalP locally for offline runs that bypass
  the queue entirely.

### `DTU job failed after <N>s`

The DTU server reported a job failure. Most often: malformed input
(non-AA characters in the protein FASTA, or a sequence > 70,000 aa).
Check the input proteins for unusual records. If it persists with
clean input, it is a DTU-side issue; report to DTU and retry later.

## Resume and partial runs

### `ssign run --resume` re-runs a step that already finished

The progress manifest at `<outdir>/.ssign/<sample-id>_progress.json`
has an entry for the step but the output file the step wrote is
missing or empty. Resume validates outputs exist before skipping; it
re-runs missing ones. Causes:

- The temp work directory was cleaned. ssign keeps `work_dir` only on
  failure; if the previous run technically succeeded but you killed it,
  intermediates may be gone. Re-run from scratch.
- The output is on a network filesystem that lagged. Wait a few seconds
  and retry.

### `--resume` does not pick up where I expected

Verify the manifest matches your config:

```bash
cat ecoli_results/.ssign/ecoli_k12_progress.json | python -m json.tool | head
```

The manifest stores the config used in the previous run. If you changed
flags between runs, ssign treats the changed config as a different run
and starts over.

## GPU and memory errors

### `MemoryError` during DeepSecE checkpoint loading

DeepSecE wraps a ~7.3 GB ESM-1b model. Free up memory or use a machine
with at least 12 GB RAM. CPU is supported but slow; on a 16 GB laptop,
close other heavy applications first.

### `RuntimeError: CUDA out of memory` during PLM-Effector or pLM-BLAST

Your GPU does not have enough VRAM for the protein-language models.
PLM-Effector wants ~16 GB VRAM; pLM-BLAST is similar. Options:

- Run on CPU (slow but works for pLM-BLAST; PLM-Effector hard-skips on
  CPU because runtime is impractical).
- Skip the step: `--skip-plm-effector` or `--skip-plmblast`.
- Run on a GPU node with more VRAM (HPC).

### `nvidia-smi` works but ssign does not see the GPU

PyTorch was installed in a CPU-only configuration. Reinstall with the
CUDA build:

```bash
pip uninstall torch
pip install torch --index-url https://download.pytorch.org/whl/cu121
```

(Replace `cu121` with the CUDA version that matches your driver; check
with `nvidia-smi`.)

## Output and parsing errors

### `master_substrates.csv` is empty (only the header)

ssign found no proteins that pass the secretion-prediction + proximity
gates. Common causes:

- The genome has no detected secretion systems (check chunk 3 of
  `<sample-id>_results.csv`). MacSyFinder requires a complete machinery
  signature; partial systems below `--wholeness-threshold` (default
  0.8) are not scored.
- T3SS-only genome and you ran with the default `--excluded-systems
  Flagellum Tad T3SS`. Override the exclusion: `--excluded-systems
  Flagellum Tad`.
- Too aggressive `--conf-threshold` (default 0.8). Lower it to 0.5 to
  see weaker DLP hits and decide if any look like real substrates.

### Cohort run: substrate counts vary across genomes despite identical input

Stochastic tools (the DTU webserver in particular) can produce
slightly different probability scores between runs. ssign's pipeline
itself is deterministic; if the variation is large, it is a
DTU-webserver effect. For reproducible runs, install DLP and SignalP
locally.

## Cohort and multi-genome errors

### `ortholog_groups.csv` is empty after a multi-genome run

The cross-genome ortholog step requires BLAST+ (for all-vs-all BLASTp).
Install it:

```bash
sudo apt install ncbi-blast+      # Debian/Ubuntu
brew install blast                 # macOS
conda install -c bioconda blast    # cross-platform
```

If BLAST+ is installed but ortholog grouping still produces no rows,
check that each genome's `<sample-id>_results.csv` has substrates to
group. Empty inputs produce empty groups.

### One genome in a cohort kills the whole batch

Each genome's pipeline is independent; a single-genome failure should
not stop the others. If it does, you are likely running with `set -e`
in a shell loop. Either drop `set -e` for the loop, or wrap each
invocation:

```bash
for g in genomes/*.gbff; do
    ssign run "$g" --outdir results/$(basename "$g" .gbff) || \
        echo "WARN: $g failed; continuing"
done
```

## Asking for help

If none of the above matches, open a GitHub issue with:

- Full ssign log (the stderr output from your `ssign run` command).
- The exact command you ran.
- Your environment: `pip freeze`, `python --version`, OS and version.
- Whether the failure is reproducible.

Issues with all four are usually solvable in one round.
