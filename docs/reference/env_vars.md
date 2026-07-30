# Environment variables

ssign reads a small set of environment variables. Most are convenience
aliases that `scripts/fetch_databases.sh` exports for shell rc files; only a
handful are read by the runtime.

## Read at run time

These are the env vars ssign actually consults during a pipeline run:

| Variable | Purpose |
|---|---|
| `SSIGN_HHSUITE_PFAM` | Path to HH-suite Pfam database. Fallback for `--hhsuite-pfam-db`. |
| `SSIGN_HHSUITE_PDB70` | Same, for PDB70. Fallback for `--hhsuite-pdb70-db`. |
| `SSIGN_HHSUITE_UNICLUST` | Same, for UniRef30 / UniClust30. Fallback for `--hhsuite-uniclust-db`. |
| `SSIGN_BLAST_SWISSPROT` | Directory holding the BLAST Swiss-Prot database (contains `swissprot.pdb`), the full-tier default BLASTp DB. When `--blastp-db` is unset, ssign resolves the `-db` prefix as `<this dir>/swissprot`. Fallback for `--blastp-db`; the flag wins. NR is opt-in: fetch it and pass `--blastp-db <nr-dir>/nr` (there is no auto-resolution for NR). |
| `SSIGN_DEEPSECE_CHECKPOINT_URL` | Replace the canonical DeepSecE checkpoint URL with an institutional mirror. Useful inside firewalled networks. |
| `SSIGN_PLMBLAST_SCRIPT` | Path to the upstream `plmblast.py` script (clone of `labstructbioinf/pLM-BLAST`). |

CLI flags always take precedence: if both `--hhsuite-pfam-db` and
`SSIGN_HHSUITE_PFAM` are set, the CLI flag wins.

## Runtime tuning / hardware

These control device selection and resource sizing. None have a CLI flag;
they exist for HPC and mixed-hardware environments.

| Variable | Purpose |
|---|---|
| `SSIGN_DEEPLOCPRO_FORCE_CPU` | Set to `1`/`true`/`yes` to force DeepLocPro onto CPU even when a CUDA GPU is visible. Default: auto-detect (GPU if torch sees one, else CPU). |
| `SSIGN_DEEPSECE_FORCE_CPU` | Same, for DeepSecE's ESM-1b step. Default: auto-detect. |
| `SSIGN_PLMBLAST_FORCE_CPU` | Same, for pLM-BLAST's ProtT5 embedding (`--cuda`). Default: auto-detect. CPU embedding is ~100x slower. |
| `SSIGN_GPU_SERIALIZE` | Force the GPU-bound steps (DeepLocPro/DeepSecE/pLM-BLAST) to run one at a time (`1`) or fully parallel (`0`). Default: auto-detect from the GPU compute mode — serialize on an exclusive-process GPU (it rejects a second CUDA context with `cudaErrorDevicesUnavailable`), stay parallel on a shareable Default-mode GPU. Set `1` if you ever see that CUDA error. |
| `SSIGN_MAX_RAM_GB` | Override the detected RAM budget (GB). ssign otherwise takes the minimum of the SLURM/PBS allocation, cgroup limit, and host total; set this on an HPC where none of those report the real cap. |
| `SSIGN_PARALLEL_GROUP_SIZE` | Divisor `N` (>= 2) each tool wrapper uses to self-limit its CPU/RAM share when the runner launches prediction/annotation tools concurrently. Set by the runner automatically; rarely set by hand. Unset or <= 1 means "standalone, use full budget". |
| `SSIGN_DEEPLOCPRO_MAX_AA` | Max protein length (aa) DeepLocPro will embed; longer sequences are withheld and marked "not predicted" to avoid O(L^2) attention OOM. Default 5000. |

## Database paths read as CLI-flag fallbacks

Each of these is read at run time as a fallback used only when the matching CLI
flag is unset, resolving its tool's database or install
directory. `scripts/fetch_databases.sh` prints a matching `set ...` line for the
databases it downloads so you can copy a one-liner into your shell rc file.
Bakta and EggNOG are read under both their **native** tool env var names
(`BAKTA_DB`, `EGGNOG_DATA_DIR`, the same names the fetcher prints) and a
`SSIGN_`-prefixed alias (`SSIGN_BAKTA_DB`, `SSIGN_EGGNOG_DB`); either works.

| Variable | CLI equivalent |
|---|---|
| `BAKTA_DB` / `SSIGN_BAKTA_DB` | `--bakta-db` |
| `EGGNOG_DATA_DIR` / `SSIGN_EGGNOG_DB` | `--eggnog-db` |
| `SSIGN_INTERPROSCAN_PATH` | `--interproscan-db` |
| `SSIGN_ECOD_DB` | `--plmblast-db` |
| `SSIGN_DEEPLOCPRO_PATH` | `--deeplocpro-path` (also used by the DLP integration test) |
| `SSIGN_SIGNALP_PATH` | `--signalp-path` (also used by the SignalP integration test) |

## Python dependency pins (extended tier)

These are not env vars but version constraints. Listed here so a maintainer
debugging a fresh install knows where the upper bounds come from. They are
captured in `[project.optional-dependencies].extended` in `pyproject.toml`,
so `pip install -e '.[extended]'` resolves them automatically.

| Package | Pin | Reason |
|---|---|---|
| `transformers` | `>=4.38,<5.0` | 5.0 removed `batch_encode_plus`, used by pLM-BLAST tokenizers. |
| `numpy` | `>=1.26,<2.0` | 2.0 removed `np.issubsctype`, used by pLM-BLAST's alignment code. |
| `protobuf` | any | Required by ProtT5's SentencePiece tokenizer at load time. |
| `mkl`, `mkl-service` | any | pLM-BLAST's `plmblast.py` imports them directly. |

The `transformers` and `numpy` upper bounds are revisited in the v1.x roadmap
once upstream pLM-BLAST publishes 5.0/2.0-compatible code.

## Test and developer-only

| Variable | Purpose |
|---|---|
| `SSIGN_TEST_FIXTURE_FULL` | Set to `1` to use the full Xanthobacter contig fixture (~1 Mb) instead of the minimal T5aSS fixture (~213 kb). Slower; closer to a real run. |
| `SSIGN_TEST_OUTDIR` | Output directory for the multi-genome integration test. |
| `SSIGN_RUN_FULL_PIPELINE` | Set to `1` to opt into the long-running fixture pipeline test. Skipped by default. |
| `SSIGN_GOLDEN_REGEN_DIR` | Where regenerated golden-output files are written when refreshing `tests/fixtures/golden/`. See `tests/fixtures/golden/REGENERATE.md`. |
