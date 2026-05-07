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
| `SSIGN_DEEPSECE_CHECKPOINT_URL` | Replace the canonical DeepSecE checkpoint URL with an institutional mirror. Useful inside firewalled networks. |
| `SSIGN_PLMBLAST_SCRIPT` | Path to the upstream `plmblast.py` script (clone of `labstructbioinf/pLM-BLAST`). |
| `SSIGN_TAXDUMP_DIR` | NCBI taxdump directory used by `resolve_taxonomy.py`. Defaults to a bundled snapshot if unset. |

CLI flags always take precedence: if both `--hhsuite-pfam-db` and
`SSIGN_HHSUITE_PFAM` are set, the CLI flag wins.

## Convenience aliases (set by the database fetcher)

`scripts/fetch_databases.sh` exports these after a successful download so you
can copy a one-liner into your shell rc file. They are not read by ssign at
run time; the matching CLI flag is the load-bearing handle.

| Variable | CLI equivalent |
|---|---|
| `SSIGN_BAKTA_DB` | `--bakta-db` |
| `SSIGN_EGGNOG_DB` | `--eggnog-db` |
| `SSIGN_INTERPROSCAN_PATH` | `--interproscan-db` |
| `SSIGN_ECOD70_DB` | `--plmblast-db` |
| `SSIGN_PLM_EFFECTOR_WEIGHTS` | (consumed by the PLM-Effector wrapper) |
| `SSIGN_DEEPSECE_CHECKPOINT` | (alternative to `SSIGN_DEEPSECE_CHECKPOINT_URL` for already-downloaded files) |
| `SSIGN_DEEPLOCPRO_PATH` | `--deeplocpro-path` (also used by the DLP integration test) |
| `SSIGN_SIGNALP_PATH` | `--signalp-path` (also used by the SignalP integration test) |

## Python dependency pins (extended tier)

These are not env vars but version constraints. Listed here so a maintainer
debugging a fresh install knows where the upper bounds come from. They are
captured in `[project.optional-dependencies].extended` in `pyproject.toml`,
so `pip install ssign[extended]` resolves them automatically.

| Package | Pin | Reason |
|---|---|---|
| `transformers` | `>=4.38,<5.0` | 5.0 removed `batch_encode_plus`, used by pLM-BLAST and PLM-Effector tokenizers. |
| `numpy` | `>=1.26,<2.0` | 2.0 removed `np.issubsctype`, used by pLM-BLAST's alignment code. |
| `protobuf` | any | Required by ProtT5's SentencePiece tokenizer at load time. |
| `mkl`, `mkl-service` | any | pLM-BLAST's `plmblast.py` imports them directly. |

The `transformers` and `numpy` upper bounds are revisited in the v1.x roadmap
once upstream pLM-BLAST and PLM-Effector publish 5.0/2.0-compatible code.

## Test and developer-only

| Variable | Purpose |
|---|---|
| `SSIGN_TEST_FIXTURE_FULL` | Set to `1` to use the full Xanthobacter contig fixture (~1 Mb) instead of the minimal T5aSS fixture (~213 kb). Slower; closer to a real run. |
| `SSIGN_TEST_OUTDIR` | Output directory for the multi-genome integration test. |
| `SSIGN_RUN_FULL_PIPELINE` | Set to `1` to opt into the long-running fixture pipeline test. Skipped by default. |
| `SSIGN_GOLDEN_REGEN_DIR` | Where regenerated golden-output files are written when refreshing `tests/fixtures/golden/`. See `tests/fixtures/golden/REGENERATE.md`. |
