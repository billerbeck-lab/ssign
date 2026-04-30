# Changelog

All notable changes to **ssign** are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Roadmap toward v1.0.0 lives in the [README](README.md#roadmap-to-v100).

## [Unreleased]

### Added

- `taxopy>=0.12` for local NCBI taxdump lookup in `resolve_taxonomy.py`,
  replacing remote E-utilities. Dump defaults to `~/.ssign/taxdump/`;
  override with `SSIGN_TAXDUMP_DIR`. Pipeline degrades gracefully if the
  dump is missing.
- `extended` install extras tier (pip extras for the ~130 GB workflow).
- HH-suite per-protein parallelism via `ThreadPoolExecutor` (~4× speedup).
- `SSIGN_DEEPSECE_CHECKPOINT_URL` env override for institutional mirrors.

### Changed

- DeepSecE checkpoint fetched from a Zenodo mirror first (URL placeholder
  until the v1.0.0 deposit), SJTU origin retained as fallback.
- Repository moved to the `billerbeck-lab` GitHub organisation. Old
  `reidmat/ssign` URLs continue to redirect.

### Removed

- Remote modes for BLASTp (NCBI), HH-suite (MPI Toolkit), and InterProScan
  (EBI web service). All three now require a local install.
- Foldseek scaffolding (never reached first-class status; dropped for v1.0.0).
- `pybiolib` dependency (unused in the codebase) and DTU diagnostic scripts.
- GUI mode toggles for tools whose remote path has been removed.
- **Nextflow "Power Mode" pipeline.** `main.nf`, `nextflow.config`,
  `workflows/`, `modules/local/`, `bin/`, and the per-tool `containers/`
  Dockerfiles are gone. The pure-Python `ssign run` CLI plus a Phase 4b
  Singularity image cover the HPC use case without the double-bookkeeping
  overhead that kept every script mirrored between `bin/` and
  `src/ssign_app/scripts/`. Migration: replace any
  `nextflow run main.nf --input X --outdir Y -profile docker` invocation
  with `ssign run X --outdir Y`. See `docs/design_decisions.md` § 6.3.

## [0.9.0-prerefactor] — 2026-04-22

Pre-publication baseline snapshot, tagged as `v0.9.0-prerefactor` on GitHub
for regression testing throughout the publication roadmap.

### Current features

- Secretion-system detection via MacSyFinder v2 + TXSScan.
- Secreted-protein prediction via DeepLocPro (BioLib), DeepSecE (local),
  SignalP 6.0 (BioLib).
- Per-component genomic proximity analysis (+/- 3 genes by default, same
  contig only).
- T5SS barrel-domain handling (PF03797).
- Optional annotation: BLASTp (local/remote), HH-suite (remote via MPI
  Toolkit), InterProScan (local/remote), ProtParam, Foldseek (scaffolded).
- Streamlit GUI with dark mode, per-genome parallelism, resume-from-checkpoint.
- Nextflow DSL2 "power mode" for HPC batch runs.
- Semaphore-based per-API job scheduling (DTU, NCBI, MPI, EBI).
- Multi-genome cross-genome summary.
- Annotation-consensus voting across tools (17 broad functional categories).

### Known limitations at baseline

- Relies on external APIs (BioLib, NCBI, MPI, EBI) for several tools — will
  break if those services change. Addressed by v1.0.0 offline-first work.
- DeepSecE checkpoint hosted on an unreliable SJTU server. Will be mirrored
  to Zenodo for v1.0.0.
- DTU academic licenses (SignalP, DeepLocPro): redistribution in a public
  Docker image requires confirmation from DTU. Pending.

---

Earlier development history is preserved in git; see `git log v0.9.0-prerefactor`
for the full pre-baseline commit record.
