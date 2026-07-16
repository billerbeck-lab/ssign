# Changelog

All notable changes to **ssign** are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Roadmap toward v1.0.0 lives in the [README](README.md#roadmap-to-v100).

## [Unreleased]

## [1.0.0] - 2026-07-16

### Added

- Run-start ASCII banner with version line, shown by `ssign run` and on GUI
  launch.
- Pre-run rough runtime estimate printed before the first step, sized from
  machine resources and genome count (it refines once the first step
  completes). Suppressed on GPU-less hosts, where the GPU-based reference would
  mislead more than help.
- Per-step remaining-time ETA is now the last line printed for each step, so
  it is easy to find.
- Aggregated multi-genome `combined_summary.txt`: pooled secreted-protein
  counts, per-type totals, per-tool annotation coverage, and pooled enrichment.
- `results.csv` opens with a short `#`-prefixed overview block naming its three
  sections (readers and the GUI parser skip it).
- **SignalP as a third enrichment predictor track** (Sec signal peptide),
  alongside DeepLocPro and DeepSecE.
- `--scratch-dir` flag and automatic scratch resolution: keep `$TMPDIR` when it
  has room, else fall back to a directory under `--outdir`. Stops tools (Bakta
  especially) from dying with "No space left on device" when a container's
  `/tmp` is a small tmpfs.
- `--skip-annotation` master switch to turn off all six annotation tools
  (BLASTp, HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam) at once.
- `taxopy>=0.12` for local NCBI taxdump lookup in `resolve_taxonomy.py`,
  replacing remote E-utilities. Dump defaults to `~/.ssign/taxdump/`;
  override with `SSIGN_TAXDUMP_DIR`. Pipeline degrades gracefully if the
  dump is missing.
- `extended` install extras tier (pip extras for the ~130 GB workflow).
- HH-suite per-protein parallelism via `ThreadPoolExecutor` (~4× speedup).
- `SSIGN_DEEPSECE_CHECKPOINT_URL` env override for institutional mirrors.

### Changed

- **Enrichment test rewritten as a per-SS-type circular-shift permutation.**
  Emits fold enrichment, permutation _p_, and BH _q_ per system type for each
  predictor track (DeepLocPro, DeepSecE, SignalP) plus a `COMBINED` union track.
  Replaces the old binomial / sampled-background test, which understated
  significance by ignoring the way secreted genes cluster along the genome. The
  exact rotation null needs whole-genome predictions, so `--enrichment-stats`
  now forces whole-genome DeepLocPro / DeepSecE / SignalP (all run locally). The
  old sampled-background machinery (`n_null_proteins`, `null_seed`,
  `sample_null_proteins.py`) is retired.
- **T5aSS/T5cSS enrichment emits two results**: a `self` result (autotransporter
  self-detection, DLP-or-SignalP) and a `window` "hitchhiker" result (secreted
  neighbours in the ±3 window that may piggyback the pore, DLP-or-DSE). The old
  standalone autotransporter self-detection figure is gone; that signal is now
  the `self` enrichment result.
- **Run figures revamped** (`figures-v2`): curated per-genome set `01`–`06` on a
  shared house style, adding size / physicochemical panels and functional-category
  breakdowns (COG, KEGG, EggNOG, curated consensus). Multi-genome runs emit the
  same set pooled over all genomes.
- **Summary report reworked**: run-metadata header (version | tier | date),
  secreted-protein count with per-type breakdown, and annotation coverage
  grouped by tool.
- **Tool subprocess timeouts are size-aware**: each modelled tool gets
  `max(4 h floor, 2 × predicted runtime)` instead of a flat 4 h cap, so pooled
  whole-genome steps no longer die mid-run.
- DeepSecE inference auto-sizes its batch to available VRAM and falls back to
  batch size 1 on CUDA out-of-memory.
- **Extended/full container is now self-contained.** The image bakes the heavy
  free toolchain (Bakta with its BLAST 2.17 / HMMER / tRNAscan-SE / diamond
  stack, EggNOG-mapper, and HH-suite) into pinned micromamba envs, plus the
  DeepLocPro, DeepSecE, and pLM-BLAST (ProtT5) model weights and their ESM
  backbones. Extended/full runs now need only the licensed DTU tools and the
  reference databases mounted from the host, instead of borrowing tools from
  fragile host conda envs.
- DeepSecE checkpoint is baked into the container image (container runs never
  download it); plain-pip installs fetch it from the authors' SJTU origin with
  a retry + size check and the `SSIGN_DEEPSECE_CHECKPOINT_URL` override. ssign
  does not re-host the checkpoint.
- Repository moved to the `billerbeck-lab` GitHub organisation. Old
  `reidmat/ssign` URLs continue to redirect.
- Bakta minimum bumped from `>=1.5` (2022) to `>=1.10` (2024), required
  to read the Bakta DB v6 schema. v1.5 cannot parse modern Bakta DBs.

### Fixed

- **Package `runtime/coefficients.json`** so it ships in the wheel. It was
  omitted from `[tool.setuptools.package-data]`, so `pip install`ed and
  containerised copies had no effort-model coefficients. That silently (the
  effort model degrades rather than crashes) disabled the live runtime ETA and
  reverted size-aware tool timeouts to the flat 4 h floor; the latter can kill
  a legitimately-long pooled/whole-genome step. Source checkouts were
  unaffected. Added a guard test that every non-`.py` data file under
  `src/ssign_app` is covered by a package-data glob.

### Notes

- **EggNOG database v7.0** (Hernández-Plaza et al., NAR 2025, 54:D402) is
  the new state-of-the-art, but ssign v1.0.0 ships against EggNOG v6.0.
  Reason: eggnog-mapper 2.1.13 does not yet read v7; see
  [eggnog-mapper#588](https://github.com/eggnogdb/eggnog-mapper/issues/588).
  When upstream adds v7 support we'll bump to track it.
- **DeepSecE error message** in `run_deepsece.py` previously pointed at
  `github.com/SijinHuang/DeepSecE` (404, dead fork). Fixed to point at the
  real upstream `github.com/zhangyumeng1sjtu/DeepSecE`.
- **The distributed container image is non-commercial** (research use only),
  because it bakes DeepLocPro (CC BY-NC-SA 4.0). The pip install is unaffected.

### Removed

- **PLM-Effector predictor, entirely.** Trialled as a fourth secreted-protein
  predictor, then removed: it over-called ~25% of a proteome (the model is
  validated only on balanced sets), which swamped the enrichment test.
  `run_plm_effector.py`, the vendored `plm_effector/`, the
  `--skip-plm-effector` / `--no-skip-plm-effector` / `--plme-batch-size` flags,
  and its PLM-Effector-only weights (ProtBert + the effector checkpoint) are
  gone. Rationale retained in `docs/explanation/design_decisions.md`.
- `scripts/fetch_weights.sh`: model weights now download automatically on
  first run (or come pre-baked in the container image), so the separate
  fetch step was dead scaffolding (it existed mainly to pre-stage PLM-Effector
  weights).
- `eggnog-mapper` from the `[extended]` / `[full]` pip extras: its
  `biopython==1.76` pin clashes with `bakta>=1.78`. Install separately
  via `conda install -c bioconda eggnog-mapper`; full recipe in
  `docs/how-to/install.md` § EggNOG-mapper. EggNOG is off by default
  (`--skip-eggnog`), so most users are unaffected.
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
  with `ssign run X --outdir Y`. See `docs/explanation/design_decisions.md` § 6.3.

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
