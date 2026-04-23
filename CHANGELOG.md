# Changelog

All notable changes to **ssign** are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Work toward v1.0.0 publication release. See the roadmap in project memory
(`project_longevity_commitment.md` + plan file) for the full phased schedule.

### Planned

- **Offline-first conversion** — remove all external API dependencies
  (BioLib for DeepLocPro/SignalP, NCBI remote BLAST, EBI InterProScan remote,
  HHpred MPI remote). Required for post-departure longevity.
- **New tools** — Bakta + EggNOG (whole-genome annotation), PLM-Effector
  (first-class prediction tool), pLM-BLAST / ECOD70 (substrate annotation).
- **Dropped tools** — Foldseek, ESM3.
- **Nextflow "power mode" deprecated** — code in `bin/`, `modules/local/`,
  `workflows/`, and `main.nf` to be removed. The Python runner covers all
  power-mode use cases (including HPC batches via Singularity + CLI). This
  removes the double-maintenance burden of mirroring every script between
  `bin/` and `src/ssign_app/scripts/`. Users who relied on Nextflow batching
  will be guided toward CLI + Singularity in the v1.0.0 migration notes.
- **Three install tiers** — `base` (~17 GB), `extended` (~130 GB), `full`
  (~630 GB) — with tier-aware database fetcher pulling from pinned Zenodo
  DOIs. The `extended` tier is the sensible default for lab researchers who
  don't want to download 390 GB of BLAST NR.
- **Re-annotate all inputs with Bakta by default** — even GenBank input gets
  re-run through Bakta (and EggNOG at the extended/full tier), since we can't
  verify what pipeline produced the incoming GenBank annotations. Provides
  uniform, reproducible annotation for consensus voting. Users with curated
  GenBanks can opt out via `--use-input-annotations`.
- **DeepSecE checkpoint mirrored to Zenodo** (URL finalised at v1.0.0 release,
  placeholder in the script until then). SJTU origin server kept as a
  secondary fallback. Env var `SSIGN_DEEPSECE_CHECKPOINT_URL` lets
  institutional mirrors override at runtime. GUI download instructions
  updated with the Zenodo mirror first.
- **Pipeline-order change** — move `enrichment_testing` before
  `filter_by_stats_and_dlp`; stats filter default ON for ≥10 genomes,
  OFF with warning otherwise.
- **Cross-validate rule change** — DLP / DSE / PLM-Effector treated as equal
  secretion predictors (any one flagging = candidate). SignalP becomes
  evidence-only, not a trigger. New `n_prediction_tools_agreeing` column.
- **Docker bundle image** — SHA-pinned, self-contained, all model weights
  baked in, reproducible for 5+ years.
- **Zenodo deposits** — separate DOIs for source, model weights, database
  bundle.
- **FAIR-compliant repository layout** — per iAMB template with extensions
  for software distribution (Docker, PyPI, CLI).
- **Diataxis documentation** structure — tutorials / how-to / reference /
  explanation.
- **bio.tools registration**.

## [0.9.0-prerefactor] — 2026-04-22

Pre-publication baseline snapshot captured before the v1.0.0 refactor. Tagged
as `v0.9.0-prerefactor` on GitHub for regression testing throughout the
publication roadmap.

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
- `CLAUDE.md` developer-instructions file previously tracked in git —
  untracked in this release.

---

Earlier development history is preserved in git; see `git log v0.9.0-prerefactor`
for the full pre-baseline commit record.
