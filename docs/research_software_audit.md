# Research-software audit (v1.0.0 readiness)

Gap analysis of ssign against published standards for research-software
publication and FAIR-ness. Snapshot: `2026-06-03` at commit `9d7f860` on `main`.

This document does not propose code changes. It is a checklist with verdicts
and pointers; pick which gaps are worth closing before v1.0.0 release.

## Frameworks audited against

1. **JOSS submission checklist** (Journal of Open Source Software)
   — <https://joss.readthedocs.io/en/latest/review_criteria.html>
2. **FAIR Principles for Research Software (FAIR4RS)** — Chue Hong et al.
   2022, _Scientific Data_, <https://doi.org/10.1038/s41597-022-01710-x>
3. **OpenSSF Best Practices Badge, passing tier** — <https://www.bestpractices.dev>
4. **Software Sustainability Institute (SSI) practical recommendations** —
   informal cross-check on long-term maintainability.

Legend: ✅ in place — ⚠ partial — ❌ missing — — n/a

---

## 1. JOSS submission readiness

JOSS reviews a code repository and a short companion paper together. The
items below are the repository half; the paper half (`paper.md` +
`paper.bib`) is a separate v1.0.0 deliverable.

### Software license, scope, significance

| Item | Verdict | Pointer |
| --- | --- | --- |
| OSI-approved license file present (not just README reference) | ✅ | `LICENSE` (GPL-3.0-or-later, 35 KB full text) |
| Scope is a substantive scholarly contribution, not a single-use script | ✅ | Integrates 11+ external tools in a documented pipeline; ~280 commits, 54 unit-test files, 1145 collected tests |
| Authorship reasonable and complete | ✅ | `CITATION.cff` lists 4 authors with ORCIDs + affiliations |

### Development history and open-source practice

| Item | Verdict | Pointer |
| --- | --- | --- |
| Public repository under version control | ✅ | `github.com/billerbeck-lab/ssign` |
| ≥ 6 months of public commit history | ❌ | First commit `2026-03-30`, latest `2026-06-03` (~2 months public). JOSS reviewers flag commits concentrated near submission. Mitigation: extensive private history pre-dates the public repo; document this in the paper |
| Documented releases / version tags | ⚠ | Only `v0.9.0-prerefactor`; no v1.0.0 yet. CHANGELOG follows Keep a Changelog format |
| Public issues and pull requests | ⚠ | Repo is open; need to confirm there are actually filed issues/PRs to demonstrate it gets used |
| Multiple developers contributed | ❌ | Git history shows a single committer (`reidmat`). CITATION.cff lists 4 authors — the gap is between commit attribution and paper authorship. JOSS flags this for "Collaborative Effort"; the paper's authorship narrative needs to explain (e.g. PI advisor + collaborators) |
| Code review via pull requests | ❌ | No PR template; default branch direct-push pattern visible in commit log. Adopt PR-based review before submission |
| External user engagement | ❓ | Unknown until v1.0.0 is publicly announced. Track adoption signals (issues, citations, forks) before submitting |

### Documentation

| Item | Verdict | Pointer |
| --- | --- | --- |
| README present with installation, quickstart, example | ✅ | `README.md` (310 lines) |
| Diátaxis-style documentation layout | ⚠ | `docs/{tutorials,how-to,reference,explanation}/` all present, but `tutorials/` has only one entry (`first_run.md`); expand before submission (task #57 already tracks this) |
| Statement of need clearly labelled | ⚠ | README opens with what ssign does, but no headed "Statement of need" or "Why ssign?" section. JOSS treats this as a required item with an explicit heading |
| State-of-the-field comparison | ❌ | README does not compare ssign to closest alternatives (e.g. T3SS-finder, BastionX, EffectiveDB, SecretEPDB). Required by JOSS for paper |
| Installation instructions: dependencies + automated procedure | ✅ | `docs/how-to/install.md`, `pip install ssign[...]`, `scripts/fetch_databases.sh` |
| Example usage with real-world demo | ✅ | `docs/tutorials/first_run.md` |
| API / CLI reference | ✅ | `docs/reference/cli.md`, `env_vars.md`, `output_files.md` |
| Community guidelines: contribute, report, support | ✅ | `CONTRIBUTING.md` (133 lines) |
| AI usage disclosure (now required by JOSS) | ❌ | No statement of which generative-AI tools were used during development or how output was verified |

### Functionality and tests

| Item | Verdict | Pointer |
| --- | --- | --- |
| Software installs and core functionality works | ✅ | Tested manually (CX3 K-12 validation runs, tutorial run) |
| Automated test suite present and documented | ✅ | 1145 collected unit tests across 54 files in `tests/unit/`; integration suite in `tests/integration/` with `@pytest.mark.integration` opt-in |
| Continuous integration | ✅ | `.github/workflows/{test,lint}.yml` — pytest across Python 3.10-3.13, ruff, mypy, dep-check |
| Coverage measured | ⚠ | `pytest-cov` is in `[dev]`; `.coverage` artefact exists locally but CI does not generate, publish, or enforce a threshold. No coverage badge in README |

### Software paper deliverables (separate from this audit)

| Item | Verdict | Pointer |
| --- | --- | --- |
| `paper.md` + `paper.bib` written | ❌ | Not present in repo |
| Statement of need section in paper | ❌ | Pending paper draft |
| State of the field comparison | ❌ | Pending paper draft |
| Software design rationale | ✅ (source material exists) | `docs/explanation/design_decisions.md` already explains many design choices; lift into paper |
| Research impact statement | ❌ | Pending |
| AI usage disclosure in paper | ❌ | Pending |

---

## 2. FAIR4RS principles

Findable, Accessible, Interoperable, Reusable — applied to software.

### Findable

| Principle | Verdict | Pointer |
| --- | --- | --- |
| F1: globally unique persistent identifier (DOI) | ❌ | No Zenodo DOI yet. README roadmap commits to one at v1.0.0. Without a DOI the software is not citable in the FAIR sense |
| F1.1: components also have identifiers (e.g. release tags) | ⚠ | One git tag only; semver releases planned but not yet exercised |
| F1.2: identifier resolves to a landing page | ❌ | Pending Zenodo deposit |
| F2: rich metadata for the software | ⚠ | `CITATION.cff` and `pyproject.toml` cover human-readable metadata. Missing machine-readable `codemeta.json` (FAIR4RS-preferred format) |
| F3: metadata explicitly includes the identifier | ❌ | Depends on F1 |
| F4: metadata indexed in a searchable resource | ❌ | `bio.tools` registration is on the roadmap (task #56). Also: PyPI listing pending; package not yet published |

### Accessible

| Principle | Verdict | Pointer |
| --- | --- | --- |
| A1: retrievable via standardised, open protocol | ✅ | `git clone` from GitHub over HTTPS; `pip install` once published |
| A1.1: open, free, universal | ✅ | GitHub + (planned) PyPI + (planned) Zenodo |
| A1.2: auth where necessary | — | No auth needed for the open-source code path |
| A2: metadata accessible even when software is no longer available | ❌ | Without a Zenodo deposit, repo deletion or org migration would orphan everything. Zenodo deposit at v1.0.0 closes this |

### Interoperable

| Principle | Verdict | Pointer |
| --- | --- | --- |
| I1: reads/writes domain-standard formats | ✅ | Inputs: GenBank, FASTA, GFF3. Outputs: CSV, TSV, FASTA (per `docs/reference/output_files.md`) |
| I2: qualified references to other research objects | ✅ | README and CITATION.cff cite upstream tools with DOIs |

### Reusable

| Principle | Verdict | Pointer |
| --- | --- | --- |
| R1: described with attributes (purpose, version, license, authors) | ✅ | CITATION.cff + pyproject.toml + README |
| R1.1: clear, accessible licence | ✅ | `LICENSE` + classifier + `docs/explanation/licensing.md` |
| R1.2: detailed provenance | ⚠ | CHANGELOG documents per-release changes; specific commits exist. Provenance for the trained model weights and DB bundles is not yet on a citable archive (Zenodo deposit pending) |
| R2: qualified references to other software | ✅ | pyproject.toml pins minimum versions with rationale comments; README cites all integrated tools |
| R3: meets domain-relevant community standards | ✅ | Pipeline outputs follow bioinformatics convention (TSV, FASTA, CSV); CLI follows POSIX conventions |

---

## 3. OpenSSF Best Practices, passing tier

Most of the criteria for the passing badge are already met. Listing the
gaps; everything not listed is in place.

| Item | Verdict | Pointer |
| --- | --- | --- |
| Basics: project actively maintained | ✅ | Daily commits, named maintainer in CONTRIBUTING.md |
| Basics: discussion mechanism | ⚠ | GitHub Issues only; no Discussions tab. Issues are sufficient for the badge |
| Change Control: unique version per release | ✅ | SemVer in pyproject + tags |
| Change Control: release notes that identify vulnerabilities (with CVEs if any) | ⚠ | CHANGELOG present and well-structured; no template for vulnerability section in releases — add when needed |
| Reporting: vulnerability reporting process documented | ❌ | No `SECURITY.md` describing where to send security reports. GitHub default ("Report a vulnerability") is enabled by default but not advertised |
| Reporting: private reporting channel | ❌ | Same fix: a `SECURITY.md` |
| Quality: build system reproducible | ✅ | `pyproject.toml` declares build-backend; `uv.lock` pins resolved dependencies |
| Quality: automated test suite documented | ✅ | `CLAUDE.md` and `CONTRIBUTING.md` document how to run tests |
| Quality: test policy for new features | ✅ | CONTRIBUTING.md § Tests |
| Quality: warning flags / linter enabled | ✅ | Ruff `E F I W`, mypy with bug-finder config |
| Quality: warnings addressed | ✅ | CI gates merges on ruff + mypy |
| Security: secrets not in repo | ✅ | `.gitignore` covers usual suspects; no scanned violations |
| Analysis: static analysis on every commit | ⚠ | Ruff + mypy in CI cover style/types. No security-oriented static analyser (e.g. `bandit`, `pip-audit`, `safety`, GitHub Dependabot alerts) — recommended for the badge but not required at passing level |

---

## 4. Software Sustainability Institute cross-check

Additional best-practice items that are not in the three frameworks above
but matter for a publication-grade research-software project.

| Item | Verdict | Pointer |
| --- | --- | --- |
| `CODE_OF_CONDUCT.md` | ❌ | CONTRIBUTING.md has a one-paragraph conduct statement; a standalone CoC (Contributor Covenant 2.1 is the usual choice) is the convention |
| `SECURITY.md` | ❌ | See OpenSSF gap above |
| Pull-request template | ❌ | `.github/PULL_REQUEST_TEMPLATE.md` would standardise the review checklist already informally documented in CONTRIBUTING.md |
| CI status badge in README | ❌ | README has License / Python / Status badges; no test or lint status badge |
| Coverage badge in README | ❌ | Same root cause as the coverage gap above |
| `codemeta.json` (FAIR-preferred machine-readable metadata) | ❌ | Generates automatically from CITATION.cff via `cffconvert` — trivial to add |
| `.zenodo.json` for automated Zenodo-on-release | ❌ | Manual deposit also fine; this file just removes a step |
| Reproducible container image with real SHA pin | ⚠ | `containers/Dockerfile` has the scaffolding but the base-image SHA at line 35 is a placeholder (`@sha256:000…000`). Image is not actually reproducible until pinned |
| Dead-code cleanup (post Nextflow removal) | ⚠ | `conf/` (Nextflow configs) and `bin/ssign_lib/` (Nextflow-bin pycache) are no longer referenced from any source file (grep finds no references in src, scripts, docs). Delete before publication |
| Code of practice for AI-generated content | ❌ | Not declared anywhere (also required by JOSS, see § 1) |
| ORCID for every author in pyproject.toml | ⚠ | CITATION.cff has ORCIDs but pyproject.toml authors-table does not (PEP 621 doesn't have an `orcid` field; acceptable to leave CITATION.cff as the canonical source) |
| CITATION.cff version matches pyproject.toml version | ⚠ | `pyproject.toml` says `0.9.0`; `CITATION.cff` says `0.9.0-prerefactor`. Re-sync at v1.0.0 release; tooling like `cffconvert --validate` can catch drift |

---

## 5. Cross-cutting strengths

Items already in good shape — keep these intact through the v1.0.0 sprint.

- LICENSE file with full text, OSI-approved, declared in pyproject + classifier + CITATION.cff.
- CITATION.cff with ORCIDs and explicit upstream-tool attribution.
- CHANGELOG.md in Keep-a-Changelog format with `[Unreleased]` working section.
- Diátaxis-shaped documentation, populated.
- CI matrix across Python 3.10-3.13; ruff + mypy + import smoke-test (`ssign doctor --imports-only`).
- Three install tiers, with rationale documented in pyproject comments and `docs/how-to/install.md`.
- Per-dep pinning with comments explaining why bounds exist (e.g. `transformers<5.0`, `numpy<2.0`).
- `ssign doctor` self-diagnostic surfaces missing externals with fix commands.
- README cites every integrated upstream tool with DOI.
- `docs/explanation/design_decisions.md` already captures most rationale a JOSS paper needs.

## 6. Cross-cutting gaps (consolidated)

The most consequential gaps, in approximate priority order for a publication
release. Detail in the framework tables above.

1. **No persistent identifier.** No Zenodo DOI, no PyPI release, no machine-
   indexed registry entry. Blocks FAIR F1-F4 and the JOSS citation flow.
2. **No software paper.** `paper.md` + `paper.bib` need writing; AI-disclosure
   section is mandatory.
3. **State-of-the-field comparison missing.** Both README and paper need
   it.
4. **Statement-of-need not labelled.** Reword opening of README so a JOSS
   reviewer can find it.
5. **No `SECURITY.md`, no `CODE_OF_CONDUCT.md`, no PR template.** Cheap to
   add.
6. **Single-author git history vs four-author citation.** Either bring
   collaborators into the public history (PRs, co-authored-by trailers) or
   address the discrepancy in the paper's author-contributions section.
7. **Dockerfile SHA pin is a placeholder.** "SHA-pinned" claim in README
   and Dockerfile preamble is currently aspirational.
8. **Dead Nextflow remnants** (`conf/`, `bin/ssign_lib/__pycache__`)
   should be deleted now that the Nextflow path is gone.
9. **Coverage not gated or reported in CI.** `pytest-cov` is installed but
   unused on the CI side.
10. **`codemeta.json` and `.zenodo.json` not present.** Cheap to generate
    from CITATION.cff once the Zenodo deposit exists.

---

## 7. Suggested follow-ups, grouped by effort

These are candidate tasks, not a commitment. Pick which to take on as
v1.0.0 work.

### Trivial (≤ 1 hour each, low-risk file additions)

- Add `SECURITY.md` (one-page: where to report, expected response time).
- Add `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1 boilerplate).
- Add `.github/PULL_REQUEST_TEMPLATE.md` mirroring the checklist already in
  CONTRIBUTING.md § Pull requests.
- Add CI status badges to README (test, lint, coverage if added).
- Reword README opening to lead with a headed "Statement of need" section.
- Add an AI-usage disclosure section to README (and later to `paper.md`).
- Generate `codemeta.json` via `cffconvert --infile CITATION.cff --format codemeta`.
- Resync `CITATION.cff` version with `pyproject.toml` at the v1.0.0 cut.
- Delete `conf/` and `bin/ssign_lib/__pycache__/` (verify no references first;
  `grep -rln "conf/" src/ scripts/ docs/` came up empty).

### Moderate (a few hours each)

- Mint a Zenodo deposit for v1.0.0 (DOI for the source code archive); add
  `.zenodo.json` so future tags auto-deposit; cite the DOI in README and
  CITATION.cff.
- Pin the Dockerfile base-image SHA to a real value (run `docker pull
  nvidia/cuda:12.4.1-runtime-ubuntu22.04` then `docker inspect`).
- Add a "State of the field" subsection in README and paper comparing
  ssign to closest alternatives (T3SS-finder, BastionX, EffectiveDB, etc.).
- Wire `pytest-cov` into the `test.yml` workflow; emit a coverage report;
  optionally fail under a threshold; add badge.
- Add `pip-audit` (or `safety` or Dependabot) to CI for dependency-CVE
  monitoring.
- Register on `bio.tools` (task #56 already tracked).
- Write a brief `AUTHORS.md` or expand CITATION.cff with author contribution
  roles (CRediT taxonomy is the convention).

### Large (multi-day, paper-track)

- Draft `paper.md` + `paper.bib` to JOSS structure: summary, statement of
  need, state of the field, software design, research-impact statement, AI
  disclosure, acknowledgements, references.
- Publish to PyPI: pre-flight `pip install ssign` (currently 404 per the
  CLAUDE memory entry); resolve the build-backend + extras matrix; tag and
  upload.
- Run a PyPI smoke install across all three tiers on a clean machine to
  catch packaging drift.
- Bring collaborators into public commit history (PRs with co-authored-by,
  or explicit author-contributions paragraph in the paper).
- Reach 6+ months of visible commit activity on the public repo before
  JOSS submission (this is calendar time, not effort — start the JOSS
  clock now if v1.0.0 ships in summer).

---

_End of audit._
