# Contributing to ssign

How to report issues, propose changes, and submit pull requests. ssign is
developed by the Billerbeck Lab at Imperial College London.

---

## Ways to contribute

- **Bug reports or unexpected behaviour:** GitHub issue.
- **Feature requests:** GitHub issue tagged `enhancement`, before writing
  code, so we can discuss scope.
- **Code or documentation fixes:** fork, branch, PR.
- **New secretion-related tool or annotation source:** open an issue with a
  brief justification and a license check for the new dependency.
- **Documentation improvements:** the lowest-friction way to contribute.

---

## Before you start

1. Check open issues and PRs first; someone may already be working on the
   same thing.
2. Discuss non-trivial changes first via an issue.
3. Large structural changes should fit the publication roadmap. Flag in
   the issue if they don't.

---

## Reporting bugs

Open an issue with:

- **ssign version** (`ssign --version` or commit SHA).
- **Platform:** OS, Python version, CUDA version if relevant.
- **Exact command(s) run** and **full error / log output** (code fences).
- **Expected vs actual behaviour**.
- A **minimal reproducer:** smallest input that triggers the issue. A small
  GenBank or FASTA attached is ideal.
- Whether you've tried the latest `main`.

Don't paste entire genome files; attach or link them.

---

## Development setup

```bash
git clone https://github.com/billerbeck-lab/ssign.git
cd ssign
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Run tests before submitting:

```bash
pytest tests/unit/ -v                  # fast unit tests
pytest tests/integration/ -v           # integration tests (some need network)
```

---

## Coding conventions

- **Python ≥ 3.11**. Pin versions in `pyproject.toml` when adding deps.
- **Formatting / linting**: `ruff format` and `ruff check` are authoritative.
- **Type hints** required for new public functions; `mypy --strict` should
  pass on `src/ssign_app/`.
- **Comments** only when the _why_ is non-obvious.
- **Critical-path code** (listed in `docs/`) has regression tests; preserve
  the comment block documenting why each fragile section is fragile.

---

## Tests

- New features **must** include unit tests (Arrange-Act-Assert).
- Keep fixtures small and real (prefer a 3-gene test genome over a mock).
  See `tests/fixtures/`.
- Integration tests that hit external services: `@pytest.mark.integration`.
- CI must pass on your PR.

---

## Pull requests

1. Fork, branch from `main`. Naming: `feature/short-description`,
   `fix/issue-NN-short`, `docs/short-description`.
2. One logical change per PR.
3. Clear PR description: what, why, how tested. Link the issue.
4. Add a `CHANGELOG.md` entry under `## [Unreleased]` for user-facing
   changes.
5. Rebase rather than merge `main` back in.
6. Signed commits preferred (GPG or SSH).

Merges are squash-merge into `main`. Expect review within ~1 week.

---

## License

ssign is distributed under **GPL-3.0-or-later** (`LICENSE`). By submitting a
PR, you agree your contribution is licensed under the same terms.

Substantial contributions (new tool integrations, new analysis modules) are
added to `CITATION.cff`. Drive-by fixes appear in the GitHub contributor
graph.

---

## Conduct

Be respectful, assume good faith. Harassment and personal attacks are not
welcome.

---

## Contact

| Need                                                                       | Channel                                                                    | Who                                             |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ----------------------------------------------- |
| Bug reports, feature requests, install troubles, code questions            | [GitHub Issues](https://github.com/billerbeck-lab/ssign/issues)            | Active maintainer                               |
| Pull request review                                                        | GitHub PR                                                                  | Active maintainer                               |
| Scientific collaboration, data-sharing, partnerships, authorship inquiries | Email: [`s.billerbeck@imperial.ac.uk`](mailto:s.billerbeck@imperial.ac.uk) | Dr. Sonja Billerbeck (PI, corresponding author) |

**Active maintainer (as of v0.9.0):** M. Teo Reid
([`@reidmat`](https://github.com/reidmat), `t.reid25@imperial.ac.uk`). A lab
successor will be named in `SYSADMIN.md` before September 2026; until then,
GitHub Issues are the fastest path to a technical answer.

**Do not email the PI with bug reports or install trouble.** Those belong
on GitHub Issues where they're searchable.
