# Contributing to ssign

How to report issues, propose changes, and submit pull requests. ssign is
developed by the Billerbeck Lab at Imperial College London.

---

## Ways to contribute

- **Bug reports or unexpected behaviour:** GitHub issue.
- **Feature requests:** GitHub issue tagged `enhancement` before writing
  code, so we can discuss scope. Check open issues first, and flag anything
  that would not fit the publication roadmap.
- **Code or documentation fixes:** fork, branch, PR.
- **New secretion-related tool or annotation source:** open an issue with a
  brief justification and a license check for the new dependency.
- **Documentation improvements:** the lowest-friction way to contribute.

---

## Reporting bugs

Open an issue with:

- **ssign version** (`ssign --version` or commit SHA), and whether the
  issue reproduces on latest `main`.
- **Platform:** OS, Python version, CUDA version if relevant.
- **Exact command(s) run** and **full error / log output** (code fences).
- **Expected vs actual behaviour.**
- A **minimal reproducer:** smallest input that triggers the issue. A
  small GenBank or FASTA attached is ideal; don't paste entire genome
  files.

---

Run tests before submitting:

```bash
pytest tests/unit/ -v                       # fast unit tests
pytest -m integration tests/integration/    # integration tests (need external tools installed)
```

---

## Coding conventions

- **Python ≥ 3.10** (matches `pyproject.toml`'s `requires-python`). Pin
  versions in `pyproject.toml` when adding deps.
- **Formatting / linting**: `ruff format` and `ruff check`.
- **Type hints** welcome on new public functions, but not gate-enforced
  (the mypy config leaves `disallow_untyped_defs` off; it runs as a
  bug-finder, checking bodies of untyped functions).
- **Type check**: `mypy` (config in `pyproject.toml`) must pass. CI runs it
  bare, with no arguments.
- **Comments** only when the _why_ is non-obvious.
- **Critical-path code** (listed in `docs/`) has regression tests; preserve
  the comment block documenting why each fragile section is fragile.

---

## Tests

- New features **must** include unit tests.
- Keep fixtures small and real. See `tests/fixtures/`.
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

---

## License

ssign is distributed under **GPL-3.0-or-later** (`LICENSE`). By submitting a
PR, you agree your contribution is licensed under the same terms.

Substantial contributions (new tool integrations, new analysis modules) are
added to `CITATION.cff`.

---

## Conduct

Be respectful, assume good faith. Harassment and personal attacks are not
welcome.

---

## Contact

- **Bugs, feature requests, install trouble, code questions, PR review:**
  [GitHub Issues](https://github.com/billerbeck-lab/ssign/issues).
- **Scientific collaboration, data sharing, authorship:** Dr. Sonja Billerbeck
  (PI), [`s.billerbeck@imperial.ac.uk`](mailto:s.billerbeck@imperial.ac.uk).
- **The active maintainer** as of v1.0.0 is M. Teo Reid ([`@reidmat`](https://github.com/reidmat)).
