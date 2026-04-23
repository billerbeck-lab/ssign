# Contributing to ssign

Thanks for your interest in contributing to **ssign**. This document covers how
to report issues, propose changes, and submit pull requests. ssign is developed
and maintained by the Billerbeck Lab at Imperial College London.

---

## Ways to contribute

- **Report bugs or unexpected behaviour** — file an issue on GitHub.
- **Suggest features or improvements** — file an issue tagged `enhancement`
  before writing code; we'll discuss scope and design before you invest effort.
- **Submit code or documentation fixes** — fork, branch, PR.
- **Contribute biology / tool-integration ideas** — propose a new secretion-
  related tool or annotation source by opening an issue; include a brief
  justification and a license check for the new dependency.
- **Improve documentation** — docs PRs are always welcome and are the lowest-
  friction way to contribute.

---

## Before you start

1. **Check open issues and PRs** — someone may already be working on the same
   thing.
2. **Discuss non-trivial changes first** — open an issue describing the change
   you want to make. This prevents wasted work on changes the maintainers
   would rather reject or scope differently.
3. **Check the plan** — large structural changes should fit into the
   publication roadmap (see project memory / plan files). If they don't, say
   so in the issue and we'll decide whether to adjust.

---

## Reporting bugs

Open an issue with:

- **ssign version** (`ssign --version` or the git commit SHA).
- **Platform** — OS and version, Python version, CUDA version if relevant.
- **Exact command(s) you ran** and the **full error / log output**. Use code
  fences.
- **Expected vs actual behaviour**.
- A **minimal reproducer** — smallest input that triggers the issue. A small
  GenBank or FASTA file attached to the issue is ideal.
- Whether you've tried the latest `main` branch.

Please don't paste entire genome files — upload them as GitHub attachments or
link them externally.

---

## Development setup

```bash
git clone https://github.com/reidmat/ssign.git
cd ssign
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Run the test suite before submitting any change:

```bash
pytest tests/unit/ -v                  # fast unit tests
pytest tests/integration/ -v           # integration tests (needs network for some)
```

For details on running the full pipeline, see the `docs/` directory and the
README.

---

## Coding conventions

- **Python ≥ 3.11**. Pin versions in `pyproject.toml` when adding new deps.
- **Formatting / linting**: `ruff format` and `ruff check` are authoritative.
  CI will reject PRs with unresolved `ruff check` errors.
- **Type hints** — required for new public functions; encouraged for internal
  ones. `mypy --strict` should pass on `src/ssign_app/`.
- **Comments** — write a comment only when the _why_ is non-obvious. Avoid
  comments that only describe _what_ a well-named function already expresses.
- **Critical-path code** (listed in the project `CLAUDE.md` / `docs/`) has
  regression tests — do not alter without understanding the bug fix that
  motivated it. Preserve the comment block documenting why each fragile
  section is fragile.

---

## Tests

- New features **must** include unit tests. Use the Arrange-Act-Assert pattern
  (set up inputs, run the code, assert the output).
- Keep fixtures small and real — prefer a 3-gene test genome over a mocked
  database response. See `tests/fixtures/`.
- Integration tests that hit external services should be marked with
  `@pytest.mark.integration` and not block the default `pytest` run.
- Run tests locally before pushing; CI must pass on your PR.

---

## Pull requests

1. **Fork** the repo and create a branch from `main`. Branch naming:
   `feature/short-description`, `fix/issue-NN-short`, `docs/short-description`.
2. **One logical change per PR.** If you're doing both a bug fix and a feature,
   submit two PRs.
3. **Write a clear PR description** — what the change does, why, how it was
   tested. Link the issue it addresses.
4. **Add a `CHANGELOG.md` entry** under `## [Unreleased]` describing user-
   facing changes (not internal refactors).
5. **Keep commits clean** — prefer one commit per PR after review feedback is
   addressed. Rebase, don't merge `main` back in.
6. **Sign your commits** (optional but preferred) with GPG or SSH signing.

A maintainer will review within roughly a week. Expect feedback on scope, tests,
and style. Merges are by squash-merge into `main`.

---

## License of contributions

ssign is distributed under the **GNU General Public License v3.0 or later**
(see `LICENSE`). By submitting a pull request, you agree that your contribution
is licensed under the same terms and that you have the right to submit it.

---

## Citations and attribution

If your contribution is substantial (a new tool integration, a new analysis
module, significant pipeline changes), you'll be added to `CITATION.cff` as a
contributor. Drive-by fixes don't need to be in the citation file but will
appear in the GitHub contributor graph.

---

## Conduct

Be respectful, be constructive, assume good faith. Disagreements on technical
direction should stay technical. Bullying, harassment, or personal attacks are
not welcome here.

---

## Contact

Different channels for different needs — please use the right one so your
question reaches someone who can actually help.

| Need                                                                       | Channel                                                                    | Who                                             |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ----------------------------------------------- |
| Bug reports, feature requests, install troubles, code questions            | [GitHub Issues](https://github.com/reidmat/ssign/issues)                   | Active maintainer                               |
| Pull request review                                                        | GitHub PR                                                                  | Active maintainer                               |
| Scientific collaboration, data-sharing, partnerships, authorship inquiries | Email: [`s.billerbeck@imperial.ac.uk`](mailto:s.billerbeck@imperial.ac.uk) | Dr. Sonja Billerbeck (PI, corresponding author) |

**Active maintainer (as of v0.9.0):** M. Teo Reid
([`@reidmat`](https://github.com/reidmat), `t.reid25@imperial.ac.uk`). A lab
successor will be named in `SYSADMIN.md` before
September 2026; until then, GitHub Issues are the fastest path to a technical
answer.

**Please do not email the PI with bug reports or installation trouble** —
those belong on GitHub Issues where they're searchable for future users.
The PI is not set up to debug software.
