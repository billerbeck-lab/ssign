"""Guard: every non-Python data file under ``src/ssign_app`` must be declared in
``[tool.setuptools.package-data]`` so it ships in the wheel.

Regression test for a real bug: ``runtime/coefficients.json`` was omitted from
the package-data globs, so ``pip install ssign`` (and the container) shipped
without it. ``load_coefficients`` then raised ``FileNotFoundError``, which the
runtime effort model swallowed — silently disabling the live ETA *and* reverting
size-aware tool timeouts to the flat 4 h floor. The failure was invisible from a
source checkout (the file is on disk there); only an installed copy hit it. This
test catches that class of omission without building a wheel.
"""

from __future__ import annotations

import glob
import subprocess
import sys
from pathlib import Path

import pytest

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - 3.10 fallback
    tomllib = pytest.importorskip("tomli")

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PKG_DIR = _REPO_ROOT / "src" / "ssign_app"


def _declared_globs() -> list[str]:
    data = tomllib.loads((_REPO_ROOT / "pyproject.toml").read_text())
    return data["tool"]["setuptools"]["package-data"]["ssign_app"]


def _tracked_data_files() -> list[str]:
    """Non-.py files under the package that git tracks (paths relative to
    _PKG_DIR). The wheel ships tracked files only, so an untracked local file
    (e.g. a gitignored dev leftover) can't be silently dropped and must not
    count against package-data coverage."""
    out = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=str(_PKG_DIR),
        capture_output=True,
        text=True,
        check=True,
    )
    return [
        f
        for f in out.stdout.split("\0")
        if f and not f.endswith(".py") and not f.endswith(".pyc") and "__pycache__" not in f.split("/")
    ]


def test_all_data_files_covered_by_package_data():
    """Fail if a shipped (git-tracked) data file is not matched by any glob."""
    covered: set[str] = set()
    for pattern in _declared_globs():
        # root_dir + non-recursive glob mirrors setuptools' segment-bounded `*`.
        covered.update(glob.glob(pattern, root_dir=str(_PKG_DIR)))

    missing = sorted(f for f in _tracked_data_files() if f not in covered)
    assert not missing, (
        "Data files under src/ssign_app not covered by "
        "[tool.setuptools.package-data]; they would be dropped from the wheel: "
        f"{missing}"
    )


def test_coefficients_json_specifically_covered():
    """Explicit anchor for the exact file whose omission caused the bug."""
    covered: set[str] = set()
    for pattern in _declared_globs():
        covered.update(glob.glob(pattern, root_dir=str(_PKG_DIR)))
    assert "runtime/coefficients.json" in covered
