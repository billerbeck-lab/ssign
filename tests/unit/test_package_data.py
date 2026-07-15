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


def test_all_data_files_covered_by_package_data():
    """Fail if a shipped data file (non-.py) is not matched by any glob."""
    covered: set[str] = set()
    for pattern in _declared_globs():
        # root_dir + non-recursive glob mirrors setuptools' segment-bounded `*`.
        covered.update(glob.glob(pattern, root_dir=str(_PKG_DIR)))

    data_files = [
        p.relative_to(_PKG_DIR).as_posix()
        for p in _PKG_DIR.rglob("*")
        if p.is_file() and p.suffix != ".py" and "__pycache__" not in p.parts and not p.name.endswith(".pyc")
    ]

    missing = sorted(f for f in data_files if f not in covered)
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
