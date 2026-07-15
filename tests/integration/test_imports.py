"""Eager-import check: every Python dep + lazy symbol ssign reaches for.

This test exists because lazy imports inside pipeline scripts (``from
DeepSecE.model import EffectorTransformer`` inside ``run_deepsece.py``)
do not fail at ``pip install`` time. They only fail when that branch of
the pipeline actually executes — sometimes hours into a run on a multi-day
HPC job. This test imports every symbol upfront, so any ``pyproject.toml``
drift fails CI within seconds instead.

Marked ``integration`` because it requires the extended-tier extras to be
installed; a clean dev env with only the base tier would skip it.
"""

from __future__ import annotations

import os
import sys

import pytest

# Make `ssign_app.scripts.ssign_lib.dependency_manifest` importable when this
# file is run directly. The package is installed in editable mode in CI, but
# tests/integration/conftest.py also inserts the scripts dir on sys.path.
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from ssign_app.scripts.doctor import check_python_dep  # noqa: E402
from ssign_app.scripts.ssign_lib.dependency_manifest import PYTHON_DEPS  # noqa: E402

pytestmark = pytest.mark.integration


@pytest.mark.parametrize("dep", PYTHON_DEPS, ids=lambda d: d.module)
def test_python_dep_importable(dep) -> None:
    """Every entry in ``PYTHON_DEPS`` imports cleanly, including lazy symbols.

    If this fails, ``pyproject.toml`` is missing the pip package named in
    ``dep.pip_name`` (or an upstream version removed the lazy symbol —
    check the failure message for ``AttributeError`` vs ``ImportError``).
    """
    result = check_python_dep(dep)
    note_line = f"\n  note: {dep.note}" if dep.note else ""
    assert result.ok, f"{dep.module} import check failed: {result.detail}\n  fix: {result.fix}{note_line}"
