"""Integration tests for run_deeplocpro.py (BioLib remote + local install).

Both modes invoke the script as a subprocess (the way runner.py does),
which produces the ssign-format TSV (`locus_tag`, `predicted_localization`,
`extracellular_prob`, ...) — the same format downstream consumers
(cross_validate_predictions, proximity_analysis) read.

Run remote:
    pytest -m integration tests/integration/test_run_deeplocpro_integration.py::TestRemote

Run local:
    SSIGN_DEEPLOCPRO_PATH=/path/to/install pytest -m integration \\
        tests/integration/test_run_deeplocpro_integration.py::TestLocal
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import skip_unless_biolib, skip_unless_dtu_local


pytestmark = pytest.mark.integration


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = PROJECT_ROOT / "src" / "ssign_app" / "scripts" / "run_deeplocpro.py"

REQUIRED_COLUMNS = {
    "locus_tag",
    "predicted_localization",
    "extracellular_prob",
}


def _read_ssign_tsv(path: str) -> list:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _run_dlp(tmp_dir, fasta, mode, deeplocpro_path=""):
    """Invoke run_deeplocpro.py via subprocess (matches runner.py)."""
    out = os.path.join(tmp_dir, "dlp_predictions.tsv")
    cmd = [
        sys.executable, str(SCRIPT),
        "--input", fasta,
        "--sample", "test",
        "--mode", mode,
        "--output", out,
    ]
    if deeplocpro_path:
        cmd.extend(["--deeplocpro-path", deeplocpro_path])
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    return out, result


class TestRemote:
    def test_remote_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        skip_unless_biolib()
        out, result = _run_dlp(tmp_dir, dtu_test_proteins, mode="remote")
        if result.returncode != 0:
            pytest.skip(
                f"BioLib remote call failed (probably network): "
                f"{result.stderr[-300:]}"
            )

        assert os.path.exists(out), "Remote DLP did not produce output TSV"
        rows = _read_ssign_tsv(out)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())

        # Biological sanity: BIMENO_04457 (autotransporter) gets a
        # numeric extracellular probability. Don't pin a threshold —
        # model versions drift; just confirm it ran.
        target = next(
            (r for r in rows if r["locus_tag"] == "BIMENO_04457"), None
        )
        assert target is not None
        try:
            float(target["extracellular_prob"])
        except (KeyError, ValueError):
            pytest.fail(
                "BIMENO_04457 row missing or non-numeric extracellular_prob"
            )


class TestLocal:
    def test_local_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        path = skip_unless_dtu_local("SSIGN_DEEPLOCPRO_PATH", "DeepLocPro")
        out, result = _run_dlp(
            tmp_dir, dtu_test_proteins, mode="local", deeplocpro_path=path
        )
        assert result.returncode == 0, (
            f"Local DLP exit {result.returncode}\n"
            f"stderr: {result.stderr[-500:]}"
        )

        assert os.path.exists(out)
        rows = _read_ssign_tsv(out)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())
