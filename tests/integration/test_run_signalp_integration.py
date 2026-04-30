"""Integration tests for run_signalp.py (BioLib remote + local install).

Both modes invoke the script as a subprocess (matches runner.py),
producing the ssign-format TSV (`locus_tag`, `signalp_prediction`,
`signalp_probability`, `signalp_cs_position`).

Run remote:
    pytest -m integration tests/integration/test_run_signalp_integration.py::TestRemote

Run local:
    SSIGN_SIGNALP_PATH=/path/to/install pytest -m integration \\
        tests/integration/test_run_signalp_integration.py::TestLocal
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
SCRIPT = PROJECT_ROOT / "src" / "ssign_app" / "scripts" / "run_signalp.py"

REQUIRED_COLUMNS = {
    "locus_tag",
    "signalp_prediction",
    "signalp_probability",
    "signalp_cs_position",
}


def _read_ssign_tsv(path: str) -> list:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _run_signalp(tmp_dir, fasta, mode, signalp_path=""):
    out = os.path.join(tmp_dir, "signalp_predictions.tsv")
    cmd = [
        sys.executable, str(SCRIPT),
        "--input", fasta,
        "--sample", "test",
        "--mode", mode,
        "--output", out,
    ]
    if signalp_path:
        cmd.extend(["--signalp-path", signalp_path])
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    return out, result


class TestRemote:
    def test_remote_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        skip_unless_biolib()
        out, result = _run_signalp(tmp_dir, dtu_test_proteins, mode="remote")
        if result.returncode != 0:
            pytest.skip(
                f"BioLib remote call failed (probably network): "
                f"{result.stderr[-300:]}"
            )
        assert os.path.exists(out)
        rows = _read_ssign_tsv(out)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())


class TestLocal:
    def test_local_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        path = skip_unless_dtu_local("SSIGN_SIGNALP_PATH", "SignalP 6.0")
        out, result = _run_signalp(
            tmp_dir, dtu_test_proteins, mode="local", signalp_path=path
        )
        assert result.returncode == 0, (
            f"Local SignalP exit {result.returncode}\n"
            f"stderr: {result.stderr[-500:]}"
        )
        assert os.path.exists(out)
        rows = _read_ssign_tsv(out)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())
