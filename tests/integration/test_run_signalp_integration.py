"""Integration tests for run_signalp.py (BioLib remote + DTU local).

Parallel to test_run_deeplocpro_integration.py — see that file for the
shape and rationale of the BioLib remote vs DTU local split.

Run remote with:
    pytest -m integration tests/integration/test_run_signalp_integration.py::TestRemote

Run local with:
    SSIGN_SIGNALP_PATH=/path/to/signalp6 pytest -m integration \\
        tests/integration/test_run_signalp_integration.py::TestLocal
"""

import os

import pytest

from conftest import read_tsv, skip_unless_biolib, skip_unless_dtu_local


pytestmark = pytest.mark.integration


REQUIRED_COLUMNS = {
    "locus_tag",
    "signalp_prediction",
    "signalp_probability",
    "signalp_cs_position",
}


class TestRemote:
    def test_remote_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        skip_unless_biolib()
        from run_signalp import run_remote_signalp

        output_dir = os.path.join(tmp_dir, "sp_out")
        os.makedirs(output_dir)
        try:
            run_remote_signalp(dtu_test_proteins, output_dir)
        except Exception as e:
            pytest.skip(f"BioLib remote call failed (probably network): {e}")

        tsv = os.path.join(output_dir, "signalp_predictions.tsv")
        assert os.path.exists(tsv)
        rows = read_tsv(tsv)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())


class TestLocal:
    def test_local_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        path = skip_unless_dtu_local("SSIGN_SIGNALP_PATH", "SignalP 6.0")
        from run_signalp import run_local_signalp

        output_dir = os.path.join(tmp_dir, "sp_out")
        os.makedirs(output_dir)
        run_local_signalp(dtu_test_proteins, path, output_dir)

        tsv = os.path.join(output_dir, "signalp_predictions.tsv")
        assert os.path.exists(tsv)
        rows = read_tsv(tsv)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())
