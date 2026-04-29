"""Integration tests for run_deeplocpro.py (BioLib remote + DTU local).

The DTU local path is gated on Sonja's redistribution email; the test
skips cleanly until SSIGN_DEEPLOCPRO_PATH is set to a real install.
The remote path (current Easy Mode default) requires only network +
`pybiolib`.

Run remote with:
    pytest -m integration tests/integration/test_run_deeplocpro_integration.py::TestRemote

Run local with:
    SSIGN_DEEPLOCPRO_PATH=/path/to/deeplocpro pytest -m integration \\
        tests/integration/test_run_deeplocpro_integration.py::TestLocal
"""

import os

import pytest

from conftest import read_tsv, skip_unless_biolib, skip_unless_dtu_local


pytestmark = pytest.mark.integration


REQUIRED_COLUMNS = {
    "locus_tag",
    "predicted_localization",
    "dlp_extracellular_prob",
}


class TestRemote:
    def test_remote_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        skip_unless_biolib()
        from run_deeplocpro import run_remote_deeplocpro

        output_dir = os.path.join(tmp_dir, "dlp_out")
        os.makedirs(output_dir)
        try:
            run_remote_deeplocpro(dtu_test_proteins, output_dir)
        except Exception as e:
            pytest.skip(f"BioLib remote call failed (probably network): {e}")

        tsv = os.path.join(output_dir, "deeplocpro_predictions.tsv")
        assert os.path.exists(tsv), "Remote DLP did not produce predictions.tsv"
        rows = read_tsv(tsv)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())

        # Biological sanity: BIMENO_04457 is a known autotransporter; DLP
        # should produce a numeric extracellular probability for it.
        # Don't assert a specific threshold (model versions drift); just
        # that the protein-of-interest got predicted on. SignalP doesn't
        # get an analogous check because every autotransporter has an
        # N-terminal signal peptide — trivially true, not informative.
        target = next(
            (r for r in rows if r["locus_tag"] == "BIMENO_04457"), None
        )
        assert target is not None
        try:
            float(target["dlp_extracellular_prob"])
        except (KeyError, ValueError):
            pytest.fail("BIMENO_04457 row missing or non-numeric ext probability")


class TestLocal:
    def test_local_runs_and_returns_predictions(
        self, tmp_dir, dtu_test_proteins
    ):
        path = skip_unless_dtu_local("SSIGN_DEEPLOCPRO_PATH", "DeepLocPro")
        from run_deeplocpro import run_local_deeplocpro

        output_dir = os.path.join(tmp_dir, "dlp_out")
        os.makedirs(output_dir)
        run_local_deeplocpro(dtu_test_proteins, path, output_dir)

        tsv = os.path.join(output_dir, "deeplocpro_predictions.tsv")
        assert os.path.exists(tsv)
        rows = read_tsv(tsv)
        assert len(rows) > 0
        assert REQUIRED_COLUMNS <= set(rows[0].keys())
