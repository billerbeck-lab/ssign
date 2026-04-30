"""Per-tool integration test for PLM-Effector.

Runs the real vendored PLM-Effector pipeline against the T1SS fixture
proteins. Skips cleanly if any prerequisite is missing:

    - `transformers`, `xgboost` must be installed (pip install ssign[extended])
    - SSIGN_PLM_EFFECTOR_WEIGHTS env var must point to a real weights
      directory with the expected layout (see run_plm_effector.py docstring)
    - a CUDA GPU must be available (PLM-Effector is impractical on CPU)

Run explicitly with:
    SSIGN_PLM_EFFECTOR_WEIGHTS=/path/to/weights pytest -m integration \\
        tests/integration/test_run_plm_effector_integration.py
"""

import csv
import importlib
import os
import shutil

import pytest

pytestmark = pytest.mark.integration


def _skip_unless_plm_effector_deps_and_weights():
    """Skip the current test unless every PLM-Effector prerequisite is set up."""
    for pkg in ("transformers", "xgboost", "torch"):
        if importlib.util.find_spec(pkg) is None:
            pytest.skip(f"{pkg} not installed; run `pip install ssign[extended]`")

    weights_dir = os.environ.get("SSIGN_PLM_EFFECTOR_WEIGHTS")
    if not weights_dir:
        pytest.skip("SSIGN_PLM_EFFECTOR_WEIGHTS env var not set")
    if not os.path.isdir(weights_dir):
        pytest.skip(f"SSIGN_PLM_EFFECTOR_WEIGHTS={weights_dir} does not exist")

    # Quick structural check; skip rather than fail on incomplete layouts
    expected = [
        os.path.join(weights_dir, "transformers_pretrained"),
        os.path.join(weights_dir, "trained_models"),
    ]
    for path in expected:
        if not os.path.isdir(path):
            pytest.skip(f"Weights dir missing sub-path: {path}")

    import torch

    if not torch.cuda.is_available():
        pytest.skip("CUDA GPU not available (PLM-Effector is impractical on CPU)")

    return weights_dir


class TestRunPlmEffectorOnFixture:
    def test_t1se_predictions_produced_for_fixture(
        self, tmp_dir, t1ss_fixture_proteins
    ):
        """End-to-end T1SE prediction run on the 179-protein fixture.

        Asserts the output TSV has the expected shape: one row per protein
        plus a header, correct columns, and the known substrate locus tag
        appears. Does not assert the threshold flag because that depends on
        the author's calibration — we just care the pipeline runs.
        """
        weights_dir = _skip_unless_plm_effector_deps_and_weights()

        from plm_effector import predict

        out_path = os.path.join(tmp_dir, "t1se_predictions.tsv")
        n_positive = predict(
            proteins_fasta=t1ss_fixture_proteins,
            weights_dir=weights_dir,
            effector_type="T1SE",
            out_path=out_path,
            device="cuda",
            batch_size=2,
        )

        assert os.path.exists(out_path), "predict() did not write output TSV"
        assert n_positive >= 0

        with open(out_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            assert reader.fieldnames[0] == "seq_id"
            assert "stacking" in reader.fieldnames
            assert "passes_threshold" in reader.fieldnames
            assert "effector_type" in reader.fieldnames

        # One row per input CDS. Minimal fixture: 9 CDS. Full fixture: 179.
        assert 5 <= len(rows) <= 300

        # Every row is T1SE and the flag is 0/1
        for row in rows:
            assert row["effector_type"] == "T1SE"
            assert row["passes_threshold"] in ("0", "1")

        # BIMENO_04457 is the expected T1SS substrate on this contig —
        # it should appear as an output row (not necessarily passing
        # threshold; we only assert it was predicted on).
        locus_tags = {row["seq_id"] for row in rows}
        assert "BIMENO_04457" in locus_tags

    def test_wrapper_script_runs_end_to_end(self, tmp_dir, t1ss_fixture_proteins):
        """Invoke run_plm_effector.py as a subprocess to mimic ssign's runner."""
        weights_dir = _skip_unless_plm_effector_deps_and_weights()

        import subprocess
        import sys

        script = shutil.which("run_plm_effector.py")
        if script is None:
            script_path = os.path.abspath(
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "..",
                    "src",
                    "ssign_app",
                    "scripts",
                    "run_plm_effector.py",
                )
            )
            cmd = [sys.executable, script_path]
        else:
            cmd = [script]

        out_path = os.path.join(tmp_dir, "subprocess_predictions.tsv")
        cmd.extend(
            [
                "--input",
                t1ss_fixture_proteins,
                "--weights-dir",
                weights_dir,
                "--effector-type",
                "T1SE",
                "--out",
                out_path,
                "--device",
                "cuda",
                "--batch-size",
                "2",
            ]
        )

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        assert result.returncode == 0, (
            f"run_plm_effector.py exited {result.returncode}\n"
            f"stdout: {result.stdout[-500:]}\n"
            f"stderr: {result.stderr[-500:]}"
        )
        assert os.path.exists(out_path)
        assert os.path.getsize(out_path) > 0
