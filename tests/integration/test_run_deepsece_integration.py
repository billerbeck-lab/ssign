"""Integration test for run_deepsece.py.

DeepSecE is a single ESM-1b-derived classifier predicting one of
{T1SE, T2SE, T3SE, T4SE, T6SE, Non-secreted} per input protein. The
wrapper auto-downloads its checkpoint to ~/.ssign/models/ on first
run from a Zenodo URL that is currently a `PLACEHOLDER` (v1.0.0
release blocker — see design_decisions.md). The SJTU origin URL
(`tool2-mml.sjtu.edu.cn/DeepSecE/checkpoint.pt`) is the working
fallback at time of writing.

Test skips cleanly unless:
  - `deepsece` package is importable (pip install ssign[extended]
    or pip install deepsece), AND
  - the checkpoint exists at SSIGN_DEEPSECE_CHECKPOINT or
    ~/.ssign/models/deepsece_checkpoint.pt

Run with:
    pytest -m integration tests/integration/test_run_deepsece_integration.py
"""

import csv
import importlib
import os

import pytest

pytestmark = pytest.mark.integration


def _skip_unless_deepsece():
    # The pip package installs as `DeepSecE` (capital D), not lowercase.
    if importlib.util.find_spec("DeepSecE") is None:
        pytest.skip("DeepSecE not installed; pip install ssign[extended]")
    ckpt = os.environ.get(
        "SSIGN_DEEPSECE_CHECKPOINT",
        os.path.expanduser("~/.ssign/models/deepsece_checkpoint.pt"),
    )
    if not os.path.exists(ckpt):
        pytest.skip(
            f"DeepSecE checkpoint not found at {ckpt}. Download with: "
            f"curl -L -o {ckpt} https://tool2-mml.sjtu.edu.cn/DeepSecE/checkpoint.pt"
        )
    return ckpt


# Schema written by run_deepsece(): raw upstream column names. The
# parse_deepsece_output() helper renames these to ssign convention
# (`locus_tag`, `dse_ss_type`, `dse_max_prob`) at consumption time.
RAW_REQUIRED_COLUMNS = {
    "protein_id",
    "deepsece_prediction",
    "deepsece_ss_type",
    "max_prob",
    "T1_prob",
    "T2_prob",
    "T3_prob",
    "T4_prob",
    "T6_prob",
    "nonsec_prob",
}


class TestRunDeepSecE:
    def test_full_pipeline_on_fixture(
        self, tmp_dir, t1ss_fixture_proteins
    ):
        """One run, three asserts — DSE on the 9-CDS fixture takes
        ~1-2 min on CPU once ESM-1b is cached; first run is dominated
        by the ESM-1b weight download (~7.5 GB from Meta CDN, 5-15 min).
        """
        ckpt = _skip_unless_deepsece()
        from run_deepsece import parse_deepsece_output, run_deepsece

        output_dir = os.path.join(tmp_dir, "dse_out")
        os.makedirs(output_dir)
        run_deepsece(
            input_fasta=t1ss_fixture_proteins,
            output_dir=output_dir,
            checkpoint_path=ckpt,
            batch_size=1,
        )

        # The wrapper writes deepsece_predictions.csv with raw upstream
        # column names. Validate both the raw schema and the ssign-format
        # produced by parse_deepsece_output (consumed by runner.py).
        out_path = os.path.join(output_dir, "deepsece_predictions.csv")
        assert os.path.exists(out_path), "DeepSecE did not produce output"

        with open(out_path) as f:
            raw_rows = list(csv.DictReader(f))
        assert len(raw_rows) > 0
        missing = RAW_REQUIRED_COLUMNS - set(raw_rows[0].keys())
        assert not missing, (
            f"Raw DSE output missing columns: {missing}. Schema may have "
            f"drifted in run_deepsece's writer."
        )

        # Invariant 1: DSE T3SS unreliability — per CLAUDE.md Critical
        # Bug Fix #5, DSE often misclassifies flagellar proteins as T3SS.
        # The fixture is a T5aSS region with no flagellum, so few or no
        # rows should be T3SS. Sanity bound at most ~1/3 of input rows.
        t3ss_count = sum(1 for r in raw_rows if r["deepsece_ss_type"] == "T3SS")
        assert t3ss_count <= 3, (
            f"DSE called {t3ss_count}/{len(raw_rows)} as T3SS — flagellar "
            f"misclassification regression?"
        )

        # Invariant 2: parse_deepsece_output renames to ssign convention
        # (locus_tag, dse_ss_type, dse_max_prob). Returns a list of
        # dicts (one per protein), not a dict-keyed-by-locus.
        parsed = parse_deepsece_output(out_path)
        assert isinstance(parsed, list)
        assert len(parsed) > 0
        ssign_required = {"locus_tag", "dse_ss_type", "dse_max_prob"}
        assert ssign_required <= set(parsed[0].keys()), (
            f"parse_deepsece_output missing ssign columns: "
            f"{ssign_required - set(parsed[0].keys())}"
        )
        target = next(
            (e for e in parsed if e["locus_tag"] == "BIMENO_04457"), None
        )
        assert target is not None, (
            "BIMENO_04457 missing from parsed output — wrapper dropped "
            "it or column-rename map drifted."
        )
