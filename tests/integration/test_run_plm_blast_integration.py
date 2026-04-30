"""Per-tool integration test for run_plm_blast.py.

Runs the real pLM-BLAST `plmblast.py` script against the T1SS fixture
proteins using an ECOD70 database. Skips cleanly if either prerequisite
is missing.

Run explicitly with:
    SSIGN_ECOD70_DB=/path/to/ecod70 pytest -m integration \\
        tests/integration/test_run_plm_blast_integration.py
"""

import os
import shutil

import pytest
from run_plm_blast import parse_plmblast_csv, run_plmblast  # noqa: E402

pytestmark = pytest.mark.integration


@pytest.fixture
def ecod70_db():
    """Path to a pLM-BLAST ECOD70 database. Skips the test if pLM-BLAST's
    `plmblast.py` (or `plmblast`) is not on PATH, SSIGN_PLMBLAST_SCRIPT
    does not point at a real file, or SSIGN_ECOD70_DB is unset / missing."""
    script_override = os.environ.get("SSIGN_PLMBLAST_SCRIPT")
    if script_override:
        if not os.path.exists(script_override):
            pytest.skip(f"SSIGN_PLMBLAST_SCRIPT={script_override} does not exist")
    elif shutil.which("plmblast.py") is None and shutil.which("plmblast") is None:
        pytest.skip(
            "pLM-BLAST not on PATH; install with "
            "`pip install git+https://github.com/labstructbioinf/pLM-BLAST.git` "
            "or set SSIGN_PLMBLAST_SCRIPT"
        )

    db = os.environ.get("SSIGN_ECOD70_DB")
    if not db:
        pytest.skip("SSIGN_ECOD70_DB env var not set")
    if not os.path.isdir(db):
        pytest.skip(f"SSIGN_ECOD70_DB={db} is not a directory")
    return db


class TestRunPlmBlastOnFixture:
    def test_full_pipeline_on_fixture(
        self, tmp_dir, t1ss_fixture_proteins, ecod70_db
    ):
        """One run, all asserts — pLM-BLAST embed+search on the
        minimal fixture takes ~10-15 min on CPU; running it once and
        asserting all invariants is much faster than the previous
        two-test split that re-embedded each time.
        """
        out_csv = os.path.join(tmp_dir, "plm_blast.csv")
        run_plmblast(
            proteins_fasta=t1ss_fixture_proteins,
            ecod_db=ecod70_db,
            out_csv=out_csv,
            cpc=70,
            threads=4,
        )

        # Invariant 1: file produced + non-empty.
        assert os.path.exists(out_csv), "pLM-BLAST did not write an output CSV"
        assert os.path.getsize(out_csv) > 0, "pLM-BLAST CSV is empty"

        # Invariant 2: parser returns at least one hit. Both fixtures
        # contain BIMENO_04457 (autotransporter, has ECOD hits) plus
        # housekeeping proteins that almost always pick up some hit.
        entries = parse_plmblast_csv(out_csv)
        assert len(entries) > 0, "Expected at least one ECOD hit on the fixture"

        # Invariant 3: every entry has the columns downstream code reads.
        required = {
            "protein_id", "target_id", "score",
            "qstart", "qend", "tstart", "tend",
        }
        for e in entries:
            assert required <= set(e.keys())
            assert e["protein_id"], "Every entry must have a non-empty protein_id"
