"""Per-tool integration test for run_eggnog.py.

Runs the real `emapper.py` binary against the T1SS fixture proteins. Skips
cleanly if EggNOG-mapper is not installed or if SSIGN_EGGNOG_DB is not
set to a valid database path.

Run explicitly with:
    SSIGN_EGGNOG_DB=/path/to/eggnog_db pytest -m integration \\
        tests/integration/test_run_eggnog_integration.py
"""

import os

import pytest
from run_eggnog import parse_eggnog_annotations, run_emapper  # noqa: E402

pytestmark = pytest.mark.integration


class TestRunEmapperOnFixture:
    def test_emapper_runs_and_produces_annotations(
        self, tmp_dir, t1ss_fixture_proteins, emapper_db
    ):
        annotations_path = run_emapper(
            proteins_fasta=t1ss_fixture_proteins,
            db_path=emapper_db,
            sample_id="t1ss_fixture",
            output_dir=tmp_dir,
            threads=4,
        )

        assert os.path.exists(annotations_path), (
            "emapper did not produce the .annotations file"
        )
        assert os.path.getsize(annotations_path) > 0, "emapper annotations is empty"

    def test_parsed_entries_have_valid_shape(
        self, tmp_dir, t1ss_fixture_proteins, emapper_db
    ):
        """At least a handful of the ~170-180 CDS on the fixture should match
        an eggNOG orthologous group; we don't pin an exact count because it
        depends on which DB subset is installed."""
        annotations_path = run_emapper(
            proteins_fasta=t1ss_fixture_proteins,
            db_path=emapper_db,
            sample_id="t1ss_fixture",
            output_dir=tmp_dir,
            threads=4,
        )

        entries = parse_eggnog_annotations(annotations_path)
        # Minimal fixture (9 proteins) typically yields 3-6 OG annotations;
        # full fixture (179 proteins) yields 100+. >2 is the lower bound that
        # validates emapper actually annotated something.
        assert len(entries) > 2, (
            f"Expected >2 eggNOG annotations, got {len(entries)}"
        )

        required = {
            "protein_id",
            "seed_ortholog",
            "evalue",
            "description",
            "preferred_name",
        }
        for e in entries:
            assert required <= set(e.keys())
            assert e["protein_id"], "Every entry must have a non-empty protein_id"
