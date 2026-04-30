"""Per-tool integration test for run_bakta.py.

Runs the real `bakta` binary against the T1SS fixture genome. Skips
cleanly if Bakta is not installed or if SSIGN_BAKTA_DB is not set to a
valid database path.

Run explicitly with:
    SSIGN_BAKTA_DB=/path/to/bakta_db pytest -m integration \\
        tests/integration/test_run_bakta_integration.py
"""

import os

import pytest
from run_bakta import parse_bakta_tsv, run_bakta  # noqa: E402

pytestmark = pytest.mark.integration


class TestRunBaktaOnFixture:
    def test_bakta_runs_and_produces_ssign_outputs(
        self, tmp_dir, t1ss_fixture_contigs, bakta_db
    ):
        proteins_faa, tsv_path = run_bakta(
            contigs_fasta=t1ss_fixture_contigs,
            db_path=bakta_db,
            sample_id="t1ss_fixture",
            output_dir=tmp_dir,
            threads=4,
        )

        assert os.path.exists(proteins_faa), "Bakta did not produce a protein FASTA"
        assert os.path.exists(tsv_path), "Bakta did not produce a TSV annotation table"
        assert os.path.getsize(proteins_faa) > 0, "Bakta protein FASTA is empty"
        assert os.path.getsize(tsv_path) > 0, "Bakta TSV annotation is empty"

    def test_parsed_entries_match_fixture_scale(
        self, tmp_dir, t1ss_fixture_contigs, bakta_db
    ):
        """Parsed CDS count should be in the right ballpark for a 213 kb contig."""
        _, tsv_path = run_bakta(
            contigs_fasta=t1ss_fixture_contigs,
            db_path=bakta_db,
            sample_id="t1ss_fixture",
            output_dir=tmp_dir,
            threads=4,
        )

        entries = parse_bakta_tsv(tsv_path)
        # Minimal fixture: 9 CDS in 20 kb. Full fixture: 179 CDS in 213 kb.
        # Wide bounds accommodate Bakta version + DB variation.
        assert 5 <= len(entries) <= 300, (
            f"Expected 5-15 CDS (minimal) or 100-300 CDS (full), got {len(entries)}"
        )

    def test_output_rows_have_required_fields(
        self, tmp_dir, t1ss_fixture_contigs, bakta_db
    ):
        _, tsv_path = run_bakta(
            contigs_fasta=t1ss_fixture_contigs,
            db_path=bakta_db,
            sample_id="t1ss_fixture",
            output_dir=tmp_dir,
            threads=4,
        )

        entries = parse_bakta_tsv(tsv_path)
        required = {"locus_tag", "product", "contig", "start", "end", "strand"}
        for e in entries:
            assert required <= set(e.keys())
            assert e["locus_tag"], "Every entry must have a locus tag"
            assert isinstance(e["start"], int)
            assert isinstance(e["end"], int)
