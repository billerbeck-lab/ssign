"""Integration test for run_interproscan.py.

InterProScan is large (~7 GB tarball, ~12-24 GB after member-DB
setup) and slow to install, so this test skips unless the user has
already set it up and pointed SSIGN_INTERPROSCAN_PATH at the
extracted directory containing interproscan.sh.

Run with:
    SSIGN_INTERPROSCAN_PATH=/path/to/interproscan-5.77-108.0 \\
    pytest -m integration tests/integration/test_run_interproscan_integration.py
"""

import csv
import os
import shutil

import pytest


pytestmark = pytest.mark.integration


def _skip_unless_interproscan():
    """Skip unless `interproscan.sh` is on PATH or
    SSIGN_INTERPROSCAN_PATH points to a directory containing it."""
    path = shutil.which("interproscan.sh")
    if not path:
        ipr_dir = os.environ.get("SSIGN_INTERPROSCAN_PATH")
        if ipr_dir and os.path.exists(os.path.join(ipr_dir, "interproscan.sh")):
            # Add to PATH for the duration of this process
            os.environ["PATH"] = ipr_dir + os.pathsep + os.environ["PATH"]
            return ipr_dir
        pytest.skip(
            "interproscan.sh not on PATH; download from "
            "https://www.ebi.ac.uk/interpro/download/ and set "
            "SSIGN_INTERPROSCAN_PATH=<install dir>"
        )
    return os.path.dirname(path)


class TestRunInterProScan:
    def test_runs_on_fixture_proteins(
        self, tmp_dir, t1ss_fixture_proteins
    ):
        _skip_unless_interproscan()

        import sys
        sys.path.insert(
            0,
            os.path.abspath(
                os.path.join(
                    os.path.dirname(__file__), "..", "..",
                    "src", "ssign_app", "scripts",
                )
            ),
        )
        from run_interproscan import run_local_interproscan

        output_dir = os.path.join(tmp_dir, "ipr_out")
        os.makedirs(output_dir)
        result_path = run_local_interproscan(
            query_fasta=t1ss_fixture_proteins,
            db_path="",  # use default DBs from the install
            output_dir=output_dir,
        )

        assert os.path.exists(result_path)
        # InterProScan TSV columns: query, md5, length, analysis, signature_acc,
        # signature_desc, start, end, evalue, status, date, ipr_acc, ipr_desc,
        # go, pathways. We don't pin all columns since member-DB output may
        # vary; just assert at least one row + that the query column matches a
        # fixture locus_tag.
        with open(result_path) as f:
            rows = list(csv.reader(f, delimiter="\t"))
        assert len(rows) > 0, "InterProScan produced no hits"
        # First column = query (locus_tag). Fixture has BIMENO_* tags.
        queries = {r[0] for r in rows if r}
        bimeno_hits = [q for q in queries if q.startswith("BIMENO_")]
        assert bimeno_hits, (
            f"No BIMENO_* protein got an InterProScan hit. Got queries: "
            f"{list(queries)[:5]}"
        )
