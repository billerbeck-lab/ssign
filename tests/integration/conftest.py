"""Shared fixtures + skip logic for integration tests.

Every test in `tests/integration/` assumes an external bioinformatics tool
may or may not be installed on the machine running pytest. Tests are
marked with `@pytest.mark.integration` and skip automatically when their
prerequisites are missing, so a default `pytest tests/` run on a clean
dev machine still produces a green suite.

Run the integration suite explicitly with:
    pytest -m integration tests/integration/
"""

import os
import shutil
import sys
import tempfile

import pytest

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)


T1SS_FIXTURE_GBFF = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "fixtures",
        "Xanthobacter_tagetidis_ATCC_700314_contig_87.gbff",
    )
)


@pytest.fixture
def tmp_dir():
    """Temp directory auto-cleaned after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def t1ss_fixture_gbff():
    """Path to the T1SS fixture GenBank (Xanthobacter contig_87, 213 kb, 179 CDS)."""
    if not os.path.exists(T1SS_FIXTURE_GBFF):
        pytest.skip(f"T1SS fixture missing: {T1SS_FIXTURE_GBFF}")
    return T1SS_FIXTURE_GBFF


@pytest.fixture
def t1ss_fixture_contigs(tmp_dir, t1ss_fixture_gbff):
    """Extract DNA from the T1SS fixture into a temporary contigs FASTA.

    Used by tools (Bakta, Prodigal) that expect raw contig input.
    """
    from Bio import SeqIO

    out_path = os.path.join(tmp_dir, "contigs.fasta")
    records = list(SeqIO.parse(t1ss_fixture_gbff, "genbank"))
    with open(out_path, "w") as f:
        for rec in records:
            f.write(f">{rec.id}\n{str(rec.seq)}\n")
    return out_path


@pytest.fixture
def t1ss_fixture_proteins(tmp_dir, t1ss_fixture_gbff):
    """Extract protein translations from the T1SS fixture into a temp FASTA.

    Used by tools (EggNOG-mapper, InterProScan, HH-suite) that expect a
    protein FASTA input.
    """
    from Bio import SeqIO

    out_path = os.path.join(tmp_dir, "proteins.faa")
    n = 0
    with open(out_path, "w") as f:
        for rec in SeqIO.parse(t1ss_fixture_gbff, "genbank"):
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                translation = feat.qualifiers.get("translation", [None])[0]
                locus_tag = feat.qualifiers.get("locus_tag", [None])[0]
                if not translation or not locus_tag:
                    continue
                f.write(f">{locus_tag}\n{translation}\n")
                n += 1
    assert n > 0, "No CDS translations extracted from fixture"
    return out_path


def _skip_unless_tool_and_db(binary_name: str, env_var_name: str) -> str:
    """Skip the current test unless *binary_name* is on PATH and
    *env_var_name* points to an existing database path. Returns the DB path."""
    if shutil.which(binary_name) is None:
        pytest.skip(f"{binary_name} not installed or not on PATH")
    db = os.environ.get(env_var_name)
    if not db:
        pytest.skip(f"{env_var_name} env var not set")
    if not os.path.exists(db):
        pytest.skip(f"{env_var_name}={db} does not exist")
    return db


@pytest.fixture
def bakta_db():
    """Path to a Bakta database. Skips the test if `bakta` is missing or
    SSIGN_BAKTA_DB is unset / points to a non-existent path."""
    return _skip_unless_tool_and_db("bakta", "SSIGN_BAKTA_DB")


@pytest.fixture
def emapper_db():
    """Path to an EggNOG database. Skips the test if `emapper.py` is missing
    or SSIGN_EGGNOG_DB is unset / points to a non-existent path."""
    return _skip_unless_tool_and_db("emapper.py", "SSIGN_EGGNOG_DB")
