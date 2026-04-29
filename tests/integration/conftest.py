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


_FIXTURES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "fixtures")
)
T1SS_FIXTURE_GBFF = os.path.join(
    _FIXTURES_DIR, "Xanthobacter_tagetidis_ATCC_700314_contig_87.gbff"
)
# Tighter ~20 kb / 9 CDS window around BIMENO_04457 (T5aSS autotransporter)
# + flanking virulence factors. Use this for fast wrapper iteration on
# CPU-only machines; the full 213 kb fixture takes ~6 hours just on
# EggNOG --sensitive --iterate.
T5ASS_MINIMAL_FIXTURE_GBFF = os.path.join(
    _FIXTURES_DIR, "Xanthobacter_T5aSS_minimal.gbff"
)


@pytest.fixture
def tmp_dir():
    """Temp directory auto-cleaned after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def t1ss_fixture_gbff():
    """Path to the GenBank fixture used by integration tests.

    Defaults to the minimal ~20 kb / 9 CDS window around the T5aSS
    autotransporter (BIMENO_04457). Set `SSIGN_TEST_FIXTURE_FULL=1` to
    use the original 213 kb / 179 CDS fixture instead — useful for
    cohort-scale and biological-coverage tests, but ~10-50x slower.
    """
    if os.environ.get("SSIGN_TEST_FIXTURE_FULL") == "1":
        path = T1SS_FIXTURE_GBFF
    else:
        path = T5ASS_MINIMAL_FIXTURE_GBFF
    if not os.path.exists(path):
        pytest.skip(f"Fixture missing: {path}")
    return path


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


# ── DTU prediction-tool tests (DeepLocPro, SignalP) ─────────────────────
# Shared helpers for the BioLib-remote and DTU-local-install integration
# tests. Both DLP and SignalP have the same shape: small protein FASTA
# in, TSV out; remote path needs `biolib`, local path needs a binary.

# A single autotransporter sequence (BIMENO_04457 from the T5aSS fixture)
# — minimal viable input. Biolib charges by compute, so the smaller the
# payload the cheaper. One real-biology record exercises the wrapper +
# parser end-to-end and lets us sanity-check the localization / signal-
# peptide call against a known autotransporter.
_DTU_TEST_PROTEINS = (
    ">BIMENO_04457\n"
    "MKKVALLAALLPALIPSAQAATTHQKVKVGKEEFNQKAQHPDAQTLAYLAEQCNLKKTSAALDIEKHLIKQDKAS\n"
    "AKWNRELDQKGQALAAAFKEKQRADIIDPDIRAKAPATERRAVQYHIARLERHLFQNALSRRADLRKSQGDQSGV\n"
)


@pytest.fixture
def dtu_test_proteins(tmp_dir):
    """Single-record protein FASTA used by DLP + SignalP integration tests."""
    path = os.path.join(tmp_dir, "proteins.faa")
    with open(path, "w") as f:
        f.write(_DTU_TEST_PROTEINS)
    return path


def skip_unless_biolib():
    """Skip a test unless `pybiolib` is importable. Used by DLP+SignalP
    remote-mode integration tests."""
    import importlib.util

    if importlib.util.find_spec("biolib") is None:
        pytest.skip(
            "pybiolib not installed; install with `pip install ssign[extended]`"
        )


def skip_unless_dtu_local(env_var: str, tool_label: str) -> str:
    """Skip unless *env_var* points to an existing DTU-binary install path."""
    path = os.environ.get(env_var)
    if not path:
        pytest.skip(f"{env_var} env var not set ({tool_label} local install)")
    if not os.path.exists(path):
        pytest.skip(f"{env_var}={path} does not exist")
    return path


def read_tsv(path: str) -> list:
    """Parse a tab-separated file with a header row into list[dict]."""
    import csv

    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))
