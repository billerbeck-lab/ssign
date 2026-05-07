"""Shared fixtures for ssign unit tests.

Conventions:
- Bare-function helpers (writers, builders, FIELD constants) live in
  `_helpers.py` so test files can import them directly. This file
  exposes pytest fixtures that wrap those helpers for the common cases.
- Pre-baked fixtures cover a 15-gene two-contig genome and the matching
  gene_order / SS-components / predictions TSVs. Tests that need
  variations override individual rows or rebuild via the writers.
- pdb_* fixtures are kept for legacy callers; no scripts currently use
  them but they pre-date the Foldseek/ESM3 removal and are cheap to keep.

Hoisted import-path setup: most pipeline scripts live in
src/ssign_app/scripts/ as plain modules (e.g. `import run_blastp`) and
also import sibling modules without the `ssign_app.scripts.` prefix.
Adding that directory to sys.path here once removes the 7-line
SCRIPTS_DIR boilerplate from ~25 individual test files.
"""

import os
import sys
import tempfile

_SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

import pytest  # noqa: E402
from _helpers import (  # noqa: E402
    GENE_INFO_FIELDS,
    GENE_ORDER_FIELDS,
    SS_COMPONENT_FIELDS,
    gene_info_to_gene_order,
    two_contig_genes,
    write_tsv,
)

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta  # noqa: E402


@pytest.fixture
def tmp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def sample_fasta(tmp_dir):
    """Two-protein FASTA, multi-line first sequence, single-line second."""
    path = os.path.join(tmp_dir, "test.fasta")
    write_fasta(
        {
            "protein_A": "MKTLLLTLLCAFSVAQAVDLPTQEPALGK",
            "protein_B": "MFVFLVLLPLVSSQ",
        },
        path,
    )
    return path


@pytest.fixture
def gene_info_rows():
    """Override-friendly: list of dict rows for the 15-gene 2-contig genome."""
    return two_contig_genes()


@pytest.fixture
def gene_info_tsv(tmp_dir, gene_info_rows):
    """gene_info.tsv path with the 15-gene 2-contig layout."""
    return write_tsv(
        os.path.join(tmp_dir, "gene_info.tsv"),
        GENE_INFO_FIELDS,
        gene_info_rows,
    )


@pytest.fixture
def gene_order_tsv(tmp_dir, gene_info_rows):
    """gene_order.tsv path. Same layout as gene_info but with per-contig gene_index."""
    return write_tsv(
        os.path.join(tmp_dir, "gene_order.tsv"),
        GENE_ORDER_FIELDS,
        gene_info_to_gene_order(gene_info_rows),
    )


@pytest.fixture
def ss_components_tsv(tmp_dir):
    """Two T2SS components on contig_A (genes 5 + 6); none excluded.
    Schema mirrors validate_macsyfinder_systems.py output."""
    base = {
        "sample_id": "test_sample",
        "sys_id": "contig_A_T2SS_1",
        "ss_type": "T2SS",
        "gene_status": "mandatory",
        "wholeness": "1.0",
        "excluded": "False",
    }
    rows = [
        {**base, "locus_tag": "GENE_0005", "gene_name": "gspD"},
        {**base, "locus_tag": "GENE_0006", "gene_name": "gspE"},
    ]
    return write_tsv(
        os.path.join(tmp_dir, "ss_components.tsv"),
        SS_COMPONENT_FIELDS,
        rows,
    )


# ---------------------------------------------------------------------------
# pLDDT scale fixtures (legacy; kept for older callers)
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_pdb_01_scale():
    """PDB content with pLDDT on 0-1 scale (HuggingFace ESMFold)."""
    return (
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.85           N\n"
        "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.72           C\n"
        "ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.90           C\n"
        "ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00  0.68           O\n"
        "END\n"
    )


@pytest.fixture
def sample_pdb_100_scale():
    """PDB content with pLDDT on 0-100 scale (AlphaFold DB)."""
    return (
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 85.00           N\n"
        "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 72.00           C\n"
        "ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 90.00           C\n"
        "ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00 68.00           O\n"
        "END\n"
    )


@pytest.fixture
def pdb_exact_66_chars():
    """PDB line that is exactly 66 characters (tests >= 66 guard, not > 66)."""
    return "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 75.00"
