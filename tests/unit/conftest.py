"""Shared fixtures for ssign unit tests."""

import os
import sys
import tempfile

import pytest

# Add bin/ to sys.path so ssign_lib is importable
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))


@pytest.fixture
def tmp_dir():
    """Provide a temporary directory that is cleaned up after test."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def sample_fasta(tmp_dir):
    """Write a small FASTA file and return its path."""
    path = os.path.join(tmp_dir, "test.fasta")
    with open(path, 'w') as f:
        f.write(">protein_A description text\n")
        f.write("MKTLLLTLLCAFSVAQA\n")
        f.write("VDLPTQEPALGK\n")
        f.write(">protein_B\n")
        f.write("MFVFLVLLPLVSSQ\n")
    return path


@pytest.fixture
def sample_pdb_01_scale():
    """PDB content with pLDDT on 0-1 scale (HuggingFace ESMFold)."""
    # B-factor column is chars 60:66
    # Exactly 66 chars triggers the >= 66 guard
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
    # This is the minimum valid ATOM line — B-factor ends at column 66
    return "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 75.00"
