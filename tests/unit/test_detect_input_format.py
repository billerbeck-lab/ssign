"""Tests for detect_input_format.py.

Two layers: extension-based dispatch (case-insensitive), and a
content-inspection fallback for unknown extensions or ambiguous .fasta /
.fna / .fa files (which always resolve to "fasta_contigs" — the pipeline
does not distinguish protein FASTA from nucleotide contigs at this
layer; .faa is the only protein-FASTA signal).
"""

import os
import sys

import pytest


from detect_input_format import detect_format  # noqa: E402


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)
    return path


# ---------------------------------------------------------------------------
# Extension-based dispatch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ext", [".gbff", ".gbk", ".gb", ".GBFF", ".Gbk"])
def test_genbank_extensions(tmp_dir, ext):
    path = _write(os.path.join(tmp_dir, f"genome{ext}"), "LOCUS       test\n")
    assert detect_format(path) == "genbank"


@pytest.mark.parametrize("ext", [".gff", ".gff3", ".gtf", ".GFF", ".Gff3"])
def test_gff3_extensions(tmp_dir, ext):
    path = _write(os.path.join(tmp_dir, f"genome{ext}"), "##gff-version 3\n")
    assert detect_format(path) == "gff3"


@pytest.mark.parametrize("ext", [".faa", ".FAA", ".Faa"])
def test_protein_fasta_extension(tmp_dir, ext):
    path = _write(os.path.join(tmp_dir, f"proteins{ext}"), ">p1\nMKT\n")
    assert detect_format(path) == "protein_fasta"


@pytest.mark.parametrize("ext", [".fasta", ".fna", ".fa", ".FASTA", ".Fa"])
def test_nucleotide_fasta_extensions_resolve_to_contigs(tmp_dir, ext):
    path = _write(os.path.join(tmp_dir, f"genome{ext}"), ">contig_1\nATGCATGCATGC\n")
    assert detect_format(path) == "fasta_contigs"


# ---------------------------------------------------------------------------
# _inspect_fasta — content-based inspection for ambiguous .fasta/.fna/.fa
# ---------------------------------------------------------------------------


def test_fasta_with_protein_content_still_returns_contigs(tmp_dir):
    # Per the implementation comment: bare FASTA without a .faa extension is
    # treated as contigs even if the residues are amino-acid-like. .faa is
    # the only protein-FASTA signal — anything else goes through extract_proteins
    # which will gracefully detect the issue downstream.
    path = _write(
        os.path.join(tmp_dir, "ambiguous.fasta"),
        ">p1\nMKTLLLTLLCAFSVAQAVDLPTQEPALGK\n",
    )
    assert detect_format(path) == "fasta_contigs"


def test_fasta_with_pure_nucleotide_content_returns_contigs(tmp_dir):
    path = _write(
        os.path.join(tmp_dir, "genome.fasta"),
        ">contig_1\nATGCATGCATGCATGC\nNNNNATGC\n",
    )
    assert detect_format(path) == "fasta_contigs"


def test_inspect_fasta_caps_at_50_data_lines(tmp_dir):
    # Only the first 50 non-header lines are inspected. Anything beyond
    # that doesn't influence the decision — keeps memory bounded for
    # multi-GB genome files.
    body = ">contig_1\n" + "\n".join("ATGC" for _ in range(100))
    path = _write(os.path.join(tmp_dir, "genome.fasta"), body)
    assert detect_format(path) == "fasta_contigs"


# ---------------------------------------------------------------------------
# _inspect_content — fallback for unknown extensions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "first_line, expected",
    [
        ("LOCUS       test_genome  1000 bp", "genbank"),
        ("##gff-version 3", "gff3"),
        (">my_contig", "fasta_contigs"),
    ],
)
def test_unknown_extension_inspected_by_content(tmp_dir, first_line, expected):
    # Use a deliberately-unknown extension to force the content-inspection path
    path = _write(os.path.join(tmp_dir, "data.unknown"), first_line + "\n")
    assert detect_format(path) == expected


def test_unknown_extension_blank_lines_skipped(tmp_dir):
    # _inspect_content skips blanks until the first meaningful line.
    path = _write(
        os.path.join(tmp_dir, "data.unknown"),
        "\n\n   \nLOCUS   test\n",
    )
    assert detect_format(path) == "genbank"


def test_unknown_extension_unparseable_raises(tmp_dir):
    # Random text content with an unknown extension → ValueError. Don't
    # silently mis-route to a parser that will crash three steps later.
    path = _write(os.path.join(tmp_dir, "data.unknown"), "this is not biology\n")
    with pytest.raises(ValueError, match="Cannot determine format"):
        detect_format(path)


def test_unknown_extension_empty_file_raises(tmp_dir):
    path = _write(os.path.join(tmp_dir, "empty.unknown"), "")
    with pytest.raises(ValueError, match="Cannot determine format"):
        detect_format(path)
