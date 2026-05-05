"""Tests for prodigal_to_gene_info.py.

Bridge between Prodigal's FASTA header convention and ssign's canonical
gene_info schema. Three things must be right or downstream proximity
arithmetic breaks:

1. Header parsing: `>contig_1_5 # start # end # strand # attrs`.
2. Coordinate conversion: Prodigal is 1-based, gene_info is 0-based-start.
3. locus_tag format: `{sample}_{gene_num.zfill(5)}` — the zero-pad is
   load-bearing because downstream tools sort locus_tags lexicographically.
"""

import os
import sys
import tempfile

from hypothesis import given, settings
from hypothesis import strategies as st

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from _helpers import read_tsv_rows, run_script_main  # noqa: E402
from prodigal_to_gene_info import main as prodigal_main  # noqa: E402


def _write_prodigal_fasta(path, entries):
    """entries: list of (raw_id, start, end, strand_code, sequence)."""
    with open(path, "w") as f:
        for raw_id, start, end, strand_code, seq in entries:
            f.write(f">{raw_id} # {start} # {end} # {strand_code} # ID=x\n")
            f.write(seq + "\n")
    return path


def _run(monkeypatch, tmp_dir, prodigal_entries, sample="genome"):
    proteins = _write_prodigal_fasta(
        os.path.join(tmp_dir, "prodigal.faa"),
        prodigal_entries,
    )
    # --gff is a required argparse arg but unused in the body; stub it
    gff_stub = os.path.join(tmp_dir, "prodigal.gff")
    open(gff_stub, "w").close()
    out_proteins = os.path.join(tmp_dir, "proteins.faa")
    out_gene_info = os.path.join(tmp_dir, "gene_info.tsv")
    run_script_main(
        monkeypatch,
        prodigal_main,
        [
            "prodigal_to_gene_info",
            "--proteins",
            proteins,
            "--gff",
            gff_stub,
            "--sample",
            sample,
            "--out-proteins",
            out_proteins,
            "--out-gene-info",
            out_gene_info,
        ],
    )
    return read_tsv_rows(out_gene_info), out_proteins


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------


def test_basic_header_parsed(monkeypatch, tmp_dir):
    rows, _ = _run(
        monkeypatch,
        tmp_dir,
        [
            ("contig_1_1", 100, 500, 1, "MKT"),
        ],
    )
    assert len(rows) == 1
    assert rows[0]["contig"] == "contig_1"
    assert rows[0]["start"] == "99"  # 1-based → 0-based: 100 - 1
    assert rows[0]["end"] == "500"
    assert rows[0]["strand"] == "+"


def test_strand_plus_one_yields_plus(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, [("contig_1_1", 1, 100, 1, "MKT")])
    assert rows[0]["strand"] == "+"


def test_strand_minus_one_yields_minus(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, [("contig_1_1", 1, 100, -1, "MKT")])
    assert rows[0]["strand"] == "-"


def test_start_converted_to_zero_based(monkeypatch, tmp_dir):
    # Prodigal 1-based start=1 → gene_info 0-based start=0
    rows, _ = _run(monkeypatch, tmp_dir, [("c_1", 1, 100, 1, "MKT")])
    assert rows[0]["start"] == "0"


# ---------------------------------------------------------------------------
# Contig extraction + locus_tag format
# ---------------------------------------------------------------------------


def test_contig_extracted_from_underscore_suffix(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, [("scaffold_42_7", 1, 100, 1, "MKT")])
    assert rows[0]["contig"] == "scaffold_42"
    # Gene num is 7, zero-padded to 5 digits, prefixed with sample
    assert rows[0]["locus_tag"] == "genome_00007"


def test_locus_tag_uses_sample_prefix(monkeypatch, tmp_dir):
    rows, _ = _run(
        monkeypatch,
        tmp_dir,
        [("c_3", 1, 100, 1, "MKT")],
        sample="EColi_K12",
    )
    assert rows[0]["locus_tag"] == "EColi_K12_00003"


def test_locus_tag_zero_padded_to_five_digits(monkeypatch, tmp_dir):
    # Lex-sort matters: "genome_00009" < "genome_00010" without the pad would invert
    rows, _ = _run(
        monkeypatch,
        tmp_dir,
        [
            ("c_9", 1, 100, 1, "MKT"),
            ("c_10", 1, 100, 1, "MKT"),
        ],
    )
    tags = [r["locus_tag"] for r in rows]
    assert tags == sorted(tags)


def test_raw_id_without_underscore_falls_back_to_contig(monkeypatch, tmp_dir):
    # Prodigal usually emits `_N`-suffixed IDs, but if a custom contig name
    # has no trailing _N the fallback assigns gene_num "0"
    rows, _ = _run(monkeypatch, tmp_dir, [("singletoncontig", 1, 100, 1, "MKT")])
    assert rows[0]["contig"] == "singletoncontig"
    assert rows[0]["locus_tag"] == "genome_00000"


# ---------------------------------------------------------------------------
# Sequence handling: trailing stop, multi-line, multi-entry
# ---------------------------------------------------------------------------


def test_trailing_stop_codon_stripped_in_output(monkeypatch, tmp_dir):
    _, out_proteins = _run(
        monkeypatch,
        tmp_dir,
        [
            ("c_1", 1, 9, 1, "MKT*"),
        ],
    )
    with open(out_proteins) as f:
        text = f.read()
    assert "MKT" in text
    assert "*" not in text


def test_multi_contig_multi_gene(monkeypatch, tmp_dir):
    # Prodigal numbers genes globally (not per-contig), so the trailing _N
    # advances 1, 2, 3 across the whole input regardless of which contig
    # the gene falls on.
    rows, _ = _run(
        monkeypatch,
        tmp_dir,
        [
            ("contig_A_1", 1, 100, 1, "MKT"),
            ("contig_A_2", 200, 300, -1, "MFV"),
            ("contig_B_3", 50, 150, 1, "MGA"),
        ],
    )
    by_locus = {r["locus_tag"]: r for r in rows}
    assert by_locus["genome_00001"]["contig"] == "contig_A"
    assert by_locus["genome_00001"]["strand"] == "+"
    assert by_locus["genome_00002"]["contig"] == "contig_A"
    assert by_locus["genome_00002"]["strand"] == "-"
    assert by_locus["genome_00003"]["contig"] == "contig_B"
    assert by_locus["genome_00003"]["strand"] == "+"


def test_multiline_protein_sequence_concatenated(monkeypatch, tmp_dir):
    proteins = os.path.join(tmp_dir, "prodigal.faa")
    with open(proteins, "w") as f:
        f.write(">contig_1_1 # 1 # 30 # 1 # ID=x\n")
        f.write("MKTLLLT\n")  # split across lines
        f.write("LLCAFSV\n")
    gff_stub = os.path.join(tmp_dir, "prodigal.gff")
    open(gff_stub, "w").close()
    out_proteins = os.path.join(tmp_dir, "proteins.faa")
    out_gene_info = os.path.join(tmp_dir, "gene_info.tsv")
    saved = sys.argv
    try:
        sys.argv = [
            "prodigal_to_gene_info",
            "--proteins",
            proteins,
            "--gff",
            gff_stub,
            "--sample",
            "genome",
            "--out-proteins",
            out_proteins,
            "--out-gene-info",
            out_gene_info,
        ]
        prodigal_main()
    finally:
        sys.argv = saved
    with open(out_proteins) as f:
        text = f.read()
    assert "MKTLLLTLLCAFSV" in text


# ---------------------------------------------------------------------------
# Property: coord conversion is consistent and reversible
# ---------------------------------------------------------------------------


@settings(max_examples=25, deadline=None)
@given(
    one_based_start=st.integers(min_value=1, max_value=10_000_000),
    length=st.integers(min_value=3, max_value=10_000),
    strand_code=st.sampled_from([1, -1]),
)
def test_coord_conversion_consistent(one_based_start, length, strand_code):
    """Prodigal start (1-based) → gene_info start (0-based) = start - 1.
    end is unchanged. Strand encoding maps {1, -1} → {+, -}."""
    end = one_based_start + length - 1
    with tempfile.TemporaryDirectory() as td:
        proteins = os.path.join(td, "prodigal.faa")
        with open(proteins, "w") as f:
            f.write(f">contig_1_1 # {one_based_start} # {end} # {strand_code} # ID=x\n")
            f.write("MKT\n")
        gff_stub = os.path.join(td, "prodigal.gff")
        open(gff_stub, "w").close()
        out_proteins = os.path.join(td, "proteins.faa")
        out_gene_info = os.path.join(td, "gene_info.tsv")
        saved = sys.argv
        try:
            sys.argv = [
                "prodigal_to_gene_info",
                "--proteins",
                proteins,
                "--gff",
                gff_stub,
                "--sample",
                "genome",
                "--out-proteins",
                out_proteins,
                "--out-gene-info",
                out_gene_info,
            ]
            prodigal_main()
        finally:
            sys.argv = saved
        rows = read_tsv_rows(out_gene_info)

    assert int(rows[0]["start"]) == one_based_start - 1
    assert int(rows[0]["end"]) == end
    assert rows[0]["strand"] == ("+" if strand_code == 1 else "-")
