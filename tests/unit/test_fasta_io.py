"""Tests for ssign_lib.fasta_io — FASTA read/write."""

import os

import pytest
from ssign_lib.fasta_io import (  # noqa: F401  # all four used in dedicated test classes below
    count_sequences,
    read_fasta,
    read_fasta_records,
    write_fasta,
)


class TestReadFasta:
    def test_basic_read(self, sample_fasta):
        seqs = read_fasta(sample_fasta)
        assert len(seqs) == 2
        assert "protein_A" in seqs
        assert "protein_B" in seqs

    def test_multiline_sequence_joined(self, sample_fasta):
        seqs = read_fasta(sample_fasta)
        assert seqs["protein_A"] == "MKTLLLTLLCAFSVAQA" + "VDLPTQEPALGK"

    def test_single_line_sequence(self, sample_fasta):
        seqs = read_fasta(sample_fasta)
        assert seqs["protein_B"] == "MFVFLVLLPLVSSQ"

    def test_header_first_token_only(self, sample_fasta):
        """Only the first token after > is used as the ID."""
        seqs = read_fasta(sample_fasta)
        assert "protein_A" in seqs
        assert "description" not in list(seqs.keys())[0]

    def test_empty_file(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.fasta")
        with open(path, "w") as f:
            f.write("")
        seqs = read_fasta(path)
        assert seqs == {}

    def test_file_not_found(self, tmp_dir):
        with pytest.raises(FileNotFoundError):
            read_fasta(os.path.join(tmp_dir, "nonexistent.fasta"))


class TestWriteFasta:
    def test_round_trip(self, tmp_dir):
        seqs = {"seq1": "ACGT", "seq2": "MFVL"}
        path = os.path.join(tmp_dir, "out.fasta")
        n = write_fasta(seqs, path)
        assert n == 2

        recovered = read_fasta(path)
        assert recovered == seqs

    def test_empty_sequence_skipped(self, tmp_dir):
        seqs = {"good": "ACGT", "empty": ""}
        path = os.path.join(tmp_dir, "out.fasta")
        n = write_fasta(seqs, path)
        assert n == 1

        recovered = read_fasta(path)
        assert "empty" not in recovered

    def test_line_wrapping(self, tmp_dir):
        seq = "A" * 200
        path = os.path.join(tmp_dir, "out.fasta")
        write_fasta({"long": seq}, path, line_width=80)

        with open(path) as f:
            lines = f.readlines()
        # Header + ceil(200/80)=3 sequence lines
        assert len(lines) == 4
        assert len(lines[1].strip()) == 80
        assert len(lines[2].strip()) == 80
        assert len(lines[3].strip()) == 40

    def test_creates_parent_dirs(self, tmp_dir):
        path = os.path.join(tmp_dir, "sub", "dir", "out.fasta")
        write_fasta({"s": "ACGT"}, path)
        assert os.path.exists(path)


class TestCountSequences:
    """count_sequences accepts a path or in-memory bytes and returns the
    number of `>`-prefixed records. Behaviour must match the legacy
    hand-rolled `sum(1 for line in ... if line.startswith(">"))` so the
    six call sites swapped over in B.4 are byte-equivalent."""

    def test_path_with_two_records(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        write_fasta({"A": "ACGT", "B": "TGCA"}, path)
        assert count_sequences(path) == 2

    def test_path_with_zero_records(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.fasta")
        open(path, "w").close()
        assert count_sequences(path) == 0

    def test_bytes_with_two_records(self):
        content = b">A\nMKT\n>B\nGGG\n"
        assert count_sequences(content) == 2

    def test_bytes_starting_with_no_leading_newline(self):
        # Old: would split-and-count >. New: must include the leading > too.
        assert count_sequences(b">only\nMKT\n") == 1

    def test_bytes_with_only_sequence_no_header(self):
        # No `>` anywhere → 0
        assert count_sequences(b"MKTLLLT\n") == 0

    def test_bytes_empty(self):
        assert count_sequences(b"") == 0

    def test_path_and_bytes_agree(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        write_fasta({"A": "ACGT", "B": "TGCA", "C": "AAA"}, path)
        with open(path, "rb") as f:
            content = f.read()
        assert count_sequences(path) == count_sequences(content) == 3

    def test_does_not_count_internal_gt_in_sequence(self):
        # `>` only counts when it starts a line. A `>` inside an annotation
        # column or as part of a description is fine — but our test is
        # paranoid: hand-rolled sum had this property too.
        # (Legitimate biological FASTA never has bare > in sequence lines.)
        # Passing this just to confirm the bytes path doesn't false-positive
        # on `XYZ>ABC` mid-line.
        content = b">A\nMKT>frag\n>B\nGGG\n"
        assert count_sequences(content) == 2


class TestReadFastaRecords:
    """read_fasta_records yields (full_header, sequence) tuples — used by
    prodigal_to_gene_info to parse Prodigal-style headers
    `>contig_N_M # start # end # strand # attrs` where read_fasta's
    first-token-only behaviour would discard the metadata."""

    def test_preserves_full_header_after_gt(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        with open(path, "w") as f:
            f.write(">contig_1_5 # 100 # 250 # 1 # ID=1_5;partial=00\nMKT\n")
        records = list(read_fasta_records(path))
        assert records == [("contig_1_5 # 100 # 250 # 1 # ID=1_5;partial=00", "MKT")]

    def test_multiple_records(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        with open(path, "w") as f:
            f.write(">A\nMKT\n>B desc\nGGG\n")
        records = list(read_fasta_records(path))
        assert records == [("A", "MKT"), ("B desc", "GGG")]

    def test_multi_line_sequence_concatenated(self, tmp_dir):
        path = os.path.join(tmp_dir, "x.fasta")
        with open(path, "w") as f:
            f.write(">A\nMKTLLL\nTLLCAFSV\n")
        records = list(read_fasta_records(path))
        assert records == [("A", "MKTLLLTLLCAFSV")]

    def test_empty_file_yields_nothing(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.fasta")
        open(path, "w").close()
        assert list(read_fasta_records(path)) == []

    def test_missing_file_raises(self, tmp_dir):
        with pytest.raises(FileNotFoundError):
            list(read_fasta_records(os.path.join(tmp_dir, "nope.fasta")))

    def test_is_lazy_generator(self, tmp_dir):
        # Memory bounded by longest single record — confirm by checking that
        # it returns an iterator, not a materialised list.
        path = os.path.join(tmp_dir, "x.fasta")
        with open(path, "w") as f:
            f.write(">A\nMKT\n")
        result = read_fasta_records(path)
        assert hasattr(result, "__next__")
