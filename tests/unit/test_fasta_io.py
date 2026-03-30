"""Tests for ssign_lib.fasta_io — FASTA read/write."""

import os

import pytest
from ssign_lib.fasta_io import read_fasta, write_fasta


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
        with open(path, 'w') as f:
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
