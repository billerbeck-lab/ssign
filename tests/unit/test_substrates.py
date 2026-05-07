"""Tests for ssign_lib/substrates.py.

Two helpers used by every annotation wrapper:

1. `load_substrate_ids` — reads locus_tags from a TSV with defensive
   tolerance for blank-tag header artefacts (older inlined copies
   used `row["locus_tag"]` and KeyError'd on the same input).
2. `write_substrates_only_fasta` — projects an all-proteins FASTA
   onto just the substrate IDs; returns the count for short-circuit
   detection by the caller.
"""

import os
import sys


from _helpers import write_tsv  # noqa: E402
from ssign_lib.fasta_io import write_fasta  # noqa: E402
from ssign_lib.substrates import (  # noqa: E402
    load_substrate_ids,
    write_substrates_only_fasta,
)

# ---------------------------------------------------------------------------
# load_substrate_ids
# ---------------------------------------------------------------------------


class TestLoadSubstrateIds:
    def test_basic_read(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["locus_tag", "extra"],
            [
                {"locus_tag": "GENE_001", "extra": "a"},
                {"locus_tag": "GENE_002", "extra": "b"},
            ],
        )
        assert load_substrate_ids(path) == {"GENE_001", "GENE_002"}

    def test_returns_set_not_list(self, tmp_dir):
        # Set semantics matter for caller's `pid in substrate_ids` checks.
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["locus_tag"],
            [{"locus_tag": "GENE_001"}, {"locus_tag": "GENE_001"}],
        )
        result = load_substrate_ids(path)
        assert isinstance(result, set)
        assert result == {"GENE_001"}

    def test_blank_locus_tag_skipped(self, tmp_dir):
        # Empty-string locus_tag in upstream output → drop, don't pollute the set.
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["locus_tag"],
            [
                {"locus_tag": ""},
                {"locus_tag": "GENE_001"},
            ],
        )
        assert load_substrate_ids(path) == {"GENE_001"}

    def test_whitespace_only_locus_tag_skipped(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["locus_tag"],
            [
                {"locus_tag": "   "},
                {"locus_tag": "\t\t"},
                {"locus_tag": "GENE_001"},
            ],
        )
        assert load_substrate_ids(path) == {"GENE_001"}

    def test_locus_tag_stripped(self, tmp_dir):
        # Surrounding whitespace removed — defensive against TSV pasting artefacts.
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["locus_tag"],
            [{"locus_tag": "  GENE_001  "}],
        )
        assert load_substrate_ids(path) == {"GENE_001"}

    def test_missing_locus_tag_column_yields_empty(self, tmp_dir):
        # `row.get("locus_tag")` returns None → `or ""` → skipped.
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["other_col"],
            [{"other_col": "GENE_001"}],
        )
        assert load_substrate_ids(path) == set()

    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.tsv")
        open(path, "w").close()
        assert load_substrate_ids(path) == set()

    def test_header_only_returns_empty(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "substrates.tsv"),
            ["locus_tag"],
            [],
        )
        assert load_substrate_ids(path) == set()


# ---------------------------------------------------------------------------
# write_substrates_only_fasta
# ---------------------------------------------------------------------------


class TestWriteSubstratesOnlyFasta:
    def test_filters_to_substrate_ids(self, tmp_dir):
        proteins = os.path.join(tmp_dir, "proteins.faa")
        write_fasta(
            {"GENE_001": "MKTLLL", "GENE_002": "GGGGG", "GENE_003": "AAAAA"},
            proteins,
        )
        out = os.path.join(tmp_dir, "subs.faa")
        n = write_substrates_only_fasta(proteins, {"GENE_001", "GENE_003"}, out)
        assert n == 2
        # Output contains only the requested substrates
        with open(out) as f:
            content = f.read()
        assert ">GENE_001" in content
        assert ">GENE_003" in content
        assert ">GENE_002" not in content

    def test_count_returned_zero_when_no_match(self, tmp_dir):
        # Caller short-circuits on n == 0 — the contract that runs the
        # downstream tool only when there's something to annotate.
        proteins = os.path.join(tmp_dir, "proteins.faa")
        write_fasta({"GENE_001": "MKT"}, proteins)
        out = os.path.join(tmp_dir, "subs.faa")
        n = write_substrates_only_fasta(proteins, {"NOT_PRESENT"}, out)
        assert n == 0

    def test_substrate_id_with_no_sequence_skipped(self, tmp_dir):
        # Substrate ID listed but no FASTA entry — silently skip, don't crash
        proteins = os.path.join(tmp_dir, "proteins.faa")
        write_fasta({"GENE_001": "MKT"}, proteins)
        out = os.path.join(tmp_dir, "subs.faa")
        n = write_substrates_only_fasta(
            proteins,
            {"GENE_001", "GENE_NONEXISTENT"},
            out,
        )
        assert n == 1

    def test_empty_substrate_set_writes_empty_fasta(self, tmp_dir):
        proteins = os.path.join(tmp_dir, "proteins.faa")
        write_fasta({"GENE_001": "MKT"}, proteins)
        out = os.path.join(tmp_dir, "subs.faa")
        n = write_substrates_only_fasta(proteins, set(), out)
        assert n == 0
        # File still exists, just empty
        assert os.path.exists(out)
        assert os.path.getsize(out) == 0
