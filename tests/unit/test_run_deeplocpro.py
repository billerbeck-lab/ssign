"""Tests for run_deeplocpro.py.

Three testable surfaces here:

1. `_split_fasta_bytes` — split FASTA into batches of ≤500 sequences
   for the DTU 5000-protein-per-batch limit.
2. `parse_deeplocpro_output` — tolerant parser handling DTU web (ACC),
   local install (protein_id / Protein_ID / ID / Name) column variants
   and comma- or tab-delimited input.
3. `find_output_file` — locates DeepLocPro CSV/TSV in an output dir
   tree (DTU emits filenames that vary across versions).

The deep subprocess + DTU HTTP paths require live network or a local
DTU install and are exercised by tests/integration/test_run_deeplocpro_integration.py.
"""

import os
import sys

import pytest


from _helpers import write_tsv  # noqa: E402
from run_deeplocpro import (  # noqa: E402
    _split_fasta_bytes,
    find_output_file,
    parse_deeplocpro_output,
)

# ---------------------------------------------------------------------------
# _split_fasta_bytes
# ---------------------------------------------------------------------------


class TestSplitFastaBytes:
    def _make_fasta(self, n_seqs):
        return b"".join(f">P{i}\nMKT{i}\n".encode() for i in range(n_seqs))

    def test_under_batch_size_yields_one_batch(self):
        batches = _split_fasta_bytes(self._make_fasta(3), batch_size=500)
        assert len(batches) == 1
        # All three sequences in the single batch
        assert batches[0].count(b">") == 3

    def test_exactly_batch_size_yields_one_batch(self):
        batches = _split_fasta_bytes(self._make_fasta(500), batch_size=500)
        assert len(batches) == 1
        assert batches[0].count(b">") == 500

    def test_over_batch_size_splits_multiple_batches(self):
        batches = _split_fasta_bytes(self._make_fasta(750), batch_size=500)
        assert len(batches) == 2
        # First batch has exactly batch_size; second has the remainder
        assert batches[0].count(b">") == 500
        assert batches[1].count(b">") == 250

    def test_three_batches(self):
        batches = _split_fasta_bytes(self._make_fasta(1100), batch_size=500)
        assert len(batches) == 3
        counts = [b.count(b">") for b in batches]
        assert counts == [500, 500, 100]

    def test_empty_fasta_yields_empty_or_singleton(self):
        # Empty input: either [] or one empty-content batch — pin behaviour
        batches = _split_fasta_bytes(b"", batch_size=500)
        # Must not crash; total seq count is 0
        assert sum(b.count(b">") for b in batches) == 0


# ---------------------------------------------------------------------------
# parse_deeplocpro_output — tolerant column-alias parser
# ---------------------------------------------------------------------------


def _write_csv(path, header_cols, rows, sep=","):
    """Wrapper over `write_tsv` that defaults to comma — DLP's native delimiter."""
    return write_tsv(path, header_cols, rows, delimiter=sep)


class TestParseDeeplocproOutput:
    """The DTU web server, the local CLI, and prior tool versions all emit
    slightly different column names. The parser must accept all known
    variants."""

    @pytest.mark.parametrize("id_col", ["ACC", "protein_id", "Protein_ID", "ID", "Name"])
    def test_accepts_each_id_column_alias(self, tmp_dir, id_col):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            [id_col, "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic"],
            [
                {
                    id_col: "GENE_001",
                    "Extracellular": "0.9",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.03",
                    "Cytoplasmic": "0.02",
                }
            ],
        )
        entries = parse_deeplocpro_output(path)
        assert len(entries) == 1
        assert entries[0]["locus_tag"] == "GENE_001"

    def test_extracts_all_four_localisation_probabilities(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            ["ACC", "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic"],
            [
                {
                    "ACC": "GENE_001",
                    "Extracellular": "0.85",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.07",
                    "Cytoplasmic": "0.03",
                }
            ],
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["extracellular_prob"] == 0.85
        assert entry["periplasmic_prob"] == 0.05
        assert entry["outer_membrane_prob"] == 0.07
        assert entry["cytoplasmic_prob"] == 0.03

    def test_predicted_localisation_is_argmax(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            ["ACC", "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic"],
            [
                {
                    "ACC": "GENE_001",
                    "Extracellular": "0.1",
                    "Periplasmic": "0.2",
                    "Outer Membrane": "0.6",
                    "Cytoplasmic": "0.1",
                }
            ],
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["predicted_localization"] == "Outer Membrane"

    def test_lowercase_column_aliases(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            ["protein_id", "extracellular", "periplasmic", "outer_membrane", "cytoplasmic"],
            [
                {
                    "protein_id": "GENE_001",
                    "extracellular": "0.9",
                    "periplasmic": "0.05",
                    "outer_membrane": "0.03",
                    "cytoplasmic": "0.02",
                }
            ],
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["locus_tag"] == "GENE_001"
        assert entry["extracellular_prob"] == 0.9

    def test_tab_separated_input(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.tsv"),
            ["ACC", "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic"],
            [
                {
                    "ACC": "GENE_001",
                    "Extracellular": "0.9",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.03",
                    "Cytoplasmic": "0.02",
                }
            ],
            sep="\t",
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["locus_tag"] == "GENE_001"

    def test_rows_without_id_skipped(self, tmp_dir):
        # A row with no recognised ID column → silently dropped
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            ["unknown_id_col", "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic"],
            [
                {
                    "unknown_id_col": "GENE_001",
                    "Extracellular": "0.9",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.03",
                    "Cytoplasmic": "0.02",
                }
            ],
        )
        assert parse_deeplocpro_output(path) == []

    def test_product_field_picked_up(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            ["ACC", "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic", "annotation"],
            [
                {
                    "ACC": "GENE_001",
                    "Extracellular": "0.9",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.03",
                    "Cytoplasmic": "0.02",
                    "annotation": "type IV pilin",
                }
            ],
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["product"] == "type IV pilin"

    def test_probabilities_rounded_to_4_decimals(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            ["ACC", "Extracellular", "Periplasmic", "Outer Membrane", "Cytoplasmic"],
            [
                {
                    "ACC": "GENE_001",
                    "Extracellular": "0.123456789",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.03",
                    "Cytoplasmic": "0.02",
                }
            ],
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["extracellular_prob"] == 0.1235  # rounded to 4 dp


# ---------------------------------------------------------------------------
# find_output_file
# ---------------------------------------------------------------------------


class TestFindOutputFile:
    def test_finds_csv_in_top_level(self, tmp_dir):
        path = os.path.join(tmp_dir, "results.csv")
        open(path, "w").close()
        assert find_output_file(tmp_dir) == path

    def test_finds_tsv_in_top_level(self, tmp_dir):
        path = os.path.join(tmp_dir, "results.tsv")
        open(path, "w").close()
        assert find_output_file(tmp_dir) == path

    def test_finds_in_subdirectory(self, tmp_dir):
        subdir = os.path.join(tmp_dir, "deeplocpro_run_1")
        os.makedirs(subdir)
        path = os.path.join(subdir, "results.csv")
        open(path, "w").close()
        assert find_output_file(tmp_dir) == path

    def test_raises_when_no_output(self, tmp_dir):
        # Empty dir → FileNotFoundError with a useful message
        with pytest.raises(FileNotFoundError, match="No output file"):
            find_output_file(tmp_dir)

    def test_ignores_non_csv_files(self, tmp_dir):
        open(os.path.join(tmp_dir, "log.txt"), "w").close()
        open(os.path.join(tmp_dir, "results.csv"), "w").close()
        assert find_output_file(tmp_dir).endswith(".csv")
