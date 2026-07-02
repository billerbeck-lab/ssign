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

import pytest
from _helpers import read_tsv_rows, run_script_main, write_tsv
from run_deeplocpro import (
    OUTPUT_FIELDNAMES,
    _skipped_localization_row,
    _split_fasta_bytes,
    _use_cuda_for_deeplocpro,
    find_output_file,
    parse_deeplocpro_output,
    partition_by_length,
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

    def test_extracts_all_five_localisation_probabilities(self, tmp_dir):
        path = _write_csv(
            os.path.join(tmp_dir, "out.csv"),
            [
                "ACC",
                "Extracellular",
                "Periplasmic",
                "Outer Membrane",
                "Cytoplasmic",
                "Cytoplasmic Membrane",
            ],
            [
                {
                    "ACC": "GENE_001",
                    "Extracellular": "0.05",
                    "Periplasmic": "0.05",
                    "Outer Membrane": "0.07",
                    "Cytoplasmic": "0.03",
                    "Cytoplasmic Membrane": "0.80",
                }
            ],
        )
        entry = parse_deeplocpro_output(path)[0]
        assert entry["extracellular_prob"] == 0.05
        assert entry["periplasmic_prob"] == 0.05
        assert entry["outer_membrane_prob"] == 0.07
        assert entry["cytoplasmic_prob"] == 0.03
        assert entry["cytoplasmic_membrane_prob"] == 0.80
        # 5-class argmax: Cytoplasmic Membrane wins (without it the wrapper
        # used to fall back to Outer Membrane at 0.07, silently misclassifying
        # T1SS/T2SS/T4SS/T6SS inner-membrane components.
        assert entry["predicted_localization"] == "Cytoplasmic Membrane"

    def test_missing_cytoplasmic_membrane_column_defaults_to_zero(self, tmp_dir):
        # DLP outputs older than the 5-class wrapper (or DTU web variants that
        # don't include the column) must still parse, with cytoplasmic_membrane_prob = 0.
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
        assert entry["cytoplasmic_membrane_prob"] == 0
        assert entry["predicted_localization"] == "Extracellular"

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


class TestUseCudaForDeeplocpro:
    """Regression coverage for the GPU auto-detect that controls `-d cuda`
    on the deeplocpro subprocess.

    Pre-fix, the wrapper never passed `-d`, so DLP fell back to its CPU
    default even on GPU nodes. These tests pin the detection contract.
    """

    def test_force_cpu_env_var_overrides_detection(self, monkeypatch):
        monkeypatch.setenv("SSIGN_DEEPLOCPRO_FORCE_CPU", "1")
        assert _use_cuda_for_deeplocpro() is False

    def test_force_cpu_accepts_alternate_truthy_values(self, monkeypatch):
        for v in ("true", "TRUE", "yes", "Yes"):
            monkeypatch.setenv("SSIGN_DEEPLOCPRO_FORCE_CPU", v)
            assert _use_cuda_for_deeplocpro() is False, f"failed for value {v!r}"

    def test_force_cpu_falsy_value_does_not_override(self, monkeypatch):
        monkeypatch.setenv("SSIGN_DEEPLOCPRO_FORCE_CPU", "0")
        try:
            import torch

            expected = bool(torch.cuda.is_available())
        except ImportError:
            expected = False
        assert _use_cuda_for_deeplocpro() is expected

    def test_missing_torch_falls_back_to_cpu(self, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "torch":
                raise ImportError("torch not installed (test stub)")
            return real_import(name, *args, **kwargs)

        monkeypatch.delenv("SSIGN_DEEPLOCPRO_FORCE_CPU", raising=False)
        monkeypatch.setattr(builtins, "__import__", fake_import)
        assert _use_cuda_for_deeplocpro() is False


# ---------------------------------------------------------------------------
# Mega-protein length guard (partition + sentinel + main wiring)
# ---------------------------------------------------------------------------


class TestPartitionByLength:
    """Withhold over-length proteins (e.g. the 16,367-aa kolossin) that OOM DLP."""

    def test_boundary_equal_is_kept(self):
        kept, skipped = partition_by_length({"a": "M" * 10}, 10)
        assert kept == {"a": "M" * 10}
        assert skipped == []

    def test_strictly_longer_is_skipped(self):
        kept, skipped = partition_by_length({"a": "M" * 11}, 10)
        assert kept == {}
        assert skipped == [("a", 11)]

    def test_mixed_input_splits(self):
        seqs = {"short": "M" * 5, "mega": "M" * 100}
        kept, skipped = partition_by_length(seqs, 10)
        assert kept == {"short": "M" * 5}
        assert skipped == [("mega", 100)]

    def test_all_short_passthrough(self):
        seqs = {"a": "MM", "b": "MMM"}
        kept, skipped = partition_by_length(seqs, 10)
        assert kept == seqs
        assert skipped == []

    def test_empty_input(self):
        assert partition_by_length({}, 10) == ({}, [])


class TestSkippedLocalizationRow:
    def test_schema_and_zero_probs(self):
        row = _skipped_localization_row("plu2670", 16367, 5000)
        assert set(row) == set(OUTPUT_FIELDNAMES)
        assert row["locus_tag"] == "plu2670"
        assert row["predicted_localization"] == "Not predicted (too long)"
        for col in (
            "extracellular_prob",
            "periplasmic_prob",
            "outer_membrane_prob",
            "cytoplasmic_prob",
            "cytoplasmic_membrane_prob",
        ):
            assert row[col] == 0.0
        assert "5000" in row["product"]


class TestMainLengthGuard:
    """main() must withhold over-length proteins from the model and re-attach
    them as not-predicted rows. The DeepLocPro binary is never invoked: we
    monkeypatch run_local_deeplocpro to emit a known native CSV for whatever
    it is handed, and assert on what it was (not) handed."""

    def _write_input(self, tmp_dir, short_len, mega_len):
        path = os.path.join(tmp_dir, "in.faa")
        with open(path, "w") as f:
            f.write(">short\n" + "M" * short_len + "\n")
            f.write(">mega\n" + "A" * mega_len + "\n")
        return path

    def _fake_local(self, seen):
        import run_deeplocpro as rdp

        def fake(input_fasta, deeplocpro_path, output_dir, organism="gram-", timeout=None):
            seen["ids"] = list(rdp.read_fasta(input_fasta).keys())
            res = os.path.join(output_dir, "results.csv")
            with open(res, "w") as fh:
                fh.write("ACC,Extracellular,Periplasmic,Outer Membrane,Cytoplasmic,Cytoplasmic Membrane\n")
                for pid in seen["ids"]:
                    fh.write(f"{pid},0.9,0.03,0.03,0.02,0.02\n")
            return res

        return fake

    def test_over_length_excluded_from_model_and_present_as_skipped(self, tmp_dir, monkeypatch):
        import run_deeplocpro as rdp

        monkeypatch.setattr(rdp, "DEEPLOCPRO_MAX_AA", 10)
        inp = self._write_input(tmp_dir, short_len=7, mega_len=50)
        out = os.path.join(tmp_dir, "dlp.tsv")
        seen = {}
        monkeypatch.setattr(rdp, "run_local_deeplocpro", self._fake_local(seen))

        run_script_main(
            monkeypatch,
            rdp.main,
            ["run_deeplocpro.py", "--input", inp, "--sample", "S", "--mode", "local", "--output", out],
        )

        # The model only ever saw the short protein.
        assert seen["ids"] == ["short"]
        rows = {r["locus_tag"]: r for r in read_tsv_rows(out)}
        assert set(rows) == {"short", "mega"}
        assert rows["short"]["predicted_localization"] == "Extracellular"
        assert rows["mega"]["predicted_localization"] == "Not predicted (too long)"
        assert float(rows["mega"]["extracellular_prob"]) == 0.0

    def test_all_skipped_never_invokes_model(self, tmp_dir, monkeypatch):
        import run_deeplocpro as rdp

        monkeypatch.setattr(rdp, "DEEPLOCPRO_MAX_AA", 5)
        inp = self._write_input(tmp_dir, short_len=50, mega_len=50)  # both > 5
        out = os.path.join(tmp_dir, "dlp.tsv")
        called = {"n": 0}

        def boom(*a, **k):
            called["n"] += 1
            raise AssertionError("DeepLocPro must not run when every protein is skipped")

        monkeypatch.setattr(rdp, "run_local_deeplocpro", boom)
        run_script_main(
            monkeypatch,
            rdp.main,
            ["run_deeplocpro.py", "--input", inp, "--sample", "S", "--mode", "local", "--output", out],
        )
        assert called["n"] == 0
        rows = {r["locus_tag"]: r for r in read_tsv_rows(out)}
        assert set(rows) == {"short", "mega"}
        assert all(r["predicted_localization"] == "Not predicted (too long)" for r in rows.values())

    def test_nothing_skipped_all_passed_through(self, tmp_dir, monkeypatch):
        import run_deeplocpro as rdp

        monkeypatch.setattr(rdp, "DEEPLOCPRO_MAX_AA", 5000)
        inp = self._write_input(tmp_dir, short_len=7, mega_len=50)  # both <= 5000
        out = os.path.join(tmp_dir, "dlp.tsv")
        seen = {}
        monkeypatch.setattr(rdp, "run_local_deeplocpro", self._fake_local(seen))
        run_script_main(
            monkeypatch,
            rdp.main,
            ["run_deeplocpro.py", "--input", inp, "--sample", "S", "--mode", "local", "--output", out],
        )
        assert sorted(seen["ids"]) == ["mega", "short"]
        rows = {r["locus_tag"]: r for r in read_tsv_rows(out)}
        assert all(r["predicted_localization"] == "Extracellular" for r in rows.values())
