"""Tests for run_signalp.py.

Two pure-Python surfaces here:

1. `parse_signalp_output` — SignalP 6.0 prediction_results.txt parser.
   Critical contract: signalp_probability is max across the five SP-type
   probabilities (SP/LIPO/TAT/TATLIPO/PILIN), explicitly excluding the
   OTHER column so non-secreted proteins report ~0 not ~1.
2. `find_output_file` — file-discovery cascade:
   prediction_results.txt → *.txt with prediction/summary in name →
   *.signalp5 → FileNotFoundError.
3. `_CS_POS_RE` — extracts the CS position range from the cell text
   "CS pos: 20-21. Pr: 0.9756".

The DTU web submission + local subprocess paths require live network or
a local DTU install and are exercised by
tests/integration/test_run_signalp_integration.py.
"""

import logging
import os
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from run_signalp import (  # noqa: E402
    _CS_POS_RE,
    _SIGNALP6_COLS,
    find_output_file,
    parse_signalp_output,
)

# ---------------------------------------------------------------------------
# parse_signalp_output
# ---------------------------------------------------------------------------


def _write_signalp_v6(path, rows, *, with_header=True, header_drift=False):
    """Build a v6 prediction_results.txt with optional header drift."""
    with open(path, "w") as f:
        f.write("# SignalP-6.0\tOrganism: Other\tTimestamp: 2026-05-05\n")
        if with_header:
            cols = list(_SIGNALP6_COLS)
            if header_drift:
                cols[1] = "DriftedColumn"  # cause schema-drift warning
            f.write("# " + "\t".join(cols) + "\n")
        for r in rows:
            cells = [r.get(c, "") for c in _SIGNALP6_COLS]
            f.write("\t".join(str(c) for c in cells) + "\n")
    return path


def _row(
    locus_tag="GENE_001",
    prediction="SP",
    other=0.05,
    sp=0.92,
    lipo=0.01,
    tat=0.01,
    tatlipo=0.0,
    pilin=0.01,
    cs_position="CS pos: 20-21. Pr: 0.9756",
):
    return {
        "ID": locus_tag,
        "Prediction": prediction,
        "OTHER": other,
        "SP(Sec/SPI)": sp,
        "LIPO(Sec/SPII)": lipo,
        "TAT(Tat/SPI)": tat,
        "TATLIPO(Tat/SPII)": tatlipo,
        "PILIN(Sec/SPIII)": pilin,
        "CS Position": cs_position,
    }


class TestSignalpProbability:
    """signalp_probability = max(SP, LIPO, TAT, TATLIPO, PILIN). The OTHER
    column is intentionally excluded so non-secreted proteins read low,
    not high — load-bearing for downstream cross_validate."""

    @pytest.mark.parametrize(
        "sp, lipo, tat, tatlipo, pilin, expected",
        [
            # Each of the five SP-types can be the max
            (0.92, 0.05, 0.05, 0.05, 0.05, 0.92),  # SP (Sec/SPI)
            (0.05, 0.85, 0.05, 0.05, 0.05, 0.85),  # LIPO (Sec/SPII)
            (0.05, 0.05, 0.85, 0.05, 0.05, 0.85),  # TAT (Tat/SPI)
            (0.05, 0.05, 0.05, 0.85, 0.05, 0.85),  # TATLIPO (Tat/SPII)
            (0.05, 0.05, 0.05, 0.05, 0.85, 0.85),  # PILIN (Sec/SPIII)
        ],
        ids=["sp", "lipo", "tat", "tatlipo", "pilin"],
    )
    def test_max_taken_over_each_sp_type(self, tmp_dir, sp, lipo, tat, tatlipo, pilin, expected):
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(sp=sp, lipo=lipo, tat=tat, tatlipo=tatlipo, pilin=pilin)],
        )
        assert parse_signalp_output(path)[0]["signalp_probability"] == expected

    def test_other_column_excluded_from_max(self, tmp_dir):
        """Critical contract: OTHER (non-secreted probability) is NOT counted.
        Without this exclusion, every non-secreted protein would read as
        ~highly secreted and pollute downstream cross_validate."""
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(other=0.99, sp=0.05, lipo=0.05, tat=0.05, tatlipo=0.05, pilin=0.05)],
        )
        assert parse_signalp_output(path)[0]["signalp_probability"] == 0.05

    def test_non_secreted_protein_reads_low(self, tmp_dir):
        # End-to-end real-world OTHER-prediction case
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(prediction="OTHER", other=0.99, sp=0.005, lipo=0.001, tat=0.001, tatlipo=0.001, pilin=0.002)],
        )
        assert parse_signalp_output(path)[0]["signalp_probability"] < 0.05


class TestCsPositionExtraction:
    @pytest.mark.parametrize(
        "cell, expected",
        [
            ("CS pos: 20-21. Pr: 0.9756", "20-21"),
            ("CS pos:1-2", "1-2"),  # no space after colon
            ("CS pos: 100-101. Pr: 0.5", "100-101"),
            ("", ""),  # empty cell
            ("no signal peptide", ""),  # no CS pos prefix
            ("CS pos: -1-5", ""),  # malformed (leading minus)
        ],
    )
    def test_cs_position_regex(self, cell, expected):
        m = _CS_POS_RE.search(cell)
        assert (m.group(1) if m else "") == expected

    def test_cs_position_extracted_in_parsed_entry(self, tmp_dir):
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(cs_position="CS pos: 25-26. Pr: 0.99")],
        )
        entry = parse_signalp_output(path)[0]
        assert entry["signalp_cs_position"] == "25-26"

    def test_no_cs_position_yields_empty(self, tmp_dir):
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(cs_position="")],
        )
        entry = parse_signalp_output(path)[0]
        assert entry["signalp_cs_position"] == ""


class TestLocusTagExtraction:
    def test_locus_tag_is_first_whitespace_token(self, tmp_dir):
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(locus_tag="GENE_001 some annotation")],
        )
        entry = parse_signalp_output(path)[0]
        assert entry["locus_tag"] == "GENE_001"

    def test_blank_id_row_skipped(self, tmp_dir):
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row(locus_tag=""), _row(locus_tag="GENE_002")],
        )
        loci = [e["locus_tag"] for e in parse_signalp_output(path)]
        assert loci == ["GENE_002"]


class TestParserResilience:
    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "out.txt")
        open(path, "w").close()
        assert parse_signalp_output(path) == []

    def test_only_comment_lines_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "out.txt")
        with open(path, "w") as f:
            f.write("# SignalP-6.0\n")
            f.write("# " + "\t".join(_SIGNALP6_COLS) + "\n")
        assert parse_signalp_output(path) == []

    def test_short_rows_skipped(self, tmp_dir):
        # A row with fewer than _COL_PILIN+1 columns is corrupt — drop it
        path = os.path.join(tmp_dir, "out.txt")
        with open(path, "w") as f:
            f.write("# SignalP-6.0\n")
            f.write("# " + "\t".join(_SIGNALP6_COLS) + "\n")
            f.write("GENE_001\tSP\t0.05\n")  # only 3 columns
        assert parse_signalp_output(path) == []

    def test_non_numeric_probability_falls_back_to_zero(self, tmp_dir, caplog):
        """A non-float in an SP-type column logs a warning and uses 0.0
        instead of crashing — defensive against format drift."""
        path = os.path.join(tmp_dir, "out.txt")
        with open(path, "w") as f:
            f.write("# SignalP-6.0\n")
            f.write("# " + "\t".join(_SIGNALP6_COLS) + "\n")
            # SP column is "n/a" — corrupt
            f.write("GENE_001\tSP\t0.05\tn/a\t0.01\t0.01\t0.0\t0.01\tCS pos: 20-21\n")
        with caplog.at_level(logging.WARNING):
            entries = parse_signalp_output(path)
        assert len(entries) == 1
        # Max is taken over the four valid columns (0.01, 0.01, 0.0, 0.01)
        # plus the substituted 0.0 for the corrupt SP column
        assert entries[0]["signalp_probability"] == 0.01
        assert any("Non-float" in rec.message for rec in caplog.records)

    def test_header_drift_logs_warning(self, tmp_dir, caplog):
        path = _write_signalp_v6(
            os.path.join(tmp_dir, "out.txt"),
            [_row()],
            header_drift=True,
        )
        with caplog.at_level(logging.WARNING):
            parse_signalp_output(path)
        assert any("SignalP header drift" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# find_output_file — file-discovery cascade
# ---------------------------------------------------------------------------


class TestFindOutputFile:
    def test_prediction_results_txt_prioritized(self, tmp_dir):
        # Both a prediction_results.txt and a fallback file present; the
        # canonical name wins.
        target = os.path.join(tmp_dir, "prediction_results.txt")
        open(target, "w").close()
        with open(os.path.join(tmp_dir, "summary.txt"), "w") as f:
            f.write("x" * 100)
        assert find_output_file(tmp_dir) == target

    def test_falls_back_to_summary_or_prediction_in_name(self, tmp_dir):
        # No prediction_results.txt, but a foo_summary.txt with content > 10 B
        path = os.path.join(tmp_dir, "foo_summary.txt")
        with open(path, "w") as f:
            f.write("x" * 50)
        assert find_output_file(tmp_dir) == path

    def test_skips_tiny_summary_files(self, tmp_dir):
        # Empty summary file (≤10 bytes) is skipped per the size guard
        with open(os.path.join(tmp_dir, "foo_summary.txt"), "w") as f:
            f.write("tiny")
        # signalp5 should be the fallback
        with open(os.path.join(tmp_dir, "x.signalp5"), "w") as f:
            f.write("x")
        assert find_output_file(tmp_dir).endswith(".signalp5")

    def test_falls_back_to_signalp5(self, tmp_dir):
        path = os.path.join(tmp_dir, "results.signalp5")
        open(path, "w").close()
        assert find_output_file(tmp_dir) == path

    def test_raises_when_no_output(self, tmp_dir):
        with pytest.raises(FileNotFoundError, match="No SignalP output"):
            find_output_file(tmp_dir)

    def test_finds_file_in_subdirectory(self, tmp_dir):
        subdir = os.path.join(tmp_dir, "signalp_run_42")
        os.makedirs(subdir)
        target = os.path.join(subdir, "prediction_results.txt")
        open(target, "w").close()
        assert find_output_file(tmp_dir) == target
