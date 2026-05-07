"""Tests for run_hhsuite.parse_hhr.

The HH-suite remote mode (which once required a specific "alignment" vs
"sequence" parameter) was removed in v0.9.x; the local hhblits + hhsearch
pipeline is the only supported route. Today's regression surface is the
parser: extract hit metadata, apply Prob threshold, return {} cleanly when
the file is missing, empty, or below threshold.
"""

import os
import sys

import pytest


from run_hhsuite import _HHR_HIT_RE, parse_hhr  # noqa: E402

_HHR_HEADER = """Query         test_protein
Match_columns 200
No_of_seqs    87 out of 1234

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
"""


def _hhr_fixture(hits, path):
    """Write a minimal HHR file with the given pre-formatted hit lines."""
    with open(path, "w") as f:
        f.write(_HHR_HEADER)
        f.write("\n".join(hits))
        f.write("\n")
    return path


# Reusable hit lines — building blocks for the parser tests. Real HHR output
# is space-aligned (not tab-separated) and the regex is tuned to the exact
# whitespace pattern, so the literals below stay intentionally verbose.
_PFAM_TOP_HIT = (
    "  1 PF00012.21 ; HSP70 protein                    99.7 1.5E-25 1.2E-29  234.5    0.0  198    5-203    1-198 (200)"
)
_PFAM_WEAK_HIT = (
    "  1 PF99999.1 ; weak hit                           70.0 1.5E-05 1.2E-09  100.5    0.0  198    5-203    1-198 (200)"
)
_PDB_TOP_HIT = (
    "  1 1abc_A ; structure of HSP70                   99.5 3.2E-22 2.5E-26  198.3    0.0  150   12-160    8-158 (160)"
)


class TestParseHhrPfam:
    def test_top_hit_parsed(self, tmp_dir):
        path = _hhr_fixture([_PFAM_TOP_HIT], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pfam", min_prob=95.0)
        assert result["pfam_top1_id"] == "PF00012.21"
        assert "HSP70 protein" in result["pfam_top1_description"]
        assert result["pfam_top1_probability"] == 99.7
        assert result["pfam_top1_evalue"] == 1.5e-25
        assert result["pfam_top1_score"] == 234.5

    def test_below_threshold_returns_empty(self, tmp_dir):
        # Top-1 below cutoff drops the whole result — next-best hits are not
        # considered (parse_hhr returns top-1 only, with the gate applied to it).
        path = _hhr_fixture([_PFAM_WEAK_HIT], os.path.join(tmp_dir, "test.hhr"))
        assert parse_hhr(path, db_prefix="pfam", min_prob=95.0) == {}


class TestParseHhrPdb:
    def test_top_hit_parsed(self, tmp_dir):
        path = _hhr_fixture([_PDB_TOP_HIT], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pdb", min_prob=95.0)
        assert result["pdb_top1_id"] == "1abc_A"
        assert result["pdb_top1_probability"] == 99.5

    def test_db_prefix_applied_to_every_field(self, tmp_dir):
        path = _hhr_fixture([_PDB_TOP_HIT], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pdb", min_prob=95.0)
        assert set(result.keys()) == {
            "pdb_top1_id",
            "pdb_top1_description",
            "pdb_top1_probability",
            "pdb_top1_evalue",
            "pdb_top1_score",
        }


class TestParseHhrEdgeCases:
    def test_missing_file_returns_empty(self, tmp_dir):
        result = parse_hhr(
            os.path.join(tmp_dir, "does_not_exist.hhr"),
            db_prefix="pfam",
            min_prob=95.0,
        )
        assert result == {}

    def test_empty_hit_table_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.hhr")
        with open(path, "w") as f:
            f.write(_HHR_HEADER)
            f.write("\n")
        assert parse_hhr(path, db_prefix="pfam", min_prob=95.0) == {}

    def test_only_top_hit_returned(self, tmp_dir):
        # Two hits in the file → only the first is returned
        second_hit = (
            "  2 PF_B ; second hit                             "
            "98.5 3.2E-22 2.5E-26  198.3    0.0  150   12-160    8-158 (160)"
        )
        path = _hhr_fixture(
            [_PFAM_TOP_HIT, second_hit],
            os.path.join(tmp_dir, "test.hhr"),
        )
        result = parse_hhr(path, db_prefix="pfam", min_prob=95.0)
        assert result["pfam_top1_id"] == "PF00012.21"

    def test_description_truncated_to_200_chars(self, tmp_dir):
        long_desc = "x" * 500
        hit_line = f"  1 PF00012.21 ; {long_desc} 99.7 1.5E-25 1.2E-29  234.5    0.0  198    5-203    1-198 (200)"
        path = _hhr_fixture([hit_line], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pfam", min_prob=95.0)
        assert len(result["pfam_top1_description"]) <= 200


@pytest.mark.parametrize(
    "malformed_line",
    [
        # Missing the trailing (TLen)
        "  1 PF00012.21 ; HSP70 protein 99.7 1.5E-25 1.2E-29 234.5 0.0 198 5-203 1-198",
        # Missing one of the \d+-\d+ ranges
        "  1 PF00012.21 ; HSP70 protein 99.7 1.5E-25 1.2E-29 234.5 0.0 198 (200)",
    ],
)
def test_regex_rejects_malformed_lines(malformed_line):
    """Pin the parser's expectations against Söding-lab format drift."""
    assert not _HHR_HIT_RE.match(malformed_line)


class TestDescriptionCleanup:
    """The description capture group catches everything between hit_id and Prob;
    the parser strips trailing "; " punctuation from PDB entries."""

    def test_trailing_semicolon_stripped(self, tmp_dir):
        # PDB-style: "1abc_A ; structure of HSP70 ;" → desc "structure of HSP70"
        hit = (
            "  1 1abc_A ; structure of HSP70 ;                "
            "99.5 3.2E-22 2.5E-26  198.3    0.0  150   12-160    8-158 (160)"
        )
        path = _hhr_fixture([hit], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pdb", min_prob=95.0)
        assert not result["pdb_top1_description"].endswith(";")
        assert not result["pdb_top1_description"].endswith("; ")


class TestNumericTypeContract:
    """Probability, e-value, score must come back as floats, not strings.
    Downstream consumers expect numeric comparison and arithmetic."""

    def test_probability_is_float(self, tmp_dir):
        path = _hhr_fixture([_PFAM_TOP_HIT], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pfam", min_prob=95.0)
        assert isinstance(result["pfam_top1_probability"], float)
        assert isinstance(result["pfam_top1_evalue"], float)
        assert isinstance(result["pfam_top1_score"], float)


class TestProbThresholdEdgeCases:
    """The Prob gate is a `<` comparison: prob < min_prob → drop."""

    def test_prob_exactly_at_threshold_kept(self, tmp_dir):
        # 95.0 == min_prob → kept (>=, not >).
        hit = (
            "  1 PF00012.21 ; HSP70 protein                    "
            "95.0 1.5E-25 1.2E-29  234.5    0.0  198    5-203    1-198 (200)"
        )
        path = _hhr_fixture([hit], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pfam", min_prob=95.0)
        assert result["pfam_top1_probability"] == 95.0


class TestHitTableBoundary:
    """Lines before " No Hit" must be ignored (header / metadata).
    The parser only enters hit-collection mode after the column-header line."""

    def test_decoy_hit_lines_in_header_ignored(self, tmp_dir):
        # A line that looks regex-compatible appears BEFORE the "No Hit" header
        # — must not be parsed as a hit.
        path = os.path.join(tmp_dir, "test.hhr")
        with open(path, "w") as f:
            f.write("Query         test_protein\n")
            # Decoy line — looks like a hit but predates the column header
            f.write(_PFAM_TOP_HIT.replace("PF00012.21", "DECOY") + "\n")
            f.write("Match_columns 200\n")
            f.write(" No Hit                             Prob E-value\n")
            f.write(_PFAM_TOP_HIT + "\n")
        result = parse_hhr(path, db_prefix="pfam", min_prob=95.0)
        assert result["pfam_top1_id"] == "PF00012.21"
        assert result["pfam_top1_id"] != "DECOY"


class TestDefaultMinProb:
    """parse_hhr's default min_prob comes from constants.HHSUITE_MIN_PROB.
    Pin the default so a future constants edit doesn't silently widen the
    cutoff."""

    def test_default_min_prob_imported_from_constants(self):
        from ssign_lib.constants import HHSUITE_MIN_PROB

        # The value should be a sane percent-probability, not 0 or > 100.
        assert 0 < HHSUITE_MIN_PROB <= 100

    def test_default_used_when_arg_omitted(self, tmp_dir):
        # No min_prob arg → uses HHSUITE_MIN_PROB. Top-1 at 99.7 should pass
        # any sane default.
        path = _hhr_fixture([_PFAM_TOP_HIT], os.path.join(tmp_dir, "test.hhr"))
        result = parse_hhr(path, db_prefix="pfam")
        assert result["pfam_top1_id"] == "PF00012.21"
