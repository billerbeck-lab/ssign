"""Unit tests for the shared parsing + TSV-io helpers under ssign_lib/.

Consolidates what used to be three near-identical TSV loaders and two
slightly-different int parsers into one shared module.
"""

from __future__ import annotations

import csv

import pytest

from ssign_app.scripts.ssign_lib.parsing import parse_int_or_none
from ssign_app.scripts.ssign_lib.tsv_io import load_tsv_by_key


class TestParseIntOrNone:
    """Tolerant int parser used across t5_passenger and t5ss_handler."""

    @pytest.mark.parametrize("value", ["", None, "   "])
    def test_empty_and_none_return_none(self, value):
        assert parse_int_or_none(value) is None

    @pytest.mark.parametrize(
        "value,expected",
        [("22", 22), ("0", 0), ("-3", -3), (" 7 ", 7)],
    )
    def test_plain_int_strings(self, value, expected):
        assert parse_int_or_none(value) == expected

    @pytest.mark.parametrize(
        "value,expected",
        [("22.0", 22), ("3.7", 3), ("-1.5", -1)],
    )
    def test_float_strings_truncate(self, value, expected):
        # int(float(x)) — matches t5_passenger/t5ss_handler legacy behaviour.
        assert parse_int_or_none(value) == expected

    @pytest.mark.parametrize("value", ["abc", "1-2-3-4-5"])
    def test_malformed_returns_none(self, value):
        # Without allow_range, '1-2-3-4-5' is a malformed integer.
        assert parse_int_or_none(value) is None

    def test_range_default_off(self):
        # '22-23' is malformed when allow_range is False.
        assert parse_int_or_none("22-23") is None

    @pytest.mark.parametrize(
        "value,expected",
        [("22-23", 22), ("22", 22), ("22.0", 22), ("22-", 22)],
    )
    def test_range_handling_when_enabled(self, value, expected):
        # SignalP CS-position formats — exercises the t5ss_handler path.
        assert parse_int_or_none(value, allow_range=True) == expected

    def test_range_empty_first_token(self):
        # '-22' splits to '' + '22'; first token is empty → None.
        assert parse_int_or_none("-22", allow_range=True) is None


class TestLoadTsvByKey:
    def _write(self, tmp_path, name, rows, fieldnames):
        path = tmp_path / name
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return path

    def test_basic_locus_tag_keyed(self, tmp_path):
        p = self._write(
            tmp_path,
            "a.tsv",
            [{"locus_tag": "A1", "x": "11"}, {"locus_tag": "A2", "x": "12"}],
            ["locus_tag", "x"],
        )
        out = load_tsv_by_key(p)
        assert set(out.keys()) == {"A1", "A2"}
        assert out["A1"]["x"] == "11"

    def test_fallback_id_columns(self, tmp_path):
        # Tool TSVs sometimes emit protein_id or seq_id instead of locus_tag —
        # mirrors the cross_validate_predictions tolerant-fallback path.
        p = self._write(
            tmp_path,
            "a.tsv",
            [{"protein_id": "P1", "x": "11"}],
            ["protein_id", "x"],
        )
        out = load_tsv_by_key(p, key_columns=("locus_tag", "protein_id", "seq_id"))
        assert set(out.keys()) == {"P1"}

    def test_skips_empty_key_rows(self, tmp_path):
        p = self._write(
            tmp_path,
            "a.tsv",
            [{"locus_tag": "A1", "x": "11"}, {"locus_tag": "", "x": "99"}],
            ["locus_tag", "x"],
        )
        out = load_tsv_by_key(p)
        assert set(out.keys()) == {"A1"}

    def test_missing_file_missing_ok_returns_empty(self, tmp_path):
        out = load_tsv_by_key(tmp_path / "does_not_exist.tsv")
        assert out == {}

    def test_missing_file_strict_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_tsv_by_key(tmp_path / "nope.tsv", missing_ok=False)

    def test_empty_path_missing_ok_returns_empty(self):
        assert load_tsv_by_key("") == {}

    def test_empty_path_strict_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            load_tsv_by_key("", missing_ok=False)

    def test_no_recognised_key_column_returns_empty(self, tmp_path):
        # Header doesn't include any of the candidate key columns.
        p = self._write(
            tmp_path,
            "a.tsv",
            [{"some_other_col": "X", "x": "11"}],
            ["some_other_col", "x"],
        )
        out = load_tsv_by_key(p, key_columns=("locus_tag", "protein_id"))
        assert out == {}

    def test_first_key_column_in_header_wins(self, tmp_path):
        # Both locus_tag and protein_id are present → use locus_tag (first in
        # key_columns) so cross_validate_predictions's fallback chain still
        # prefers the canonical column when both are emitted.
        p = self._write(
            tmp_path,
            "a.tsv",
            [{"locus_tag": "L1", "protein_id": "P_OTHER", "x": "11"}],
            ["locus_tag", "protein_id", "x"],
        )
        out = load_tsv_by_key(p, key_columns=("locus_tag", "protein_id"))
        assert set(out.keys()) == {"L1"}
