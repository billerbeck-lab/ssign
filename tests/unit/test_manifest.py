"""Tests for ssign_lib/manifest.py — `Manifest` class.

A TSV-backed per-protein status tracker. Used by web-API wrappers
(SignalP, DeepLocPro via BioLib) for resume-on-restart: the per-protein
status TSV records "success" / "failed: X" / "skipped" so a re-run
picks up where the previous one died.

Surfaces:

- `__init__` — ensure `protein_id` is the first column; load existing
  TSV if present.
- `set` / `get` — round-trip a record.
- `get_successful` / `get_pending` — selectors used by the resume path.
- `save` — write back to TSV (with `extrasaction="ignore"` so unknown
  kwargs in `set` don't crash on save).
- `summary` — `{status: count}` for human-readable end-of-run logging.
- `__len__` — size after the in-memory dict.
"""

import csv
import os
import sys

import pytest


from _helpers import write_tsv  # noqa: E402
from ssign_lib.manifest import Manifest  # noqa: E402

# ---------------------------------------------------------------------------
# Construction + column ordering
# ---------------------------------------------------------------------------


class TestColumnOrdering:
    def test_protein_id_is_always_first(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status", "result"])
        assert m.columns[0] == "protein_id"

    def test_protein_id_added_when_omitted(self, tmp_dir):
        # Caller often forgets — the constructor must add it implicitly.
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        assert "protein_id" in m.columns

    def test_protein_id_not_duplicated_when_provided(self, tmp_dir):
        # If the caller explicitly includes protein_id, it must not appear twice.
        m = Manifest(
            os.path.join(tmp_dir, "m.tsv"),
            columns=["protein_id", "status"],
        )
        assert m.columns.count("protein_id") == 1

    def test_caller_columns_preserved_in_order(self, tmp_dir):
        m = Manifest(
            os.path.join(tmp_dir, "m.tsv"),
            columns=["status", "duration", "result"],
        )
        assert m.columns == ["protein_id", "status", "duration", "result"]


# ---------------------------------------------------------------------------
# set / get round-trip
# ---------------------------------------------------------------------------


class TestSetGet:
    def test_set_then_get_returns_full_record(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        m.set("GENE_001", status="success")
        assert m.get("GENE_001") == {"protein_id": "GENE_001", "status": "success"}

    def test_get_missing_returns_none(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        assert m.get("MISSING") is None

    def test_set_overwrites_existing(self, tmp_dir):
        # Re-`set` replaces the prior record entirely (no merge).
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status", "extra"])
        m.set("GENE_001", status="success", extra="first")
        m.set("GENE_001", status="failed")
        assert m.get("GENE_001") == {"protein_id": "GENE_001", "status": "failed"}

    def test_kwargs_with_unknown_columns_kept_in_memory(self, tmp_dir):
        # In-memory `set` stores everything; the column allow-list applies
        # only at save time.
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        m.set("GENE_001", status="success", debug_field="extra info")
        record = m.get("GENE_001")
        assert record["debug_field"] == "extra info"


# ---------------------------------------------------------------------------
# Selectors
# ---------------------------------------------------------------------------


class TestSelectors:
    def test_get_successful_returns_only_success_status(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        m.set("OK_1", status="success")
        m.set("FAIL_1", status="failed: timeout")
        m.set("OK_2", status="success")
        assert sorted(m.get_successful()) == ["OK_1", "OK_2"]

    @pytest.mark.parametrize(
        "stored_status, expected_pending",
        [
            (None, ["X"]),  # not in manifest at all → pending
            ("failed: API down", ["X"]),  # failed proteins retry on next run
            ("success", []),  # terminal: success
            ("skipped", []),  # terminal: skipped (caller chose not to retry)
        ],
        ids=["unprocessed", "failed", "success", "skipped"],
    )
    def test_get_pending_status_table(self, tmp_dir, stored_status, expected_pending):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        if stored_status is not None:
            m.set("X", status=stored_status)
        assert m.get_pending(["X"]) == expected_pending


# ---------------------------------------------------------------------------
# save / load round-trip
# ---------------------------------------------------------------------------


class TestPersistence:
    def test_save_then_reload_round_trips(self, tmp_dir):
        path = os.path.join(tmp_dir, "m.tsv")
        m = Manifest(path, columns=["status", "duration"])
        m.set("GENE_001", status="success", duration="1.2")
        m.set("GENE_002", status="failed: timeout")
        m.save()

        m2 = Manifest(path, columns=["status", "duration"])
        assert m2.get("GENE_001")["status"] == "success"
        assert m2.get("GENE_001")["duration"] == "1.2"
        assert m2.get("GENE_002")["status"] == "failed: timeout"

    def test_save_drops_columns_not_in_schema(self, tmp_dir):
        # `extrasaction="ignore"` — unknown in-memory keys must not appear
        # on disk.
        path = os.path.join(tmp_dir, "m.tsv")
        m = Manifest(path, columns=["status"])
        m.set("GENE_001", status="success", debug_field="leak")
        m.save()

        with open(path) as f:
            content = f.read()
        assert "leak" not in content

    def test_save_creates_parent_dirs(self, tmp_dir):
        path = os.path.join(tmp_dir, "subdir", "m.tsv")
        m = Manifest(path, columns=["status"])
        m.set("GENE_001", status="success")
        m.save()
        assert os.path.exists(path)

    def test_load_skips_blank_protein_id(self, tmp_dir):
        # A row whose `protein_id` cell is empty must not pollute the in-
        # memory dict.
        path = write_tsv(
            os.path.join(tmp_dir, "m.tsv"),
            ["protein_id", "status"],
            [
                {"protein_id": "GOOD", "status": "success"},
                {"protein_id": "", "status": "ghost"},
            ],
        )
        m = Manifest(path, columns=["status"])
        assert m.get("GOOD") is not None
        assert m.get("") is None

    def test_load_missing_file_starts_empty(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "nonexistent.tsv"), columns=["status"])
        assert len(m) == 0

    def test_save_is_atomic_under_write_failure(self, tmp_dir, monkeypatch):
        # Existing manifest on disk; subsequent save() that fails mid-write
        # MUST leave the original file intact, not a truncated half-write.
        path = os.path.join(tmp_dir, "m.tsv")
        m = Manifest(path, columns=["status"])
        m.set("ORIG_001", status="success")
        m.save()

        # Sanity: original is on disk.
        with open(path) as f:
            original_content = f.read()
        assert "ORIG_001" in original_content

        # Now stage a save that explodes mid-way through.
        m.set("NEW_001", status="success")
        original_writerow = csv.DictWriter.writerow

        def explode(self, row):
            if row.get("protein_id") == "NEW_001":
                raise OSError("simulated disk full")
            return original_writerow(self, row)

        monkeypatch.setattr(csv.DictWriter, "writerow", explode)

        with pytest.raises(OSError, match="simulated disk full"):
            m.save()

        # Original file must be unchanged (atomicity guarantee).
        with open(path) as f:
            assert f.read() == original_content
        # And the .tmp partial-write must have been cleaned up.
        assert not os.path.exists(path + ".tmp")


# ---------------------------------------------------------------------------
# summary + __len__
# ---------------------------------------------------------------------------


class TestSummaryAndLen:
    def test_summary_counts_per_status(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        m.set("a", status="success")
        m.set("b", status="success")
        m.set("c", status="failed: timeout")
        m.set("d", status="skipped")
        s = m.summary()
        assert s["success"] == 2
        assert s["failed: timeout"] == 1
        assert s["skipped"] == 1

    def test_summary_treats_missing_status_as_unknown(self, tmp_dir):
        # `set` without a status= kwarg → entry has no status key.
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        m.set("X")  # no status
        s = m.summary()
        assert s.get("unknown") == 1

    def test_len_matches_entry_count(self, tmp_dir):
        m = Manifest(os.path.join(tmp_dir, "m.tsv"), columns=["status"])
        m.set("a", status="success")
        m.set("b", status="success")
        m.set("c", status="success")
        assert len(m) == 3
