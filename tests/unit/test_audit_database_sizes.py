"""Tests for the database size audit script.

Covers the pure-logic helpers (byte summing, cumulative tier rollup,
markdown rendering). The end-to-end `audit()` function is exercised via
a small fixture tree that mimics one entry from each tier.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

_SCRIPT_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPT_DIR))

import audit_database_sizes as audit  # noqa: E402


class TestPathSizeBytes:
    def test_empty_dir_is_zero(self, tmp_path):
        assert audit.path_size_bytes(str(tmp_path)) == 0

    def test_missing_path_is_zero(self, tmp_path):
        assert audit.path_size_bytes(str(tmp_path / "nope")) == 0

    def test_single_file(self, tmp_path):
        f = tmp_path / "model.pt"
        f.write_bytes(b"x" * 4096)
        assert audit.path_size_bytes(str(f)) == 4096

    def test_sums_files_recursively(self, tmp_path):
        (tmp_path / "a.txt").write_bytes(b"x" * 100)
        sub = tmp_path / "nested"
        sub.mkdir()
        (sub / "b.txt").write_bytes(b"y" * 250)
        (sub / "c.bin").write_bytes(b"\0" * 1)
        assert audit.path_size_bytes(str(tmp_path)) == 351

    def test_symlinked_file_counted_as_link_not_target(self, tmp_path):
        target = tmp_path.parent / "big_target.bin"
        target.write_bytes(b"q" * 999)
        (tmp_path / "real.txt").write_bytes(b"z" * 50)
        (tmp_path / "link_to_big").symlink_to(target)
        # 50 bytes for real.txt + link's own (tiny) size, NOT 999
        total = audit.path_size_bytes(str(tmp_path))
        assert total < 200, f"symlink target leaked into count: got {total}"
        assert total >= 50

    def test_hardlinks_deduped(self, tmp_path):
        f1 = tmp_path / "a.bin"
        f2 = tmp_path / "b.bin"
        f1.write_bytes(b"x" * 1000)
        os.link(f1, f2)
        # Two dirent entries, one inode → count it once
        assert audit.path_size_bytes(str(tmp_path)) == 1000

    def test_handles_unreadable_subdir(self, tmp_path):
        (tmp_path / "ok.txt").write_bytes(b"x" * 10)
        secret = tmp_path / "secret"
        secret.mkdir()
        (secret / "hidden.bin").write_bytes(b"y" * 99)
        os.chmod(secret, 0o000)
        try:
            # Unreadable subdir is skipped, the rest still counts.
            total = audit.path_size_bytes(str(tmp_path))
        finally:
            os.chmod(secret, 0o700)
        assert total == 10


class TestDirSizeBytesAlias:
    def test_old_name_still_works(self, tmp_path):
        (tmp_path / "a.txt").write_bytes(b"x" * 7)
        assert audit.dir_size_bytes(str(tmp_path)) == 7


class TestCumulativeTierTotals:
    def test_each_tier_sums_self_plus_lower(self):
        rows = [
            {"name": "A", "tier": "base", "bytes": 10},
            {"name": "B", "tier": "extended", "bytes": 100},
            {"name": "C", "tier": "full", "bytes": 1000},
        ]
        totals = audit.cumulative_tier_totals(rows)
        assert totals == {"base": 10, "extended": 110, "full": 1110}

    def test_missing_bytes_treated_as_zero(self):
        rows = [
            {"name": "A", "tier": "base", "bytes": None},
            {"name": "B", "tier": "extended", "bytes": 50},
        ]
        totals = audit.cumulative_tier_totals(rows)
        assert totals["base"] == 0
        assert totals["extended"] == 50
        assert totals["full"] == 50

    def test_empty_input(self):
        assert audit.cumulative_tier_totals([]) == {"base": 0, "extended": 0, "full": 0}


class TestFormatGb:
    @pytest.mark.parametrize(
        "n_bytes,expected",
        [
            (None, "—"),
            (0, "0.0 GB"),
            (1024**3, "1.0 GB"),
            (2 * 1024**3 + 512 * 1024**2, "2.5 GB"),
        ],
    )
    def test_rounds_to_one_decimal(self, n_bytes, expected):
        assert audit._format_gb(n_bytes) == expected


class TestAudit:
    """End-to-end: synthesize a fake install root, confirm audit picks each entry."""

    def test_finds_present_and_missing(self, tmp_path):
        db_root = tmp_path / "db"
        weights_root = tmp_path / "weights"
        db_root.mkdir()
        weights_root.mkdir()

        # Bakta DB sentinel: db*/version.json under <db_root>/bakta
        bakta_dir = db_root / "bakta" / "db-light"
        bakta_dir.mkdir(parents=True)
        (bakta_dir / "version.json").write_text("{}")
        (bakta_dir / "blob").write_bytes(b"x" * 2048)

        rows = audit.audit(str(db_root), str(weights_root))
        by_name = {r["name"]: r for r in rows}

        bakta = by_name["Bakta DB"]
        assert bakta["path"] is not None
        assert bakta["bytes"] is not None and bakta["bytes"] >= 2048

        # Things not present should report path=None, bytes=None
        eggnog = by_name["EggNOG DB"]
        assert eggnog["path"] is None
        assert eggnog["bytes"] is None

    def test_file_weights_entry_counted(self, tmp_path):
        """DeepSecE checkpoint is a single .pt file (not a directory).

        Earlier code checked os.path.isdir and reported it missing even
        when the file was present. Audit must accept both files and dirs.
        """
        db_root = tmp_path / "db"
        weights_root = tmp_path / "weights"
        db_root.mkdir()
        (weights_root / "models").mkdir(parents=True)
        ckpt = weights_root / "models" / "deepsece_checkpoint.pt"
        ckpt.write_bytes(b"w" * 7777)

        rows = audit.audit(str(db_root), str(weights_root))
        deepsece = next(r for r in rows if r["name"] == "DeepSecE checkpoint")
        assert deepsece["path"] == str(ckpt)
        assert deepsece["bytes"] == 7777


class TestRenderMarkdown:
    def test_includes_tier_summary_section(self):
        rows = [
            {"kind": "database", "name": "Bakta DB", "tier": "extended", "path": "/x/bakta", "bytes": 2 * audit._GB},
        ]
        totals = audit.cumulative_tier_totals(rows)
        out = audit.render_markdown(rows, totals)
        assert "Bakta DB" in out
        assert "/x/bakta" in out
        assert "2.0 GB" in out
        assert "Cumulative per tier" in out

    def test_missing_path_shown_as_dash(self):
        rows = [
            {"kind": "database", "name": "BLAST NR", "tier": "full", "path": None, "bytes": None},
        ]
        totals = audit.cumulative_tier_totals(rows)
        out = audit.render_markdown(rows, totals)
        assert "BLAST NR" in out
        assert "| — |" in out
