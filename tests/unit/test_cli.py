"""Tests for ssign_app/cli.py — the `ssign` entry point.

The CLI dispatches subcommands (`run`, `doctor`, `fetch-databases`) and, with
no subcommand, prints a usage hint. Testable surfaces:

- `--version` short-circuit (prints version, never runs the pipeline).
- argparse defaults (`--help` exits 0, unknown flag exits 2).
- Bare `ssign` prints the usage hint and returns 0.
- `fetch-databases` shells out to the bundled script with translated args.
"""

import sys

import pytest

# cli imports as a sub-module of the installed package; no SCRIPTS_DIR shim.
from ssign_app import cli

# ---------------------------------------------------------------------------
# --version
# ---------------------------------------------------------------------------


class TestVersion:
    def test_version_flag_prints_and_returns(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["ssign", "--version"])
        rc = cli.main()
        assert rc == 0
        assert "ssign" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# argparse defaults
# ---------------------------------------------------------------------------


class TestArgparseDefaults:
    def test_help_exits_cleanly(self, monkeypatch):
        # `--help` triggers SystemExit(0) — argparse standard behaviour.
        monkeypatch.setattr(sys, "argv", ["ssign", "--help"])
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 0

    def test_unknown_flag_exits_with_error(self, monkeypatch):
        # argparse default: unknown flag → SystemExit(2)
        monkeypatch.setattr(sys, "argv", ["ssign", "--bogus-flag"])
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 2


# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------


class TestBanner:
    def test_banner_is_present(self):
        # Pure ASCII-art wordmark spelling "ssign"; assert it is non-empty and
        # built from the underscore/pipe strokes so an accidental blanking of
        # the constant is caught.
        assert cli.BANNER.strip()
        assert "_|" in cli.BANNER


# ---------------------------------------------------------------------------
# Bare `ssign` (no subcommand) — usage hint
# ---------------------------------------------------------------------------


class TestUsageHint:
    def test_bare_ssign_prints_hint_and_returns_zero(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["ssign"])
        rc = cli.main()
        assert rc == 0
        out = capsys.readouterr().out
        assert "ssign run" in out
        assert cli.GITHUB_URL in out


# ---------------------------------------------------------------------------
# fetch-databases
# ---------------------------------------------------------------------------


class TestFetchDatabases:
    def test_find_fetch_script_resolves_in_repo(self):
        # In the dev/editable checkout the repo-root fallback must resolve.
        p = cli._find_fetch_script()
        assert p is not None and p.endswith("scripts/fetch_databases.sh")

    def test_dispatch_shells_out_with_translated_args(self, monkeypatch):
        # `ssign fetch-databases` must exec the bundled script with the flags
        # mapped through; --dry-run passes so nothing downloads.
        called = {}

        def fake_call(cmd):
            called["cmd"] = cmd
            return 0

        monkeypatch.setattr(cli.subprocess, "call", fake_call)
        monkeypatch.setattr(sys, "argv", ["ssign", "fetch-databases", "--tier", "base", "--dry-run"])
        rc = cli.main()
        assert rc == 0
        assert called["cmd"][0] == "bash"
        assert called["cmd"][1].endswith("scripts/fetch_databases.sh")
        assert called["cmd"][2:] == ["--tier", "base", "--dry-run"]

    def test_missing_script_returns_1(self, monkeypatch, capsys):
        monkeypatch.setattr(cli, "_find_fetch_script", lambda: None)
        monkeypatch.setattr(sys, "argv", ["ssign", "fetch-databases", "--tier", "base"])
        rc = cli.main()
        assert rc == 1
        assert "fetch_databases.sh not found" in capsys.readouterr().err
