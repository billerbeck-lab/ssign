"""Pre-flight tool checks in scripts/fetch_databases.sh.

The script requires each external tool only when it reaches that tool's DB
step, and `update_blastdb.pl` (BLAST Swiss-Prot, the last DB fetched) is the
LAST step — so before the pre-flight, a full-tier fetch with a missing
update_blastdb.pl would download Bakta + all HH-suite (hours) and only then
fail. `_preflight_tier` checks every tool the tier needs up front. These tests
drive the real script via subprocess with stubbed tools.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "fetch_databases.sh"


def _stub_bin(tmp_path, names):
    """Create a dir of no-op executables for each name; return the dir."""
    d = tmp_path / "stubbin"
    d.mkdir()
    for n in names:
        f = d / n
        f.write_text("#!/bin/sh\nexit 0\n")
        f.chmod(0o755)
    return d


def _run(tmp_path, stub_dir, *args):
    import os

    env = dict(os.environ)
    env["PATH"] = f"{stub_dir}:{env['PATH']}"
    # A clean env must not leak an external Bakta DB into the bypass check.
    env.pop("BAKTA_DB", None)
    return subprocess.run(
        ["bash", str(_SCRIPT), *args, "--target", str(tmp_path / "dbs")],
        capture_output=True,
        text=True,
        env=env,
    )


class TestPreflight:
    def test_full_tier_missing_nr_tool_fails_before_any_fetch(self, tmp_path):
        # update_blastdb.pl (NR, the last fetch step) missing → must abort at
        # the pre-flight, before a single '==> ' download step runs.
        if shutil.which("update_blastdb.pl"):
            pytest.skip("update_blastdb.pl present on host; cannot test its absence")
        stub = _stub_bin(tmp_path, ["bakta_db", "amrfinder", "hf", "unzip"])
        r = _run(tmp_path, stub, "--tier", "full")
        assert r.returncode != 0
        assert "update_blastdb.pl" in r.stderr
        combined = r.stdout + r.stderr
        assert "==> " not in combined, "fetch started before the pre-flight caught the missing tool"

    def test_extended_tier_does_not_require_nr_tool(self, tmp_path):
        # NR is full-only; extended's pre-flight must not demand update_blastdb.pl.
        # --dry-run so no network fetch happens (exit code is not asserted: a
        # bakta-fetching dry-run exits 1 on a pre-existing compgen/pipefail
        # quirk, unrelated to the pre-flight under test).
        stub = _stub_bin(tmp_path, [])
        r = _run(tmp_path, stub, "--tier", "extended", "--dry-run")
        assert "update_blastdb.pl" not in (r.stdout + r.stderr)

    def test_full_tier_checks_nr_tool_up_front(self, tmp_path):
        # The whole point: the NR tool is checked in the pre-flight, BEFORE the
        # first fetch step — not lazily at the last (NR) step.
        stub = _stub_bin(tmp_path, [])
        r = _run(tmp_path, stub, "--tier", "full", "--dry-run")
        out = r.stdout + r.stderr
        assert "update_blastdb.pl" in out
        assert out.index("update_blastdb.pl") < out.index("==> NCBI taxdump")


class TestBaktaVariantSkipGuard:
    """fetch_bakta keys its skip on the EXACT variant subdir (db/ for full,
    db-light/ for light), not db*. The old db* glob matched db-light during a
    --tier full fetch, so an extended→full upgrade silently kept the light DB."""

    def _seed_bakta(self, tmp_path, subdir):
        d = tmp_path / "dbs" / "bakta" / subdir
        d.mkdir(parents=True)
        (d / "version.json").write_text("{}")

    def _full_stubs(self, tmp_path):
        return _stub_bin(tmp_path, ["bakta_db", "amrfinder", "hf", "unzip", "update_blastdb.pl"])

    def test_full_not_skipped_when_only_db_light_present(self, tmp_path):
        self._seed_bakta(tmp_path, "db-light")
        r = _run(tmp_path, self._full_stubs(tmp_path), "--tier", "full", "--dry-run")
        out = r.stdout + r.stderr
        assert "bakta_db download" in out and "--type full" in out, (
            "extended→full upgrade wrongly skipped the full Bakta download"
        )
        assert "Skipping (Bakta full DB" not in out

    def test_full_skipped_when_db_full_already_present(self, tmp_path):
        self._seed_bakta(tmp_path, "db")
        r = _run(tmp_path, self._full_stubs(tmp_path), "--tier", "full", "--dry-run")
        out = r.stdout + r.stderr
        assert "Skipping (Bakta full DB" in out
