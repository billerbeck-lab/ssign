"""Unit tests for ssign_app/scripts/doctor.py — Resources section.

Other doctor sections (Python packages, External binaries, Databases,
Model weights) are exercised by tests/integration/test_doctor_*.py
against a real install. Here we only pin the new informational
Resources section so HPC users can spot scheduler-throttle bugs at
a glance.
"""

import io


class TestReportResources:
    def test_emits_cpu_ram_gpu_and_plme_lines(self, monkeypatch):
        from ssign_app.scripts import doctor

        # Pin the underlying detectors so the test is host-independent.
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_cpu_count", lambda: 8)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_ram_gb", lambda: 32.0)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.auto_batch_size_from_vram", lambda: 16)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.probe_cuda_device", lambda: ("FAKE-GPU", 24.0))
        monkeypatch.setattr("os.cpu_count", lambda: 16)

        buf = io.StringIO()
        doctor.report_resources(buf)
        out = buf.getvalue()

        assert "Resources" in out
        assert "effective: 8" in out and "host: 16" in out
        assert "32.0 GB" in out
        assert "FAKE-GPU, 24.0 GiB VRAM" in out
        assert "auto batch size: 16" in out

    def test_flags_scheduler_throttle(self, monkeypatch):
        # When effective_cpu_count < host cpu_count, doctor should call
        # out the gap — that's the bug we most want HPC users to notice.
        from ssign_app.scripts import doctor

        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_cpu_count", lambda: 4)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_ram_gb", lambda: 16.0)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.auto_batch_size_from_vram", lambda: 4)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.probe_cuda_device", lambda: (None, None))
        monkeypatch.setattr("os.cpu_count", lambda: 32)

        buf = io.StringIO()
        doctor.report_resources(buf)
        assert "scheduler is restricting" in buf.getvalue()

    def test_omits_throttle_note_when_full_allocation(self, monkeypatch):
        from ssign_app.scripts import doctor

        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_cpu_count", lambda: 16)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.effective_ram_gb", lambda: 32.0)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.auto_batch_size_from_vram", lambda: 4)
        monkeypatch.setattr("ssign_app.scripts.ssign_lib.resources.probe_cuda_device", lambda: (None, None))
        monkeypatch.setattr("os.cpu_count", lambda: 16)

        buf = io.StringIO()
        doctor.report_resources(buf)
        assert "scheduler is restricting" not in buf.getvalue()


class TestBinaryTierMapping:
    """HH-suite + BLAST+ are only needed at the full install tier (full
    HH-suite databases and BLAST nr ship there). Doctor used to flag
    them as missing at extended-tier, causing false FAILs for the most
    common ssign install on HPC. Pin the tier mapping."""

    def test_hhsuite_and_blast_are_full_tier_only(self):
        from ssign_app.scripts.ssign_lib.dependency_manifest import (
            EXTERNAL_BINARIES,
            binaries_for_tier,
        )

        by_name = {b.name: b for b in EXTERNAL_BINARIES}
        # The two HH-suite binaries and BLAST+ all require their full-tier
        # databases to do anything useful, so they live at tier="full".
        for name in ("HH-suite — hhsearch", "HH-suite — hhblits", "BLAST+"):
            assert by_name[name].tier == "full", f"{name} should be full-tier"

        # At extended tier, none of the three appear in the required list.
        extended_names = {b.name for b in binaries_for_tier("extended")}
        assert "HH-suite — hhsearch" not in extended_names
        assert "HH-suite — hhblits" not in extended_names
        assert "BLAST+" not in extended_names

        # At full tier they DO appear (otherwise we've broken full-tier doctor).
        full_names = {b.name for b in binaries_for_tier("full")}
        assert "HH-suite — hhsearch" in full_names
        assert "HH-suite — hhblits" in full_names
        assert "BLAST+" in full_names


class TestManifestEnvVarPairing:
    """When an ExternalBinary's install_dir_env matches a DatabasePath's
    env_var, the doctor's DB-root fallback uses the DB to locate the
    binary. A silent rename of either side breaks that fallback with no
    test failure. Pin every paired env var so a rename has to update
    both sides (or explicitly drop the pairing in this test)."""

    def test_known_paired_env_vars_stay_paired(self):
        from ssign_app.scripts.ssign_lib.dependency_manifest import (
            EXTERNAL_BINARIES,
            find_db_by_env_var,
        )

        # The pairings we actually rely on. Add to this list when a new
        # binary↔DB pairing is introduced; the test then guards both
        # sides.
        EXPECTED_PAIRINGS = {
            "InterProScan": "SSIGN_INTERPROSCAN_PATH",
        }
        for binary_name, env_var in EXPECTED_PAIRINGS.items():
            binary = next(b for b in EXTERNAL_BINARIES if b.name == binary_name)
            assert binary.install_dir_env == env_var, (
                f"{binary_name} expected install_dir_env={env_var!r}, got {binary.install_dir_env!r}"
            )
            dbp = find_db_by_env_var(env_var)
            assert dbp is not None, (
                f"no DatabasePath has env_var={env_var!r}; doctor's DB-root "
                f"fallback for {binary_name} silently won't fire"
            )


class TestCheckExternalBinaryDbRootFallback:
    """When a binary's `install_dir_env` is unset at doctor-invocation
    time but the matching DatabasePath resolves under db_root (e.g.
    SSIGN_INTERPROSCAN_PATH not exported in the current shell, but the
    install lives under <db_root>/interproscan/interproscan-*/),
    doctor should still find the binary instead of false-flagging it
    as missing. Mirrors what the runner does at execution time."""

    def _ips_binary(self):
        from ssign_app.scripts.ssign_lib.dependency_manifest import EXTERNAL_BINARIES

        return next(b for b in EXTERNAL_BINARIES if b.name == "InterProScan")

    def test_resolves_via_db_root_when_env_var_unset(self, tmp_path, monkeypatch):
        from ssign_app.scripts.doctor import check_external_binary

        # Build an IPS-like layout under tmp_path acting as db_root:
        #   <db_root>/interproscan/interproscan-5.77-108.0/
        #       interproscan.sh             (executable stub)
        #       interproscan.properties     (sentinel for resolve_path)
        ips_dir = tmp_path / "interproscan" / "interproscan-5.77-108.0"
        ips_dir.mkdir(parents=True)
        (ips_dir / "interproscan.properties").write_text("")
        bin_path = ips_dir / "interproscan.sh"
        bin_path.write_text("#!/bin/sh\nexit 0\n")
        bin_path.chmod(0o755)

        monkeypatch.delenv("SSIGN_INTERPROSCAN_PATH", raising=False)
        monkeypatch.setenv("PATH", "")  # ensure shutil.which misses

        result = check_external_binary(self._ips_binary(), db_root=str(tmp_path))
        assert result.ok, f"expected DB-root fallback to find IPS; got {result.detail}"
        assert str(bin_path) in result.detail

    def test_still_fails_when_neither_env_nor_db_root_has_it(self, tmp_path, monkeypatch):
        from ssign_app.scripts.doctor import check_external_binary

        monkeypatch.delenv("SSIGN_INTERPROSCAN_PATH", raising=False)
        monkeypatch.setenv("PATH", "")

        result = check_external_binary(self._ips_binary(), db_root=str(tmp_path))
        assert not result.ok
        assert "not on PATH" in result.detail

    def test_db_root_unset_doesnt_crash(self, monkeypatch):
        # Empty db_root (default) must not raise when looking up the
        # manifest — exercises the `if b.install_dir_env and db_root`
        # guard in check_external_binary.
        from ssign_app.scripts.doctor import check_external_binary

        monkeypatch.delenv("SSIGN_INTERPROSCAN_PATH", raising=False)
        monkeypatch.setenv("PATH", "")
        result = check_external_binary(self._ips_binary())
        assert not result.ok


class TestDatabasePathResolvePrefersFull:
    """resolve_path with a `db*/version.json` sentinel (Bakta) must return the
    full `db/` over `db-light/` deterministically. glob order is filesystem-
    arbitrary, so a naive matches[0] could hand a full-tier run the light DB
    when a user upgraded extended→full and both subdirs are present."""

    def _bakta_db(self):
        from ssign_app.scripts.ssign_lib.dependency_manifest import find_db_by_env_var

        return find_db_by_env_var("SSIGN_BAKTA_DB")

    def _clear_env(self, monkeypatch):
        monkeypatch.delenv("SSIGN_BAKTA_DB", raising=False)
        monkeypatch.delenv("BAKTA_DB", raising=False)

    def test_prefers_db_full_over_db_light(self, tmp_path, monkeypatch):
        self._clear_env(monkeypatch)
        bakta = tmp_path / "bakta"
        for sub in ("db", "db-light"):
            (bakta / sub).mkdir(parents=True)
            (bakta / sub / "version.json").write_text("{}")
        resolved = self._bakta_db().resolve_path(str(tmp_path))
        assert resolved == str(bakta / "db"), f"expected full db/, got {resolved}"

    def test_light_only_still_resolves(self, tmp_path, monkeypatch):
        self._clear_env(monkeypatch)
        light = tmp_path / "bakta" / "db-light"
        light.mkdir(parents=True)
        (light / "version.json").write_text("{}")
        resolved = self._bakta_db().resolve_path(str(tmp_path))
        assert resolved == str(light), f"expected db-light/, got {resolved}"
