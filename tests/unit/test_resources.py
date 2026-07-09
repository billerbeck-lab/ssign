"""Tests for ssign_lib.resources — cgroup/scheduler-aware CPU and RAM detection.

These tests pin down the lookup priorities and the parse rules. The actual
production behavior on PBS/SLURM is hard to mock end-to-end (qstat must be
on PATH), so we mock the integration points and verify each branch in
isolation.
"""

from __future__ import annotations

import os

import pytest
from ssign_lib.resources import (
    _CGROUP_UNLIMITED_BYTES,
    _parse_size_to_gb,
    effective_cpu_count,
    effective_ram_gb,
    is_remote_filesystem,
    stage_prefix_files_to_local_ssd_if_remote,
    stage_to_local_ssd_if_remote,
)


class TestParseSizeToGb:
    @pytest.mark.parametrize(
        "raw,expected_gb",
        [
            ("32gb", 32.0),
            ("32GB", 32.0),
            ("32Gib", 32.0),
            ("4096mb", 4.0),  # 4096 MiB = 4 GiB
            ("128", 128 / 2**30),  # bare bytes
            ("0.5tb", 512.0),
        ],
    )
    def test_known_units(self, raw, expected_gb):
        assert _parse_size_to_gb(raw) == pytest.approx(expected_gb, rel=1e-3)

    def test_garbage_returns_none(self):
        assert _parse_size_to_gb("hello") is None
        assert _parse_size_to_gb("") is None


class TestEffectiveCpuCount:
    def test_returns_affinity_when_available(self, monkeypatch):
        monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: {0, 1, 2, 3}, raising=False)
        assert effective_cpu_count() == 4

    def test_falls_back_when_affinity_oserror(self, monkeypatch):
        def boom(_pid):
            raise OSError("no")

        monkeypatch.setattr(os, "sched_getaffinity", boom, raising=False)
        monkeypatch.setattr(os, "cpu_count", lambda: 16)
        assert effective_cpu_count() == 16

    def test_never_returns_zero(self, monkeypatch):
        monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: set(), raising=False)
        # Empty affinity (shouldn't happen but just in case) clamps to 1.
        assert effective_cpu_count() >= 1


class TestResolveThreads:
    """Standard "None → effective_cpu_count(), else passthrough" shape used
    by every tool wrapper that accepts a `threads` argument."""

    def test_none_resolves_to_effective_cpu_count(self, monkeypatch):
        from ssign_lib.resources import resolve_threads

        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 12)
        monkeypatch.delenv("SSIGN_PARALLEL_GROUP_SIZE", raising=False)
        assert resolve_threads(None) == 12

    def test_explicit_int_passes_through(self):
        from ssign_lib.resources import resolve_threads

        assert resolve_threads(7) == 7

    def test_zero_or_negative_clamps_to_one(self):
        from ssign_lib.resources import resolve_threads

        assert resolve_threads(0) == 1
        assert resolve_threads(-3) == 1

    def test_none_honors_parallel_group_share(self, monkeypatch):
        # The bug this guards against: K-12 RTX6000 annotation stage on CX3
        # showed cput/walltime = 6.6x because IPS/EggNOG/pLM-BLAST each
        # grabbed 16 cores on a 16-core box. resolve_threads(None) must
        # divide by SSIGN_PARALLEL_GROUP_SIZE when set, not return the
        # whole effective_cpu_count.
        from ssign_lib.resources import resolve_threads

        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 16)
        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "4")
        assert resolve_threads(None) == 4

    def test_explicit_int_overrides_parallel_share(self, monkeypatch):
        # User passing --threads N opts out of the auto-share; their value
        # wins even inside a parallel group.
        from ssign_lib.resources import resolve_threads

        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 16)
        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "4")
        assert resolve_threads(8) == 8


class TestParallelShareCpus:
    """When SSIGN_PARALLEL_GROUP_SIZE is set (by the runner before
    launching the DLP+DSE+SignalP parallel group), each wrapper should
    take 1/N of the effective CPU budget instead of fighting over the
    full allocation. Unset = standalone wrapper invocation; full budget."""

    def test_unset_returns_full_cpu_count(self, monkeypatch):
        from ssign_lib.resources import parallel_share_cpus

        monkeypatch.delenv("SSIGN_PARALLEL_GROUP_SIZE", raising=False)
        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 16)
        assert parallel_share_cpus() == 16

    def test_group_of_3_splits_evenly(self, monkeypatch):
        from ssign_lib.resources import parallel_share_cpus

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "3")
        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 16)
        # 16 // 3 == 5; each of DLP/DSE/SignalP gets 5 cores
        assert parallel_share_cpus() == 5

    def test_group_of_1_returns_full_count(self, monkeypatch):
        from ssign_lib.resources import parallel_share_cpus

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "1")
        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 16)
        assert parallel_share_cpus() == 16

    def test_low_cpu_count_clamps_to_one(self, monkeypatch):
        from ssign_lib.resources import parallel_share_cpus

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "3")
        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 2)
        # 2 // 3 == 0, must clamp up to 1
        assert parallel_share_cpus() == 1

    def test_garbage_value_returns_full_count(self, monkeypatch):
        from ssign_lib.resources import parallel_share_cpus

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "not-an-int")
        monkeypatch.setattr("ssign_lib.resources.effective_cpu_count", lambda: 16)
        assert parallel_share_cpus() == 16


class TestParallelShareRamGb:
    """RAM analogue of parallel_share_cpus: a wrapper that sizes memory to
    fill the node (EggNOG --dbmem / DIAMOND --block_size) must take 1/N of
    the RAM when it runs inside the annotation group, else the concurrent
    tools collectively overcommit and OOM-kill the biggest allocator."""

    def test_unset_returns_full_ram(self, monkeypatch):
        from ssign_lib.resources import parallel_share_ram_gb

        monkeypatch.delenv("SSIGN_PARALLEL_GROUP_SIZE", raising=False)
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 120.0)
        assert parallel_share_ram_gb() == 120.0

    def test_group_of_6_splits_evenly(self, monkeypatch):
        from ssign_lib.resources import parallel_share_ram_gb

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "6")
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 120.0)
        # 120 / 6 == 20 GB — below EggNOG's 50 GB --dbmem threshold, so the
        # full-tier 6-tool group no longer overcommits the node.
        assert parallel_share_ram_gb() == 20.0

    def test_group_of_1_returns_full_ram(self, monkeypatch):
        from ssign_lib.resources import parallel_share_ram_gb

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "1")
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 120.0)
        assert parallel_share_ram_gb() == 120.0

    def test_garbage_value_returns_full_ram(self, monkeypatch):
        from ssign_lib.resources import parallel_share_ram_gb

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "not-an-int")
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 88.0)
        assert parallel_share_ram_gb() == 88.0


class TestParallelGroupContext:
    """parallel_group() context-manager-ifies the SSIGN_PARALLEL_GROUP_SIZE
    env-var save/set/restore dance the runner used to do inline. It must
    restore the prior value (even None → unset) on exit, including when
    the body raises."""

    def test_sets_var_inside_restores_on_exit(self, monkeypatch):
        from ssign_lib.resources import parallel_group

        monkeypatch.delenv("SSIGN_PARALLEL_GROUP_SIZE", raising=False)
        with parallel_group(3):
            assert os.environ["SSIGN_PARALLEL_GROUP_SIZE"] == "3"
        assert "SSIGN_PARALLEL_GROUP_SIZE" not in os.environ

    def test_restores_prior_value(self, monkeypatch):
        from ssign_lib.resources import parallel_group

        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "OLD")
        with parallel_group(5):
            assert os.environ["SSIGN_PARALLEL_GROUP_SIZE"] == "5"
        assert os.environ["SSIGN_PARALLEL_GROUP_SIZE"] == "OLD"

    def test_restores_on_exception(self, monkeypatch):
        from ssign_lib.resources import parallel_group

        monkeypatch.delenv("SSIGN_PARALLEL_GROUP_SIZE", raising=False)
        with pytest.raises(RuntimeError):
            with parallel_group(2):
                raise RuntimeError("simulated step failure")
        assert "SSIGN_PARALLEL_GROUP_SIZE" not in os.environ


class TestProbeCudaDevice:
    """Single source of truth for "is there a CUDA GPU and how big is its VRAM"."""

    def _patch_torch(self, monkeypatch, *, has_torch=True, available=True, total_bytes=48 * 1024**3, name="FAKE"):
        import sys
        import types

        if not has_torch:
            monkeypatch.setitem(sys.modules, "torch", None)
            return
        fake = types.SimpleNamespace()
        fake.cuda = types.SimpleNamespace(
            is_available=lambda: available,
            get_device_name=lambda i: name,
            get_device_properties=lambda i: types.SimpleNamespace(total_memory=total_bytes),
        )
        monkeypatch.setitem(sys.modules, "torch", fake)

    def test_returns_name_and_gib_when_cuda_visible(self, monkeypatch):
        from ssign_lib.resources import probe_cuda_device

        self._patch_torch(monkeypatch, available=True, total_bytes=24 * 1024**3, name="RTX 4090")
        name, gib = probe_cuda_device()
        assert name == "RTX 4090"
        assert abs(gib - 24.0) < 0.01

    def test_returns_none_pair_when_cuda_unavailable(self, monkeypatch):
        from ssign_lib.resources import probe_cuda_device

        self._patch_torch(monkeypatch, available=False)
        assert probe_cuda_device() == (None, None)

    def test_returns_none_pair_when_properties_raise(self, monkeypatch):
        # Driver/runtime mismatch can throw on get_device_properties.
        import sys
        import types

        def _boom(i):
            raise RuntimeError("CUDA driver too old")

        fake = types.SimpleNamespace()
        fake.cuda = types.SimpleNamespace(
            is_available=lambda: True,
            get_device_name=lambda i: "X",
            get_device_properties=_boom,
        )
        monkeypatch.setitem(sys.modules, "torch", fake)
        from ssign_lib.resources import probe_cuda_device

        assert probe_cuda_device() == (None, None)


class TestEffectiveRamGb:
    """Lookup priority: env override < SLURM < PBS qstat < cgroup < psutil. Returns MIN."""

    def test_env_override_wins_when_smallest(self, monkeypatch):
        # 1 GB is guaranteed smaller than any CI runner or dev-host total
        # RAM. Earlier this test set 16 GB and asserted exact equality —
        # GitHub runners report ~15.61 GB total (OS overhead), so psutil
        # won the min and the test failed.
        monkeypatch.setenv("SSIGN_MAX_RAM_GB", "1")
        monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
        monkeypatch.delenv("PBS_JOBID", raising=False)
        assert effective_ram_gb() == pytest.approx(1.0)

    def test_slurm_alloc_is_parsed_in_mb(self, monkeypatch):
        monkeypatch.delenv("SSIGN_MAX_RAM_GB", raising=False)
        # 1024 MiB = 1 GB; smaller than any plausible host total, so MIN
        # logic should return the SLURM-derived value.
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "1024")
        monkeypatch.delenv("PBS_JOBID", raising=False)
        assert effective_ram_gb() == pytest.approx(1.0)

    def test_garbage_env_is_ignored(self, monkeypatch, caplog):
        monkeypatch.setenv("SSIGN_MAX_RAM_GB", "not-a-number")
        monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
        monkeypatch.delenv("PBS_JOBID", raising=False)
        # Doesn't crash; falls through to other sources.
        assert effective_ram_gb() > 0

    def test_cgroup_unlimited_sentinel_is_ignored(self):
        # The sentinel is 2**63 - 4096 (close enough to 2**63 - 1 modulo
        # page rounding). Confirms the production guard catches both.
        assert _CGROUP_UNLIMITED_BYTES < 2**63

    def test_returns_zero_when_nothing_available(self, monkeypatch):
        # Strip everything: env, scheduler, cgroup, psutil.
        for var in ("SSIGN_MAX_RAM_GB", "SLURM_MEM_PER_NODE", "PBS_JOBID"):
            monkeypatch.delenv(var, raising=False)
        # Block cgroup paths and psutil so we exercise the empty-candidates branch.
        monkeypatch.setattr("builtins.open", lambda *a, **kw: (_ for _ in ()).throw(OSError("no")))
        import sys

        monkeypatch.setitem(sys.modules, "psutil", None)
        # psutil=None makes the `import psutil` line raise TypeError, which the
        # bare except catches. Final result is 0.0.
        assert effective_ram_gb() == 0.0


class TestIsRemoteFilesystem:
    """is_remote_filesystem should detect networked FS via /proc/mounts."""

    def _mock_mounts(self, monkeypatch, mounts_text):
        from io import StringIO

        def fake_open(path, *a, **kw):
            if path == "/proc/mounts":
                return StringIO(mounts_text)
            raise OSError("no")

        monkeypatch.setattr("builtins.open", fake_open)

    def test_gpfs_path_is_remote(self, monkeypatch):
        self._mock_mounts(monkeypatch, "rds /rds gpfs rw 0 0\n")
        monkeypatch.setattr("os.path.realpath", lambda p: "/rds/general/user/x/db")
        assert is_remote_filesystem("/rds/general/user/x/db") is True

    def test_local_xfs_is_not_remote(self, monkeypatch):
        self._mock_mounts(monkeypatch, "/dev/sda1 / xfs rw 0 0\n")
        monkeypatch.setattr("os.path.realpath", lambda p: "/home/user/db")
        assert is_remote_filesystem("/home/user/db") is False

    def test_nested_mount_picks_longest_prefix(self, monkeypatch):
        # /home is xfs but /home/user/scratch is mounted as gpfs underneath.
        self._mock_mounts(monkeypatch, "/dev/sda1 /home xfs rw 0 0\nrds /home/user/scratch gpfs rw 0 0\n")
        monkeypatch.setattr("os.path.realpath", lambda p: "/home/user/scratch/db")
        assert is_remote_filesystem("/home/user/scratch/db") is True

    def test_unreadable_mounts_returns_false(self, monkeypatch):
        monkeypatch.setattr("builtins.open", lambda *a, **kw: (_ for _ in ()).throw(OSError("no")))
        # No info → conservative: assume local, skip cache.
        assert is_remote_filesystem("/any/path") is False


class TestStageToLocalSsdIfRemote:
    """stage_to_local_ssd_if_remote: copy only when src is networked."""

    def _make_db(self, src, files, nbytes=128):
        os.makedirs(src, exist_ok=True)
        for name in files:
            with open(os.path.join(src, name), "wb") as f:
                f.write(b"\0" * nbytes)

    def test_skips_copy_when_src_is_local(self, tmp_dir, monkeypatch):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        self._make_db(src, ("eggnog.db",))
        os.makedirs(cache)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: False)
        out = stage_to_local_ssd_if_remote(src, cache, required=("eggnog.db",))
        assert out == src  # unchanged
        assert not os.listdir(cache)  # nothing copied

    def test_copies_when_src_is_remote(self, tmp_dir, monkeypatch):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        self._make_db(src, ("eggnog.db", "extra.dmnd"))
        os.makedirs(cache)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: True)
        out = stage_to_local_ssd_if_remote(
            src,
            cache,
            required=("eggnog.db",),
            optional=("extra.dmnd",),
            min_free_gb=0.0,  # tmpfs in tests
        )
        assert out != src
        assert os.path.exists(os.path.join(out, "eggnog.db"))
        assert os.path.exists(os.path.join(out, "extra.dmnd"))

    def test_raises_on_missing_required_file(self, tmp_dir, monkeypatch):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        self._make_db(src, ("eggnog.db",))
        os.makedirs(cache)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: True)
        with pytest.raises(FileNotFoundError, match="missing_required_file.fasta"):
            stage_to_local_ssd_if_remote(
                src,
                cache,
                required=("eggnog.db", "missing_required_file.fasta"),
                min_free_gb=0.0,
            )

    def test_falls_back_when_cache_too_small(self, tmp_dir, monkeypatch, caplog):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        self._make_db(src, ("eggnog.db",))
        os.makedirs(cache)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: True)
        # 1 PB free required — definitely won't have it locally.
        out = stage_to_local_ssd_if_remote(src, cache, required=("eggnog.db",), min_free_gb=1_000_000.0)
        assert out == src  # falls back to remote path


class TestStagePrefixFiles:
    """stage_prefix_files_to_local_ssd_if_remote: ffindex DBs (HH-suite)."""

    def _write(self, path, nbytes=64):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as f:
            f.write(b"\0" * nbytes)

    def test_skips_when_local(self, tmp_dir, monkeypatch):
        prefix = os.path.join(tmp_dir, "uniclust", "uniclust30_2018_08")
        cache = os.path.join(tmp_dir, "cache")
        os.makedirs(cache)
        for suf in ("_a3m.ffdata", "_a3m.ffindex"):
            self._write(prefix + suf)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: False)
        assert stage_prefix_files_to_local_ssd_if_remote(prefix, cache) == prefix

    def test_copies_only_matching_files(self, tmp_dir, monkeypatch):
        # Two DBs sharing the same parent dir. Only the prefix's files should copy.
        parent = os.path.join(tmp_dir, "hhdb")
        cache = os.path.join(tmp_dir, "cache")
        os.makedirs(cache)
        prefix_a = os.path.join(parent, "uniclust30_2018_08")
        prefix_b = os.path.join(parent, "pfamA_32.0")
        for suf in ("_a3m.ffdata", "_a3m.ffindex", "_hhm.ffdata", "_hhm.ffindex"):
            self._write(prefix_a + suf)
            self._write(prefix_b + suf)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: True)

        new_a = stage_prefix_files_to_local_ssd_if_remote(prefix_a, cache, min_free_gb=0.0)
        # New prefix points to cache; uniclust files staged, pfam files NOT.
        assert os.path.dirname(new_a).startswith(cache)
        assert os.path.exists(new_a + "_a3m.ffdata")
        assert not os.path.exists(os.path.join(os.path.dirname(new_a), "pfamA_32.0_a3m.ffdata"))

    def test_empty_prefix_returns_unchanged(self, tmp_dir):
        # Optional pfam/pdb70 args default to "" — must be no-op.
        assert stage_prefix_files_to_local_ssd_if_remote("", "/tmp") == ""

    def test_no_matching_files_returns_original(self, tmp_dir, monkeypatch, caplog):
        cache = os.path.join(tmp_dir, "cache")
        os.makedirs(cache)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: True)
        prefix = os.path.join(tmp_dir, "missing", "nonexistent_db")
        # No files at prefix*; helper warns and returns the original.
        assert stage_prefix_files_to_local_ssd_if_remote(prefix, cache) == prefix

    def test_two_prefixes_dont_collide_in_cache(self, tmp_dir, monkeypatch):
        # Each prefix gets its own target subdir under cache/db_prefixes/.
        parent = os.path.join(tmp_dir, "hhdb")
        cache = os.path.join(tmp_dir, "cache")
        os.makedirs(cache)
        prefix_a = os.path.join(parent, "uniclust30")
        prefix_b = os.path.join(parent, "pfam")
        for suf in ("_a3m.ffdata",):
            self._write(prefix_a + suf)
            self._write(prefix_b + suf)
        monkeypatch.setattr("ssign_lib.resources.is_remote_filesystem", lambda p: True)
        new_a = stage_prefix_files_to_local_ssd_if_remote(prefix_a, cache, min_free_gb=0.0)
        new_b = stage_prefix_files_to_local_ssd_if_remote(prefix_b, cache, min_free_gb=0.0)
        assert os.path.dirname(new_a) != os.path.dirname(new_b)


# ---------------------------------------------------------------------------
# auto_batch_size_from_vram — VRAM-tier table for PLM-E
# ---------------------------------------------------------------------------


class TestAutoBatchSizeFromVram:
    """Pin the VRAM → batch size tier mapping so a future tier-table tweak
    can't silently shrink throughput on A40 / A100 / L40S."""

    def _patch_torch(self, monkeypatch, available: bool, total_bytes: int = 0):
        """Stand in for torch.cuda so the test doesn't need a real GPU."""
        import sys
        import types

        fake_torch = types.SimpleNamespace()
        fake_torch.cuda = types.SimpleNamespace(
            is_available=lambda: available,
            get_device_name=lambda i: "FAKE-GPU",
            get_device_properties=lambda i: types.SimpleNamespace(total_memory=total_bytes),
        )
        monkeypatch.setitem(sys.modules, "torch", fake_torch)

    def test_no_cuda_returns_safe_default(self, monkeypatch):
        from ssign_lib.resources import auto_batch_size_from_vram

        self._patch_torch(monkeypatch, available=False)
        assert auto_batch_size_from_vram() == 4

    def test_48gb_a40_picks_batch_64(self, monkeypatch):
        from ssign_lib.resources import auto_batch_size_from_vram

        self._patch_torch(monkeypatch, available=True, total_bytes=48 * 1024**3)
        assert auto_batch_size_from_vram() == 64

    def test_24gb_4090_picks_batch_16(self, monkeypatch):
        from ssign_lib.resources import auto_batch_size_from_vram

        self._patch_torch(monkeypatch, available=True, total_bytes=24 * 1024**3)
        assert auto_batch_size_from_vram() == 16

    def test_12gb_picks_batch_8(self, monkeypatch):
        from ssign_lib.resources import auto_batch_size_from_vram

        self._patch_torch(monkeypatch, available=True, total_bytes=12 * 1024**3)
        assert auto_batch_size_from_vram() == 8

    def test_8gb_picks_batch_4(self, monkeypatch):
        from ssign_lib.resources import auto_batch_size_from_vram

        self._patch_torch(monkeypatch, available=True, total_bytes=8 * 1024**3)
        assert auto_batch_size_from_vram() == 4

    def test_real_a40_reports_44gib_still_picks_batch_64(self, monkeypatch):
        # An A40 advertises 48 GB but torch.cuda.get_device_properties
        # reports ~44.4 GiB (firmware reserve + GB->GiB drop). Real-world
        # threshold has to sit below that or the A40 falls to tier-2.
        from ssign_lib.resources import auto_batch_size_from_vram

        a40_real_bytes = int(44.4 * 1024**3)
        self._patch_torch(monkeypatch, available=True, total_bytes=a40_real_bytes)
        assert auto_batch_size_from_vram() == 64

    def test_real_rtx_4090_reports_22gib_still_picks_batch_16(self, monkeypatch):
        # Same firmware-reserve gap on the 24 GB consumer card (~22 GiB).
        from ssign_lib.resources import auto_batch_size_from_vram

        rtx_4090_real_bytes = int(22 * 1024**3)
        self._patch_torch(monkeypatch, available=True, total_bytes=rtx_4090_real_bytes)
        assert auto_batch_size_from_vram() == 16
