"""Tests for the GPU-serialization guard.

Bug: DeepLocPro and DeepSecE run as a parallel group, so on an exclusive-process
GPU (the HPC norm) the second process to reach `model.to(cuda)` dies with
`cudaErrorDevicesUnavailable`. Fix: `cuda_compute_mode()` detects the mode and the
runner serialises GPU-bound steps behind a single lock only when the device can't
be shared. These tests pin the mode parsing and the serialise/parallel decision.
"""

from __future__ import annotations

from contextlib import nullcontext
from unittest import mock

from ssign_lib import resources

from ssign_app.core.runner import PipelineConfig, PipelineRunner

_MODE_FN = "ssign_app.scripts.ssign_lib.resources.cuda_compute_mode"


def _fake_run(stdout: str = "", returncode: int = 0):
    return mock.Mock(stdout=stdout, returncode=returncode)


class TestCudaComputeMode:
    def test_default_mode(self, monkeypatch):
        monkeypatch.setattr(resources.subprocess, "run", lambda *a, **k: _fake_run("Default\n"))
        assert resources.cuda_compute_mode() == "Default"

    def test_exclusive_mode(self, monkeypatch):
        monkeypatch.setattr(resources.subprocess, "run", lambda *a, **k: _fake_run("Exclusive_Process\n"))
        assert resources.cuda_compute_mode() == "Exclusive_Process"

    def test_nonzero_exit_is_none(self, monkeypatch):
        monkeypatch.setattr(resources.subprocess, "run", lambda *a, **k: _fake_run("err", returncode=1))
        assert resources.cuda_compute_mode() is None

    def test_missing_binary_is_none(self, monkeypatch):
        def _raise(*a, **k):
            raise FileNotFoundError("nvidia-smi")

        monkeypatch.setattr(resources.subprocess, "run", _raise)
        assert resources.cuda_compute_mode() is None


class TestGpuSerializeDecision:
    def _runner(self, tmp_dir):
        return PipelineRunner(PipelineConfig(outdir=tmp_dir, sample_id="t"))

    def test_no_gpu_stays_parallel(self, tmp_dir):
        r = self._runner(tmp_dir)
        r._gpu_present_cache = False
        assert r._gpu_needs_serial() is False
        assert isinstance(r._gpu_serialize(), nullcontext)

    def test_default_mode_stays_parallel(self, tmp_dir, monkeypatch):
        r = self._runner(tmp_dir)
        r._gpu_present_cache = True
        monkeypatch.setattr(_MODE_FN, lambda: "Default")
        assert r._gpu_needs_serial() is False
        assert isinstance(r._gpu_serialize(), nullcontext)

    def test_exclusive_mode_serializes_on_the_shared_lock(self, tmp_dir, monkeypatch):
        r = self._runner(tmp_dir)
        r._gpu_present_cache = True
        monkeypatch.setattr(_MODE_FN, lambda: "Exclusive_Process")
        assert r._gpu_needs_serial() is True
        # the returned context manager IS the one shared lock, so DLP/DSE contend
        assert r._gpu_serialize() is r._gpu_lock

    def test_unreadable_mode_serializes_safely(self, tmp_dir, monkeypatch):
        r = self._runner(tmp_dir)
        r._gpu_present_cache = True
        monkeypatch.setattr(_MODE_FN, lambda: None)
        assert r._gpu_needs_serial() is True
