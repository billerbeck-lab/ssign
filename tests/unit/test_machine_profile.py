"""Unit tests for the reference/current machine profile + resource-ratio prior.

The ratio only seeds the estimator's t=0 rates, so the tests pin the coarse
behaviour (more cores -> faster; no GPU -> heavy penalty; missing data -> parity)
rather than exact multipliers beyond the documented formulas.
"""

from ssign_app.runtime.effort_model import load_coefficients
from ssign_app.runtime.machine import (
    _CPU_DIMINISHING_EXP,
    _NO_GPU_PENALTY,
    MachineProfile,
    reference_profile,
    resource_ratio,
)

REF = MachineProfile(cpu_cores=24, gpu="A40", gpu_vram_gb=48)


def test_reference_profile_read_from_bundled_coeffs():
    prof = reference_profile(load_coefficients())
    assert prof is not None
    assert prof.cpu_cores == 24 and prof.gpu == "A40"


def test_reference_profile_absent_returns_none():
    assert reference_profile({"_meta": {}}) is None
    assert reference_profile({}) is None


def test_more_cores_gives_cpu_rate_above_one():
    # 96-core node vs the 24-core reference: faster, but sub-linear.
    r = resource_ratio(REF, MachineProfile(96, "A40", 48))
    assert r["cpu"] == (96 / 24) ** _CPU_DIMINISHING_EXP
    assert 1.0 < r["cpu"] < (96 / 24)  # diminishing returns, not linear


def test_fewer_cores_gives_cpu_rate_below_one():
    r = resource_ratio(REF, MachineProfile(6, "A40", 48))
    assert r["cpu"] < 1.0


def test_no_gpu_applies_cpu_fallback_penalty():
    # ThinkPad: 12 effective cores, no CUDA device.
    r = resource_ratio(REF, MachineProfile(12, None, None))
    assert r["gpu"] == _NO_GPU_PENALTY
    assert r["gpu"] < 0.1  # GPU tools ~50x slower on CPU


def test_gpu_present_both_sides_is_parity():
    r = resource_ratio(REF, MachineProfile(24, "RTX6000", 24))
    assert r["gpu"] == 1.0  # coarse: model FLOPs unknown, learned online


def test_io_prior_is_always_parity():
    assert resource_ratio(REF, MachineProfile(96, "A40", 48))["io"] == 1.0
    assert resource_ratio(None, MachineProfile(96, "A40", 48))["io"] == 1.0


def test_missing_reference_degrades_to_parity():
    r = resource_ratio(None, MachineProfile(96, None, None))
    assert r == {"cpu": 1.0, "gpu": 1.0, "io": 1.0}


def test_missing_core_counts_degrade_cpu_to_parity():
    r = resource_ratio(MachineProfile(None, "A40", 48), MachineProfile(None, "A40", 48))
    assert r["cpu"] == 1.0
