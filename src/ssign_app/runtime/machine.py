"""Current-vs-reference machine profile + a coarse resource-ratio prior.

The estimator needs a first ETA *before any tool has completed*, on a machine it
has never seen. That prior comes from a coarse ratio of the current machine's
resources to the calibration reference machine's (encoded in ``coefficients.json``
``_meta.reference_profile``): more CPU cores -> CPU-bound tools run faster; a GPU
present (or absent) -> GPU-bound tools much faster (or ~50x slower on the CPU
fallback). The ratio is deliberately rough. The online per-factor rate inference
(``estimator.py``) supersedes it within the first couple of completed tools, so
the prior only has to get the run start into the right ballpark.
"""

from __future__ import annotations

from typing import NamedTuple

# CPU-bound tools rarely parallelise perfectly, so wall-clock does not fall
# linearly with core count. A square-root law is a conservative diminishing-
# returns default: 4x the cores -> ~2x the rate.
_CPU_DIMINISHING_EXP = 0.5
# A GPU-bound tool (ESM / ProtT5 transformer) forced onto CPU runs ~50-100x
# slower. One conservative number for the prior; the online rate correction
# replaces it as soon as a GPU tool completes, and the CI stays wide until then.
_NO_GPU_PENALTY = 0.02


class MachineProfile(NamedTuple):
    cpu_cores: int | None
    gpu: str | None  # coarse CUDA device name, or None when no device is visible
    gpu_vram_gb: float | None


def machine_profile() -> MachineProfile:
    """Detect the current machine's effective cores + GPU, reusing doctor's probes."""
    from ssign_app.scripts.ssign_lib.resources import effective_cpu_count, probe_cuda_device

    name, vram = probe_cuda_device()
    return MachineProfile(cpu_cores=effective_cpu_count(), gpu=name, gpu_vram_gb=vram)


def reference_profile(coeffs: dict) -> MachineProfile | None:
    """The calibration reference machine's profile from ``_meta.reference_profile``.

    Returns None when the coefficients predate this field, so the caller falls
    back to a parity prior (ratio 1.0) rather than crashing.
    """
    rp = coeffs.get("_meta", {}).get("reference_profile")
    if not rp:
        return None
    return MachineProfile(
        cpu_cores=rp.get("cpu_cores"),
        gpu=rp.get("gpu"),
        gpu_vram_gb=rp.get("gpu_vram_gb"),
    )


def resource_ratio(reference: MachineProfile | None, current: MachineProfile) -> dict[str, float]:
    """Per-limiting-factor prior rate multiplier: >1 = current beats the reference.

    Keys match ``limiting_factor()``: ``cpu`` / ``gpu`` / ``io``. Used only to
    seed the estimator's rates at t=0; each factor is overwritten by the inferred
    rate once a tool of that class finishes. Missing data degrades to 1.0 (parity)
    rather than a wrong guess.
    """
    if reference is None:
        return {"cpu": 1.0, "gpu": 1.0, "io": 1.0}

    if reference.cpu_cores and current.cpu_cores:
        cpu = (current.cpu_cores / reference.cpu_cores) ** _CPU_DIMINISHING_EXP
    else:
        cpu = 1.0

    ref_has_gpu = reference.gpu is not None
    cur_has_gpu = current.gpu is not None
    if ref_has_gpu and not cur_has_gpu:
        gpu = _NO_GPU_PENALTY  # reference measured on GPU, we're on the CPU fallback
    elif cur_has_gpu and not ref_has_gpu:
        gpu = 1.0 / _NO_GPU_PENALTY
    else:
        gpu = 1.0  # both GPU (coarse: model FLOPs unknown) or both CPU

    # No cheap proxy for storage speed (the io factor); parity prior, learned online.
    return {"cpu": cpu, "gpu": gpu, "io": 1.0}
