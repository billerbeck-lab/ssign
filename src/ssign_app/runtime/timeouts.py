"""Size-aware per-tool subprocess timeouts derived from the effort model.

A flat cap doesn't scale with input size: the whole-genome prediction pool of a
52-genome batch is ~160k proteins and needs ~5-6h, well past the flat 4h
`TOOL_TIMEOUT_S`. `scaled_timeout` sizes each tool's cap from the fitted effort
model (`effort_model` / `coefficients.json`, reference-machine seconds) times a
margin, floored at `TOOL_TIMEOUT_S` so small/neighborhood runs are unchanged.
Unmodelled tools (blastp, hh_suite, orthologs-blast) return the floor untouched.

The margin absorbs machine variance (reference is CX3-A40; a run GPU within ~1.5x
is covered by 2x). Low-confidence fits get a wider margin so a thin/noisy fit
can't produce a cap that kills a still-running tool; the floor is the backstop
for the degenerate case where the fit predicts far too low.
"""

from __future__ import annotations

import math

from ssign_app.scripts.ssign_lib.constants import TOOL_TIMEOUT_S

from .effort_model import effort, load_coefficients

# Margin on the reference-machine effort estimate. 2x for well-calibrated fits
# (DLP/DSE whole_genome are 6.5-7% LOO error); wider for low-confidence fits
# (thin/noisy, e.g. eggnog LOO 137%) so an under-estimate can't false-kill.
DEFAULT_MARGIN = 2.0
LOW_CONFIDENCE_MARGIN = 3.0

_coeffs_cache: dict | None = None


def _coefficients() -> dict:
    """Load and cache coefficients.json once; degrade to floor-only if unreadable."""
    global _coeffs_cache
    if _coeffs_cache is None:
        try:
            _coeffs_cache = load_coefficients()
        except (OSError, ValueError):
            # Missing/corrupt coefficients must degrade to the floor, not crash a run.
            _coeffs_cache = {}
    return _coeffs_cache


def scaled_timeout(
    tool: str,
    size: int,
    regime: str,
    *,
    margin: float = DEFAULT_MARGIN,
    low_confidence_margin: float = LOW_CONFIDENCE_MARGIN,
    floor: int = TOOL_TIMEOUT_S,
    machine_rate: float | None = None,
) -> int:
    """Subprocess timeout (seconds) for `tool` processing `size` items in `regime`.

    Returns `max(floor, ceil(m * predicted_seconds))`, where `predicted_seconds`
    is the effort model's reference-machine estimate and `m` is `margin` for a
    well-calibrated fit or `low_confidence_margin` for a low-confidence one.
    Returns `floor` unchanged when the tool/regime has no fit.

    `machine_rate` is the estimator's inferred rate for this tool's limiting
    factor (>1 = this machine beats the reference, <1 = slower). When given, the
    reference-machine estimate is divided by it before the margin, so a slow
    machine gets a proportionally longer kill-timeout. Omitted (None) leaves the
    raw reference estimate untouched, preserving the pre-estimator behaviour.
    """
    e = effort(tool, size, regime, _coefficients())
    if e is None:
        return int(floor)
    predicted = e.seconds / machine_rate if machine_rate and machine_rate > 0 else e.seconds
    m = low_confidence_margin if e.low_confidence else margin
    return int(max(floor, math.ceil(m * predicted)))
