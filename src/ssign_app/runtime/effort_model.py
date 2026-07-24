"""Per-tool effort model: machine-agnostic runtime in reference-machine seconds.

`effort(tool, size, regime)` returns how long a tool *would* take on the
reference machine (CX3-A40, rate 1.0). The online estimator divides this by the
inferred per-machine rate to get a real ETA (see `estimator.py`).

Coefficients are fit offline from cleaned calibration data and bundled as
`coefficients.json` (regenerate with `calibration/fit.py`). This module only
*reads* them; it has no dependency on the calibration data or numpy.
"""

from __future__ import annotations

import json
import os
from typing import NamedTuple

_COEFFS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "coefficients.json")

# Which physical resource gates each tool. The online estimator infers one rate
# per factor (cpu_rate / gpu_rate / io_factor) and applies it to every tool that
# shares the factor. Mirrors the README "per-tool limiting factors" table; keep
# the two in sync. `calibration/fit.py` mirrors the gpu set for its CPU-run filter.
LIMITING_FACTOR: dict[str, str] = {
    "bakta": "cpu",
    "macsyfinder": "cpu",
    "interproscan": "cpu",
    "protparam": "cpu",
    "extract_proteins": "cpu",
    "blastp": "cpu",
    "deeplocpro": "gpu",
    "signalp": "gpu",
    "deepsece": "gpu",
    "plm_blast_embed": "gpu",
    "plm_blast": "gpu",  # ProtT5 embedding dominates; ~100x slower on CPU
    "plm_blast_search": "cpu",
    "hh_suite": "cpu",
    "eggnog": "io",  # 39 GB DB: NFS-IO bound, or CPU bound with --dbmem
}

# Tool -> size axis / regime class. A tool's regime selects which coefficient
# block applies (see resolve_regime). Mirrors calibration/clean.py's sets.
_WHOLE = {"bakta", "macsyfinder", "extract_proteins"}  # always whole proteome -> "fixed"
_PREDICT = {"deeplocpro", "signalp", "deepsece"}  # neighborhood OR whole_genome
_ANNOT = {
    "eggnog",
    "interproscan",
    "plm_blast",
    "protparam",
    "plm_blast_embed",
    "plm_blast_search",
    "hh_suite",
    "blastp",  # BLASTs the substrate/passenger set, not the whole proteome
}  # run on the substrate set -> "substrates"

# Annotation tools whose wall-clock is driven by DB cache state / domain complexity,
# NOT by substrate count: a linear a*size+b effort is misleading for them (it produced
# the ~1693 min pre-run over-estimate when the substrate count was still a rough prior).
# For these, the estimator uses a fixed historical p10/p50/p90 wallclock band as-is (no
# size, no machine-rate division: they are DB-load/startup-bound, not CPU-bound) instead
# of a*size+b. eggnog/interproscan/plm_blast are the three with enough real calibration points to fit
# a percentile band; blastp/hh_suite are hand-set priors with no data (they stay on the
# effort path). See coefficients.json _meta and calibration fit.py.
FIXED_RANGE_TOOLS = frozenset({"eggnog", "interproscan", "plm_blast"})


class Effort(NamedTuple):
    seconds: float  # reference-machine seconds
    low_confidence: bool  # fit was thin / high-error -> widen the CI
    method: str  # "linear" | "mean" | "negligible"
    n: int  # calibration points behind the fit
    loo_pct: float | None  # leave-one-out median % error of the fit (None if unmeasured)


def limiting_factor(tool: str) -> str | None:
    """'cpu' | 'gpu' | 'io', or None if the tool isn't modelled."""
    return LIMITING_FACTOR.get(tool)


def resolve_regime(tool: str, whole_genome: bool = False) -> str | None:
    """The coefficient block to use for this tool given how it will run.

    Predictors run on the +/-3 neighborhood by default, or the whole proteome
    when their --*-whole-genome flag is set. Whole-proteome tools are "fixed";
    annotation tools are always "substrates".
    """
    if tool in _WHOLE:
        return "fixed"
    if tool in _PREDICT:
        return "whole_genome" if whole_genome else "neighborhood"
    if tool in _ANNOT:
        return "substrates"
    return None


def load_coefficients(path: str | None = None) -> dict:
    """Load the bundled (or given) coefficients.json."""
    with open(path or _COEFFS_PATH) as f:
        return json.load(f)


def effort(tool: str, size: int, regime: str, coeffs: dict) -> Effort | None:
    """Reference-machine seconds for `tool` processing `size` items in `regime`.

    Returns None when no fit exists for (tool, regime) — the caller decides how
    to handle an unmodelled tool (e.g. omit from the ETA, or use a wide prior).
    """
    block = coeffs.get("models", {}).get(tool, {}).get(regime)
    if block is None:
        return None
    seconds = max(block["a"] * size + block["b"], 0.0)
    return Effort(
        seconds=seconds,
        # Conservative default: a block missing the flag (e.g. an older coefficients.json)
        # is treated as low-confidence so a stale fit widens the CI rather than overstating it.
        low_confidence=bool(block.get("low_confidence", True)),
        method=block.get("method", "mean"),
        n=int(block.get("n", 0)),
        loo_pct=block.get("loo_med_pct"),
    )


def fixed_range(tool: str, regime: str, coeffs: dict) -> tuple[float, float, float] | None:
    """Historical (p10, p50, p90) wall-clock seconds for a fixed-range annotation tool,
    or None if the tool isn't fixed-range or the block carries no ``range_s``.

    Pooled RAW seconds across calibration machines and cache states. The estimator uses
    them as-is (no size, no machine-rate division): these tools are DB-load/startup-bound,
    so the cross-machine spread already captures their variance.
    """
    if tool not in FIXED_RANGE_TOOLS:
        return None
    block = coeffs.get("models", {}).get(tool, {}).get(regime)
    if not block:
        return None
    rng = block.get("range_s")
    if not rng:
        return None
    try:
        return float(rng["p10"]), float(rng["p50"]), float(rng["p90"])
    except (KeyError, TypeError, ValueError):
        return None
