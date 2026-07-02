"""Runtime ETA estimation for ssign pipeline runs.

`effort_model` holds the machine-agnostic per-tool effort fit; `estimator`
turns it into a live ETA by inferring this machine's speed from completed steps.
"""

from .effort_model import Effort, effort, limiting_factor, load_coefficients, resolve_regime
from .timeouts import scaled_timeout

__all__ = [
    "Effort",
    "effort",
    "limiting_factor",
    "load_coefficients",
    "resolve_regime",
    "scaled_timeout",
]
