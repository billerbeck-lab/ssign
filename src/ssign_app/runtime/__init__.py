"""Runtime ETA estimation for ssign pipeline runs.

`effort_model` holds the machine-agnostic per-tool effort fit; `estimator`
turns it into a live ETA by inferring this machine's speed from completed steps.
"""

from .effort_model import Effort, effort, limiting_factor, load_coefficients, resolve_regime
from .estimator import Estimator, EtaResult, Step
from .machine import MachineProfile, machine_profile, reference_profile, resource_ratio
from .timeouts import scaled_timeout

__all__ = [
    "Effort",
    "Estimator",
    "EtaResult",
    "MachineProfile",
    "Step",
    "effort",
    "limiting_factor",
    "load_coefficients",
    "machine_profile",
    "reference_profile",
    "resolve_regime",
    "resource_ratio",
    "scaled_timeout",
]
