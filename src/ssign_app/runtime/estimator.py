"""Online ETA estimator: project a run's remaining time, tightening as it goes.

The machine's speed is unknown at run start, so the estimate begins as a wide
prior built from machine-agnostic effort (effort_model) at the reference rate,
then infers this machine's cpu/gpu/io rate from each completed tool and
re-projects the rest.

A *plan* is the ordered list of stages a run will execute. A stage is a list of
Steps that run concurrently (e.g. the DLP/DSE/SignalP predictor block); stages
run serially. A stage's wall-clock is gated by its slowest member, so the ETA
sums the per-stage maxima.

This module is pipeline-agnostic: it takes a plan and completion events. Building
the plan from a live PipelineRunner (or replaying a finished run) is the caller's
job — see calibration/replay.py for the replay path.
"""

from __future__ import annotations

from collections import defaultdict
from typing import NamedTuple

from .effort_model import effort, fixed_range, limiting_factor
from .machine import machine_profile, reference_profile, resource_ratio

# Relative uncertainty of a fit with no measured LOO error (thin whole-genome regimes).
U_THIN = 0.5
# Machine-speed uncertainty before any rate has been inferred for a limiting factor.
PRIOR_MACHINE_U = 0.5
# CI half-width in relative-sigma units (~90%).
Z = 1.64


class Step(NamedTuple):
    tool: str
    regime: str
    size: int


# A stage is a list of Steps that run concurrently; a plan is a list of stages.
Stage = list  # list[Step]
Plan = list  # list[Stage]


class EtaResult(NamedTuple):
    total_s: float  # point estimate of total run wall-clock
    lo_s: float  # total CI lower bound
    hi_s: float  # total CI upper bound
    remaining_s: float  # point estimate of wall-clock left (un-observed stages only)
    remaining_lo_s: float  # remaining CI lower bound
    remaining_hi_s: float  # remaining CI upper bound
    rel_uncertainty: float  # relative half-width of the total band, (hi-lo)/(2*total)
    rates: dict  # inferred {factor: rate} (>1 = faster than reference); empty until a tool completes
    n_unmodeled: int  # planned tools with no fit -> omitted from total_s (so the ETA is a lower bound)


class Estimator:
    """Holds inferred machine rates and projects a plan's ETA from them."""

    def __init__(self, coeffs: dict, prior_rates: dict[str, float] | None = None):
        self.coeffs = coeffs
        # Coarse per-factor rate prior (from resource_ratio) used before any tool
        # of that factor has completed. Overwritten by the inferred rate once one does.
        self._prior_rates: dict[str, float] = dict(prior_rates or {})
        self._rate_samples: dict[str, list[float]] = defaultdict(list)
        self._actual: dict[tuple[int, str], float] = {}

    @classmethod
    def from_profile(cls, coeffs: dict, current=None) -> "Estimator":
        """Build an estimator whose t=0 prior rates come from the current-vs-reference
        resource ratio, so the first ETA is machine-adjusted even with no history.

        ``current`` defaults to the detected `machine_profile()`; pass one to test.
        """
        cur = current if current is not None else machine_profile()
        prior = resource_ratio(reference_profile(coeffs), cur)
        return cls(coeffs, prior_rates=prior)

    def observe(self, stage_idx: int, tool: str, regime: str, size: int, wallclock: float) -> None:
        """Record a completed tool's real wall-clock and update its limiting-factor rate."""
        self._actual[(stage_idx, tool)] = wallclock
        # Don't infer a machine rate from a tool whose measured wall-clock doesn't match
        # its effort model, or the rate is garbage and poisons every tool sharing the factor:
        #   - extract_proteins internally runs Bakta (minutes) but is modelled as a ~second
        #     parse -> effort/wallclock reads as a ~90x-slow CPU (the 661 min post-Bakta spike).
        #   - fixed-range tools are DB/startup-bound, not count- or CPU-bound.
        if tool == "extract_proteins" or fixed_range(tool, regime, self.coeffs) is not None:
            return
        e = effort(tool, size, regime, self.coeffs)
        # Skip negligible/zero-effort tools: effort/wallclock is unstable near zero.
        if e is not None and e.method != "negligible" and e.seconds > 0 and wallclock > 0:
            self._rate_samples[limiting_factor(tool)].append(e.seconds / wallclock)

    def rate(self, factor: str | None) -> float | None:
        """Mean *inferred* rate for a limiting factor (observations only; None until
        a tool of that factor completes). Excludes the prior so callers can see
        exactly what has been learned."""
        samples = self._rate_samples.get(factor)
        return sum(samples) / len(samples) if samples else None

    def _effective_rate(self, factor: str | None) -> tuple[float | None, bool]:
        """(rate, is_observed): the inferred rate if any tool of this factor has
        completed, else the t=0 prior, else (None, False) for a factor with no prior."""
        observed = self.rate(factor)
        if observed is not None:
            return observed, True
        return self._prior_rates.get(factor), False

    def rates(self) -> dict:
        return {f: self.rate(f) for f in self._rate_samples}

    def _project(self, step: Step) -> tuple[float, float, float] | None:
        """(point_s, lo_s, hi_s) for a not-yet-run step, or None if unmodelled.

        Fixed-range annotation tools (eggnog/interproscan/plm_blast) project a raw
        historical p10/p50/p90 wall-clock band instead of a*size+b. Their runtime is
        DB-load/startup-bound, not substrate-count- or CPU-speed-bound (a 1-substrate
        run still took ~28 min), so neither the linear effort NOR a machine-rate
        division applies: the cross-machine percentile spread already IS the band.
        """
        rng = fixed_range(step.tool, step.regime, self.coeffs)
        if rng is not None:
            p10, p50, p90 = rng
            return p50, p10, p90
        e = effort(step.tool, step.size, step.regime, self.coeffs)
        if e is None:
            return None
        r, is_observed = self._effective_rate(limiting_factor(step.tool))
        # Divide reference-machine effort by this machine's rate (any r>0: r>1 faster
        # than reference, r<1 slower). r is None only when no prior AND no observation.
        seconds = e.seconds / r if r else e.seconds
        u_eff = (e.loo_pct / 100.0) if e.loo_pct is not None else U_THIN
        # An inferred rate carries no machine uncertainty; a prior (or no rate) does.
        u_machine = 0.0 if is_observed else PRIOR_MACHINE_U
        u = (u_eff**2 + u_machine**2) ** 0.5
        return seconds, max(0.0, seconds * (1 - Z * u)), seconds * (1 + Z * u)

    @staticmethod
    def _stage_bounds(triples: list[tuple[float, float, float]]) -> tuple[float, float, float]:
        """(point, lo, hi) for a concurrent stage: each bound gated by its slowest member,
        since the stage finishes only when its slowest member does."""
        return (
            max(p for p, _, _ in triples),
            max(lo for _, lo, _ in triples),
            max(hi for _, _, hi in triples),
        )

    def eta(self, plan: Plan) -> EtaResult:
        """Whole-run and remaining ETA: observed wall-clock for completed steps,
        projections for the rest.

        Each stage contributes a (point, lo, hi) triple gated per-bound by its slowest
        member (a concurrent stage finishes when its slowest member does). ``total_*``
        sums every stage; ``remaining_*`` sums only stages with an un-observed step, so
        it never underflows to ~0 while work is still running (the old total-minus-
        elapsed did, because unmodelled steps count in elapsed but not in the total).

        Tools with no fit are skipped (counted in ``n_unmodeled``), so a plan with
        unmodelled tools has a lower-bound ``total_s``. An empty plan returns zeros.
        """
        total = total_lo = total_hi = 0.0
        rem = rem_lo = rem_hi = 0.0
        n_unmodeled = 0
        for stage_idx, stage in enumerate(plan):
            done: list[tuple[float, float, float]] = []  # observed members (point=lo=hi)
            todo: list[tuple[float, float, float]] = []  # projected, not-yet-run members
            for step in stage:
                key = (stage_idx, step.tool)
                if key in self._actual:
                    a = self._actual[key]
                    done.append((a, a, a))
                else:
                    proj = self._project(step)
                    if proj is not None:
                        todo.append(proj)
                    else:
                        n_unmodeled += 1  # no fit for this tool -> ETA is a lower bound
            members = done + todo
            if members:
                p, lo, hi = self._stage_bounds(members)
                total, total_lo, total_hi = total + p, total_lo + lo, total_hi + hi
            if todo:  # a stage with un-run members still owes wall-clock
                p, lo, hi = self._stage_bounds(todo)
                rem, rem_lo, rem_hi = rem + p, rem_lo + lo, rem_hi + hi
        u_rel = ((total_hi - total_lo) / (2 * total)) if total > 0 else 0.0
        return EtaResult(
            total_s=total,
            lo_s=total_lo,
            hi_s=total_hi,
            remaining_s=rem,
            remaining_lo_s=rem_lo,
            remaining_hi_s=rem_hi,
            rel_uncertainty=u_rel,
            rates=self.rates(),
            n_unmodeled=n_unmodeled,
        )
