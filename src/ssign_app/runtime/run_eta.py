"""Translate a live PipelineRunner's stage list into an estimator ``Plan``.

Pure, runner-free translation so it is unit-testable on its own: given the
ordered per-stage step-ids, the sizes known so far, and whether the predictors
run whole-genome, produce the list-of-concurrent-stages ``Plan`` the
``Estimator`` consumes. Steps that are unmodelled (no effort fit) or whose size
is not yet known are dropped from their stage; the stage itself is kept (even if
it becomes empty) so a plan stage index always equals the runner's stage index,
which is what ``Estimator.observe(stage_idx, ...)`` relies on.
"""

from __future__ import annotations

from typing import NamedTuple

from .effort_model import resolve_regime
from .estimator import Estimator, Step

# Runner step_id -> effort-model tool name (identity except where they differ:
# the runner's `_step_hhsuite` is the model's `hh_suite`). Steps absent from this
# map (detect_format, cross_validate, proximity, filtering, integrate, ...) are
# unmodelled and contribute nothing to the ETA.
TOOL_BY_STEP: dict[str, str] = {
    "extract_proteins": "extract_proteins",
    "macsyfinder": "macsyfinder",
    "deeplocpro": "deeplocpro",
    "deepsece": "deepsece",
    "signalp": "signalp",
    "blastp": "blastp",
    "hhsuite": "hh_suite",
    "interproscan": "interproscan",
    "protparam": "protparam",
    "eggnog": "eggnog",
    "plm_blast": "plm_blast",
}
_PREDICTORS = frozenset({"deeplocpro", "deepsece", "signalp"})


def step_to_plan(step_id: str, sizes: dict[str, int], whole_genome_tools: frozenset[str] | set[str]) -> Step | None:
    """(tool, regime, size) Step for one runner step, or None if unmodelled / size unknown.

    ``sizes`` carries what is known so far: ``proteins`` (after input parsing),
    ``neighborhood`` (after neighborhood extraction), ``substrates`` (after
    filtering). ``whole_genome_tools`` is the set of predictor tool names running
    whole-genome (the three predictors have independent ``*_whole_genome`` flags);
    it only affects the predictor regime.
    """
    tool = TOOL_BY_STEP.get(step_id)
    if tool is None:
        return None
    regime = resolve_regime(tool, whole_genome=tool in whole_genome_tools)
    if regime == "fixed" or regime == "whole_genome":
        size = sizes.get("proteins")
    elif regime == "neighborhood":
        # Neighborhood count once known; the full proteome is a safe upper-bound proxy before then.
        size = sizes.get("neighborhood", sizes.get("proteins"))
    elif regime == "substrates":
        size = sizes.get("substrates")
    else:
        size = None
    if size is None:
        return None
    return Step(tool, regime, int(size))


def build_plan(
    stage_step_ids: list[list[str]], sizes: dict[str, int], whole_genome_tools: frozenset[str] | set[str]
) -> list[list[Step]]:
    """Map the runner's ordered stages (each a list of step-ids) to an estimator Plan.

    One plan stage per runner stage (empty stages preserved) so plan stage
    indices line up with the runner's stage iteration order.
    """
    return [
        [s for s in (step_to_plan(sid, sizes, whole_genome_tools) for sid in stage) if s is not None]
        for stage in stage_step_ids
    ]


class ReplayStep(NamedTuple):
    step: int  # 1-based completion index
    tool: str
    total_s: float  # projected TOTAL run wall-clock after this completion
    rel_uncertainty: float
    error_pct: float  # |projected total - true total| / true total, in %


class ReplayResult(NamedTuple):
    true_total_s: float  # the run's actual total (all steps observed)
    steps: list[ReplayStep]


def replay(coeffs: dict, plan: list[list[Step]], completions: list[tuple]) -> ReplayResult:
    """Offline convergence harness: feed a finished run's per-tool wall-clocks
    through a fresh estimator in completion order, recording the projected total
    (and its error vs the true total) after each step, using only data available
    up to that step.

    ``completions`` is an ordered list of ``(stage_idx, tool, regime, size,
    wallclock_s)`` — the order the tools finished. Starts from a bare estimator
    (reference-rate prior, no machine detection) so the result is deterministic
    and isolates the online rate inference. The true total is the estimator's ETA
    once every step is observed (each stage contributes the max of its members'
    real wall-clocks).
    """
    truth = Estimator(coeffs)
    for stage_idx, tool, regime, size, wc in completions:
        truth.observe(stage_idx, tool, regime, size, wc)
    true_total = truth.eta(plan).total_s

    est = Estimator(coeffs)
    steps: list[ReplayStep] = []
    for k, (stage_idx, tool, regime, size, wc) in enumerate(completions, start=1):
        est.observe(stage_idx, tool, regime, size, wc)
        res = est.eta(plan)
        err = abs(res.total_s - true_total) / true_total * 100 if true_total > 0 else 0.0
        steps.append(ReplayStep(k, tool, res.total_s, res.rel_uncertainty, err))
    return ReplayResult(true_total_s=true_total, steps=steps)
