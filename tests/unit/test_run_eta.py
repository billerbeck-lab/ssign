"""Tests for the runner<->estimator glue: the pure plan builder (run_eta) and
the PipelineRunner's defensive ETA hooks (_eta_sizes / _eta_step).

The ETA is informational, so the key guarantees are: it maps steps to the right
effort (tool, regime, size), it refines the substrate size once filtering runs,
and a failure inside it never propagates out of a run.
"""

from ssign_app.core.runner import PipelineConfig, PipelineRunner
from ssign_app.runtime.effort_model import load_coefficients
from ssign_app.runtime.estimator import Estimator, Step
from ssign_app.runtime.run_eta import build_plan, replay, step_to_plan

# --- pure plan builder ----------------------------------------------------


def test_step_to_plan_maps_predictor_to_neighborhood():
    s = step_to_plan("deeplocpro", {"proteins": 4000, "neighborhood": 120}, set())
    assert (s.tool, s.regime, s.size) == ("deeplocpro", "neighborhood", 120)


def test_step_to_plan_whole_genome_uses_proteins():
    s = step_to_plan("signalp", {"proteins": 4000, "neighborhood": 120}, {"signalp"})
    assert (s.tool, s.regime, s.size) == ("signalp", "whole_genome", 4000)


def test_step_to_plan_neighborhood_falls_back_to_proteins():
    # Before neighborhood extraction, the full proteome is the upper-bound proxy.
    s = step_to_plan("deepsece", {"proteins": 4000}, set())
    assert s.regime == "neighborhood" and s.size == 4000


def test_step_to_plan_whole_genome_is_per_predictor():
    # The three predictors have independent whole-genome flags: with only DSE
    # whole-genome, DSE is sized on the full proteome while DLP stays neighborhood.
    sizes = {"proteins": 4000, "neighborhood": 120}
    wg = {"deepsece"}
    assert step_to_plan("deepsece", sizes, wg).regime == "whole_genome"
    assert step_to_plan("deeplocpro", sizes, wg).regime == "neighborhood"


def test_step_to_plan_hhsuite_step_id_maps_to_hh_suite_tool():
    s = step_to_plan("hhsuite", {"substrates": 50}, set())
    assert s.tool == "hh_suite" and s.regime == "substrates" and s.size == 50


def test_step_to_plan_unmodelled_step_is_none():
    assert step_to_plan("cross_validate", {"proteins": 4000}, set()) is None
    assert step_to_plan("integrate", {"proteins": 4000}, set()) is None


def test_step_to_plan_unknown_size_is_none():
    # eggnog needs a substrate count; without it the step drops out of the plan.
    assert step_to_plan("eggnog", {"proteins": 4000}, set()) is None


def test_build_plan_preserves_stage_indices():
    stages = [["extract_proteins"], ["deeplocpro", "deepsece", "signalp"], ["cross_validate"], ["eggnog"]]
    plan = build_plan(stages, {"proteins": 4000, "substrates": 40}, set())
    assert len(plan) == 4  # one plan stage per runner stage, even the empty one
    assert [s.tool for s in plan[1]] == ["deeplocpro", "deepsece", "signalp"]
    assert plan[2] == []  # cross_validate is unmodelled -> empty stage, index kept
    assert plan[3][0].tool == "eggnog"


# --- runner ETA hooks -----------------------------------------------------


def _write_fasta(path, n):
    with open(path, "w") as f:
        for i in range(n):
            f.write(f">rec{i}\nMKT\n")
    return str(path)


def test_eta_sizes_uses_substrate_prior_before_filtering(tmp_path):
    r = PipelineRunner(PipelineConfig(outdir=str(tmp_path), sample_id="x"))
    r.files["proteins"] = _write_fasta(tmp_path / "p.faa", 1000)
    r.files["neighborhood_proteins"] = _write_fasta(tmp_path / "n.faa", 20)
    sizes = r._eta_sizes()
    assert sizes["proteins"] == 1000 and sizes["neighborhood"] == 20
    assert sizes["substrates"] == 6  # round(0.006 * 1000), the pre-filtering prior


def test_eta_sizes_uses_measured_substrates_after_filtering(tmp_path, monkeypatch):
    r = PipelineRunner(PipelineConfig(outdir=str(tmp_path), sample_id="x"))
    r.files["proteins"] = _write_fasta(tmp_path / "p.faa", 100)
    r.files["substrates_filtered"] = "whatever.tsv"  # presence flips to the measured path
    monkeypatch.setattr(r, "_substrate_count", lambda: 7)
    assert r._eta_sizes()["substrates"] == 7  # measured count, not the fraction prior


def test_eta_step_observes_and_returns_line(tmp_path):
    r = PipelineRunner(PipelineConfig(outdir=str(tmp_path), sample_id="x"))
    r.files["proteins"] = _write_fasta(tmp_path / "p.faa", 4000)
    r._estimator = Estimator(load_coefficients())
    r._eta_stage_ids = [["deeplocpro", "deepsece", "signalp"], ["eggnog"]]
    msg = r._eta_step(0, "deeplocpro", 120.0)  # a real GPU completion
    # The gpu rate is now inferred from the observation...
    assert r._estimator.rate("gpu") is not None
    # ...and a one-line ETA was returned for the caller to print last.
    assert msg is not None and "estimated" in msg


def test_eta_step_returns_none_when_estimator_absent(tmp_path):
    r = PipelineRunner(PipelineConfig(outdir=str(tmp_path), sample_id="x"))
    r._estimator = None
    assert r._eta_step(0, "deeplocpro", 120.0) is None  # no-op, no crash


def test_eta_step_never_raises_on_bad_state(tmp_path):
    # No proteins file, garbage stage ids: the hook must swallow everything.
    r = PipelineRunner(PipelineConfig(outdir=str(tmp_path), sample_id="x"))
    r._estimator = Estimator(load_coefficients())
    r._eta_stage_ids = [["not_a_real_step"]]
    r._eta_step(0, "not_a_real_step", 5.0)  # should not raise


# --- replay convergence harness (§5) --------------------------------------

# Synthetic fit: one cpu tool, a concurrent gpu predictor block, one io tool.
_REPLAY_COEFFS = {
    "models": {
        "macsyfinder": {
            "fixed": {"a": 0.0, "b": 200.0, "method": "mean", "n": 9, "low_confidence": False, "loo_med_pct": 10.0}
        },
        "deeplocpro": {
            "neighborhood": {
                "a": 1.0,
                "b": 0.0,
                "method": "linear",
                "n": 9,
                "low_confidence": False,
                "loo_med_pct": 10.0,
            }
        },
        "deepsece": {
            "neighborhood": {
                "a": 2.0,
                "b": 0.0,
                "method": "linear",
                "n": 9,
                "low_confidence": False,
                "loo_med_pct": 10.0,
            }
        },
        "eggnog": {
            "substrates": {
                "a": 10.0,
                "b": 100.0,
                "method": "linear",
                "n": 9,
                "low_confidence": True,
                "loo_med_pct": 40.0,
            }
        },
    }
}
_REPLAY_PLAN = [
    [Step("macsyfinder", "fixed", 4000)],
    [Step("deeplocpro", "neighborhood", 200), Step("deepsece", "neighborhood", 200)],
    [Step("eggnog", "substrates", 30)],
]
# This machine runs at 0.5x the reference (every tool takes 2x its reference effort).
_REPLAY_COMPLETIONS = [
    (0, "macsyfinder", "fixed", 4000, 400.0),  # effort 200 -> rate 0.5 (cpu)
    (1, "deeplocpro", "neighborhood", 200, 400.0),  # effort 200 -> rate 0.5 (gpu)
    (1, "deepsece", "neighborhood", 200, 800.0),  # effort 400 -> rate 0.5 (gpu)
    (2, "eggnog", "substrates", 30, 800.0),  # effort 400 -> rate 0.5 (io)
]


def test_replay_converges_to_true_total():
    res = replay(_REPLAY_COEFFS, _REPLAY_PLAN, _REPLAY_COMPLETIONS)
    # True total on a 0.5x machine: 400 (cpu) + max(400,800) (gpu) + 800 (io) = 2000.
    assert res.true_total_s == 2000.0
    # Final step observed everything -> exact.
    assert res.steps[-1].error_pct == 0.0
    # Convergence: the first estimate (reference-rate prior for the unobserved
    # tools) is further off than the last.
    assert res.steps[0].error_pct > res.steps[-1].error_pct
    # The CI tightens as classes get observed (last step has no unknown machine rate).
    assert res.steps[-1].rel_uncertainty <= res.steps[0].rel_uncertainty


def test_replay_error_is_zero_when_model_is_perfect():
    # If the machine runs exactly at reference rate, even the first projection
    # nails the total (no machine correction needed).
    perfect = [
        (0, "macsyfinder", "fixed", 4000, 200.0),
        (1, "deeplocpro", "neighborhood", 200, 200.0),
        (1, "deepsece", "neighborhood", 200, 400.0),
        (2, "eggnog", "substrates", 30, 400.0),
    ]
    res = replay(_REPLAY_COEFFS, _REPLAY_PLAN, perfect)
    assert res.true_total_s == 1000.0  # 200 + max(200,400) + 400
    assert all(s.error_pct == 0.0 for s in res.steps)
