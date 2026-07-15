"""Unit tests for the runtime ETA estimator (src/ssign_app/runtime/).

Uses a small synthetic coefficients dict so the assertions are deterministic and
independent of the evolving bundled coefficients.json. One test loads the real
bundled artifact to confirm it stays valid.
"""

from ssign_app.runtime import effort_model as em
from ssign_app.runtime.estimator import Estimator, Step

# Synthetic fit: a cpu tool (macsyfinder), a gpu predictor (deeplocpro, both regimes),
# an io tool (eggnog), and a negligible tool (protparam).
COEFFS = {
    "_meta": {"reference_machine": "CX3-A40"},
    "models": {
        "macsyfinder": {
            "fixed": {"a": 0.0, "b": 100.0, "method": "mean", "n": 20, "low_confidence": False, "loo_med_pct": 20.0}
        },
        "deeplocpro": {
            "neighborhood": {
                "a": 1.0,
                "b": 50.0,
                "method": "linear",
                "n": 17,
                "low_confidence": False,
                "loo_med_pct": 10.0,
            },
            "whole_genome": {
                "a": 0.0,
                "b": 600.0,
                "method": "mean",
                "n": 2,
                "low_confidence": True,
                "loo_med_pct": None,
            },
        },
        "signalp": {
            "neighborhood": {
                "a": 0.4,
                "b": 46.0,
                "method": "linear",
                "n": 17,
                "low_confidence": False,
                "loo_med_pct": 12.0,
            }
        },
        "eggnog": {
            "substrates": {
                "a": 10.0,
                "b": 1000.0,
                "method": "linear",
                "n": 24,
                "low_confidence": True,
                "loo_med_pct": 85.0,
            }
        },
        "protparam": {
            "substrates": {
                "a": 0.0,
                "b": 0.5,
                "method": "negligible",
                "n": 8,
                "low_confidence": False,
                "loo_med_pct": None,
            }
        },
    },
}


# --- effort_model ---------------------------------------------------------


def test_effort_linear_and_clamp():
    e = em.effort("deeplocpro", 200, "neighborhood", COEFFS)
    assert e.seconds == 1.0 * 200 + 50.0
    assert e.low_confidence is False and e.method == "linear"
    # b can be negative in a fit; effort must never go below zero
    neg = {"models": {"x": {"substrates": {"a": 1.0, "b": -100.0}}}}
    assert em.effort("x", 10, "substrates", neg).seconds == 0.0


def test_effort_unknown_returns_none():
    assert em.effort("nope", 10, "substrates", COEFFS) is None
    assert em.effort("deeplocpro", 10, "no_such_regime", COEFFS) is None


def test_resolve_regime():
    assert em.resolve_regime("macsyfinder") == "fixed"  # whole-proteome tool
    assert em.resolve_regime("deeplocpro", whole_genome=False) == "neighborhood"
    assert em.resolve_regime("deeplocpro", whole_genome=True) == "whole_genome"
    assert em.resolve_regime("eggnog") == "substrates"
    assert em.resolve_regime("totally_unknown") is None


def test_limiting_factor():
    assert em.limiting_factor("macsyfinder") == "cpu"
    assert em.limiting_factor("deeplocpro") == "gpu"
    assert em.limiting_factor("eggnog") == "io"
    assert em.limiting_factor("unknown") is None


def test_bundled_coefficients_load():
    c = em.load_coefficients()
    assert "models" in c and c["_meta"]["n_clean_points"] > 0
    # every block exposes the fields the estimator reads
    for tool, regimes in c["models"].items():
        for regime, block in regimes.items():
            assert {"a", "b", "method", "low_confidence"} <= block.keys()


# --- estimator ------------------------------------------------------------


def _plan():
    return [
        [Step("macsyfinder", "fixed", 4000)],
        [Step("deeplocpro", "neighborhood", 200)],
        [Step("eggnog", "substrates", 30)],
    ]


def test_prior_is_wide_with_no_rates():
    est = Estimator(COEFFS)
    r = est.eta(_plan())
    assert r.rates == {}
    assert r.rel_uncertainty > 0.3  # machine unknown -> wide
    assert r.lo_s < r.total_s < r.hi_s


def test_cpu_completion_infers_only_cpu_rate():
    est = Estimator(COEFFS)
    est.observe(0, "macsyfinder", "fixed", 4000, 50.0)  # effort 100 / 50 = rate 2.0
    assert est.rate("cpu") == 2.0
    assert est.rate("gpu") is None  # untouched


def test_gpu_completion_tightens_remaining_gpu_tools():
    # Two GPU tools: inferring the rate from the first must narrow the second's projection.
    est = Estimator(COEFFS)
    plan = [[Step("deeplocpro", "neighborhood", 200)], [Step("signalp", "neighborhood", 200)]]
    before = est.eta(plan).rel_uncertainty
    est.observe(0, "deeplocpro", "neighborhood", 200, 125.0)  # effort 250 / 125 = rate 2.0
    after = est.eta(plan)
    assert after.rates.get("gpu") == 2.0
    assert after.rel_uncertainty < before  # signalp (also GPU) now projected with a known rate


def test_negligible_tool_does_not_pollute_rates():
    est = Estimator(COEFFS)
    est.observe(0, "protparam", "substrates", 30, 1.0)  # negligible -> no rate sample
    assert est.rate("cpu") is None


def test_parallel_stage_gated_by_slowest():
    # Two concurrent predictors: stage time should equal the slower one's projection, not the sum.
    coeffs = {
        "models": {
            "deeplocpro": {
                "neighborhood": {
                    "a": 0.0,
                    "b": 100.0,
                    "method": "mean",
                    "n": 9,
                    "low_confidence": False,
                    "loo_med_pct": 10.0,
                }
            },
            "deepsece": {
                "neighborhood": {
                    "a": 0.0,
                    "b": 300.0,
                    "method": "mean",
                    "n": 9,
                    "low_confidence": False,
                    "loo_med_pct": 10.0,
                }
            },
        }
    }
    est = Estimator(coeffs)
    parallel = [[Step("deeplocpro", "neighborhood", 200), Step("deepsece", "neighborhood", 200)]]
    serial = [[Step("deeplocpro", "neighborhood", 200)], [Step("deepsece", "neighborhood", 200)]]
    assert est.eta(parallel).total_s == 300.0  # max(100, 300)
    assert est.eta(serial).total_s == 400.0  # 100 + 300


def test_unmodelled_tool_omitted_not_crash():
    est = Estimator(COEFFS)
    plan = [[Step("frobnicate", "substrates", 10)], [Step("macsyfinder", "fixed", 4000)]]
    r = est.eta(plan)
    assert r.total_s == 100.0  # only the modelled tool contributes
    assert r.n_unmodeled == 1  # the unmodelled tool is surfaced, not silently dropped


def test_cherrypicked_single_tool_subset():
    est = Estimator(COEFFS)
    r = est.eta([[Step("eggnog", "substrates", 30)]])
    assert r.total_s == 10.0 * 30 + 1000.0  # composes from the one tool's model


def test_observed_step_uses_actual_not_projection():
    est = Estimator(COEFFS)
    est.observe(0, "macsyfinder", "fixed", 4000, 999.0)  # absurd actual
    r = est.eta([[Step("macsyfinder", "fixed", 4000)]])
    assert r.total_s == 999.0  # the completed step contributes its real wall-clock


# --- t=0 resource prior ---------------------------------------------------

from ssign_app.runtime.machine import MachineProfile  # noqa: E402

COEFFS_REF = {
    **COEFFS,
    "_meta": {**COEFFS["_meta"], "reference_profile": {"cpu_cores": 24, "gpu": "A40", "gpu_vram_gb": 48}},
}


def test_prior_rate_scales_projection_before_any_observation():
    # A seeded gpu prior of 2.0 halves the reference-machine projection at t=0,
    # but is NOT reported as an inferred rate and keeps the wide machine CI.
    est = Estimator(COEFFS, prior_rates={"gpu": 2.0})
    r = est.eta([[Step("deeplocpro", "neighborhood", 200)]])  # effort 250
    assert r.total_s == 125.0  # 250 / 2.0
    assert r.rates == {}  # prior is not a learned rate
    assert r.rel_uncertainty > 0.3  # still wide until a real observation lands


def test_observation_overrides_prior():
    est = Estimator(COEFFS, prior_rates={"gpu": 2.0})
    plan = [[Step("deeplocpro", "neighborhood", 200)], [Step("signalp", "neighborhood", 200)]]
    est.observe(0, "deeplocpro", "neighborhood", 200, 62.5)  # effort 250 / 62.5 = rate 4.0
    r = est.eta(plan)
    # signalp (also gpu) now projects with the inferred 4.0, not the 2.0 prior.
    sp_effort = 0.4 * 200 + 46.0  # 126
    assert r.rates["gpu"] == 4.0
    # stage 0 observed (62.5) + stage 1 projected at rate 4.0
    assert r.total_s == 62.5 + sp_effort / 4.0


def test_from_profile_slow_machine_widens_prior_eta():
    # 12-core, no-GPU laptop vs the 24-core A40 reference: cpu prior 0.707, gpu 0.02.
    est = Estimator.from_profile(COEFFS_REF, current=MachineProfile(12, None, None))
    cpu_eta = est.eta([[Step("macsyfinder", "fixed", 4000)]]).total_s  # effort 100
    assert cpu_eta == 100.0 / (12 / 24) ** 0.5  # slower than the reference 100 s
    gpu_eta = est.eta([[Step("deeplocpro", "neighborhood", 200)]]).total_s  # effort 250
    assert gpu_eta == 250.0 / 0.02  # CPU-fallback penalty -> far slower


def test_from_profile_no_reference_profile_is_parity():
    # COEFFS lacks _meta.reference_profile -> prior is parity, same as bare Estimator.
    est = Estimator.from_profile(COEFFS, current=MachineProfile(12, None, None))
    assert est.eta([[Step("macsyfinder", "fixed", 4000)]]).total_s == 100.0
