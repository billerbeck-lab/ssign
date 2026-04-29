"""Unit tests for cross_validate_predictions.py.

Exercises the pure-Python `cross_validate()` generator directly with
in-memory dicts — no filesystem writes, no subprocess. Covers each
trigger tool independently (DLP, DSE, PLM-Effector), the SignalP
evidence-only rule, the DSE T3SS flagging guard, and the
`n_prediction_tools_agreeing` count.
"""

import os
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from cross_validate_predictions import cross_validate  # noqa: E402, F401


def _dlp_row(locus, ext_prob, om_prob=0.0):
    return {
        "locus_tag": locus,
        "extracellular_prob": str(ext_prob),
        "outer_membrane_prob": str(om_prob),
        "predicted_localization": "Extracellular" if ext_prob >= 0.5 else "Cytoplasm",
        "max_localization": "Extracellular" if ext_prob >= 0.5 else "Cytoplasm",
        "max_probability": str(max(ext_prob, om_prob)),
    }


def _dse_row(locus, ss_type, max_prob=0.9):
    return {
        "locus_tag": locus,
        "dse_ss_type": ss_type,
        "dse_max_prob": str(max_prob),
    }


def _plm_e_row(locus, passes, effector_type="T1SE"):
    return {
        "locus_tag": locus,
        "passes_threshold": "1" if passes else "0",
        "effector_type": effector_type,
    }


def _sp_row(locus, prediction, probability=0.95):
    return {
        "locus_tag": locus,
        "signalp_prediction": prediction,
        "signalp_probability": str(probability),
        "signalp_cs_position": "",
    }


def _run(
    dlp=None, dse=None, plm_e=None, sp=None,
    has_t3ss=False, conf_threshold=0.8, ss_component_types=None,
):
    return list(
        cross_validate(
            dlp_data=dlp or {},
            dse_data=dse or {},
            plm_e_data=plm_e or {},
            sp_data=sp or {},
            sample_id="sample1",
            conf_threshold=conf_threshold,
            has_t3ss=has_t3ss,
            ss_component_types=ss_component_types,
        )
    )


class TestEqualTriggers:
    """DLP, DSE, and PLM-Effector are equal secretion predictors."""

    def test_dlp_alone_marks_secreted(self):
        rows = _run(dlp={"G1": _dlp_row("G1", 0.95)})
        assert len(rows) == 1
        r = rows[0]
        assert r["is_secreted"] is True
        assert r["n_prediction_tools_agreeing"] == 1
        assert r["secretion_evidence"] == "DeepLocPro"

    def test_dse_alone_marks_secreted(self):
        rows = _run(dse={"G1": _dse_row("G1", "T1SS")})
        assert rows[0]["is_secreted"] is True
        assert rows[0]["n_prediction_tools_agreeing"] == 1
        assert rows[0]["secretion_evidence"] == "DeepSecE"

    def test_plm_effector_alone_marks_secreted(self):
        rows = _run(plm_e={"G1": _plm_e_row("G1", passes=True)})
        assert rows[0]["is_secreted"] is True
        assert rows[0]["n_prediction_tools_agreeing"] == 1
        assert rows[0]["secretion_evidence"] == "PLM-Effector"
        assert rows[0]["plm_effector_secreted"] is True

    def test_below_dlp_threshold_not_secreted(self):
        rows = _run(dlp={"G1": _dlp_row("G1", 0.5)}, conf_threshold=0.8)
        assert rows[0]["is_secreted"] is False
        assert rows[0]["n_prediction_tools_agreeing"] == 0
        assert rows[0]["secretion_evidence"] == ""

    def test_dse_non_secreted_not_counted(self):
        rows = _run(dse={"G1": _dse_row("G1", "Non-secreted", max_prob=0.95)})
        assert rows[0]["is_secreted"] is False

    def test_plm_effector_not_passing_not_counted(self):
        rows = _run(plm_e={"G1": _plm_e_row("G1", passes=False)})
        assert rows[0]["is_secreted"] is False
        assert rows[0]["plm_effector_secreted"] is False


class TestAgreeingCount:
    """`n_prediction_tools_agreeing` counts only DLP+DSE+PLM-E (0-3)."""

    def test_two_tools_agree(self):
        rows = _run(
            dlp={"G1": _dlp_row("G1", 0.95)},
            dse={"G1": _dse_row("G1", "T1SS")},
        )
        assert rows[0]["n_prediction_tools_agreeing"] == 2
        assert set(rows[0]["secretion_evidence"].split(",")) == {
            "DeepLocPro",
            "DeepSecE",
        }

    def test_all_three_tools_agree(self):
        rows = _run(
            dlp={"G1": _dlp_row("G1", 0.95)},
            dse={"G1": _dse_row("G1", "T1SS")},
            plm_e={"G1": _plm_e_row("G1", passes=True)},
        )
        assert rows[0]["n_prediction_tools_agreeing"] == 3

    def test_signalp_does_not_increment_count(self):
        rows = _run(sp={"G1": _sp_row("G1", "SP(Sec/SPI)")})
        assert rows[0]["n_prediction_tools_agreeing"] == 0
        assert rows[0]["is_secreted"] is False
        assert rows[0]["signalp_supports_secretion"] is True
        assert "SignalP" not in rows[0]["secretion_evidence"]


class TestSignalPEvidenceOnly:
    """SignalP is recorded but does not trigger `is_secreted`."""

    def test_signalp_alone_not_secreted(self):
        rows = _run(sp={"G1": _sp_row("G1", "LIPO(Sec/SPII)")})
        r = rows[0]
        assert r["is_secreted"] is False
        assert r["signalp_supports_secretion"] is True
        assert r["signalp_prediction"] == "LIPO(Sec/SPII)"

    def test_signalp_negative_recorded(self):
        rows = _run(sp={"G1": _sp_row("G1", "OTHER")})
        assert rows[0]["signalp_supports_secretion"] is False

    def test_signalp_plus_dlp_still_counts_one(self):
        """SignalP doesn't contribute to the agreeing count — even
        together with DLP, the count is 1, not 2."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", 0.95)},
            sp={"G1": _sp_row("G1", "SP(Sec/SPI)")},
        )
        assert rows[0]["n_prediction_tools_agreeing"] == 1
        assert rows[0]["secretion_evidence"] == "DeepLocPro"
        assert rows[0]["signalp_supports_secretion"] is True


class TestDseT3ssGuard:
    """DSE T3SS calls are flagged and excluded when MacSyFinder found
    no T3SS in the genome."""

    def test_dse_t3ss_flagged_without_macsy_t3ss(self):
        rows = _run(dse={"G1": _dse_row("G1", "T3SS")}, has_t3ss=False)
        r = rows[0]
        assert r["dse_T3SS_flagged"] is True
        assert r["is_secreted"] is False
        assert r["n_prediction_tools_agreeing"] == 0

    def test_dse_t3ss_accepted_with_macsy_t3ss(self):
        rows = _run(dse={"G1": _dse_row("G1", "T3SS")}, has_t3ss=True)
        r = rows[0]
        assert r["dse_T3SS_flagged"] is False
        assert r["is_secreted"] is True
        assert r["n_prediction_tools_agreeing"] == 1

    def test_non_t3ss_dse_calls_never_flagged(self):
        rows = _run(dse={"G1": _dse_row("G1", "T1SS")}, has_t3ss=False)
        assert rows[0]["dse_T3SS_flagged"] is False


class TestT5SSLocalisationRule:
    """T5SS substrates: DLP triggers on Extracellular OR Outer membrane.

    Biology: T5aSS autotransporter passenger can be cleaved (extracellular)
    or remain tethered (outer membrane). T5b/c/d/e all have similar duality
    or surface display. The standard extracellular-only rule under-calls
    these — relax to max(ext, om) >= conf_threshold for T5*SS components.
    """

    def test_t5a_om_protein_triggers_dlp(self):
        """A T5aSS component with OM=0.95, ext=0.05 should be flagged."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_types={"G1": "T5aSS"},
        )
        assert rows[0]["is_secreted"] is True
        assert rows[0]["secretion_evidence"] == "DeepLocPro"

    @pytest.mark.parametrize(
        "subtype", ["T5SS", "T5aSS", "T5bSS", "T5cSS", "T5dSS", "T5eSS"]
    )
    def test_all_t5_subtypes_apply_rule(self, subtype):
        """Every TXSScan T5SS subtype gets the relaxed rule."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_types={"G1": subtype},
        )
        assert rows[0]["is_secreted"] is True

    def test_t5a_extracellular_protein_still_triggers(self):
        """Standard extracellular case still works for T5SS."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.95, om_prob=0.0)},
            ss_component_types={"G1": "T5aSS"},
        )
        assert rows[0]["is_secreted"] is True

    def test_non_t5_om_protein_does_not_trigger(self):
        """A T1SS component with OM=0.95, ext=0.05 should NOT be flagged
        — T1SS substrates are extracellular, not OM-tethered."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_types={"G1": "T1SS"},
        )
        assert rows[0]["is_secreted"] is False
        assert rows[0]["n_prediction_tools_agreeing"] == 0

    def test_unmapped_protein_uses_standard_rule(self):
        """Proteins absent from ss_component_types (neighborhood, not a
        component) get the standard extracellular-only rule."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_types={},  # G1 not a component of any system
        )
        assert rows[0]["is_secreted"] is False

    def test_t5a_below_threshold_in_both_does_not_trigger(self):
        """If neither ext nor OM crosses threshold, even T5SS shouldn't pass."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.3, om_prob=0.4)},
            ss_component_types={"G1": "T5aSS"},
            conf_threshold=0.8,
        )
        assert rows[0]["is_secreted"] is False

    def test_no_ss_components_arg_preserves_old_behaviour(self):
        """Calling cross_validate without ss_component_types kwarg = same
        behaviour as before this change."""
        rows = _run(dlp={"G1": _dlp_row("G1", ext_prob=0.95)})  # no ss_component_types
        assert rows[0]["is_secreted"] is True


class TestOutputShape:
    def test_union_of_loci(self):
        """Output should have one row per protein in the union of inputs."""
        rows = _run(
            dlp={"A": _dlp_row("A", 0.5)},
            dse={"B": _dse_row("B", "Non-secreted")},
            plm_e={"C": _plm_e_row("C", passes=False)},
            sp={"D": _sp_row("D", "OTHER")},
        )
        assert [r["locus_tag"] for r in rows] == ["A", "B", "C", "D"]

    def test_all_required_fields_present(self):
        rows = _run(dlp={"G1": _dlp_row("G1", 0.95)})
        required = {
            "locus_tag",
            "sample_id",
            "is_secreted",
            "n_prediction_tools_agreeing",
            "secretion_evidence",
            "dse_T3SS_flagged",
            "plm_effector_secreted",
            "plm_effector_type",
            "signalp_supports_secretion",
        }
        assert required <= set(rows[0].keys())

    def test_sample_id_propagated(self):
        rows = _run(dlp={"G1": _dlp_row("G1", 0.95)})
        assert rows[0]["sample_id"] == "sample1"
