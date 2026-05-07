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
    has_t3ss=False, conf_threshold=0.8, ss_component_info=None,
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
            ss_component_info=ss_component_info,
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
    """Per-component T5SS DLP rules — match real TXSScan biology.

    TXSScan v2 models exactly three T5 subtypes (T5a/b/c), each with one
    HMM-modelled component (the β-barrel/translocator). The DLP rule is
    per (ss_type, gene_name), not per subtype.
    """

    def test_t5a_autotransporter_om_passes(self):
        """T5aSS_PF03797 with OM=0.95 (untethered passenger) passes."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_info={"G1": ("T5aSS", "T5aSS_PF03797")},
        )
        assert rows[0]["is_secreted"] is True
        assert rows[0]["secretion_evidence"] == "DeepLocPro"

    def test_t5a_autotransporter_extracellular_passes(self):
        """T5aSS_PF03797 with ext=0.95 (cleaved passenger) also passes."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.95, om_prob=0.0)},
            ss_component_info={"G1": ("T5aSS", "T5aSS_PF03797")},
        )
        assert rows[0]["is_secreted"] is True

    def test_t5b_translocator_om_passes(self):
        """T5bSS_translocator (TpsB pore) with OM=0.95 passes — biology."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_info={"G1": ("T5bSS", "T5bSS_translocator")},
        )
        assert rows[0]["is_secreted"] is True

    def test_t5b_translocator_extracellular_does_NOT_pass(self):
        """The headline correctness fix: a TpsB pore predicted as
        extracellular is biologically wrong — pores live in the outer
        membrane. Don't accept extracellular as evidence for TpsB."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.95, om_prob=0.05)},
            ss_component_info={"G1": ("T5bSS", "T5bSS_translocator")},
        )
        assert rows[0]["is_secreted"] is False
        assert rows[0]["n_prediction_tools_agreeing"] == 0

    @pytest.mark.parametrize(
        "ss_type,gene_name",
        [
            ("T5aSS", "T5aSS_PF03797"),
            ("T5cSS", "T5cSS_PF03895"),
        ],
    )
    def test_t5a_t5c_both_localizations_valid(self, ss_type, gene_name):
        """T5a (cleaved-or-tethered passenger) and T5c (surface-displayed
        trimeric AT) both legitimately localize to either Extracellular
        or Outer membrane. T5b is excluded — its translocator is OM-only."""
        rows_om = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_info={"G1": (ss_type, gene_name)},
        )
        rows_ext = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.95, om_prob=0.05)},
            ss_component_info={"G1": (ss_type, gene_name)},
        )
        assert rows_om[0]["is_secreted"] is True
        assert rows_ext[0]["is_secreted"] is True

    def test_non_t5_om_protein_does_not_trigger(self):
        """A T1SS component with OM=0.95, ext=0.05 should NOT be flagged
        — T1SS substrates are extracellular, not OM-tethered."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_info={"G1": ("T1SS", "T1SS_abc")},
        )
        assert rows[0]["is_secreted"] is False
        assert rows[0]["n_prediction_tools_agreeing"] == 0

    def test_unmapped_protein_uses_standard_rule(self):
        """Proteins absent from ss_component_info (neighborhood, not a
        component) get the standard extracellular-only rule. The TpsA
        passenger of T5bSS lands here because TXSScan doesn't model it
        as a system component — and the standard ext-only rule is
        correct biology for the secreted partner."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_info={},
        )
        assert rows[0]["is_secreted"] is False

    def test_t5a_below_threshold_in_both_does_not_trigger(self):
        """If neither ext nor OM crosses threshold, even T5SS shouldn't pass."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.3, om_prob=0.4)},
            ss_component_info={"G1": ("T5aSS", "T5aSS_PF03797")},
            conf_threshold=0.8,
        )
        assert rows[0]["is_secreted"] is False

    def test_no_ss_components_arg_preserves_old_behaviour(self):
        """Calling cross_validate without ss_component_info kwarg = same
        behaviour as before this change."""
        rows = _run(dlp={"G1": _dlp_row("G1", ext_prob=0.95)})
        assert rows[0]["is_secreted"] is True

    def test_t5_subtype_with_unknown_gene_name_falls_back(self):
        """If MacSyFinder/TXSScan ever emits a new T5 component we
        haven't mapped (e.g. a new T5 subtype model), fall back to the
        safe extracellular-only rule."""
        rows = _run(
            dlp={"G1": _dlp_row("G1", ext_prob=0.05, om_prob=0.95)},
            ss_component_info={"G1": ("T5aSS", "T5aSS_some_future_gene")},
        )
        assert rows[0]["is_secreted"] is False


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
