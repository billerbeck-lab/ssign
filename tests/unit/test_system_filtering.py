"""Tests for system_filtering.py.

Covers proximity ⊕ T5SS merging, the substrate_source tag, exclusion
of substrates whose only nearby SS type is in the excluded list, the
T5SS-self carve-out, and the optional --filter-dse-type-mismatch flag
(DSE-only substrates whose predicted type doesn't match a nearby
MacSyFinder system).
"""

import os

import pytest
from _helpers import read_tsv_rows, run_script_main, write_tsv
from system_filtering import main as system_filtering_main

# Substrate-table schema produced by proximity_analysis + t5ss_handler. The
# script reads whatever columns are present and only relies on:
#   locus_tag, nearby_ss_types, tool, dse_type_match (if --filter-dse-type-mismatch)
SUBSTRATE_FIELDS = [
    "locus_tag",
    "sample_id",
    "tool",
    "nearby_ss_types",
    "dse_type_match",
    "dlp_extracellular_prob",
    "outer_membrane_prob",
    "signalp_prediction",
]

# valid_systems and predictions are required CLI args but the script doesn't
# actually consume them — pass empty TSVs to satisfy argparse.
_VALID_SYSTEMS_FIELDS = ["sys_id", "ss_type", "wholeness"]
_PREDICTIONS_FIELDS = ["locus_tag"]


def _make_substrate(
    locus_tag,
    *,
    nearby_ss_types="T2SS",
    tool="DLP",
    dse_type_match="True",
    dlp_extracellular_prob="",
    outer_membrane_prob="",
    signalp_prediction="",
):
    return {
        "locus_tag": locus_tag,
        "sample_id": "test_sample",
        "tool": tool,
        "nearby_ss_types": nearby_ss_types,
        "dse_type_match": dse_type_match,
        "dlp_extracellular_prob": dlp_extracellular_prob,
        "outer_membrane_prob": outer_membrane_prob,
        "signalp_prediction": signalp_prediction,
    }


def _run_filter(
    monkeypatch,
    tmp_dir,
    proximity_rows,
    t5ss_rows=None,
    *,
    excluded_systems="Flagellum,Tad,T3SS",
    filter_dse_mismatch=False,
):
    """Run system_filtering.main() and return (filtered_rows, all_rows)."""
    proximity = write_tsv(
        os.path.join(tmp_dir, "proximity.tsv"),
        SUBSTRATE_FIELDS,
        proximity_rows,
    )
    t5ss = write_tsv(
        os.path.join(tmp_dir, "t5ss.tsv"),
        SUBSTRATE_FIELDS,
        t5ss_rows or [],
    )
    valid = write_tsv(
        os.path.join(tmp_dir, "valid.tsv"),
        _VALID_SYSTEMS_FIELDS,
        [],
    )
    predictions = write_tsv(
        os.path.join(tmp_dir, "predictions.tsv"),
        _PREDICTIONS_FIELDS,
        [],
    )
    out_filtered = os.path.join(tmp_dir, "substrates_filtered.tsv")
    out_all = os.path.join(tmp_dir, "substrates_all.tsv")

    argv = [
        "system_filtering",
        "--proximity-substrates",
        proximity,
        "--t5ss-substrates",
        t5ss,
        "--valid-systems",
        valid,
        "--predictions",
        predictions,
        "--sample",
        "test_sample",
        "--excluded-systems",
        excluded_systems,
        "--out-filtered",
        out_filtered,
        "--out-all",
        out_all,
    ]
    if filter_dse_mismatch:
        argv.append("--filter-dse-type-mismatch")
    run_script_main(monkeypatch, system_filtering_main, argv)
    return read_tsv_rows(out_filtered), read_tsv_rows(out_all)


class TestSourceMerging:
    def test_substrates_from_both_sources_appear_in_output(self, monkeypatch, tmp_dir):
        _, all_rows = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("PROX_1")],
            [_make_substrate("T5SS_1")],
        )
        assert {r["locus_tag"] for r in all_rows} == {"PROX_1", "T5SS_1"}

    def test_substrate_source_tag_distinguishes_origin(self, monkeypatch, tmp_dir):
        _, all_rows = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("PROX_1")],
            [_make_substrate("T5SS_1")],
        )
        by_locus = {r["locus_tag"]: r for r in all_rows}
        assert by_locus["PROX_1"]["substrate_source"] == "proximity"
        assert by_locus["T5SS_1"]["substrate_source"] == "T5SS-self"


class TestExclusionFilter:
    def test_substrate_only_near_excluded_system_dropped(self, monkeypatch, tmp_dir):
        # Genome contains a flagellum component (excluded by default).
        # The substrate has no other nearby system → drop.
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", nearby_ss_types="Flagellum")],
        )
        assert filtered == []

    def test_substrate_near_mixed_systems_kept(self, monkeypatch, tmp_dir):
        # Substrate is near Flagellum + T2SS — only T2SS survives the filter,
        # nearby_ss_types is rewritten accordingly.
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", nearby_ss_types="Flagellum,T2SS")],
        )
        assert len(filtered) == 1
        assert filtered[0]["nearby_ss_types"] == "T2SS"

    def test_unfiltered_output_keeps_excluded_substrates(self, monkeypatch, tmp_dir):
        # The audit trail (out_all) preserves the substrate even if the
        # filtered view drops it.
        _, all_rows = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", nearby_ss_types="Flagellum")],
        )
        assert {r["locus_tag"] for r in all_rows} == {"X"}

    def test_t5ss_self_with_evidence_kept_despite_excluded_only(self, monkeypatch, tmp_dir):
        # T5SS-self carve-out: even if nearby_ss_types is empty / fully excluded, a
        # T5SS substrate WITH evidence (here a Sec signal peptide) is kept (it IS
        # its own system).
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [],
            [_make_substrate("T5SS_1", nearby_ss_types="", signalp_prediction="SP")],
        )
        assert {r["locus_tag"] for r in filtered} == {"T5SS_1"}

    def test_t5ss_self_without_evidence_dropped(self, monkeypatch, tmp_dir):
        # New evidence gate: a T5SS component with neither DeepLocPro localization
        # nor a Sec signal peptide is dropped (previously kept unconditionally).
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [],
            [_make_substrate("T5SS_1", nearby_ss_types="T5aSS")],  # no DLP, no SignalP
        )
        assert filtered == []

    def test_t5ss_self_kept_via_dlp_localization(self, monkeypatch, tmp_dir):
        # The DeepLocPro path still calls a T5 component (extracellular/OM) with no
        # signal peptide.
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [],
            [_make_substrate("T5SS_1", nearby_ss_types="T5aSS", outer_membrane_prob="0.95")],
        )
        assert {r["locus_tag"] for r in filtered} == {"T5SS_1"}

    def test_t5bss_translocator_om_only_rule_preserved(self, monkeypatch, tmp_dir):
        # Bug-fix #6: T5bSS translocator passes on outer-membrane only, never
        # extracellular. A T5bSS row that is extracellular-high but OM-low and has
        # no signal peptide must be dropped.
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [],
            [
                _make_substrate(
                    "T5B_1", nearby_ss_types="T5bSS", dlp_extracellular_prob="0.99", outer_membrane_prob="0.1"
                )
            ],
        )
        assert filtered == []

    def test_custom_excluded_systems_override(self, monkeypatch, tmp_dir):
        # User passes --excluded-systems="" → nothing excluded → Flagellum kept
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", nearby_ss_types="Flagellum")],
            excluded_systems="",
        )
        assert {r["locus_tag"] for r in filtered} == {"X"}

    def test_appendage_substrate_dropped_under_default_exclusions(self, monkeypatch, tmp_dir):
        # The production default excludes non-secretion appendages: a protein whose
        # only nearby system is a T4aP pilus must not be substrate-called.
        from ssign_app.scripts.ssign_lib.constants import DEFAULT_EXCLUDED_SYSTEMS

        default = ",".join(DEFAULT_EXCLUDED_SYSTEMS)
        for appendage in ("T4aP", "T4bP", "MSH", "ComM", "Archaeal-T4P"):
            filtered, _ = _run_filter(
                monkeypatch,
                tmp_dir,
                [_make_substrate("X", nearby_ss_types=appendage)],
                excluded_systems=default,
            )
            assert filtered == [], f"{appendage} substrate should be dropped by default"


class TestDseTypeMismatchFilter:
    """--filter-dse-type-mismatch is opt-in; only affects DSE-only substrates."""

    def test_disabled_by_default(self, monkeypatch, tmp_dir):
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", tool="DSE", dse_type_match="False")],
        )
        # Default: DSE-only with mismatched type still appears in filtered
        assert {r["locus_tag"] for r in filtered} == {"X"}

    def test_dse_only_mismatch_dropped_when_enabled(self, monkeypatch, tmp_dir):
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", tool="DSE", dse_type_match="False")],
            filter_dse_mismatch=True,
        )
        assert filtered == []

    def test_dse_only_match_kept_when_enabled(self, monkeypatch, tmp_dir):
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", tool="DSE", dse_type_match="True")],
            filter_dse_mismatch=True,
        )
        assert {r["locus_tag"] for r in filtered} == {"X"}

    @pytest.mark.parametrize("tool", ["DLP", "DLP+DSE"])
    def test_dlp_substrates_immune_to_dse_filter(self, monkeypatch, tmp_dir, tool):
        # DSE type-match filter only acts on DSE-only substrates;
        # DLP and DLP+DSE substrates are kept regardless of match status.
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", tool=tool, dse_type_match="False")],
            filter_dse_mismatch=True,
        )
        assert {r["locus_tag"] for r in filtered} == {"X"}


class TestFieldnamesUnion:
    """T5SS-self substrates carry t5_quality_flag; proximity substrates don't.

    The merged output must include columns from BOTH sources — otherwise
    whichever schema appears first wins and the other side's columns are
    silently dropped.
    """

    def test_union_of_keys_preserves_t5_only_column(self, monkeypatch, tmp_dir):
        # Proximity row uses the standard schema; T5 row has an extra column.
        proximity_path = os.path.join(tmp_dir, "proximity.tsv")
        t5ss_path = os.path.join(tmp_dir, "t5ss.tsv")
        write_tsv(proximity_path, SUBSTRATE_FIELDS, [_make_substrate("PROX_1")])
        write_tsv(
            t5ss_path,
            SUBSTRATE_FIELDS + ["t5_quality_flag"],
            [{**_make_substrate("T5SS_1"), "t5_quality_flag": "barrel_only"}],
        )
        valid = write_tsv(os.path.join(tmp_dir, "valid.tsv"), _VALID_SYSTEMS_FIELDS, [])
        predictions = write_tsv(os.path.join(tmp_dir, "preds.tsv"), _PREDICTIONS_FIELDS, [])
        out_filtered = os.path.join(tmp_dir, "filtered.tsv")
        out_all = os.path.join(tmp_dir, "all.tsv")
        run_script_main(
            monkeypatch,
            system_filtering_main,
            [
                "system_filtering",
                "--proximity-substrates",
                proximity_path,
                "--t5ss-substrates",
                t5ss_path,
                "--valid-systems",
                valid,
                "--predictions",
                predictions,
                "--sample",
                "test_sample",
                "--out-filtered",
                out_filtered,
                "--out-all",
                out_all,
            ],
        )
        all_rows = read_tsv_rows(out_all)
        by_locus = {r["locus_tag"]: r for r in all_rows}
        assert "t5_quality_flag" in by_locus["T5SS_1"]
        assert by_locus["T5SS_1"]["t5_quality_flag"] == "barrel_only"
        # Proximity row doesn't have the column populated but the column exists
        assert by_locus["PROX_1"].get("t5_quality_flag", "") == ""


class TestEmptyInputs:
    def test_no_substrates_writes_header_only(self, monkeypatch, tmp_dir):
        filtered, all_rows = _run_filter(monkeypatch, tmp_dir, [])
        assert filtered == []
        assert all_rows == []
        # Header still present in both output files
        with open(os.path.join(tmp_dir, "substrates_filtered.tsv")) as f:
            assert f.readline().strip()
        with open(os.path.join(tmp_dir, "substrates_all.tsv")) as f:
            assert f.readline().strip()


class TestLocalizationGateEndToEnd:
    """Smoke-test the gate's real code path with non-empty ss_components +
    predictions inputs.

    Earlier tests left both inputs empty, which triggered the no-systems
    short-circuit in _run_localization_gate and skipped over the
    load_tsv_by_key call. That hid a missing import (caught on CX3
    Jun 6 2026, jobs 2939885 + 2939886, crash at step 13 of 24). This
    test exercises the same code path the runner does so a NameError of
    that kind shows up in CI.
    """

    def test_gate_runs_with_real_components_and_predictions(self, monkeypatch, tmp_dir):
        # One system with one component; tag it with a real ss_type so the
        # shipped literature table has rules to score against. The exact
        # verdict (pass/fail) doesn't matter — we just need the gate to
        # complete without raising.
        ss_components = write_tsv(
            os.path.join(tmp_dir, "ss_components.tsv"),
            ["sys_id", "ss_type", "locus_tag", "gene_name", "excluded"],
            [
                {
                    "sys_id": "sys_1",
                    "ss_type": "T2SS",
                    "locus_tag": "COMP_1",
                    "gene_name": "T2SS_gspD",
                    "excluded": "False",
                }
            ],
        )
        predictions = write_tsv(
            os.path.join(tmp_dir, "predictions.tsv"),
            [
                "locus_tag",
                "dlp_extracellular_prob",
                "outer_membrane_prob",
                "periplasm_prob",
                "cytoplasm_prob",
                "cytoplasmic_membrane_prob",
            ],
            [
                {
                    "locus_tag": "COMP_1",
                    "dlp_extracellular_prob": "0.05",
                    "outer_membrane_prob": "0.90",
                    "periplasm_prob": "0.02",
                    "cytoplasm_prob": "0.02",
                    "cytoplasmic_membrane_prob": "0.01",
                }
            ],
        )
        proximity = write_tsv(
            os.path.join(tmp_dir, "proximity.tsv"),
            SUBSTRATE_FIELDS,
            [_make_substrate("PROX_1", nearby_ss_types="T2SS")],
        )
        t5ss = write_tsv(os.path.join(tmp_dir, "t5ss.tsv"), SUBSTRATE_FIELDS, [])
        valid = write_tsv(
            os.path.join(tmp_dir, "valid.tsv"),
            _VALID_SYSTEMS_FIELDS,
            [{"sys_id": "sys_1", "ss_type": "T2SS", "wholeness": "1.0"}],
        )
        out_filtered = os.path.join(tmp_dir, "substrates_filtered.tsv")
        out_all = os.path.join(tmp_dir, "substrates_all.tsv")
        out_verdicts = os.path.join(tmp_dir, "gate_verdicts.tsv")
        run_script_main(
            monkeypatch,
            system_filtering_main,
            [
                "system_filtering",
                "--proximity-substrates",
                proximity,
                "--t5ss-substrates",
                t5ss,
                "--valid-systems",
                valid,
                "--predictions",
                predictions,
                "--ss-components",
                ss_components,
                "--sample",
                "test_sample",
                "--out-filtered",
                out_filtered,
                "--out-all",
                out_all,
                "--out-gate-verdicts",
                out_verdicts,
            ],
        )
        # Gate ran end-to-end; verdicts file exists with at least the header row.
        with open(out_verdicts) as f:
            assert f.readline().strip(), "gate verdicts file is empty"
