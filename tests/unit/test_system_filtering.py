"""Tests for system_filtering.py.

Covers proximity ⊕ T5SS merging, the substrate_source tag, exclusion
of substrates whose only nearby SS type is in the excluded list, the
T5SS-self carve-out, and the optional --filter-dse-type-mismatch flag
(DSE-only substrates whose predicted type doesn't match a nearby
MacSyFinder system).
"""

import os
import sys

import pytest


from _helpers import read_tsv_rows, run_script_main, write_tsv  # noqa: E402
from system_filtering import main as system_filtering_main  # noqa: E402

# Substrate-table schema produced by proximity_analysis + t5ss_handler. The
# script reads whatever columns are present and only relies on:
#   locus_tag, nearby_ss_types, tool, dse_type_match (if --filter-dse-type-mismatch)
SUBSTRATE_FIELDS = [
    "locus_tag",
    "sample_id",
    "tool",
    "nearby_ss_types",
    "dse_type_match",
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
):
    return {
        "locus_tag": locus_tag,
        "sample_id": "test_sample",
        "tool": tool,
        "nearby_ss_types": nearby_ss_types,
        "dse_type_match": dse_type_match,
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

    def test_t5ss_self_kept_despite_excluded_only(self, monkeypatch, tmp_dir):
        # T5SS-self carve-out: even if nearby_ss_types is empty / fully
        # excluded, T5SS substrates are kept (they ARE their own system).
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [],
            [_make_substrate("T5SS_1", nearby_ss_types="")],
        )
        assert {r["locus_tag"] for r in filtered} == {"T5SS_1"}

    def test_custom_excluded_systems_override(self, monkeypatch, tmp_dir):
        # User passes --excluded-systems="" → nothing excluded → Flagellum kept
        filtered, _ = _run_filter(
            monkeypatch,
            tmp_dir,
            [_make_substrate("X", nearby_ss_types="Flagellum")],
            excluded_systems="",
        )
        assert {r["locus_tag"] for r in filtered} == {"X"}


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
