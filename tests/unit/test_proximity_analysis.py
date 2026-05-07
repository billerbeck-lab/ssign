"""Regression tests for proximity_analysis.py.

Replaces the prior sham implementation that re-implemented the proximity
logic inline instead of testing the production module. These tests:

- Import dse_type_in_genome and main directly from the production script.
- Cover Critical Bug Fix #1 (per-component ±N window, NOT full-system span).
- Cover Critical Bug Fix #2 (DSE cross-genome leakage guard).

Inputs are real TSV files written to tmp_dir; main() runs end-to-end via
sys.argv monkey-patching, then the output TSV is parsed and asserted on.
"""

import os
import sys

import pytest

# Production module

from _helpers import (  # noqa: E402
    PREDICTIONS_FIELDS,
    SS_COMPONENT_FIELDS,
    make_prediction_row,
    make_ss_component_row,
    read_tsv_rows,
    run_script_main,
    write_tsv,
)
from proximity_analysis import dse_type_in_genome  # noqa: E402
from proximity_analysis import main as proximity_main  # noqa: E402

# ---------------------------------------------------------------------------
# dse_type_in_genome — Critical Bug Fix #2
# ---------------------------------------------------------------------------
# DSE_TO_MACSYFINDER expansions:
#   T1SS→[T1SS], T2SS→[T2SS], T3SS→[T3SS],
#   T4SS→[pT4SSt, T4SS], T6SS→[T6SSi, T6SS]
# The check is a substring match (`mf_name in genome_type`), so T6SSi
# covers T6SSii / T6SSiii naturally.


@pytest.mark.parametrize(
    "dse_type, genome, expected",
    [
        ("T1SS", {"T1SS"}, True),
        ("T1SS", {"T2SS", "T6SSi"}, False),
        ("T6SS", {"T6SSi"}, True),
        ("T6SS", {"T6SSii"}, True),
        ("T6SS", {"T6SSiii"}, True),
        ("T4SS", {"pT4SSt"}, True),
        ("T4SS", {"T4SS"}, True),
        ("T1SS", set(), False),
        # Unknown DSE type → fallback list is [literal]; substring match applies
        ("Foo", {"Foobar"}, True),
        ("Foo", {"Bar"}, False),
    ],
)
def test_dse_type_in_genome(dse_type, genome, expected):
    assert dse_type_in_genome(dse_type, genome) is expected


# ---------------------------------------------------------------------------
# proximity_analysis.main() — Critical Bug Fix #1 (per-component window)
# ---------------------------------------------------------------------------


def _run_proximity(
    monkeypatch,
    tmp_dir,
    gene_order,
    ss_components,
    predictions,
    *,
    window=3,
    conf=0.8,
):
    """Invoke proximity_analysis.main() in-process and return parsed output rows."""
    output_path = os.path.join(tmp_dir, "substrates.tsv")
    run_script_main(
        monkeypatch,
        proximity_main,
        [
            "proximity_analysis",
            "--gene-order",
            gene_order,
            "--ss-components",
            ss_components,
            "--predictions",
            predictions,
            "--sample",
            "test_sample",
            "--window",
            str(window),
            "--conf-threshold",
            str(conf),
            "--output",
            output_path,
        ],
    )
    return read_tsv_rows(output_path)


def _all_loci():
    """Match the conftest two_contig_genes layout."""
    return [f"GENE_{i:04d}" for i in range(10)] + [f"GENEB_{i:04d}" for i in range(5)]


def _predictions_with_dlp(tmp_dir, dlp_positives):
    """Predictions TSV: every locus present, listed loci flagged DLP+ (≥0.8)."""
    rows = [make_prediction_row(locus, dlp_ext=0.95 if locus in dlp_positives else 0.05) for locus in _all_loci()]
    return write_tsv(os.path.join(tmp_dir, "predictions.tsv"), PREDICTIONS_FIELDS, rows)


def _predictions_with_dse(tmp_dir, dse_calls):
    """Predictions TSV: dse_calls is {locus: dse_ss_type}; DLP=0 for everyone."""
    rows = [
        make_prediction_row(
            locus,
            dlp_ext=0.0,
            dse_type=dse_calls[locus] if locus in dse_calls else "Non-secreted",
            dse_prob=0.95 if locus in dse_calls else 0.0,
        )
        for locus in _all_loci()
    ]
    return write_tsv(os.path.join(tmp_dir, "predictions.tsv"), PREDICTIONS_FIELDS, rows)


class TestPerComponentWindow:
    """Components live at GENE_0005 + GENE_0006 on contig_A. Window ±3 around either
    covers genes 2-9 (union), excluding the components themselves."""

    def test_per_component_window_returns_expected_neighbors(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        dlp_positives = {f"GENE_{i:04d}" for i in [2, 3, 4, 7, 8, 9]}
        predictions = _predictions_with_dlp(tmp_dir, dlp_positives)
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert {r["locus_tag"] for r in rows} == dlp_positives

    def test_window_does_not_cross_contigs(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        # All DLP-positive loci sit on contig_B; components are on contig_A.
        predictions = _predictions_with_dlp(tmp_dir, {f"GENEB_{i:04d}" for i in range(5)})
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert rows == []

    def test_components_themselves_are_not_substrates(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        predictions = _predictions_with_dlp(tmp_dir, {"GENE_0005", "GENE_0006"})
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert rows == []

    def test_excluded_components_skipped(self, monkeypatch, tmp_dir, gene_order_tsv):
        ss_components = write_tsv(
            os.path.join(tmp_dir, "ss_components.tsv"),
            SS_COMPONENT_FIELDS,
            [make_ss_component_row("GENE_0005", "Flagellum", "fliC", excluded="True")],
        )
        predictions = _predictions_with_dlp(
            tmp_dir,
            {"GENE_0003", "GENE_0004", "GENE_0006"},
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components,
            predictions,
        )
        assert rows == []

    def test_window_at_contig_start(self, monkeypatch, tmp_dir, gene_order_tsv):
        # Component at GENE_0000 → window ±3 stops at gene index 0 on the low side
        ss_components = write_tsv(
            os.path.join(tmp_dir, "ss_components.tsv"),
            SS_COMPONENT_FIELDS,
            [make_ss_component_row("GENE_0000", "T1SS", "tolC")],
        )
        # Genes 1-3 (in window) and 4 (out of window) flagged positive
        predictions = _predictions_with_dlp(
            tmp_dir,
            {"GENE_0001", "GENE_0002", "GENE_0003", "GENE_0004"},
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components,
            predictions,
        )
        assert {r["locus_tag"] for r in rows} == {"GENE_0001", "GENE_0002", "GENE_0003"}


class TestDseLeakageGuard:
    """Critical Bug Fix #2: DSE-only calls (DLP negative) require the predicted
    SS type to exist in the genome. Otherwise we treat it as cross-genome leakage."""

    def test_dse_t2ss_kept_when_genome_has_t2ss(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        predictions = _predictions_with_dse(
            tmp_dir,
            {"GENE_0003": "T2SS", "GENE_0007": "T2SS"},
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert {r["locus_tag"] for r in rows} == {"GENE_0003", "GENE_0007"}

    def test_dse_t1ss_dropped_when_genome_lacks_t1ss(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        predictions = _predictions_with_dse(
            tmp_dir,
            {"GENE_0003": "T1SS", "GENE_0007": "T1SS"},
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert rows == []

    def test_dse_t3ss_dropped_unconditionally(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        # T3SS is excluded outright in the is_dse predicate (Fix #4 territory)
        predictions = _predictions_with_dse(
            tmp_dir,
            {"GENE_0003": "T3SS", "GENE_0007": "T3SS"},
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert rows == []

    def test_dlp_positive_skips_dse_leakage_check(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        # When DLP also flags, the leakage guard is bypassed (DLP carries the call).
        rows_pred = [
            make_prediction_row(
                locus,
                dlp_ext=0.95 if locus == "GENE_0003" else 0.0,
                dse_type="T1SS" if locus == "GENE_0003" else "Non-secreted",
                dse_prob=0.95 if locus == "GENE_0003" else 0.0,
            )
            for locus in _all_loci()
        ]
        predictions = write_tsv(
            os.path.join(tmp_dir, "predictions.tsv"),
            PREDICTIONS_FIELDS,
            rows_pred,
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert {r["locus_tag"] for r in rows} == {"GENE_0003"}


class TestToolAttribution:
    """The `tool` column reports DLP, DSE, or DLP+DSE depending on which fired."""

    def test_dlp_only(self, monkeypatch, tmp_dir, gene_order_tsv, ss_components_tsv):
        predictions = _predictions_with_dlp(tmp_dir, {"GENE_0003"})
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert len(rows) == 1
        assert rows[0]["tool"] == "DLP"

    def test_dse_only(self, monkeypatch, tmp_dir, gene_order_tsv, ss_components_tsv):
        predictions = _predictions_with_dse(tmp_dir, {"GENE_0003": "T2SS"})
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert len(rows) == 1
        assert rows[0]["tool"] == "DSE"

    def test_both_tools_fire(
        self,
        monkeypatch,
        tmp_dir,
        gene_order_tsv,
        ss_components_tsv,
    ):
        rows_pred = [
            make_prediction_row(
                locus,
                dlp_ext=0.95 if locus == "GENE_0003" else 0.0,
                dse_type="T2SS" if locus == "GENE_0003" else "Non-secreted",
                dse_prob=0.95 if locus == "GENE_0003" else 0.0,
            )
            for locus in _all_loci()
        ]
        predictions = write_tsv(
            os.path.join(tmp_dir, "predictions.tsv"),
            PREDICTIONS_FIELDS,
            rows_pred,
        )
        rows = _run_proximity(
            monkeypatch,
            tmp_dir,
            gene_order_tsv,
            ss_components_tsv,
            predictions,
        )
        assert len(rows) == 1
        assert rows[0]["tool"] == "DLP+DSE"
