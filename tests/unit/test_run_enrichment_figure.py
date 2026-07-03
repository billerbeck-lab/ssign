"""Unit tests for the enrichment fold/significance bar chart (headless, Agg).

The figure previously had only opt-in integration coverage, so a missing-import
crash shipped undetected. These render-level tests guard the default suite.
"""

from _helpers import run_script_main, write_tsv

from ssign_app.scripts.run_enrichment_figure import build_columns
from ssign_app.scripts.run_enrichment_figure import main as enr_main

_STATS_FIELDS = ["ss_type", "tool", "observed", "fold", "qvalue", "mode"]

# Covers the tricky rows: significant + non-significant, a finite and an infinite
# fold (null mean 0), a self-mode autotransporter, a skipped DSE-T3SS row (empty
# observed -> no bar), the SignalP third predictor (incl. an inf-fold self row and a
# low-fold n.s. T3SS contrast), and a T5aSS hitchhiker(window) set so the type
# renders two adjacent x-groups.
_STATS_ROWS = [
    {"ss_type": "T1SS", "tool": "DLP", "observed": "5", "fold": "3.2", "qvalue": "0.0005", "mode": "window"},
    {"ss_type": "T1SS", "tool": "DSE", "observed": "4", "fold": "2.1", "qvalue": "0.04", "mode": "window"},
    {"ss_type": "T1SS", "tool": "SignalP", "observed": "3", "fold": "1.8", "qvalue": "0.03", "mode": "window"},
    {"ss_type": "T3SS", "tool": "DLP", "observed": "2", "fold": "inf", "qvalue": "0.02", "mode": "window"},
    {"ss_type": "T3SS", "tool": "DSE", "observed": "", "fold": "", "qvalue": "", "mode": "window"},
    {"ss_type": "T3SS", "tool": "SignalP", "observed": "1", "fold": "0.9", "qvalue": "0.6", "mode": "window"},
    {"ss_type": "T5aSS", "tool": "DLP", "observed": "3", "fold": "1.5", "qvalue": "0.2", "mode": "self"},
    {"ss_type": "T5aSS", "tool": "DSE", "observed": "1", "fold": "0.8", "qvalue": "0.9", "mode": "self"},
    {"ss_type": "T5aSS", "tool": "SignalP", "observed": "4", "fold": "inf", "qvalue": "0.01", "mode": "self"},
    # T5aSS hitchhiker (window) rows: DLP-or-DSE scored, distinct from the self set.
    {"ss_type": "T5aSS", "tool": "DLP", "observed": "2", "fold": "2.5", "qvalue": "0.03", "mode": "window"},
    {"ss_type": "T5aSS", "tool": "DSE", "observed": "2", "fold": "2.7", "qvalue": "0.02", "mode": "window"},
    {"ss_type": "T5aSS", "tool": "SignalP", "observed": "1", "fold": "0.9", "qvalue": "0.8", "mode": "window"},
    # Combined rows (one per type per mode) — only plotted by --combined. Window
    # types carry the DLP-or-DSE score; the T5aSS self row carries SignalP; the
    # T5aSS hitchhiker window row is DLP-or-DSE again.
    {"ss_type": "T1SS", "tool": "COMBINED", "observed": "6", "fold": "3.6", "qvalue": "0.0003", "mode": "window"},
    {"ss_type": "T3SS", "tool": "COMBINED", "observed": "2", "fold": "8.0", "qvalue": "0.02", "mode": "window"},
    {"ss_type": "T5aSS", "tool": "COMBINED", "observed": "4", "fold": "inf", "qvalue": "0.01", "mode": "self"},
    {"ss_type": "T5aSS", "tool": "COMBINED", "observed": "3", "fold": "3.1", "qvalue": "0.02", "mode": "window"},
]


def test_renders_fold_chart(monkeypatch, tmp_path):
    stats = write_tsv(str(tmp_path / "enrichment_stats.tsv"), _STATS_FIELDS, _STATS_ROWS)
    out = tmp_path / "enrichment_fold.png"
    run_script_main(
        monkeypatch, enr_main, ["run_enrichment_figure.py", "--stats", stats, "--out", str(out), "--dpi", "80"]
    )
    assert out.exists() and out.stat().st_size > 0


def test_renders_combined_chart(monkeypatch, tmp_path):
    stats = write_tsv(str(tmp_path / "enrichment_stats.tsv"), _STATS_FIELDS, _STATS_ROWS)
    out = tmp_path / "enrichment_fold_combined.png"
    run_script_main(
        monkeypatch,
        enr_main,
        ["run_enrichment_figure.py", "--stats", stats, "--out", str(out), "--combined", "--dpi", "80"],
    )
    assert out.exists() and out.stat().st_size > 0


def test_no_detected_systems_still_writes_placeholder(monkeypatch, tmp_path):
    # All rows skipped (empty observed) -> the empty-branch placeholder figure.
    rows = [{"ss_type": "T3SS", "tool": "DSE", "observed": "", "fold": "", "qvalue": "", "mode": "window"}]
    stats = write_tsv(str(tmp_path / "empty_stats.tsv"), _STATS_FIELDS, rows)
    out = tmp_path / "enrichment_fold.png"
    run_script_main(
        monkeypatch, enr_main, ["run_enrichment_figure.py", "--stats", stats, "--out", str(out), "--dpi", "80"]
    )
    assert out.exists() and out.stat().st_size > 0


def test_build_columns_splits_autotransporter_into_two_groups():
    # T5aSS carries both modes -> two adjacent columns (self then hitch); a plain
    # window type stays one column; a self-only type keeps its (self) tag.
    ss_types = ["T1SS", "T5aSS", "T5cSS"]
    col_modes = {"T1SS": {"window"}, "T5aSS": {"self", "window"}, "T5cSS": {"self"}}
    cols = build_columns(ss_types, col_modes)
    assert cols == [
        ("T1SS", "window", "T1SS"),
        ("T5aSS", "self", "T5aSS\n(self)"),
        ("T5aSS", "window", "T5aSS\n(hitchhiker)"),
        ("T5cSS", "self", "T5cSS\n(self)"),
    ]
