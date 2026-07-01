"""Unit tests for the shared plot_style module."""

import os

from ssign_app.scripts.ssign_lib import plot_style as ps


def test_ss_type_palette_is_stable_and_known():
    """Known SS types map to their fixed canonical colours, deterministically."""
    p1 = ps.ss_type_palette(["T1SS", "T6SS", "T3SS"])
    p2 = ps.ss_type_palette(["T3SS", "T1SS", "T6SS"])
    # Same type -> same colour regardless of input order or subset.
    assert p1["T1SS"] == p2["T1SS"]
    assert p1["T6SS"] == p2["T6SS"]
    assert p1["T1SS"] != p1["T6SS"]


def test_ss_type_palette_variants_inherit_parent():
    pal = ps.ss_type_palette(["T4SS", "pT4SSt", "T6SS", "T6SSi", "T6SSii"])
    assert pal["pT4SSt"] == pal["T4SS"]
    assert pal["T6SSi"] == pal["T6SS"] == pal["T6SSii"]


def test_ss_type_palette_unknown_is_deterministic():
    a = ps.ss_type_palette(["WeirdSS", "OtherSS"])
    b = ps.ss_type_palette(["OtherSS", "WeirdSS"])
    assert a == b  # fallback assigned by sorted order, order-independent


def test_numbered_and_pooled_paths_zero_pad():
    assert ps.numbered_path("/out", 1, "foo").endswith("/out/01_foo.png")
    assert ps.numbered_path("/out", 12, "bar").endswith("/out/12_bar.png")
    assert ps.pooled_path("/out", 2, "baz").endswith("/out/P02_baz.png")


def test_clear_figure_set_removes_owned_keeps_others(tmp_path):
    owned = ["fig1_old.png", "fig7_legacy.png", "01_new.png", "P02_pooled.png"]
    kept = ["pao1_enrichment_fold.png", "report.png", "notes.txt"]
    for fn in owned + kept:
        (tmp_path / fn).write_text("x")

    ps.clear_figure_set(str(tmp_path))

    remaining = set(os.listdir(tmp_path))
    for fn in owned:
        assert fn not in remaining, f"{fn} should have been cleared"
    for fn in kept:
        assert fn in remaining, f"{fn} should have been kept"


def test_clear_figure_set_missing_dir_is_noop():
    ps.clear_figure_set("/nonexistent/path/xyz")  # must not raise


def test_print_figure_index(capsys):
    ps.print_figure_index([("01", "01_foo.png", "context"), ("P02", "P02_bar.png", "pooled")])
    out = capsys.readouterr().out
    assert "Figure index:" in out
    assert "01_foo.png" in out
    assert "P02_bar.png" in out
