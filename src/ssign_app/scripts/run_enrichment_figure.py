#!/usr/bin/env python3
"""Render the per-SS-type circular-shift null-distribution figure.

One panel per (SS type, predictor) cell: the circular-shift null histogram with
the observed value, the null mean, and the fold + permutation p annotated. Rows
are SS types, columns are DeepLocPro and DeepSecE. Window types show "secreted
proteins near the system"; autotransporters (T5aSS/T5cSS) show "self-detection
of the component". Consumes the stats TSV + null-array dump from
enrichment_testing.py.

Null counts are small integers, so each panel uses integer-aligned bins (one bar
per achievable count). Non-significant panels (q >= 0.05) are visually muted so
the eye lands on the significant cells; their statistics are still shown.

    run_enrichment_figure.py --stats <sid>_enrichment_stats.tsv \
        --nulls <sid>_enrichment_nulls.npz --out <sid>_enrichment_null_distributions.png
"""

import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import ENRICH_TOOLS, enrich_null_key  # noqa: E402
from ssign_lib.plot_style import THEME, apply_house_style  # noqa: E402

TOOL_LABEL = {"DLP": "DeepLocPro (extracellular)", "DSE": "DeepSecE (secreted type)"}
_MUTED_RED = "#D9A6A0"  # observed line on a non-significant panel

apply_house_style()


def sig_stars(q: float) -> str:
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "n.s."


def load_stats(path: str) -> list:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fnum(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _integer_bins(null):
    """Bin edges placing one bar per achievable integer count."""
    lo = int(np.floor(null.min()))
    hi = int(np.ceil(null.max()))
    if hi <= lo:
        hi = lo + 1
    return np.arange(lo, hi + 2) - 0.5


def draw_panel(ax, row, null, tool):
    """One null-histogram panel; `row` may be None / skipped (renders an n/a note)."""
    if row is None or row.get("observed", "").strip() == "" or null is None or null.size == 0:
        ax.text(0.5, 0.5, "n/a", ha="center", va="center", transform=ax.transAxes, color="#999999")
        ax.set_xticks([])
        ax.set_yticks([])
        return
    observed = int(float(row["observed"]))
    null_mean = _fnum(row.get("null_mean"))
    fold = _fnum(row.get("fold"))
    q = _fnum(row.get("qvalue"), 1.0)
    sig = q < 0.05

    hist_col = THEME[tool] if sig else THEME["muted"]
    obs_col = THEME["observed"] if sig else _MUTED_RED
    box_edge = THEME["observed"] if sig else THEME["muted"]
    box_txt = THEME["observed"] if sig else "#8A8A8A"

    ax.hist(
        null, bins=_integer_bins(null), color=hist_col, alpha=0.65 if sig else 0.45, edgecolor="white", linewidth=0.4
    )
    ax.axvline(observed, color=obs_col, lw=2.2, label=f"observed = {observed}")
    ax.axvline(null_mean, color=THEME["null_mean"], lw=1.3, ls="--", label=f"null mean = {null_mean:.1f}")
    # Legend top-left, fold/significance box top-right: the observed line sits at
    # the right (enrichment), so keeping the two corners apart avoids collisions.
    ax.text(
        0.97,
        0.93,
        f"{fold:.1f}x  {sig_stars(q)}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=box_txt,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=box_edge, alpha=0.9),
    )
    ax.legend(frameon=False, fontsize=7, loc="upper left")


def main():
    ap = argparse.ArgumentParser(description="Per-SS-type circular-shift null figure")
    ap.add_argument("--stats", required=True)
    ap.add_argument("--nulls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    rows = load_stats(args.stats)
    with np.load(args.nulls) as npz:
        nulls = {k: npz[k] for k in npz.files}
    by_key = {(r["ss_type"], r["tool"]): r for r in rows}
    ss_types = sorted({r["ss_type"] for r in rows})
    tools = list(ENRICH_TOOLS)

    if not ss_types:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No secretion systems detected — nothing to plot", ha="center", va="center")
        ax.axis("off")
        fig.savefig(args.out, dpi=args.dpi)
        plt.close(fig)
        print(f"Saved (empty): {args.out}")
        return

    nrow, ncol = len(ss_types), len(tools)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.6 * nrow), squeeze=False)
    for i, st in enumerate(ss_types):
        for j, tool in enumerate(tools):
            ax = axes[i][j]
            row = by_key.get((st, tool))
            null = nulls.get(enrich_null_key(st, tool))
            draw_panel(ax, row, null, tool)
            if i == 0:
                ax.set_title(TOOL_LABEL[tool])
            if j == 0:
                mode = (row or {}).get("mode", "")
                tag = " (self)" if mode == "self" else ""
                ax.set_ylabel(f"{st}{tag}\nfrequency", fontsize=9)
            if i == nrow - 1:
                ax.set_xlabel("count near SS components (per rotation)")
    fig.suptitle(
        "Circular-shift permutation enrichment (genome-structure-preserving null)",
        fontsize=13,
        fontweight="bold",
        y=1.0,
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
