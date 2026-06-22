#!/usr/bin/env python3
"""Render the per-SS-type circular-shift null-distribution figure.

One panel per (SS type, predictor) cell: the circular-shift null histogram with
the observed value, the null mean, and the fold + permutation p annotated. Rows
are SS types, columns are DeepLocPro and DeepSecE. Window types show "secreted
proteins near the system"; autotransporters (T5aSS/T5cSS) show "self-detection
of the component". Consumes the stats TSV + null-array dump from
enrichment_testing.py.

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

TOOL_LABEL = {"DLP": "DeepLocPro (extracellular)", "DSE": "DeepSecE (secreted type)"}
THEME = {
    "DLP": "#3F8E8C",
    "DSE": "#E0884B",
    "observed": "#C0392B",
    "null_mean": "#333333",
}
plt.rcParams.update(
    {
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "xtick.color": "#444444",
        "ytick.color": "#444444",
    }
)


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
    col = THEME[tool]
    ax.hist(null, bins=40, color=col, alpha=0.6, edgecolor="white", linewidth=0.4)
    ax.axvline(observed, color=THEME["observed"], lw=2.2, label=f"observed = {observed}")
    ax.axvline(null_mean, color=THEME["null_mean"], lw=1.3, ls="--", label=f"null mean = {null_mean:.1f}")
    top = ax.get_ylim()[1]
    ax.annotate(
        f"{fold:.1f}x  {sig_stars(q)}",
        xy=(observed, top * 0.82),
        ha="center",
        fontsize=9,
        fontweight="bold",
        color=THEME["observed"],
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=THEME["observed"], alpha=0.9),
    )
    ax.legend(frameon=False, fontsize=7, loc="upper right")


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
