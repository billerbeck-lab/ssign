#!/usr/bin/env python3
"""Render the enrichment fold/significance bar chart.

Two figures share this script. The per-tool figure gives each secretion-system
type three bars (DeepLocPro, DeepSecE, SignalP) of fold enrichment (observed /
expected), annotated with the fold value and a BH q-value significance star. The
combined figure (``--combined``) gives one bar per type: DLP-or-DSE for ordinary
window types and DLP-or-SignalP for the Sec-dependent T5 results (autotransporter
self-detection + T5bSS, drawn in the SignalP colour). Window types
(T1/T2/T3/T4/T5b/T6) measure secreted-predicted proteins clustering near the
system's components. Autotransporters (T5aSS/T5cSS) render TWO adjacent x-groups:
"(self)" (the component detecting itself) and "(hitchhiker)" (the hitchhiker window,
secreted-predicted neighbours that may piggyback through the T5 pore, scored like
any window type: DLP-or-DSE combined). DeepSecE is not tested for T3SS.
Non-significant bars (q >= 0.05) keep their colour but are faded so the eye lands
on the real enrichments.

Reads only the per-type stats TSV from enrichment_testing.py (fold + qvalue per
row); the null arrays are no longer needed for the figure. Identical figure at
single-genome and pooled (multi-genome) scale.

    run_enrichment_figure.py --stats <sid>_enrichment_stats.tsv --out <sid>_enrichment_fold.png
"""

import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import (  # noqa: E402
    ENRICH_COMBINED_TOOL,
    ENRICH_MODE_SELF,
    ENRICH_MODE_WINDOW,
    ENRICH_TOOLS,
    T5_HITCH_TAG,
    T5_SELF_TAG,
    enrich_combined_uses_signalp,
)
from ssign_lib.plot_style import THEME, apply_house_style, muted, ordered_ss_types  # noqa: E402

apply_house_style()


def sig_stars(q: float) -> str:
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "n.s."


def build_columns(ss_types, col_modes):
    """Ordered x-columns as (ss_type, mode, label). A type carrying both modes
    (an autotransporter with a self and a hitchhiker window result) yields two
    adjacent columns, self first then hitchhiker; every other type yields one."""
    columns = []
    for st in ss_types:
        ms = col_modes.get(st, {ENRICH_MODE_WINDOW})
        if ENRICH_MODE_SELF in ms and ENRICH_MODE_WINDOW in ms:
            columns.append((st, ENRICH_MODE_SELF, f"{st}\n{T5_SELF_TAG}"))
            columns.append((st, ENRICH_MODE_WINDOW, f"{st}\n{T5_HITCH_TAG}"))
        elif ENRICH_MODE_SELF in ms:
            columns.append((st, ENRICH_MODE_SELF, f"{st}\n{T5_SELF_TAG}"))
        else:
            columns.append((st, ENRICH_MODE_WINDOW, st))
    return columns


def load_stats(path: str) -> list:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fnum(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main():
    ap = argparse.ArgumentParser(description="Enrichment fold/significance bar chart")
    ap.add_argument("--stats", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--title", default="")
    ap.add_argument(
        "--combined",
        action="store_true",
        help="Plot the single combined bar per type (DLP-or-DSE; DLP-or-SignalP for T5SS) instead of "
        "the per-tool DLP/DSE/SignalP bars.",
    )
    args = ap.parse_args()

    tools = [ENRICH_COMBINED_TOOL] if args.combined else list(ENRICH_TOOLS)
    title = args.title or (
        "Combined enrichment significance by secretion-system type (DLP or DSE; DLP or SignalP for T5SS)"
        if args.combined
        else "Enrichment significance by secretion-system type"
    )

    rows = load_stats(args.stats)
    # Index by (ss_type, tool, mode); keep only rows that actually carry an
    # observed/fold (skipped rows like DSE-T3SS have an empty observed and get no
    # bar) and only the tool(s) this figure plots. Autotransporters (T5a/c) carry
    # two modes (self + hitchhiker window), each its own x-group.
    by_key = {}
    col_modes: dict[str, set[str]] = {}  # ss_type -> set of modes present
    for r in rows:
        if str(r.get("observed", "")).strip() == "" or r.get("tool") not in tools:
            continue
        m = r.get("mode", ENRICH_MODE_WINDOW) or ENRICH_MODE_WINDOW
        by_key[(r["ss_type"], r["tool"], m)] = r
        col_modes.setdefault(r["ss_type"], set()).add(m)

    ss_types = ordered_ss_types(set(col_modes))
    # One x-column per (ss_type, mode); T5a/c with both modes get two adjacent ones.
    columns = build_columns(ss_types, col_modes)

    if not ss_types:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No secretion systems detected — nothing to plot", ha="center", va="center")
        ax.axis("off")
        fig.savefig(args.out, dpi=args.dpi)
        plt.close(fig)
        print(f"Saved (empty): {args.out}")
        return

    # Cap the display height of an infinite fold (null mean 0, observed > 0) at a
    # little above the tallest finite bar, so it doesn't blow the axis; the exact
    # fold is in the annotation. inf only happens at single-genome scale.
    folds = [_fnum(r.get("fold")) for r in by_key.values()]
    max_finite = max((f for f in folds if np.isfinite(f) and f > 0), default=1.0)
    inf_cap = max_finite * 1.15

    x = np.arange(len(columns))
    # Fit all bars of a group within ~0.85 of the unit x-spacing so adjacent
    # groups never overlap (the 3-bar per-tool figure needs the narrower width;
    # a single combined bar is capped so it doesn't get comically wide).
    width = min(0.55, 0.85 / len(tools))
    fig, ax = plt.subplots(figsize=(max(9, 1.7 * len(columns) + 2), 5.8))

    for i, tool in enumerate(tools):
        offset = (i - (len(tools) - 1) / 2) * width
        for xi, (st, mode, _label) in zip(x, columns):
            row = by_key.get((st, tool, mode))
            if row is None:
                continue
            fold = _fnum(row.get("fold"))
            if np.isnan(fold):
                continue
            q = _fnum(row.get("qvalue"), 1.0)
            sig = q < 0.05
            height = inf_cap if np.isinf(fold) else fold
            # In the combined figure the Sec-dependent T5 bars (autotransporter self
            # + T5bSS) are DLP-or-SignalP, so colour them in the SignalP hue; the
            # T5a/c hitchhiker window is DLP-or-DSE and keeps the COMBINED hue.
            base = THEME["SignalP"] if (args.combined and enrich_combined_uses_signalp(st, mode)) else THEME[tool]
            # Non-significant bars keep their hue but are desaturated/faded (so you
            # can still tell which tool it is), rather than flattened to grey.
            color = base if sig else muted(base)
            ax.bar(xi + offset, height, width, color=color, zorder=3)
            label = "∞" if np.isinf(fold) else f"{fold:.0f}x" if fold >= 10 else f"{fold:.1f}x"
            ax.text(
                xi + offset,
                height,
                f"{label}\n{sig_stars(q)}",
                ha="center",
                va="bottom",
                fontsize=7,
                fontweight="bold",
                color=base if sig else muted(base, sat_scale=0.55, val_boost=0.9),
            )

    ax.axhline(1, color=THEME["ref_line"], ls="--", lw=1, alpha=0.7, zorder=1, label="no enrichment")
    ax.set_xticks(x)
    ax.set_xticklabels([label for _st, _mode, label in columns])
    ax.set_ylabel("fold enrichment (observed / expected)")
    ax.set_ylim(0, inf_cap * 1.18)
    ax.set_title(title)
    if args.combined:
        legend_handles = [
            Patch(color=THEME["COMBINED"], label="DLP or DSE"),
            Patch(color=THEME["SignalP"], label="DLP or SignalP (T5SS)"),
            Patch(color=muted(THEME["COMBINED"]), label="faded = not significant (q ≥ 0.05)"),
        ]
    else:
        legend_handles = [
            Patch(color=THEME["DLP"], label="DeepLocPro"),
            Patch(color=THEME["DSE"], label="DeepSecE"),
            Patch(color=THEME["SignalP"], label="SignalP"),
            Patch(color=muted(THEME["DLP"]), label="faded = not significant (q ≥ 0.05)"),
        ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
