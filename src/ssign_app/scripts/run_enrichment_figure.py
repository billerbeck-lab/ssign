#!/usr/bin/env python3
"""Render the enrichment fold/significance bar chart.

Two figures share this script. The per-tool figure gives each secretion-system
type three bars (DeepLocPro, DeepSecE, SignalP) of fold enrichment (observed /
expected), annotated with the fold value and a BH q-value significance star. The
combined figure (``--combined``) gives one bar per type: DLP-or-DSE for non-T5
types and DLP-or-SignalP for all T5SS (T5a/b/c, drawn in its own colour, since T5
substrates are Sec-dependent and DSE is unreliable for T5). Window types
(T1/T2/T3/T4/T5b/T6) measure secreted-predicted proteins clustering near the
system's components; autotransporters (T5aSS/T5cSS) measure self-detection of the
component itself (tagged "(self)"). DeepSecE is not tested for T3SS.
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
from ssign_lib.constants import ENRICH_COMBINED_TOOL, ENRICH_T5SS_TYPES, ENRICH_TOOLS  # noqa: E402
from ssign_lib.plot_style import THEME, apply_house_style, muted, ordered_ss_types  # noqa: E402

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
    # Index by (ss_type, tool); keep only rows that actually carry an observed/fold
    # (skipped rows like DSE-T3SS have an empty observed and get no bar) and only
    # the tool(s) this figure plots.
    by_key = {}
    modes = {}
    for r in rows:
        if str(r.get("observed", "")).strip() == "" or r.get("tool") not in tools:
            continue
        by_key[(r["ss_type"], r["tool"])] = r
        modes[r["ss_type"]] = r.get("mode", "")

    ss_types = ordered_ss_types({k[0] for k in by_key})

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

    x = np.arange(len(ss_types))
    # Fit all bars of a group within ~0.85 of the unit x-spacing so adjacent
    # groups never overlap (the 3-bar per-tool figure needs the narrower width;
    # a single combined bar is capped so it doesn't get comically wide).
    width = min(0.55, 0.85 / len(tools))
    fig, ax = plt.subplots(figsize=(max(9, 1.7 * len(ss_types) + 2), 5.8))

    for i, tool in enumerate(tools):
        offset = (i - (len(tools) - 1) / 2) * width
        for xi, st in zip(x, ss_types):
            row = by_key.get((st, tool))
            if row is None:
                continue
            fold = _fnum(row.get("fold"))
            if np.isnan(fold):
                continue
            q = _fnum(row.get("qvalue"), 1.0)
            sig = q < 0.05
            height = inf_cap if np.isinf(fold) else fold
            # In the combined figure the T5SS bar is the DLP-or-SignalP score, not
            # DLP-or-DSE; colour all T5 (a/b/c) in the SignalP hue so it reads true.
            base = THEME["SignalP"] if (args.combined and st in ENRICH_T5SS_TYPES) else THEME[tool]
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
    ax.set_xticklabels([f"{st}\n(self)" if modes.get(st) == "self" else st for st in ss_types])
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
