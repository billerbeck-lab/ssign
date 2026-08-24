#!/usr/bin/env python3
"""Draw the recall-by-SS-type figure on docs/benchmark/benchmarks.md.

The figure predates this script and was originally produced ad hoc, which meant
every edit to the benchmarking list left it silently stale. It is regenerated
from the list itself, so the bars, the per-bar counts and the title's total all
move together.

Three stacked bands per type, matching the list's own columns:
  predicted                 found_by_ssign
  reachable, not predicted  reachable_within_3 in {yes, self}, but not found
  unreachable               own machinery further than the proximity window

`self` counts as reachable: for the T5aSS/T5cSS autotransporters the passenger
and the translocator are one protein, so the component is its own neighbour and
proximity has nothing to reach for.

Usage: recall_figure.py [--out docs/benchmark/secreted_protein_prediction_by_type.png]
"""

import argparse
import collections
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from _paths import BENCHMARK_CSV, ROOT  # noqa: E402

PREDICTED = "#2e7d6f"
REACHABLE = "#e0a33b"
UNREACHABLE = "#c7cdd4"
REACHABLE_VALUES = {"yes", "self"}
TYPE_ORDER = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS"]


def truthy(v):
    return v.strip().lower() in ("yes", "true", "1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "benchmark", "secreted_protein_prediction_by_type.png"))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(BENCHMARK_CSV)))
    counts = collections.defaultdict(collections.Counter)
    for r in rows:
        c = counts[r["ss_type"]]
        c["n"] += 1
        if truthy(r["found_by_ssign"]):
            c["predicted"] += 1
        elif r["reachable_within_3"].strip().lower() in REACHABLE_VALUES:
            c["reachable"] += 1
        else:
            c["unreachable"] += 1

    types = [t for t in TYPE_ORDER if t in counts] + sorted(set(counts) - set(TYPE_ORDER))
    n_total = sum(counts[t]["n"] for t in types)
    n_pred = sum(counts[t]["predicted"] for t in types)

    fig, ax = plt.subplots(figsize=(11.2, 5.6))
    x = range(len(types))
    bottoms = [0] * len(types)
    for key, colour, label in (
        ("predicted", PREDICTED, "predicted"),
        ("reachable", REACHABLE, "reachable, not predicted"),
        ("unreachable", UNREACHABLE, "unreachable (own machinery >3 genes away)"),
    ):
        vals = [counts[t][key] for t in types]
        ax.bar(x, vals, 0.62, bottom=bottoms, color=colour, label=label)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v:
                ax.text(
                    xi,
                    b + v / 2,
                    str(v),
                    ha="center",
                    va="center",
                    fontsize=13,
                    fontweight="bold",
                    color="white" if key != "unreachable" else "#444444",
                )
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    for xi, t in enumerate(types):
        ax.text(
            xi,
            counts[t]["n"] + 0.45,
            f"{counts[t]['predicted']}/{counts[t]['n']}",
            ha="center",
            va="bottom",
            fontsize=12,
            color="#444444",
        )

    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{t}\n(n={counts[t]['n']})" for t in types], fontsize=13)
    ax.set_ylabel("secreted proteins", fontsize=13)
    ax.set_title(
        f"Secreted-protein prediction by secretion-system type: {n_pred}/{n_total} predicted",
        fontsize=15,
        fontweight="bold",
        pad=16,
    )
    ax.set_ylim(0, max(counts[t]["n"] for t in types) + 3)
    ax.tick_params(axis="y", labelsize=12)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False, fontsize=12)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote {args.out}  ({n_pred}/{n_total} predicted)")
    for t in types:
        c = counts[t]
        print(f"  {t}: n={c['n']} predicted={c['predicted']} reachable={c['reachable']} unreachable={c['unreachable']}")


if __name__ == "__main__":
    main()
