#!/usr/bin/env python3
"""Circular-shift enrichment on the 67-genome fleet, with actual fold values.

Reproduces the Xanthobacter-era figure4/figure5 style (genome-structure-
preserving circular-shift null + per-type FOLD enrichment) on the fleet, and
fixes the T5SS handling the per-system binomial got wrong:

  - T5aSS / T5cSS are AUTOTRANSPORTERS: the protein is both machinery and
    substrate (self-secreting). A +/-W neighbour window is meaningless (mostly
    empty). Instead we report SELF-DETECTION: is the autotransporter component
    itself predicted Outer-Membrane-or-Extracellular (DLP) / a secreted type
    (DSE)? (t5ss_aware_analysis.py SECTION A.)
  - T5bSS (two-partner) + T1/T2/T3/T4/T6: window types. The +/-W window finds
    the TpsA hitchhiker (T5b) / nearby effector substrates.

Method (per tool, per scope = genome-wide or one SS type):
  observed = # tool-positive proteins within +/-W of any in-scope SS component,
             summed over genomes.
  null     = circularly rotate each genome's tool-positive positions by a random
             offset (preserves their clustering), recount, sum over genomes.
             The exact all-rotations count per genome is the circular cross-
             correlation of the positive-indicator and the window mask (FFT);
             we then pool random rotations across genomes.
  fold     = observed / null-mean ;  p = (#null >= observed + 1)/(reps + 1).
BH-correct the per-type p-values across types (within tool).

    .venv/bin/python 04_circular_shift_enrichment.py     # needs neighborhoods/ (with components)

Approximation: multi-contig genomes are ordered contig-then-start and rotated as
one circular replicon (a rotation can wrap across a contig junction).
"""

from __future__ import annotations

import csv
import json
import os
import re
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ssign_app.scripts.enrichment_testing import bh_fdr, is_dlp_positive, is_dse_positive

FLEET = "/tmp/ssign_fleet_67"
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figures")
NEIGH = os.path.join(HERE, "neighborhoods")
CONF = 0.8
WINDOW = 3
REPS = 10_000
SEED = 42

# substrate-secreting systems where the +/-W window finds nearby substrates
WINDOW_TYPES = {"T1SS", "T2SS", "T3SS", "T4SS", "pT4SSt", "T5bSS", "T6SSi", "T6SSii"}
# autotransporters: the component IS the substrate -> self-detection, not window
AUTOTRANSPORTER_TYPES = {"T5aSS", "T5cSS"}
# DSE cannot call T3SS (excluded in ssign), so it gets no T3SS window
DSE_NO_WINDOW = {"T3SS"}

THEME = {
    "DLP": "#3F8E8C",
    "DSE": "#E0884B",
    "observed": "#C0392B",
    "null_mean": "#333333",
    "ns": "#BDBDBD",
    "self_hit": "#2E7D6F",
    "self_miss": "#E0E0E0",
}
plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlepad": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "xtick.color": "#444444",
        "ytick.color": "#444444",
    }
)


def sig_stars(p: float) -> str:
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def display_type(ss: str) -> str:
    """Keep T5 subtypes distinct; collapse T6SSi/ii -> T6SS, pT4SSt -> T4SS."""
    if ss.startswith("T5"):
        return ss
    m = re.match(r"p?(T\d+)[a-z]*SS", ss)
    return f"{m.group(1)}SS" if m else ss


def load_genome(genome: str):
    """Ordered per-protein arrays + component (index, type), or None if no cache."""
    cache = os.path.join(NEIGH, f"{genome}.json")
    if not os.path.exists(cache):
        return None
    comps = json.load(open(cache)).get("components")
    if comps is None:
        return None
    recs = []
    with open(os.path.join(FLEET, genome, "results", f"{genome}_results_raw.csv")) as fh:
        for row in csv.DictReader(fh):
            try:
                start = int(float(row.get("start") or 0))
            except ValueError:
                start = 0
            recs.append((row.get("contig", ""), start, row))
    recs.sort(key=lambda r: (r[0], r[1]))
    n = len(recs)
    idx = {r[2].get("locus_tag", ""): i for i, r in enumerate(recs)}
    dlp = np.array([is_dlp_positive(r[2], CONF) for r in recs], dtype=float)
    dse = np.array([is_dse_positive(r[2], CONF) for r in recs], dtype=float)

    # autotransporter self-detection signal: OM>=0.8 OR Extracellular>=0.8
    def om_or_extra(row):
        try:
            return (
                float(row.get("outer_membrane_prob") or 0) >= CONF
                or float(row.get("dlp_extracellular_prob") or 0) >= CONF
            )
        except ValueError:
            return False

    dlp_self = np.array([om_or_extra(r[2]) for r in recs], dtype=float)
    comp_by_type = defaultdict(list)
    for c in comps:
        if c["locus_tag"] in idx:
            comp_by_type[c["ss_type"]].append(idx[c["locus_tag"]])
    return {"genome": genome, "n": n, "dlp": dlp, "dse": dse, "dlp_self": dlp_self, "comp_by_type": comp_by_type}


def window_mask(comp_idx, n, w=WINDOW):
    """Boolean length-n mask: positions within +/-w of any component (circular)."""
    mask = np.zeros(n, dtype=float)
    for p in comp_idx:
        for d in range(-w, w + 1):
            mask[(p + d) % n] = 1.0
    return mask


def rotation_counts(pos_vec, win_mask):
    """All-rotations count array c[k] = # positives landing in window after shift k.

    Circular cross-correlation via FFT; c[0] == the un-rotated (observed) count.
    """
    n = len(pos_vec)
    c = np.fft.irfft(np.fft.rfft(win_mask) * np.conj(np.fft.rfft(pos_vec)), n)
    return np.rint(c).astype(int)


def pooled_test(per_genome_counts, rng):
    """observed (sum c[0]) + pooled null over REPS random rotations; fold, p, null."""
    observed = int(sum(c[0] for c in per_genome_counts))
    null = np.zeros(REPS, dtype=int)
    for c in per_genome_counts:
        n = len(c)
        if n <= 1:
            continue
        null += c[rng.integers(1, n, size=REPS)]  # shift in 1..n-1 (exclude identity)
    mean = null.mean() if len(null) else 0.0
    fold = observed / mean if mean > 0 else (float("inf") if observed > 0 else 0.0)
    p = (np.sum(null >= observed) + 1) / (REPS + 1)
    return {"observed": observed, "null": null, "null_mean": float(mean), "fold": float(fold), "pvalue": float(p)}


def scope_counts(genomes, tool, types):
    """Per-genome rotation-count arrays for `tool` positives near components of `types`."""
    out = []
    for g in genomes:
        comp_idx = [i for t in types for i in g["comp_by_type"].get(t, [])]
        if not comp_idx:
            continue
        out.append(rotation_counts(g[tool], window_mask(comp_idx, g["n"])))
    return out


def main():
    os.makedirs(FIGS, exist_ok=True)
    rng = np.random.default_rng(SEED)
    genomes = [x for x in (load_genome(g) for g in sorted(os.listdir(FLEET))) if x]
    have_comp = sum(1 for g in genomes if g["comp_by_type"])
    print(f"loaded {len(genomes)} genomes ({have_comp} with >=1 component)\n")

    # ---- genome-wide circular shift (window types; DSE drops T3SS) ----
    dlp_types = WINDOW_TYPES
    dse_types = WINDOW_TYPES - DSE_NO_WINDOW
    gw_dlp = pooled_test(scope_counts(genomes, "dlp", dlp_types), rng)
    gw_dse = pooled_test(scope_counts(genomes, "dse", dse_types), rng)
    print("=== genome-wide enrichment near window-type SS components ===")
    for lab, r in (("DLP extracellular", gw_dlp), ("DSE substrate (excl T3SS)", gw_dse)):
        print(
            f"  {lab:26s}: observed {r['observed']:4d}  null {r['null_mean']:6.1f}  "
            f"fold {r['fold']:5.1f}x  p {r['pvalue']:.2e}"
        )

    # ---- per-type fold ----
    types = sorted(
        {
            display_type(t)
            for g in genomes
            for t in g["comp_by_type"]
            if display_type(t) in {display_type(x) for x in WINDOW_TYPES}
        }
    )
    per_type = {}
    for tool in ("dlp", "dse"):
        rows = []
        for dt in types:
            if tool == "dse" and dt in DSE_NO_WINDOW:
                rows.append(
                    {
                        "display": dt,
                        "tool": tool,
                        "observed": 0,
                        "null_mean": 0.0,
                        "fold": 0.0,
                        "pvalue": 1.0,
                        "n_sys": 0,
                        "skip": True,
                    }
                )
                continue
            raw_types = {t for g in genomes for t in g["comp_by_type"] if display_type(t) == dt}
            counts = scope_counts(genomes, tool, raw_types)
            n_sys = sum(len(g["comp_by_type"].get(t, [])) > 0 for g in genomes for t in raw_types)
            r = pooled_test(counts, rng)
            r.pop("null")
            rows.append({"display": dt, "tool": tool, "n_sys": n_sys, "skip": False, **r})
        # BH only over real tests; skipped (T3SS-DSE) rows must not pad the denominator
        bh_fdr([r for r in rows if not r["skip"]], pvalue_key="pvalue", q_key="qvalue", sig_key="significant")
        for r in rows:
            if r["skip"]:
                r["qvalue"], r["significant"] = 1.0, False
        per_type[tool] = rows

    print("\n=== per-type fold enrichment (circular shift, BH-corrected) ===")
    print(f"  {'type':7s} {'tool':4s} {'obs':>4s} {'null':>6s} {'fold':>6s} {'q':>9s}  sig")
    for tool in ("dlp", "dse"):
        for r in per_type[tool]:
            if r.get("skip"):
                continue
            print(
                f"  {r['display']:7s} {tool.upper():4s} {r['observed']:4d} {r['null_mean']:6.1f} "
                f"{r['fold']:5.1f}x {r['qvalue']:9.2e}  {sig_stars(r['qvalue'])}"
            )

    # ---- autotransporter self-detection (T5aSS / T5cSS) ----
    self_rows = []
    for at in sorted(AUTOTRANSPORTER_TYPES):
        n_comp = dlp_hit = dse_hit = both = 0
        for g in genomes:
            for i in g["comp_by_type"].get(at, []):
                n_comp += 1
                d = bool(g["dlp_self"][i])
                s = bool(g["dse"][i])
                dlp_hit += d
                dse_hit += s
                both += d and s
        self_rows.append({"type": at, "n": n_comp, "dlp": dlp_hit, "dse": dse_hit, "both": both})
    print("\n=== autotransporter self-detection (component itself predicted secreted) ===")
    for r in self_rows:
        if r["n"]:
            print(
                f"  {r['type']}: {r['n']} components | DLP OM/extra {r['dlp']} ({r['dlp'] / r['n']:.0%})  "
                f"DSE {r['dse']} ({r['dse'] / r['n']:.0%})  both {r['both']} ({r['both'] / r['n']:.0%})"
            )

    fig_genomewide_null(gw_dlp, gw_dse)
    fig_per_type(per_type, types)
    fig_self_detection(self_rows)
    print("\nFigure index:")
    print("  09  09_circular_shift_null.png        — genome-wide circular-shift null (DLP, DSE)")
    print("  10  10_per_type_fold_enrichment.png   — per-type fold (window types, T5 subtypes split)")
    print("  11  11_autotransporter_self_detect.png— T5aSS/T5cSS self-detection rates")


def fig_genomewide_null(gw_dlp, gw_dse):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    panels = [
        (axes[0], gw_dlp, THEME["DLP"], "DeepLocPro — extracellular near SS"),
        (axes[1], gw_dse, THEME["DSE"], "DeepSecE (excl. T3SS) — substrate near SS"),
    ]
    for ax, r, col, title in panels:
        null = r["null"]
        ax.hist(
            null,
            bins=40,
            color=col,
            alpha=0.6,
            edgecolor="white",
            linewidth=0.5,
            label=f"null ({REPS:,} circular shifts)",
        )
        ax.axvline(r["observed"], color=THEME["observed"], lw=2.5, label=f"observed = {r['observed']}")
        ax.axvline(r["null_mean"], color=THEME["null_mean"], lw=1.5, ls="--", label=f"null mean = {r['null_mean']:.1f}")
        top = ax.get_ylim()[1]
        ax.fill_betweenx(
            [0, top * 1.1],
            r["observed"],
            max(null.max(), r["observed"]) + 5,
            alpha=0.13,
            color=THEME["observed"],
            label=f"p < {max(r['pvalue'], 1e-4):.0e}",
        )
        ax.annotate(
            f"{r['fold']:.1f}x enrichment",
            xy=(r["observed"], top * 0.85),
            ha="center",
            fontsize=12,
            fontweight="bold",
            color=THEME["observed"],
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=THEME["observed"], alpha=0.9),
        )
        ax.set_xlabel("count of positive proteins near SS components")
        ax.set_ylabel("frequency (permutations)")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle(
        "Circular-Shift Permutation (genome-structure-preserving), 67-genome fleet",
        fontsize=15,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "09_circular_shift_null.png"))
    plt.close(fig)


def fig_per_type(per_type, types):
    fig, ax = plt.subplots(figsize=(max(9, 1.3 * len(types) + 2), 5.6))
    x = np.arange(len(types))
    w = 0.4
    for i, tool in enumerate(("dlp", "dse")):
        by = {r["display"]: r for r in per_type[tool]}
        folds = [by[t]["fold"] if not by[t].get("skip") else 0 for t in types]
        for xi, t in zip(x, types):
            r = by[t]
            if r.get("skip"):
                continue
            col = THEME[tool.upper()] if r["qvalue"] < 0.05 else THEME["ns"]
            bx = xi + (i - 0.5) * w
            ax.bar(bx, r["fold"], w, color=col)
            ax.text(
                bx,
                r["fold"],
                f"{r['fold']:.1f}x\n{sig_stars(r['qvalue'])}",
                ha="center",
                va="bottom",
                fontsize=8,
                fontweight="bold",
                color=THEME[tool.upper()] if r["qvalue"] < 0.05 else "#999999",
            )
    ax.axhline(1, color="#999999", ls="--", lw=1, alpha=0.7, label="no enrichment")
    ax.set_xticks(x)
    ax.set_xticklabels(types)
    ax.set_ylabel("fold enrichment (observed / circular-shift null mean)")
    ax.set_title(
        "Per-type fold enrichment near SS components (window types)\n"
        "left bar DLP, right bar DSE; T5a/T5c are self-secreting (see self-detection fig)"
    )
    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(color=THEME["DLP"], label="DLP (sig)"),
            Patch(color=THEME["DSE"], label="DSE (sig)"),
            Patch(color=THEME["ns"], label="n.s."),
        ],
        frameon=False,
        fontsize=9,
        loc="upper right",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "10_per_type_fold_enrichment.png"))
    plt.close(fig)


def fig_self_detection(self_rows):
    rows = [r for r in self_rows if r["n"]]
    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(len(rows))
    w = 0.38
    for i, (key, lab, col) in enumerate(
        (("dlp", "DLP (OM or extracellular)", THEME["DLP"]), ("dse", "DSE (secreted type)", THEME["DSE"]))
    ):
        rates = [r[key] / r["n"] * 100 for r in rows]
        bx = x + (i - 0.5) * w
        ax.bar(bx, rates, w, color=col, label=lab)
        for xi, r, rate in zip(bx, rows, rates):
            ax.text(xi, rate, f"{rate:.0f}%\n{r[key]}/{r['n']}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['type']}\n({r['n']} components)" for r in rows])
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of autotransporter components self-detected")
    ax.set_title("T5aSS / T5cSS autotransporter self-detection\n(the component IS the substrate)")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "11_autotransporter_self_detect.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
