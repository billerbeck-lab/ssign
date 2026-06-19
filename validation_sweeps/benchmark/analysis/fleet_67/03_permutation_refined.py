#!/usr/bin/env python3
"""Refined circular-permutation validation (three refinements over script 02).

Script 02 showed the binomial test is anti-conservative vs a spatial-permutation
null, but its permutation had a known weakness: the null windows could land on
*other* systems' substrate-rich neighborhoods, inflating the background secreted
density and making it a conservative floor. This addresses three asks:

(a) MASKED null — exclude every system's +/-3 neighborhood from the null window
    pool (loci from 03a_regen_neighborhoods.py). A system is now compared only
    against genuinely non-system background, the fair spatial test.

(b) PER SS TYPE — break the significant-call counts down by broad SS type, so we
    can see which system types the binomial over-calls relative to the spatial null.

(c) EXACT vs n=1000 permutation — the masked null computed over ALL clean windows
    (exact) vs estimated from 1000 randomly-sampled clean window placements. Same
    window width M, so counts are directly comparable; this is the spatial analogue
    of the n_null=200-vs-1000 background-sampling question. (NOT genome-thinning:
    thinned windows would span wider regions and make k incomparable.)

    .venv/bin/python 03_permutation_refined.py     # needs neighborhoods/ cached first

Approximations (same as 02): neighborhood treated as one contiguous M-window;
multi-contig genomes ordered contig-then-position and wrapped circularly.
"""

from __future__ import annotations

import csv
import json
import os
import random
import re
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ssign_app.scripts.enrichment_testing import (
    bh_fdr,
    binom_pvalue,
    broad_type,
    is_dlp_positive,
    is_dse_positive,
)

FLEET = "/tmp/ssign_fleet_67"
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figures")
NEIGH = os.path.join(HERE, "neighborhoods")
CONF = 0.8
SEED = 42
ALPHA = 0.05
N_SAMPLE = 1000  # sampled clean-window placements for the estimated permutation

THEME = {
    "perm_masked": "#2E7D6F",
    "perm_unmasked": "#9E9E9E",
    "binom_exact": "#3F8E8C",
    "binom_n1000": "#4C72B0",
    "perm_n1000": "#C58A2E",
    "ref": "#444444",
    "both": "#2E7D6F",
    "binom_only": "#C44E52",
    "perm_only": "#7A5C9E",
    "neither": "#C9C9C9",
}
plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlepad": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "xtick.color": "#444444",
        "ytick.color": "#444444",
    }
)


def ordered_genome(genome: str):
    """Return (dlp_vec, dse_vec, loci) in circular gene order; bools as int arrays."""
    raw = os.path.join(FLEET, genome, "results", f"{genome}_results_raw.csv")
    recs = []
    with open(raw) as fh:
        for row in csv.DictReader(fh):
            try:
                start = int(float(row.get("start") or 0))
            except ValueError:
                start = 0
            recs.append(
                (
                    row.get("contig", ""),
                    start,
                    row.get("locus_tag", ""),
                    is_dlp_positive(row, CONF),
                    is_dse_positive(row, CONF),
                )
            )
    recs.sort(key=lambda r: (r[0], r[1]))
    dlp = np.array([r[3] for r in recs], dtype=int)
    dse = np.array([r[4] for r in recs], dtype=int)
    loci = [r[2] for r in recs]
    return dlp, dse, loci


def neighborhood_mask(genome: str, loci: list[str]) -> np.ndarray:
    """Bool mask over the ordered vector: True where a gene is in ANY system's neighborhood."""
    path = os.path.join(NEIGH, f"{genome}.json")
    if not os.path.exists(path):
        return np.zeros(len(loci), dtype=bool)
    neigh = set(json.load(open(path))["neigh_ids"])
    return np.array([l in neigh for l in loci], dtype=bool)


def _clean_window_counts(vec: np.ndarray, mask: np.ndarray, M: int):
    """Secreted-counts of all circular M-windows that don't overlap the mask."""
    n = len(vec)
    if M <= 0 or M >= n:
        return np.array([], dtype=int)
    sec = np.cumsum(np.concatenate([[0], vec, vec[:M]]))  # prefix sums, len n+M+1
    msk = np.cumsum(np.concatenate([[0], mask.astype(int), mask[:M].astype(int)]))
    win_sec = sec[M : M + n] - sec[:n]
    win_msk = msk[M : M + n] - msk[:n]
    return win_sec[win_msk == 0]


def masked_perm(vec, mask, M, k, rng):
    """(p_exact, p_sampled, n_clean) for the masked spatial null, +1 smoothed.

    p_exact uses all clean (non-neighborhood) M-windows; p_sampled estimates the
    same from N_SAMPLE drawn without replacement. One shared window-count pass.
    """
    counts = _clean_window_counts(vec, mask, M)
    nc = len(counts)
    if nc == 0 or k <= 0:
        return 1.0, 1.0, nc
    p_exact = float((np.sum(counts >= k) + 1) / (nc + 1))
    draw = min(N_SAMPLE, nc)
    sub = counts[rng.sample(range(nc), draw)]
    p_sampled = float((np.sum(sub >= k) + 1) / (draw + 1))
    return p_exact, p_sampled, nc


def perm_unmasked(vec, M, k):
    """Original script-02 permutation (no mask), for the masking-effect comparison."""
    n = len(vec)
    if M <= 0 or M >= n or k <= 0:
        return 1.0
    cs = np.concatenate([[0], np.cumsum(np.concatenate([vec, vec[:M]]))])
    counts = cs[M : M + n] - cs[:n]
    return float((np.sum(counts >= k) + 1) / (n + 1))


def parse_systems(genome: str):
    """(ss_type, tool, M, k, sig_binom_exact) per system row in the genome summary."""
    summ = os.path.join(FLEET, genome, "results", f"{genome}_summary.txt")
    out = []
    if not os.path.exists(summ):
        return out
    with open(summ) as fh:
        for line in fh:
            if not re.match(r"\s*" + re.escape(genome) + r"\s+system\s", line):
                continue
            p = line.split()
            try:
                out.append(
                    {"ss_type": p[3], "tool": p[4], "M": int(p[5]), "k": int(p[6]), "sig_binom_exact": p[11] == "True"}
                )
            except (ValueError, IndexError):
                continue
    return out


def main():
    os.makedirs(FIGS, exist_ok=True)
    if not os.path.isdir(NEIGH) or not os.listdir(NEIGH):
        raise SystemExit(f"no neighborhood cache at {NEIGH}; run 03a_regen_neighborhoods.py first")
    rng = random.Random(SEED)
    rows = []
    small_clean = 0
    for g in sorted(os.listdir(FLEET)):
        systems = parse_systems(g)
        if not systems:
            continue
        dlp, dse, loci = ordered_genome(g)
        mask = neighborhood_mask(g, loci)
        vecs = {"DLP": dlp, "DSE": dse}
        grows = []
        for s in systems:
            vec = vecs.get(s["tool"])
            if vec is None or len(vec) == 0:
                continue
            n = len(vec)
            p_me, p_ms, nc = masked_perm(vec, mask, s["M"], s["k"], rng)
            if 0 < nc < 50:
                small_clean += 1
            p_pu = perm_unmasked(vec, s["M"], s["k"])
            # binomial with an n=1000 sampled background (default users get this)
            n_draw = min(1000, n)
            p_bg = vec[rng.sample(range(n), n_draw)].mean()
            p_bn = binom_pvalue(s["k"], s["M"], float(p_bg))
            grows.append(
                {
                    **s,
                    "genome": g,
                    "broad": broad_type(s["ss_type"]),
                    "n_clean": nc,
                    "p_perm_masked": p_me,
                    "p_perm_sampled": p_ms,
                    "p_perm_unmasked": p_pu,
                    "p_binom_n1000": p_bn,
                }
            )
        # BH within genome (matches production), one pass per p-source
        bh_fdr(grows, pvalue_key="p_perm_masked", q_key="q_pm", sig_key="sig_perm_masked")
        bh_fdr(grows, pvalue_key="p_perm_sampled", q_key="q_ps", sig_key="sig_perm_sampled")
        bh_fdr(grows, pvalue_key="p_perm_unmasked", q_key="q_pu", sig_key="sig_perm_unmasked")
        bh_fdr(grows, pvalue_key="p_binom_n1000", q_key="q_bn", sig_key="sig_binom_n1000")
        rows.extend(grows)

    report(rows, small_clean)
    fig_masked_vs_binomial(rows)
    fig_per_sstype(rows)
    fig_exact_vs_sampled(rows)
    print("\nFigure index:")
    print("  06  06_permutation_masked_vs_binomial.png  — refined spatial null vs binomial")
    print("  07  07_permutation_per_sstype.png          — significant calls per SS type + per tool (DLP/DSE)")
    print("  08  08_permutation_exact_vs_sampled.png    — exact vs n=1000-sampled permutation")


def _sig(rows, key):
    return sum(r[key] for r in rows)


def report(rows, small_clean):
    n = len(rows)
    ng = len({r["genome"] for r in rows})
    print(f"{n} per-system tests (DLP+DSE) across {ng} genomes\n")
    print("significant systems (q < 0.05) by method:")
    print(f"  permutation masked-exact   : {_sig(rows, 'sig_perm_masked')}")
    print(f"  permutation masked-n1000   : {_sig(rows, 'sig_perm_sampled')}")
    print(f"  permutation unmasked (02)  : {_sig(rows, 'sig_perm_unmasked')}")
    print(f"  binomial exact (production): {sum(r['sig_binom_exact'] for r in rows)}")
    print(f"  binomial n=1000 sampled    : {_sig(rows, 'sig_binom_n1000')}")
    if small_clean:
        print(
            f"\n  note: {small_clean} systems had <50 clean (non-neighborhood) windows — masked null is low-power there"
        )

    # per tool (DLP vs DSE): does the binomial over-call more for one predictor?
    print("\nsignificant by tool (masked-perm / binom-exact / binom-n1000):")
    for tool in ("DLP", "DSE"):
        tr = [r for r in rows if r["tool"] == tool]
        print(
            f"  {tool}: tested {len(tr):3d}  |  perm {_sig(tr, 'sig_perm_masked'):3d}  "
            f"binom-exact {sum(r['sig_binom_exact'] for r in tr):3d}  "
            f"binom-n1000 {_sig(tr, 'sig_binom_n1000'):3d}"
        )

    # per (SS type x tool)
    print("\nsignificant by SS type x tool (perm-masked / binom-exact):")
    pairs = sorted({(r["broad"], r["tool"]) for r in rows})
    for bt, tool in pairs:
        sub = [r for r in rows if r["broad"] == bt and r["tool"] == tool]
        pm, be = _sig(sub, "sig_perm_masked"), sum(r["sig_binom_exact"] for r in sub)
        if len(sub):
            print(f"  {bt:6s} {tool}: tested {len(sub):3d}  perm {pm:3d}  binom-exact {be:3d}")

    # masking effect: unmasked vs masked
    gained = sum(1 for r in rows if r["sig_perm_masked"] and not r["sig_perm_unmasked"])
    lost = sum(1 for r in rows if r["sig_perm_unmasked"] and not r["sig_perm_masked"])
    print("\nmasking effect (excluding other systems' neighborhoods from the null):")
    print(f"  newly significant after masking: {gained}   no-longer significant: {lost}")

    # exact vs sampled agreement
    both = sum(1 for r in rows if r["sig_perm_masked"] and r["sig_perm_sampled"])
    only_e = sum(1 for r in rows if r["sig_perm_masked"] and not r["sig_perm_sampled"])
    only_s = sum(1 for r in rows if r["sig_perm_sampled"] and not r["sig_perm_masked"])
    rec = both / (both + only_e) if (both + only_e) else float("nan")
    print(
        f"\nexact vs n=1000-sampled permutation: both {both}, exact-only {only_e}, "
        f"sampled-only {only_s} -> sampled recovers {rec:.0%} of exact-significant"
    )

    # binomial vs refined spatial null
    confirm = sum(1 for r in rows if r["sig_binom_exact"] and r["sig_perm_masked"])
    bex = sum(r["sig_binom_exact"] for r in rows)
    print(
        f"\nof {bex} production-significant systems, the refined (masked) permutation "
        f"confirms {confirm} ({confirm / bex:.0%})"
        if bex
        else ""
    )


def fig_masked_vs_binomial(rows):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))

    def cat(r):
        if r["sig_perm_masked"] and r["sig_binom_n1000"]:
            return "both"
        if r["sig_binom_n1000"]:
            return "binom_only"
        if r["sig_perm_masked"]:
            return "perm_only"
        return "neither"

    for c in ("neither", "both", "binom_only", "perm_only"):
        pts = [r for r in rows if cat(r) == c]
        if not pts:
            continue
        a1.scatter(
            [max(r["p_binom_n1000"], 1e-6) for r in pts],
            [max(r["p_perm_masked"], 1e-6) for r in pts],
            s=18,
            color=THEME[c],
            alpha=0.7,
            label=f"{c.replace('_', '-')} ({len(pts)})",
        )
    a1.plot([1e-6, 1], [1e-6, 1], ls=":", color=THEME["ref"], lw=0.8)
    a1.axhline(ALPHA, ls="--", color=THEME["ref"], lw=0.6)
    a1.axvline(ALPHA, ls="--", color=THEME["ref"], lw=0.6)
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlabel("binomial (n=1000) p-value")
    a1.set_ylabel("masked circular-permutation p-value")
    a1.set_title("Per-system: binomial vs refined spatial null")
    a1.legend(frameon=False, fontsize=8, loc="lower right")

    labels = ["permutation\nmasked-exact", "permutation\nunmasked (02)", "binomial\nexact bg", "binomial\nn=1000 bg"]
    vals = [
        _sig(rows, "sig_perm_masked"),
        _sig(rows, "sig_perm_unmasked"),
        sum(r["sig_binom_exact"] for r in rows),
        _sig(rows, "sig_binom_n1000"),
    ]
    colors = [THEME["perm_masked"], THEME["perm_unmasked"], THEME["binom_exact"], THEME["binom_n1000"]]
    a2.bar(labels, vals, color=colors)
    for i, v in enumerate(vals):
        a2.text(i, v, str(v), ha="center", va="bottom", weight="bold")
    a2.set_ylabel("significant systems (q < 0.05)")
    a2.set_title("Significant systems by test method")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "06_permutation_masked_vs_binomial.png"))
    plt.close(fig)


def _grouped_sig_bars(ax, groups, label_of, title):
    """Grouped bars: one group per key, 3 method bars (perm-masked, binom-exact, binom-n1000)."""
    methods = [
        ("sig_perm_masked", "perm masked-exact", THEME["perm_masked"]),
        ("sig_binom_exact", "binomial exact", THEME["binom_exact"]),
        ("sig_binom_n1000", "binomial n=1000", THEME["binom_n1000"]),
    ]
    keys = list(groups)
    x = np.arange(len(keys))
    w = 0.26
    for i, (key, lab, col) in enumerate(methods):
        vals = [sum(r[key] for r in groups[k]) for k in keys]
        ax.bar(x + (i - 1) * w, vals, w, color=col, label=lab)
        for xi, v in zip(x + (i - 1) * w, vals):
            if v:
                ax.text(xi, v, str(v), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([label_of(k, groups[k]) for k in keys])
    ax.set_ylabel("significant systems (q < 0.05)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)


def fig_per_sstype(rows):
    by_type = defaultdict(list)
    by_tool = defaultdict(list)
    for r in rows:
        by_type[r["broad"]].append(r)
        by_tool[r["tool"]].append(r)
    by_type = dict(sorted(by_type.items(), key=lambda kv: -len(kv[1])))
    by_tool = {t: by_tool[t] for t in ("DLP", "DSE") if t in by_tool}

    fig, (a1, a2) = plt.subplots(
        1,
        2,
        figsize=(max(11, 1.2 * len(by_type) + 5), 5.2),
        gridspec_kw={"width_ratios": [max(2.2, len(by_type) * 0.7), 1.2]},
    )
    _grouped_sig_bars(
        a1, by_type, lambda k, v: f"{k}\n(n={len(v)})", "Significant systems per SS type: spatial null vs binomial"
    )
    _grouped_sig_bars(a2, by_tool, lambda k, v: f"{k}\n(n={len(v)})", "Per predictor (DLP vs DSE)")
    a2.get_legend().remove()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "07_permutation_per_sstype.png"))
    plt.close(fig)


def fig_exact_vs_sampled(rows):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))

    def cat(r):
        if r["sig_perm_masked"] and r["sig_perm_sampled"]:
            return "both"
        if r["sig_perm_masked"]:
            return "exact_only"
        if r["sig_perm_sampled"]:
            return "sampled_only"
        return "neither"

    colors = {
        "both": THEME["perm_masked"],
        "exact_only": THEME["perm_only"],
        "sampled_only": THEME["perm_n1000"],
        "neither": THEME["neither"],
    }
    a1.plot([1e-4, 1], [1e-4, 1], ls=":", color=THEME["ref"], lw=0.8)
    for c, lab in (
        ("neither", "neither"),
        ("both", "both"),
        ("exact_only", "exact-only"),
        ("sampled_only", "sampled-only"),
    ):
        pts = [r for r in rows if cat(r) == c]
        if not pts:
            continue
        a1.scatter(
            [max(r["p_perm_masked"], 1e-4) for r in pts],
            [max(r["p_perm_sampled"], 1e-4) for r in pts],
            s=18,
            color=colors[c],
            alpha=0.7,
            label=f"{lab} ({len(pts)})",
        )
    a1.axhline(ALPHA, ls="--", color=THEME["ref"], lw=0.6)
    a1.axvline(ALPHA, ls="--", color=THEME["ref"], lw=0.6)
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlabel("exact masked-permutation p-value")
    a1.set_ylabel("n=1000-sampled permutation p-value")
    a1.set_title("Exact vs sampled permutation (masked null)")
    a1.legend(frameon=False, fontsize=8, loc="upper left")

    labels = ["exact\n(all clean windows)", "n=1000\nsampled windows"]
    vals = [_sig(rows, "sig_perm_masked"), _sig(rows, "sig_perm_sampled")]
    a2.bar(labels, vals, color=[THEME["perm_masked"], THEME["perm_n1000"]])
    for i, v in enumerate(vals):
        a2.text(i, v, str(v), ha="center", va="bottom", weight="bold")
    a2.set_ylabel("significant systems (q < 0.05)")
    a2.set_title("Significant systems: exact vs sampled")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "08_permutation_exact_vs_sampled.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
