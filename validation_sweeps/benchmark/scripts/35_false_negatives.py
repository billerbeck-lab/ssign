#!/usr/bin/env python3
"""Phase 2: WHY does ssign miss effectors it could have reached? (false-negative failure modes.)

The "reachable @±3 but missed" effectors are ssign's true false negatives: the literature machinery is
within ±3 genes, yet ssign did not emit them. This script classifies each by failure mode, using the
per-protein tool signals in the actual-call table:

  detection_miss     - ssign detected NO secretion system near the effector and never ran the
                       secreted-protein predictors on it (all tool signals blank). The dominant mode.
  processed_not_emit - the predictors DID run (signals present) but cross-validation / filtering did
                       not emit it.

T1SS now has 0 false negatives: its RTX-toxin misses (HlyA, ApxIA, LtxA, LktA) were a benchmark
staging artifact, each effector was staged on the single plasmid/contig carrying its operon, which
lacks the chromosomal TolC (OMF). The TXSScan T1SS model marks OMF as a loner, so a distant TolC is
fine; staging the complete assembly (script 50) restores it and all four detect a T1SS. The clean_dataset
filter folds those recoveries in. Remaining misses are genuine detection failures (mostly T3SS/T4SS).

Inputs : data/phase2/actual_per_effector.panel_genbank_{default,t3ss}.tsv  (T3SS uses the t3ss tag)
Output : data/phase2/figures/summary/05_false_negative_modes.png
Run:     <repo>/.venv/bin/python scripts/35_false_negatives.py
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import clean_dataset  # noqa: E402,F401

BENCH = Path(__file__).resolve().parents[1]
FIGDIR = BENCH / "data" / "phase2" / "figures" / "summary"
TYPES = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]
SIG_COLS = (
    "dlp_extracellular_prob",
    "dse_ss_type",
    "signalp_prediction",
    "plm_effector_secreted",
    "predicted_localization",
)
C = {"detection": "#A93232", "processed": "#E0A33B"}
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
    }
)


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    ad = clean_dataset.load_clean_actual(BENCH / "data" / "phase2" / "actual_per_effector.panel_genbank_default.tsv")
    a3 = clean_dataset.load_clean_actual(BENCH / "data" / "phase2" / "actual_per_effector.panel_genbank_t3ss.tsv")
    per = {}
    for ss in TYPES:
        src = a3 if ss == "T3SS" else ad
        miss = [
            r
            for r in src
            if r["ss_type"] == ss
            and r["testable"] == "yes"
            and r["reachable_n3"] == "true"
            and r["ssign_call"] != "emitted_secreted"
        ]
        c = Counter()
        for r in miss:
            processed = any(r.get(col, "").strip() for col in SIG_COLS)
            c["processed" if processed else "detection"] += 1
        per[ss] = c

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    x = np.arange(len(TYPES))
    det = np.array([per[ss]["detection"] for ss in TYPES])
    pro = np.array([per[ss]["processed"] for ss in TYPES])
    ax.bar(x, det, color=C["detection"], width=0.6, label="detection miss — ssign found no system, never evaluated it")
    ax.bar(
        x,
        pro,
        bottom=det,
        color=C["processed"],
        width=0.6,
        label="processed but not emitted — predictors ran, filter rejected",
    )
    for xi, d, p in zip(x, det, pro):
        if d:
            ax.text(xi, d / 2, str(d), ha="center", va="center", color="white", fontsize=9.5, fontweight="bold")
        if p:
            ax.text(xi, d + p / 2, str(p), ha="center", va="center", color="white", fontsize=9.5, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{ss}\n(n={det[i] + pro[i]})" for i, ss in enumerate(TYPES)])
    ax.set_ylabel("reachable @±3 but missed")
    ax.set_ylim(top=max(det + pro) * 1.18 if max(det + pro) else 1)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1, fontsize=8.5)
    tot_d, tot_p = int(det.sum()), int(pro.sum())
    ax.set_title(
        f"ssign's false negatives are mostly DETECTION failures: {tot_d} of {tot_d + tot_p}\n"
        "the secretion system was never detected, so the effector was never evaluated"
    )
    ax.text(
        0.02,
        0.96,
        "T1SS: 4 RTX toxins (HlyA, ApxIA, LtxA, LktA) missed as detection\n"
        "failures because each draft assembly fragments the locus away from the\n"
        "chromosomal TolC, so MacSyFinder calls no T1SS. Full-assembly MacSyFinder\n"
        "(script 50) finds a complete T1SS for all four, so the whole-genome CX3\n"
        "rerun is expected to recover them. Other misses are genuine detection\n"
        "failures (T3SS/T4SS).",
        transform=ax.transAxes,
        fontsize=8,
        color="#555",
        va="top",
    )
    fig.tight_layout()
    fig.savefig(FIGDIR / "05_false_negative_modes.png")
    plt.close(fig)
    print("wrote 05_false_negative_modes.png")
    for ss in TYPES:
        print(f"  {ss}: detection={per[ss]['detection']} processed={per[ss]['processed']}")


if __name__ == "__main__":
    main()
