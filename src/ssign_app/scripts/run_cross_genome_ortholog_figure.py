#!/usr/bin/env python3
"""Render the cross-genome ortholog conservation figure (two panels).

Panel 1 (left): how many genomes share each ortholog group -- a histogram of
``n_genomes`` over all groups, with a "N unique / M shared" callout. This is the
core/accessory/unique shape of the secreted-protein pangenome.

Panel 2 (right): the 15 most-conserved groups (highest ``n_genomes``) as
horizontal stacked bars. Each bar's total width is the number of genomes carrying
the group; the stack segments split that width by secretion-system type (paralogs
deduped: a segment counts UNIQUE genomes carrying that type in the group). Each
bar is labelled with the group's most-common functional annotation and mean
within-group %identity.

Reads the augmented ``cross_genome_ortholog_groups.csv`` (ortholog_group,
n_members, n_genomes, genomes, members, mean_pident) plus every genome's
integrated CSV (for the per-protein annotation + SS-type + genome map). Uses the
shared house style so colours match the run figures 01-06.

    run_cross_genome_ortholog_figure.py \\
        --groups-csv cross_genome_ortholog_groups.csv \\
        --integrated-csvs a_integrated.csv b_integrated.csv ... \\
        --outfile figures/07_cross_genome_orthologs.png
"""

import argparse
import os
import sys
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.plot_style import THEME, apply_house_style, ordered_ss_types, ss_type_palette  # noqa: E402

from ssign_app.core._pool_utils import split_prefixed_id  # noqa: E402

# Annotation columns to try, in order of preference (curated consensus first).
_ANN_COLS = ("broad_consensus_annotation", "broad_annotation", "product")
_TOP_N = 15


def _primary_ss_type(raw: object) -> str:
    """First SS type in a (possibly multi-valued) nearby_ss_types cell."""
    s = str(raw or "").strip()
    if not s or s.lower() == "nan":
        return ""
    return s.replace(",", ";").split(";")[0].strip()


def build_protein_maps(integrated_csvs: list[str]) -> tuple[dict, dict, dict]:
    """Return (locus_tag -> annotation, locus_tag -> ss_type, locus_tag -> genome)
    pooled across every genome's integrated CSV."""
    pid_to_ann: dict[str, str] = {}
    pid_to_ss: dict[str, str] = {}
    pid_to_genome: dict[str, str] = {}
    for path in integrated_csvs:
        if not path or not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if "locus_tag" not in df.columns:
            continue
        ann_col = next((c for c in _ANN_COLS if c in df.columns), None)
        for _, row in df.iterrows():
            pid = str(row["locus_tag"])
            if ann_col:
                val = str(row.get(ann_col, "") or "").strip()
                if val and val.lower() != "nan":
                    pid_to_ann[pid] = val
            ss = _primary_ss_type(row.get("nearby_ss_types", ""))
            if ss:
                pid_to_ss[pid] = ss
            genome = str(row.get("sample_id", "") or "").strip()
            pid_to_genome[pid] = genome or pid
    return pid_to_ann, pid_to_ss, pid_to_genome


def _group_bar(members: str, pid_to_ann: dict, pid_to_ss: dict, pid_to_genome: dict, mean_pident: float):
    """Per-group: (label, {ss_type: n_unique_genomes}) for the panel-2 stacked bar.

    Members are ``sample_id__locus_tag`` (make_prefixed_id), so a group's genome
    comes straight from each member's prefix -- this keeps the per-SS-type genome
    count right even when a locus_tag is shared across genomes. Bare (un-prefixed)
    members fall back to the annotation-CSV genome map.
    """
    anns: list[str] = []
    ss_genomes: dict[str, set] = {}
    for raw in str(members).split(";"):
        m = raw.strip()
        if not m:
            continue
        try:
            genome, locus = split_prefixed_id(m)
        except ValueError:
            locus, genome = m, pid_to_genome.get(m, m)
        a = pid_to_ann.get(locus)
        if a:
            anns.append(a)
        ss = pid_to_ss.get(locus, "")
        if ss:
            ss_genomes.setdefault(ss, set()).add(genome)
    ann = Counter(anns).most_common(1)[0][0] if anns else "Unknown"
    ss_counts = {ss: len(gs) for ss, gs in ss_genomes.items()}
    return f"{ann} ({mean_pident:.0f}%id)", ss_counts


def render(groups_csv: str, integrated_csvs: list[str], outfile: str, dpi: int = 200) -> bool:
    og_df = pd.read_csv(groups_csv)
    if og_df.empty or "n_genomes" not in og_df.columns:
        return False

    pid_to_ann, pid_to_ss, pid_to_genome = build_protein_maps(integrated_csvs)

    genome_counts = og_df["n_genomes"].to_numpy()
    max_g = int(genome_counts.max())
    n_unique = int((og_df["n_genomes"] == 1).sum())
    n_shared = int((og_df["n_genomes"] > 1).sum())

    top = og_df.nlargest(_TOP_N, "n_genomes")
    labels, ss_count_rows = [], []
    for _, row in top.iterrows():
        label, ss_counts = _group_bar(row["members"], pid_to_ann, pid_to_ss, pid_to_genome, float(row["mean_pident"]))
        labels.append(label)
        ss_count_rows.append(ss_counts)

    ss_present = ordered_ss_types({ss for sc in ss_count_rows for ss in sc})
    palette = ss_type_palette(ss_present)

    apply_house_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

    # Panel 1: how many genomes share each group.
    ax1.hist(genome_counts, bins=range(1, max_g + 2), color=THEME["neutral_bar"], edgecolor="white")
    ax1.set_xlabel("Number of genomes sharing ortholog group")
    ax1.set_ylabel("Number of ortholog groups")
    ax1.set_title("Ortholog group size distribution")
    ax1.text(
        0.95,
        0.95,
        f"{n_unique} unique\n{n_shared} shared",
        transform=ax1.transAxes,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round", fc="#FFF8E1", ec=THEME["caption"]),
    )

    # Panel 2: most-conserved groups, stacked by SS type (width = unique genomes).
    y_pos = np.arange(len(top))
    left = np.zeros(len(top))
    for ss in ss_present:
        widths = np.array([sc.get(ss, 0) for sc in ss_count_rows], dtype=float)
        if not widths.any():
            continue
        ax2.barh(y_pos, widths, left=left, color=palette.get(ss, THEME["neutral_bar"]), edgecolor="white", label=ss)
        left += widths
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=8)
    ax2.set_xlabel("Number of genomes")
    ax2.set_title(f"{len(top)} most-conserved secreted-protein groups")
    ax2.invert_yaxis()
    handles = [mpatches.Patch(facecolor=palette.get(ss, THEME["neutral_bar"]), label=ss) for ss in ss_present]
    if handles:
        ax2.legend(handles=handles, title="SS type", loc="lower right", fontsize=8, title_fontsize=8, frameon=False)

    fig.suptitle(
        "Cross-genome secreted-protein ortholog conservation\nBLASTp-based ortholog groups",
        fontweight="bold",
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description="Render the cross-genome ortholog conservation figure")
    parser.add_argument("--groups-csv", required=True, help="cross_genome_ortholog_groups.csv (with n_genomes column)")
    parser.add_argument(
        "--integrated-csvs", nargs="+", required=True, help="Per-genome integrated CSVs (annotation + SS-type map)"
    )
    parser.add_argument("--outfile", required=True, help="Output PNG path")
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()

    if not os.path.exists(args.groups_csv):
        print(f"groups CSV not found: {args.groups_csv}", file=sys.stderr)
        return 1
    if render(args.groups_csv, args.integrated_csvs, args.outfile, args.dpi):
        print(f"Wrote {args.outfile}")
        return 0
    print("No ortholog groups to plot -- skipping figure", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
