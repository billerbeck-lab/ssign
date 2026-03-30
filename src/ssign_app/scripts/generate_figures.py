#!/usr/bin/env python3
"""Generate publication-quality figures for ssign pipeline output.

Adapted from:
- pipeline/scripts/generate_visualizations.py (VIZ-01 through VIZ-04)
- analyze_substrate_annotations.py (fig1-fig8)
- generate_characterization_figures.py (figA-figF)
"""

import argparse
import logging
import os
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def load_data(master_csvs):
    """Load and combine all master CSVs."""
    dfs = []
    for f in master_csvs:
        try:
            dfs.append(pd.read_csv(f))
        except Exception as e:
            logger.warning(f"Could not read {f}: {e}")
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def fig_ss_type_distribution(df, outdir, dpi):
    """Bar chart of substrate counts per SS type."""
    if 'nearby_ss_types' not in df.columns:
        return

    ss_counts = Counter()
    for val in df['nearby_ss_types'].dropna():
        for ss in str(val).split(','):
            ss = ss.strip()
            if ss:
                ss_counts[ss] += 1

    if not ss_counts:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    types = sorted(ss_counts.keys())
    counts = [ss_counts[t] for t in types]
    colors = sns.color_palette("Set2", len(types))

    ax.bar(types, counts, color=colors)
    ax.set_xlabel("Secretion System Type")
    ax.set_ylabel("Number of Substrates")
    ax.set_title("Substrate Distribution by Secretion System Type")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig1_ss_type_distribution.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig1_ss_type_distribution.png")


def fig_tool_coverage(df, outdir, dpi):
    """Heatmap of annotation tool coverage per substrate."""
    tool_prefixes = {
        'BLASTp': 'blastp_hit',
        'Foldseek': 'foldseek_best',
        'InterProScan': 'interpro_domains',
        'HHpred Pfam': 'pfam_top1',
        'HHpred PDB': 'pdb_top1',
        'pLM-BLAST': 'ecod70_top1',
        'SignalP': 'signalp_prediction',
    }

    coverage = {}
    for tool_name, prefix in tool_prefixes.items():
        matching_cols = [c for c in df.columns if c.startswith(prefix)]
        if matching_cols:
            coverage[tool_name] = df[matching_cols[0]].notna().sum()

    if not coverage:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    tools = list(coverage.keys())
    counts = [coverage[t] for t in tools]
    total = len(df)

    bars = ax.barh(tools, counts, color=sns.color_palette("Blues_d", len(tools)))
    ax.set_xlabel(f"Proteins with hits (out of {total})")
    ax.set_title("Annotation Tool Coverage")

    for bar, count in zip(bars, counts):
        pct = 100 * count / max(total, 1)
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                f"{pct:.0f}%", va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig2_tool_coverage.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig2_tool_coverage.png")


def fig_protein_lengths(df, outdir, dpi):
    """Violin plot of protein lengths by SS type."""
    if 'aa_length' not in df.columns and 'nearby_ss_types' not in df.columns:
        # Try to compute from sequence if available
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    if 'nearby_ss_types' in df.columns:
        # Explode multi-type substrates
        rows = []
        for _, row in df.iterrows():
            length = row.get('aa_length', 0)
            if pd.isna(length) or length == 0:
                continue
            for ss in str(row.get('nearby_ss_types', '')).split(','):
                ss = ss.strip()
                if ss:
                    rows.append({'ss_type': ss, 'length': length})

        if rows:
            plot_df = pd.DataFrame(rows)
            sns.violinplot(data=plot_df, x='ss_type', y='length', ax=ax,
                           inner='box', palette='Set2')
            ax.set_xlabel("Secretion System Type")
            ax.set_ylabel("Protein Length (aa)")
            ax.set_title("Substrate Protein Length by SS Type")
            plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig3_protein_lengths.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig3_protein_lengths.png")


def fig_physicochemical(df, outdir, dpi):
    """Violin plots of physicochemical properties."""
    props = ['gravy', 'mw_da', 'isoelectric_point', 'instability_index']
    available = [p for p in props if p in df.columns]

    if not available:
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4 * len(available), 6))
    if len(available) == 1:
        axes = [axes]

    for ax, prop in zip(axes, available):
        data = df[prop].dropna()
        if len(data) > 0:
            sns.violinplot(y=data, ax=ax, color='steelblue', inner='box')
            ax.set_title(prop.replace('_', ' ').title())

    plt.suptitle("Physicochemical Properties of Substrates", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig4_physicochemical.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig4_physicochemical.png")


def fig_category_distribution(df, outdir, dpi):
    """Bar chart of functional category distribution across substrates."""
    # Look for annotation columns that indicate functional categories
    cat_col = None
    for candidate in ['broad_annotation', 'functional_category', 'blastp_hit_description']:
        if candidate in df.columns:
            cat_col = candidate
            break

    if cat_col is None:
        logger.info("Skipping category distribution — no category column found")
        return

    cats = df[cat_col].dropna()
    if cats.empty:
        return

    counts = cats.value_counts().head(20)
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = sns.color_palette("Set3", len(counts))
    ax.barh(range(len(counts)), counts.values, color=colors)
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=9)
    ax.set_xlabel("Number of Substrates")
    ax.set_title("Functional Category Distribution")
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig5_category_distribution.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig5_category_distribution.png")


def fig_substrate_count_per_genome(df, outdir, dpi):
    """Bar chart of substrate counts per genome/sample."""
    sample_col = None
    for candidate in ['sample_id', 'genome', 'genome_id']:
        if candidate in df.columns:
            sample_col = candidate
            break

    if sample_col is None:
        logger.info("Skipping substrate count per genome — no sample column found")
        return

    counts = df[sample_col].value_counts().sort_index()
    if len(counts) < 2:
        return  # Only useful with multiple genomes

    fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.4), 6))
    ax.bar(range(len(counts)), counts.values, color='steelblue')
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=90, fontsize=7)
    ax.set_xlabel("Genome")
    ax.set_ylabel("Number of Substrates")
    ax.set_title("Substrate Count per Genome")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig6_substrate_count.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig6_substrate_count.png")


def fig_functional_summary(df, outdir, dpi):
    """Stacked bar chart: functional categories per SS type."""
    if 'nearby_ss_types' not in df.columns:
        return

    cat_col = None
    for candidate in ['broad_annotation', 'functional_category', 'blastp_hit_description']:
        if candidate in df.columns:
            cat_col = candidate
            break

    if cat_col is None:
        return

    rows = []
    for _, row in df.iterrows():
        cat = row.get(cat_col, '')
        if pd.isna(cat) or not cat:
            cat = 'Unknown'
        for ss in str(row.get('nearby_ss_types', '')).split(','):
            ss = ss.strip()
            if ss:
                rows.append({'ss_type': ss, 'category': cat})

    if not rows:
        return

    plot_df = pd.DataFrame(rows)
    ct = pd.crosstab(plot_df['ss_type'], plot_df['category'])

    fig, ax = plt.subplots(figsize=(12, 7))
    ct.plot(kind='bar', stacked=True, ax=ax, colormap='Set3')
    ax.set_xlabel("Secretion System Type")
    ax.set_ylabel("Number of Substrates")
    ax.set_title("Functional Categories per SS Type")
    ax.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig7_functional_summary.png"), dpi=dpi)
    plt.close()
    logger.info("Generated fig7_functional_summary.png")


def main():
    parser = argparse.ArgumentParser(description="Generate ssign figures")
    parser.add_argument("--master-csvs", nargs='+', required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    # Figure toggles
    parser.add_argument("--no-category", action="store_true",
                        help="Skip category distribution figure")
    parser.add_argument("--no-ss-comp", action="store_true",
                        help="Skip SS type distribution figure")
    parser.add_argument("--no-tool-heatmap", action="store_true",
                        help="Skip tool coverage heatmap figure")
    parser.add_argument("--no-substrate-count", action="store_true",
                        help="Skip substrate count per genome figure")
    parser.add_argument("--no-func-summary", action="store_true",
                        help="Skip functional summary figure")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = load_data(args.master_csvs)
    if df.empty:
        logger.warning("No data to plot")
        return

    logger.info(f"Generating figures for {len(df)} substrates...")

    # Always generate core figures (SS distribution, tool coverage, lengths, physico)
    if not args.no_ss_comp:
        fig_ss_type_distribution(df, args.outdir, args.dpi)
    if not args.no_tool_heatmap:
        fig_tool_coverage(df, args.outdir, args.dpi)
    fig_protein_lengths(df, args.outdir, args.dpi)
    fig_physicochemical(df, args.outdir, args.dpi)

    # Optional figures controlled by toggles
    if not args.no_category:
        fig_category_distribution(df, args.outdir, args.dpi)
    if not args.no_substrate_count:
        fig_substrate_count_per_genome(df, args.outdir, args.dpi)
    if not args.no_func_summary:
        fig_functional_summary(df, args.outdir, args.dpi)

    logger.info(f"Figures saved to {args.outdir}")


if __name__ == '__main__':
    main()
