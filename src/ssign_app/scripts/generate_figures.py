#!/usr/bin/env python3
"""Generate publication-quality figures for ssign pipeline output.

Two modes:
- ``per_genome`` (default): the curated per-run figure set (``01_*`` ... ``07_*``),
  drawn from one genome's integrated substrate CSV.
- ``pooled``: cross-genome figures (``P01_*`` ...) drawn from several genomes'
  integrated CSVs; used by the multi-genome runner. No-op for a single genome.

Style (theme, palette, numbering, figure index) comes from ``ssign_lib.plot_style``
so the look never drifts from the enrichment figure. Every figure guards the
columns it needs and is skipped with a logged note if they are absent.
"""

import argparse
import logging
import os
import sys
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import CONF_THRESHOLD  # noqa: E402
from ssign_lib.plot_style import (  # noqa: E402
    SEQUENTIAL_CMAP,
    SS_TYPE_ORDER,
    THEME,
    apply_house_style,
    clear_figure_set,
    numbered_path,
    pooled_path,
    print_figure_index,
    ss_type_palette,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_SIGNALP_NEGATIVE = {"", "OTHER", "NO_SP", "NONE", "-", "NAN"}


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


def _ordered_types(types):
    """Canonical SS-type order first, then any extras alphabetically."""
    present = set(types)
    known = [t for t in SS_TYPE_ORDER if t in present]
    extra = sorted(t for t in present if t not in SS_TYPE_ORDER)
    return known + extra


def _explode_ss_types(df, cols):
    """Long frame with one row per (substrate, SS type) carrying `cols`."""
    keep = [c for c in cols if c in df.columns]
    if "nearby_ss_types" not in df.columns:
        return pd.DataFrame(columns=["ss_type", *keep])
    long = df[keep].copy()
    long["ss_type"] = df["nearby_ss_types"].fillna("").astype(str).str.split(",")
    long = long.explode("ss_type")
    long["ss_type"] = long["ss_type"].str.strip()
    return long[long["ss_type"] != ""].reset_index(drop=True)


def _signalp_positive(v) -> bool:
    """True when SignalP assigned a signal-peptide class (not OTHER / missing)."""
    return str(v).strip().upper() not in _SIGNALP_NEGATIVE


def _bar_counts(ax, labels, counts, colors, annotate=True):
    bars = ax.bar(range(len(labels)), counts, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    if annotate:
        for b, c in zip(bars, counts):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(), str(int(c)), ha="center", va="bottom", fontsize=8)
    return bars


def _hbar_counts(ax, counts, xlabel, title):
    """Horizontal count bars (counts = a value_counts Series), most-frequent on top."""
    ax.barh(range(len(counts)), counts.values, color=THEME["neutral_bar"])
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(counts.index, fontsize=9)
    ax.invert_yaxis()
    for i, c in enumerate(counts.values):
        ax.text(c, i, f" {int(c)}", va="center", fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title)


# --- per-genome figures --------------------------------------------------------


def fig01_substrates_per_type(df, outdir, dpi):
    if "nearby_ss_types" not in df.columns:
        logger.info("Skip 01 substrates-per-type: no nearby_ss_types")
        return None
    counts: Counter = Counter()
    for val in df["nearby_ss_types"].dropna():
        for ss in str(val).split(","):
            ss = ss.strip()
            if ss:
                counts[ss] += 1
    if not counts:
        return None
    types = _ordered_types(counts)
    pal = ss_type_palette(types)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    _bar_counts(ax, types, [counts[t] for t in types], [pal[t] for t in types])
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("Secreted proteins")
    ax.set_title("Secreted proteins by secretion-system type")
    fig.tight_layout()
    out = numbered_path(outdir, 1, "substrates_per_type")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("01", os.path.basename(out), "secreted proteins per SS type")


def fig02_secretion_evidence(df, outdir, dpi):
    """Which trusted predictor(s) flagged each protein (the `tool` column)."""
    if "tool" not in df.columns or df["tool"].dropna().empty:
        logger.info("Skip 02 secretion-evidence: no tool column")
        return None
    counts = df["tool"].fillna("(none)").value_counts()
    fig, ax = plt.subplots(figsize=(8, 5))
    _hbar_counts(ax, counts, "Secreted proteins", "Secretion-call support (predictor / combination per protein)")
    fig.tight_layout()
    out = numbered_path(outdir, 2, "secretion_evidence")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("02", os.path.basename(out), "predictor support per secreted protein")


def fig03_localization_confidence(df, outdir, dpi):
    long = _explode_ss_types(df, ["dlp_extracellular_prob"])
    if long.empty or "dlp_extracellular_prob" not in long.columns:
        logger.info("Skip 03 localization-confidence: no dlp_extracellular_prob")
        return None
    long["dlp_extracellular_prob"] = pd.to_numeric(long["dlp_extracellular_prob"], errors="coerce")
    long = long.dropna(subset=["dlp_extracellular_prob"])
    if long.empty:
        return None
    types = _ordered_types(long["ss_type"].unique())
    pal = ss_type_palette(types)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.stripplot(
        data=long,
        x="ss_type",
        y="dlp_extracellular_prob",
        hue="ss_type",
        order=types,
        palette=pal,
        legend=False,
        ax=ax,
        size=5,
        jitter=0.25,
        alpha=0.8,
    )
    ax.axhline(
        CONF_THRESHOLD, color=THEME["ref_line"], ls="--", lw=0.8, alpha=0.8, label=f"call threshold = {CONF_THRESHOLD}"
    )
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("DeepLocPro extracellular probability")
    ax.set_title("Localization confidence of secreted proteins by SS type")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    out = numbered_path(outdir, 3, "localization_confidence")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("03", os.path.basename(out), "DLP extracellular probability by SS type")


def fig04_signalp_by_type(df, outdir, dpi):
    long = _explode_ss_types(df, ["signalp_prediction"])
    if long.empty or "signalp_prediction" not in long.columns:
        logger.info("Skip 04 signalp-by-type: no signalp_prediction")
        return None
    long["sp_pos"] = long["signalp_prediction"].apply(_signalp_positive)
    types = _ordered_types(long["ss_type"].unique())
    pal = ss_type_palette(types)
    frac, ann = [], []
    for t in types:
        sub = long[long["ss_type"] == t]
        f = sub["sp_pos"].mean() if len(sub) else 0.0
        frac.append(f)
        ann.append(f"{int(sub['sp_pos'].sum())}/{len(sub)}")
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(range(len(types)), frac, color=[pal[t] for t in types])
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(types, rotation=45, ha="right")
    for b, a in zip(bars, ann):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), a, ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("Fraction with SignalP signal peptide")
    ax.set_title("SignalP-positive fraction per SS type (predictor call)")
    ax.text(
        0.5,
        -0.28,
        "SignalP-negative is expected for injected types (T3SS/T6SS) and T1SS; positive for T2SS/T5SS.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7.5,
        color=THEME["caption"],
    )
    fig.tight_layout()
    out = numbered_path(outdir, 4, "signalp_by_type")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("04", os.path.basename(out), "SignalP-positive fraction per SS type")


def fig05_tool_coverage(df, outdir, dpi):
    tool_prefixes = {
        "BLASTp": "blastp_hit",
        "InterProScan": "interpro_domains",
        "HHpred Pfam": "pfam_top1",
        "HHpred PDB": "pdb_top1",
        "pLM-BLAST": "ecod_top1",
        "SignalP": "signalp_prediction",
    }
    coverage = {}
    for name, prefix in tool_prefixes.items():
        cols = [c for c in df.columns if c.startswith(prefix)]
        if cols:
            coverage[name] = int(df[cols[0]].notna().sum())
    if not coverage:
        logger.info("Skip 05 tool-coverage: no annotation columns")
        return None
    tools = list(coverage.keys())
    counts = [coverage[t] for t in tools]
    total = len(df)
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.barh(tools, counts, color=THEME["neutral_bar"])
    ax.set_xlabel(f"Proteins with a hit (of {total})")
    ax.set_title("Annotation-tool coverage")
    for b, c in zip(bars, counts):
        ax.text(
            b.get_width(), b.get_y() + b.get_height() / 2, f" {100 * c / max(total, 1):.0f}%", va="center", fontsize=9
        )
    fig.tight_layout()
    out = numbered_path(outdir, 5, "tool_coverage")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("05", os.path.basename(out), "annotation-tool hit coverage")


def fig06_protein_length(df, outdir, dpi):
    long = _explode_ss_types(df, ["aa_length"])
    if long.empty or "aa_length" not in long.columns:
        logger.info("Skip 06 protein-length: no aa_length")
        return None
    long["aa_length"] = pd.to_numeric(long["aa_length"], errors="coerce")
    long = long[long["aa_length"] > 0].dropna(subset=["aa_length"])
    if long.empty:
        return None
    types = _ordered_types(long["ss_type"].unique())
    pal = ss_type_palette(types)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.violinplot(
        data=long,
        x="ss_type",
        y="aa_length",
        hue="ss_type",
        order=types,
        palette=pal,
        legend=False,
        ax=ax,
        inner="box",
        cut=0,
    )
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("Protein length (aa)")
    ax.set_title("Secreted-protein length by SS type")
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    out = numbered_path(outdir, 6, "protein_length")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("06", os.path.basename(out), "protein length by SS type")


def _category_column(df):
    for c in ["broad_annotation", "broad_consensus_annotation", "functional_category"]:
        if c in df.columns:
            return c
    return None


def fig07_functional_categories(df, outdir, dpi, top_n=8):
    cat_col = _category_column(df)
    if cat_col is None or "nearby_ss_types" not in df.columns:
        logger.info("Skip 07 functional-categories: no annotation/ss column")
        return None
    long = _explode_ss_types(df, [cat_col])
    if long.empty:
        return None
    long[cat_col] = long[cat_col].fillna("Unknown").replace("", "Unknown")
    top = long[cat_col].value_counts().head(top_n).index
    long["cat"] = long[cat_col].where(long[cat_col].isin(top), "Other")
    types = _ordered_types(long["ss_type"].unique())
    ct = pd.crosstab(long["ss_type"], long["cat"]).reindex(types).fillna(0)
    cats = list(ct.columns)
    palette = sns.color_palette("muted", len(cats))
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ct.plot(kind="bar", stacked=True, ax=ax, color=palette, width=0.8)
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("Secreted proteins")
    ax.set_title("Functional categories per SS type")
    ax.legend(title="Category", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    out = numbered_path(outdir, 7, "functional_categories")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("07", os.path.basename(out), "functional categories per SS type")


def fig_physicochemical(df, outdir, dpi):
    """Opt-in only: needs ProtParam columns the standard CSV does not carry."""
    props = ["gravy", "mw_da", "isoelectric_point", "instability_index"]
    available = [p for p in props if p in df.columns and df[p].notna().any()]
    if not available:
        logger.info("Skip physicochemical: ProtParam columns absent")
        return None
    fig, axes = plt.subplots(1, len(available), figsize=(4 * len(available), 5.5))
    if len(available) == 1:
        axes = [axes]
    for ax, prop in zip(axes, available):
        data = pd.to_numeric(df[prop], errors="coerce").dropna()
        if len(data):
            sns.violinplot(y=data, ax=ax, color=THEME["neutral_bar"], inner="box", cut=0)
            ax.set_title(prop.replace("_", " ").title())
            ax.set_ylabel("")
    fig.suptitle("Physicochemical properties of secreted proteins", y=1.02, fontweight="bold")
    fig.tight_layout()
    out = numbered_path(outdir, 8, "physicochemical")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("08", os.path.basename(out), "physicochemical properties (opt-in)")


# --- pooled cross-genome figures ----------------------------------------------


def _sample_column(df):
    for c in ["sample_id", "genome", "genome_id"]:
        if c in df.columns:
            return c
    return None


def pooled01_substrates_per_genome(df, outdir, dpi):
    sample_col = _sample_column(df)
    if sample_col is None:
        logger.info("Skip P01 substrates-per-genome: no sample column")
        return None
    counts = df[sample_col].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.5), 5.5))
    _bar_counts(ax, list(counts.index), list(counts.values), [THEME["neutral_bar"]] * len(counts))
    ax.set_xlabel("Genome")
    ax.set_ylabel("Secreted proteins")
    ax.set_title("Secreted proteins per genome")
    fig.tight_layout()
    out = pooled_path(outdir, 1, "substrates_per_genome")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("P01", os.path.basename(out), "secreted proteins per genome")


def pooled02_sstype_by_genome(df, outdir, dpi):
    sample_col = _sample_column(df)
    if sample_col is None or "nearby_ss_types" not in df.columns:
        logger.info("Skip P02 sstype-by-genome: no sample/ss column")
        return None
    long = _explode_ss_types(df.assign(_g=df[sample_col]), ["_g"])
    if long.empty:
        return None
    types = _ordered_types(long["ss_type"].unique())
    ct = pd.crosstab(long["_g"], long["ss_type"]).reindex(columns=types).fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(max(7, len(types) * 0.9), max(4, len(ct) * 0.6) + 1))
    sns.heatmap(ct, annot=True, fmt="d", cmap=SEQUENTIAL_CMAP, ax=ax, cbar_kws={"label": "Secreted proteins"})
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("Genome")
    ax.set_title("Secreted proteins: SS type x genome")
    fig.tight_layout()
    out = pooled_path(outdir, 2, "sstype_by_genome")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("P02", os.path.basename(out), "SS-type x genome substrate heatmap")


def pooled03_evidence_basis(df, outdir, dpi):
    if "tool" not in df.columns or df["tool"].dropna().empty:
        logger.info("Skip P03 evidence-basis: no tool column")
        return None
    counts = df["tool"].fillna("(none)").value_counts()
    fig, ax = plt.subplots(figsize=(8, 5))
    _hbar_counts(ax, counts, "Secreted proteins (all genomes)", "Secretion-call support, pooled across genomes")
    fig.tight_layout()
    out = pooled_path(outdir, 3, "evidence_basis")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("P03", os.path.basename(out), "pooled predictor support")


PER_GENOME_FIGS = [
    ("ss_comp", fig01_substrates_per_type),
    ("evidence", fig02_secretion_evidence),
    ("localization", fig03_localization_confidence),
    ("signalp", fig04_signalp_by_type),
    ("tool_heatmap", fig05_tool_coverage),
    ("length", fig06_protein_length),
    ("func_summary", fig07_functional_categories),
]
POOLED_FIGS = [pooled01_substrates_per_genome, pooled02_sstype_by_genome, pooled03_evidence_basis]


def main():
    parser = argparse.ArgumentParser(description="Generate ssign figures")
    parser.add_argument("--master-csvs", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--mode", choices=["per_genome", "pooled"], default="per_genome")
    # Per-figure toggles (default on except physicochemical).
    parser.add_argument("--no-ss-comp", action="store_true")
    parser.add_argument("--no-evidence", action="store_true")
    parser.add_argument("--no-localization", action="store_true")
    parser.add_argument("--no-signalp", action="store_true")
    parser.add_argument("--no-tool-heatmap", action="store_true")
    parser.add_argument("--no-length", action="store_true")
    parser.add_argument("--no-func-summary", action="store_true")
    parser.add_argument("--physicochemical", action="store_true", help="Opt-in ProtParam figure (off by default).")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load_data(args.master_csvs)
    if df.empty:
        logger.warning("No data to plot")
        return

    apply_house_style()
    clear_figure_set(args.outdir)
    entries = []

    if args.mode == "pooled":
        sample_col = _sample_column(df)
        n_genomes = df[sample_col].nunique() if sample_col else 0
        if n_genomes < 2:
            logger.info("Pooled mode needs >=2 genomes; got %d. Nothing to plot.", n_genomes)
            return
        logger.info("Pooling %d substrates across %d genomes...", len(df), n_genomes)
        for fn in POOLED_FIGS:
            e = fn(df, args.outdir, args.dpi)
            if e:
                entries.append(e)
    else:
        logger.info("Generating per-genome figures for %d secreted proteins...", len(df))
        # Each PER_GENOME_FIGS key has a matching --no-<key> flag (dest no_<key>).
        for key, fn in PER_GENOME_FIGS:
            if getattr(args, f"no_{key}", False):
                continue
            e = fn(df, args.outdir, args.dpi)
            if e:
                entries.append(e)
        if args.physicochemical:
            e = fig_physicochemical(df, args.outdir, args.dpi)
            if e:
                entries.append(e)

    print_figure_index(entries, logger=logger)
    logger.info("Figures saved to %s", args.outdir)


if __name__ == "__main__":
    main()
