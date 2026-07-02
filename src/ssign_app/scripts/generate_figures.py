#!/usr/bin/env python3
"""Generate publication-quality figures for ssign pipeline output.

Two modes:
- ``per_genome`` (default): the curated per-run figure set (``01_*``-``06_*``),
  drawn from one genome's integrated substrate CSV.
- ``pooled``: the same curated set computed over several genomes' integrated CSVs
  and written as ``0N_pooled_*`` (used by the multi-genome runner). The ``01``
  figure (secreted proteins per genome, stacked by SS type) is itself the
  cross-genome overview. No-op for a single genome.

Curated set:
  01 secreted proteins by SS type (one bar per SS type for a single genome; one
     stacked bar per genome for a group)
  02 size & physicochemical properties by SS type (length, MW, GRAVY, pI, ...)
  03 COG functional category by SS type
  04 KEGG function by SS type
  05 EggNOG description by SS type
  06 curated consensus function by SS type

Autotransporter (T5aSS/T5cSS) self-detection is not a standalone figure; the
signal is carried by the enrichment self-mode result. T5a/c substrates are split
into ``(self)`` (the component) and ``(hitch)`` (proximity "hitchhiker"
neighbours) categories in the count and functional figures.

Style (theme, palette, numbering, figure index) comes from ``ssign_lib.plot_style``
so the look never drifts from the enrichment figure. Every figure guards the
columns it needs and is skipped with a logged note if they are absent.
"""

import argparse
import logging
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import (  # noqa: E402  # display_type: collapse SS subtypes
    ENRICH_AUTOTRANSPORTER_TYPES,
    SUBSTRATE_SOURCE_T5_SELF,
    T5_HITCH_TAG,
    T5_SELF_TAG,
    display_type,
)
from ssign_lib.functional_vocab import (  # noqa: E402
    cog_category_names,
    consensus_bucket,
    is_missing,
    kegg_descriptions,
)
from ssign_lib.plot_style import (  # noqa: E402
    THEME,
    apply_house_style,
    clear_figure_set,
    numbered_path,
    ordered_ss_types,
    print_figure_index,
    ss_type_palette,
)

# Autotransporters (T5aSS/T5cSS) are the substrate themselves; their rows are split
# into (self) vs (hitch) categories in the curated figures (see _t5_split_label).
AUTOTRANSPORTER_TYPES = ENRICH_AUTOTRANSPORTER_TYPES

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Size + ProtParam properties (column, axis label), in reading order. Length leads
# (it lives with the physicochemical panels rather than as its own figure).
# Column-guarded per panel.
_PHYSCHEM_PROPS = [
    ("aa_length", "Length (aa)"),
    ("mw_da", "Molecular weight (Da)"),
    ("gravy", "GRAVY (hydropathy)"),
    ("isoelectric_point", "Isoelectric point (pI)"),
    ("instability_index", "Instability index"),
    ("aromaticity", "Aromaticity"),
    ("charge_ph7", "Charge at pH 7"),
]


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


def _sample_column(df):
    for c in ["sample_id", "genome", "genome_id"]:
        if c in df.columns:
            return c
    return None


def _t5_split_label(ss_type: str, substrate_source) -> str:
    """Split autotransporter (T5aSS/T5cSS) rows into self vs hitchhiker categories:
    the component itself (`substrate_source == "T5SS-self"`) is `T5aSS (self)`, a
    secreted-predicted proximity neighbour is `T5aSS (hitch)`. Other types pass
    through unchanged."""
    if ss_type in AUTOTRANSPORTER_TYPES:
        is_self = str(substrate_source) == SUBSTRATE_SOURCE_T5_SELF
        return f"{ss_type} {T5_SELF_TAG}" if is_self else f"{ss_type} {T5_HITCH_TAG}"
    return ss_type


def _explode_ss_types(df, cols):
    """Long frame with one row per (substrate, SS type) carrying `cols`. T5a/c rows
    are split into `(self)` / `(hitch)` categories via `substrate_source`."""
    keep = [c for c in cols if c in df.columns]
    if "nearby_ss_types" not in df.columns:
        return pd.DataFrame(columns=["ss_type", *keep])
    src = df["substrate_source"] if "substrate_source" in df.columns else None
    long = df[keep].copy()
    long["ss_type"] = df["nearby_ss_types"].fillna("").astype(str).str.split(",")
    if src is not None:
        long["_src"] = src.values
    long = long.explode("ss_type")
    # Collapse subtype labels to canonical display types (T6SSi/ii -> T6SS,
    # pT4SSt -> T4SS; T5 subtypes kept distinct) so figures match the enrichment.
    long["ss_type"] = long["ss_type"].str.strip().map(lambda s: display_type(s) if s else s)
    long = long[long["ss_type"] != ""].reset_index(drop=True)
    if "_src" in long.columns:
        long["ss_type"] = [_t5_split_label(t, s) for t, s in zip(long["ss_type"], long["_src"])]
        long = long.drop(columns=["_src"])
    return long


def _stacked_by_type(ct, ax, types):
    """Stacked bar of a (rows x ss_type) crosstab, using the stable SS-type palette."""
    pal = ss_type_palette(types)
    ct = ct.reindex(columns=types).fillna(0)
    ct.plot(kind="bar", stacked=True, ax=ax, color=[pal[t] for t in ct.columns], width=0.8)
    ax.legend(title="SS type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    ax.tick_params(axis="x", rotation=45)
    for lbl in ax.get_xticklabels():
        lbl.set_ha("right")


# --- curated figures -----------------------------------------------------------


def fig01_secreted_by_genome(df, outdir, dpi):
    """Secreted proteins by SS type. A single genome gets one bar per SS type; a
    group gets one stacked bar per genome (the cross-genome overview)."""
    if "nearby_ss_types" not in df.columns:
        logger.info("Skip 01 secreted-by-genome: no nearby_ss_types")
        return None
    sample_col = _sample_column(df)
    g = df[sample_col].astype(str) if sample_col else "genome"
    long = _explode_ss_types(df.assign(_g=g), ["_g"])
    if long.empty:
        return None
    types = ordered_ss_types(long["ss_type"].unique())
    pal = ss_type_palette(types)

    if long["_g"].nunique() < 2:
        # Single genome: a plain count histogram, one bar per SS type.
        counts = long["ss_type"].value_counts().reindex(types).fillna(0)
        fig, ax = plt.subplots(figsize=(max(7, len(types) * 0.9 + 1), 5.5))
        ax.bar(range(len(types)), counts.values, color=[pal[t] for t in types])
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels(types, rotation=45, ha="right")
    else:
        # Group: one stacked bar per genome, coloured by SS type.
        ct = pd.crosstab(long["_g"], long["ss_type"])
        fig, ax = plt.subplots(figsize=(max(6, len(ct) * 0.6 + 2.5), 5.8))
        _stacked_by_type(ct, ax, types)
        ax.set_xlabel("Genome")
    ax.set_ylabel("Secreted proteins")
    ax.set_title("Secreted proteins by secretion-system type")
    fig.tight_layout()
    out = numbered_path(outdir, 1, "secreted_by_genome")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("01", os.path.basename(out), "secreted proteins by SS type (per genome for a group)")


def fig02_physicochemical(df, outdir, dpi):
    """One violin-by-SS-type panel per size/physicochemical property (length, MW,
    GRAVY, pI, instability, aromaticity, charge). Column-guarded per panel; emitted
    when at least one property is present (length is always; the ProtParam columns
    when ProtParam ran)."""
    avail = [
        (c, lbl)
        for c, lbl in _PHYSCHEM_PROPS
        if c in df.columns and pd.to_numeric(df[c], errors="coerce").notna().any()
    ]
    if not avail:
        logger.info("Skip 02 physicochemical: no size/ProtParam columns")
        return None
    long = _explode_ss_types(df, [c for c, _ in avail])
    if long.empty:
        return None
    types = ordered_ss_types(long["ss_type"].unique())
    pal = ss_type_palette(types)
    ncol = 2
    nrow = (len(avail) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(7 * ncol, 4.2 * nrow), squeeze=False)
    for idx, (col, lbl) in enumerate(avail):
        ax = axes[idx // ncol][idx % ncol]
        sub = long[["ss_type", col]].copy()
        sub[col] = pd.to_numeric(sub[col], errors="coerce")
        sub = sub.dropna(subset=[col])
        if sub.empty:
            ax.axis("off")
            continue
        sns.violinplot(
            data=sub,
            x="ss_type",
            y=col,
            hue="ss_type",
            order=types,
            palette=pal,
            legend=False,
            ax=ax,
            inner="box",
            cut=0,
        )
        ax.set_title(lbl)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=45)
        for t in ax.get_xticklabels():
            t.set_ha("right")
    for j in range(len(avail), nrow * ncol):  # hide unused panels
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Size & physicochemical properties of secreted proteins by SS type", y=1.0, fontweight="bold")
    fig.tight_layout()
    out = numbered_path(outdir, 2, "physicochemical")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return ("02", os.path.basename(out), "size & physicochemical properties by SS type")


# --- functional-category figures (by SS type) ----------------------------------


def _consensus_column(df):
    for c in ("broad_annotation", "broad_consensus_annotation"):
        if c in df.columns:
            return c
    return None


def _eggnog_desc(value) -> list:
    """One `eggnog_description` cell -> [description] (truncated for the axis), or
    [] when blank / `-`. High-cardinality free text; the top-N keeps it readable."""
    if is_missing(value):
        return []
    s = str(value).strip()
    return [s if len(s) <= 50 else s[:47] + "..."]


# Cap how many categories each functional figure shows before the rest collapse
# to a single grey "Other" bar, so the stacks stay readable and colours distinct.
_FUNCTIONAL_TOP_N = 10


def _functional_sources(df):
    """Available (column, value_fn, slug, label, fig_number). `value_fn` maps one
    cell to a list of category strings (multi-valued sources explode)."""
    sources = []
    if "cog_category" in df.columns:
        sources.append(("cog_category", cog_category_names, "cog_category", "COG functional category", 3))
    if "kegg_ko" in df.columns:
        sources.append(("kegg_ko", kegg_descriptions, "kegg_function", "KEGG function", 4))
    if "eggnog_description" in df.columns:
        sources.append(("eggnog_description", _eggnog_desc, "eggnog_description", "EggNOG description", 5))
    cons = _consensus_column(df)
    if cons:
        sources.append((cons, lambda v: [consensus_bucket(v)], "consensus_function", "Consensus function (curated)", 6))
    return sources


def _functional_by_type(df, col, value_fn):
    """Long (ss_type, cat) rows for the stacked-by-SS-type scope."""
    long = _explode_ss_types(df, [col])
    if long.empty or col not in long.columns:
        return pd.DataFrame(columns=["ss_type", "cat"])
    long["cat"] = long[col].apply(value_fn)
    long = long.explode("cat").dropna(subset=["cat"])
    return long[["ss_type", "cat"]].reset_index(drop=True)  # explode dup'd the index; crosstab needs it unique


def _plot_functional_by_type(long, outdir, n, slug, label, dpi):
    if long.empty:
        return None
    top = long["cat"].value_counts().head(_FUNCTIONAL_TOP_N).index
    long = long.assign(catg=long["cat"].where(long["cat"].isin(top), "Other"))
    types = ordered_ss_types(long["ss_type"].unique())
    ct = pd.crosstab(long["ss_type"], long["catg"]).reindex(types).fillna(0)
    if ct.values.sum() == 0:
        return None
    # Columns most-common-first, "Other" last; distinct colours from tab20 with
    # "Other" greyed out (so colours never repeat within the legend).
    cats = [c for c in ct.sum().sort_values(ascending=False).index if c != "Other"]
    if "Other" in ct.columns:
        cats.append("Other")
    ct = ct[cats]
    base = sns.color_palette("tab20", max(1, sum(c != "Other" for c in cats)))
    colors, bi = [], 0
    for c in cats:
        colors.append(THEME["muted"] if c == "Other" else base[bi])
        bi += 0 if c == "Other" else 1
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ct.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.8)
    ax.set_xlabel("Secretion system type")
    ax.set_ylabel("Secreted proteins")
    ax.set_title(f"{label} by SS type")
    ax.legend(title=label, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    out = numbered_path(outdir, n, f"{slug}_by_sstype")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return (f"{n:02d}", os.path.basename(out), f"{label}: by SS type")


def fig_functional(df, outdir, dpi):
    """Functional-category figures from controlled vocabularies (COG / KEGG /
    EggNOG / curated consensus), stacked by SS type. Column-guarded per source;
    returns a list of index entries."""
    if "nearby_ss_types" not in df.columns:
        logger.info("Skip functional figures: no nearby_ss_types")
        return []
    sources = _functional_sources(df)
    if not sources:
        logger.info(
            "Skip functional figures: no functional source column "
            "(cog_category / kegg_ko / eggnog_description / broad_annotation)"
        )
        return []
    entries = []
    for col, value_fn, slug, label, n in sources:
        e = _plot_functional_by_type(_functional_by_type(df, col, value_fn), outdir, n, slug, label, dpi)
        if e:
            entries.append(e)
    return entries


# `func_summary` is the functional dispatcher: it emits up to four numbered
# figures (03-06) and returns a list of index entries, so the runner handles a
# list or a single tuple.
PER_GENOME_FIGS = [
    ("ss_comp", fig01_secreted_by_genome),
    ("physicochemical", fig02_physicochemical),
    ("func_summary", fig_functional),
]


def _run_per_genome_figs(df, outdir, dpi, skip):
    """Run the curated per-genome figures, flattening list-returning fns (the
    functional dispatcher). `skip` is a set of PER_GENOME_FIGS keys to omit."""
    entries = []
    for key, fn in PER_GENOME_FIGS:
        if key in skip:
            continue
        e = fn(df, outdir, dpi)
        if isinstance(e, list):
            entries.extend(e)
        elif e:
            entries.append(e)
    return entries


def _emit_pooled_per_genome(df, outdir, dpi, skip):
    """The curated set over the combined frame, renamed ``0N_pooled_<name>.png``
    (the genome-group view; the 01 figure becomes the per-genome stacked overview)."""
    entries = []
    for label, fname, desc in _run_per_genome_figs(df, outdir, dpi, skip):
        renamed = f"{label}_pooled_{fname.split('_', 1)[1]}"  # 01_x.png -> 01_pooled_x.png
        os.replace(os.path.join(outdir, fname), os.path.join(outdir, renamed))
        entries.append((label, renamed, f"pooled: {desc}"))
    return entries


def main():
    parser = argparse.ArgumentParser(description="Generate ssign figures")
    parser.add_argument("--master-csvs", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--mode", choices=["per_genome", "pooled"], default="per_genome")
    # Per-figure toggles (all default on). Each PER_GENOME_FIGS key has a --no-<key>.
    parser.add_argument("--no-ss-comp", action="store_true")
    parser.add_argument("--no-physicochemical", action="store_true")
    parser.add_argument("--no-func-summary", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load_data(args.master_csvs)
    if df.empty:
        logger.warning("No data to plot")
        return

    apply_house_style()
    clear_figure_set(args.outdir)
    skip = {key for key, _ in PER_GENOME_FIGS if getattr(args, f"no_{key}", False)}

    if args.mode == "pooled":
        sample_col = _sample_column(df)
        n_genomes = df[sample_col].nunique() if sample_col else 0
        if n_genomes < 2:
            logger.info("Pooled mode needs >=2 genomes; got %d. Nothing to plot.", n_genomes)
            return
        logger.info("Pooling %d substrates across %d genomes...", len(df), n_genomes)
        # The curated set over all genomes combined (01 = stacked per-genome overview).
        entries = _emit_pooled_per_genome(df, args.outdir, args.dpi, skip)
    else:
        logger.info("Generating per-genome figures for %d secreted proteins...", len(df))
        entries = _run_per_genome_figs(df, args.outdir, args.dpi, skip)

    print_figure_index(entries, logger=logger)
    logger.info("Figures saved to %s", args.outdir)


if __name__ == "__main__":
    main()
