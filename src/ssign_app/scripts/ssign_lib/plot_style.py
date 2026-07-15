#!/usr/bin/env python3
"""Shared plotting style for ssign figures (publication-plots house rules).

One THEME, one rcParams setup, one stable SS-type palette, and small file
helpers (numbered filenames, stale-figure cleanup, figure index). Imported by
both generate_figures.py (regular + pooled figures) and run_enrichment_figure.py
(enrichment null grid) so colours and style never drift between the two scripts.
"""

from __future__ import annotations

import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .constants import SS_TYPE_DISPLAY_ORDER, T5_HITCH_TAG, T5_SELF_TAG, display_type  # noqa: E402

# Semantic theme keys. Figure scripts read colours from here, never inline.
THEME = {
    "DLP": "#3F8E8C",  # DeepLocPro (teal)
    "DSE": "#E0884B",  # DeepSecE (amber)
    "SignalP": "#9B5DA0",  # SignalP (magenta-purple)
    "COMBINED": "#7E6BA8",  # DLP-or-DSE combined predictor (muted purple)
    "observed": "#C0392B",  # observed-value line
    "null_mean": "#333333",  # null-mean reference line
    "neutral_bar": "#6C8EAD",  # single-hue bars (no categorical meaning)
    "ref_line": "#444444",  # cutoff / zero reference lines
    "muted": "#B5B5B8",  # de-emphasised (non-significant) bars/elements
    "caption": "#666666",  # sub-axis caption / annotation text
}

# Canonical SS-type order comes from constants (single source of truth). A type
# keeps its colour across every figure and every run, whichever subset is present;
# variant labels (pT4SSt, T6SSi/T6SSii/...) inherit their parent family's colour
# via display_type().
SS_TYPE_ORDER = SS_TYPE_DISPLAY_ORDER
_SS_TYPE_COLORS = {
    "T1SS": "#4E79A7",
    "T2SS": "#59A14F",
    "T3SS": "#E15759",
    "T4SS": "#B07AA1",
    "T5aSS": "#F28E2B",
    "T5bSS": "#EDC948",
    "T5cSS": "#FF9DA7",
    "T6SS": "#76B7B2",
}
_FALLBACK_CYCLE = ["#9C755F", "#BAB0AC", "#86BCB6", "#D37295", "#A0CBE8"]

# Regular-figure PNGs this pipeline owns: legacy fig<N>_ plus current 0N_ / P0N_.
# Used by clear_figure_set to wipe a stale set without touching the enrichment
# figure (*_enrichment_fold.png) or any foreign PNG.
_OWNED_FIG_RE = re.compile(r"^(fig\d+_|\d{2}_|P\d{2}_).*\.png$")


def apply_house_style() -> None:
    """Set matplotlib rcParams once for a uniform look. Call at script start."""
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
            "axes.labelcolor": "#222222",
            "xtick.color": "#444444",
            "ytick.color": "#444444",
        }
    )


def muted(hex_color: str, sat_scale: float = 0.4, val_boost: float = 1.1) -> tuple:
    """Desaturated, slightly lightened version of a colour for de-emphasising a
    non-significant bar while keeping its hue (so you can still tell which tool it
    is, unlike a flat grey). Returns an RGB tuple for matplotlib.
    """
    import colorsys

    import matplotlib.colors as mcolors

    r, g, b = mcolors.to_rgb(hex_color)
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    return colorsys.hsv_to_rgb(h, s * sat_scale, min(1.0, v * val_boost))


def _split_base(t: str) -> str:
    """Strip a `(self)`/`(hitchhiker)` autotransporter-split suffix, if present. Robust to
    either separator: producers use a space (`T5aSS (self)`, annotation figures) or a
    newline (`T5aSS\n(self)`, enrichment tick labels). No SS-type name contains `(`."""
    return t.split("(", 1)[0].rstrip()


def ordered_ss_types(present) -> list:
    """Canonical SS-type display order first, then any extras alphabetically.

    Autotransporter self/hitchhiker split labels (`T5aSS (self)` / `T5aSS (hitchhiker)`)
    are grouped under their base type's slot, self before hitch, so the two columns
    sit adjacent.
    """
    present = list(dict.fromkeys(present))

    def rank(t: str) -> int:
        if t.endswith(T5_SELF_TAG):
            return 1
        if t.endswith(T5_HITCH_TAG):
            return 2
        return 0

    ordered = []
    for canon in SS_TYPE_ORDER:
        variants = [t for t in present if _split_base(t) == canon]
        ordered.extend(sorted(variants, key=lambda t: (rank(t), t)))
    extra = sorted(t for t in present if _split_base(t) not in SS_TYPE_ORDER)
    return ordered + extra


def ss_type_palette(types) -> dict:
    """Stable {ss_type: hex} for the given types.

    Same type maps to the same colour everywhere; unknown types get deterministic
    fallback colours assigned in sorted order (so a given input always yields the
    same mapping).
    """
    palette: dict[str, str | tuple[float, ...]] = {}
    unknown: list[str] = []
    for t in types:
        base = display_type(_split_base(t))  # strip self/hitch, collapse pT4SSt/T6SSi.. -> parent
        if base in _SS_TYPE_COLORS:
            color = _SS_TYPE_COLORS[base]
            # A `(hitchhiker)` autotransporter column shares its type's hue but muted, so
            # self vs hitchhiker are distinguishable while reading as the same type.
            palette[t] = muted(color) if t.endswith(T5_HITCH_TAG) else color
        elif t not in palette:
            unknown.append(t)
    for i, t in enumerate(sorted(dict.fromkeys(unknown))):
        palette[t] = _FALLBACK_CYCLE[i % len(_FALLBACK_CYCLE)]
    return palette


def numbered_path(outdir: str, n: int, name: str) -> str:
    """Zero-padded numbered figure path, e.g. n=1 -> ``<outdir>/01_<name>.png``."""
    return os.path.join(outdir, f"{n:02d}_{name}.png")


def pooled_path(outdir: str, n: int, name: str) -> str:
    """Pooled cross-genome figure path: ``<outdir>/P0N_<name>.png``."""
    return os.path.join(outdir, f"P{n:02d}_{name}.png")


def clear_figure_set(outdir: str) -> None:
    """Remove regular-figure PNGs this pipeline owns (legacy ``fig<N>_`` and
    numbered ``0N_`` / ``P0N_``) so a regenerated set never mixes with a stale
    one. Leaves the enrichment figure and any non-figure PNG untouched.
    """
    if not os.path.isdir(outdir):
        return
    for fn in os.listdir(outdir):
        if _OWNED_FIG_RE.match(fn):
            try:
                os.remove(os.path.join(outdir, fn))
            except OSError:
                pass


def print_figure_index(entries, logger=None) -> None:
    """Print a numbered 'Figure index' so figures can be referred to by number.

    ``entries``: iterable of ``(label, filename, description)`` where ``label`` is
    the figure's prefix (e.g. ``"01"`` or ``"P02"``).
    """
    lines = ["Figure index:"]
    for label, fname, desc in entries:
        lines.append(f"  {label}  {fname}  - {desc}")
    msg = "\n".join(lines)
    if logger is not None:
        logger.info(msg)
    else:
        print(msg)
