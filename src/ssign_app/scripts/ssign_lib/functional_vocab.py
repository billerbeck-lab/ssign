#!/usr/bin/env python3
"""Controlled functional vocabularies for the figure generator.

Two reference tables, kept out of the plotting code so they can be tested and
reused:

- ``COG_CATEGORY_NAMES`` maps the single-letter COG functional-category codes
  that EggNOG-mapper emits (``cog_category``) to their standard human-readable
  names, so figures show "Cell wall/membrane/envelope biogenesis" instead of
  "M". A protein may carry several letters (e.g. ``EGP``); ``cog_category_names``
  splits them.
- ``consensus_bucket`` collapses a free-text functional annotation (broad
  consensus / EggNOG description) to a small curated vocabulary, routing
  secretion-machinery hits to a labelled "Apparatus-associated" bucket rather
  than dropping them.
"""

from __future__ import annotations

import functools
import gzip
import os

# Reuse the live consensus keyword vocabulary (single source of truth) so the
# figure's curated buckets track any category the consensus voter gains/renames.
try:  # ssign_lib is top-level when the script runs (scripts dir on sys.path)...
    from annotation_consensus import CATEGORY_PATTERNS
except ImportError:  # ...and a subpackage when imported as ssign_app.scripts.*
    from ssign_app.scripts.annotation_consensus import CATEGORY_PATTERNS

# Tokens that mean "no value" across the integrated CSV (EggNOG uses "-").
_MISSING = ("", "-", "nan", "none")


def is_missing(value) -> bool:
    """True for blank / ``-`` / NaN-like cells (one rule for every vocabulary)."""
    return ("" if value is None else str(value).strip()).lower() in _MISSING


# Standard COG functional categories (the 1-letter code -> name table EggNOG and
# NCBI share). X (Mobilome) is an EggNOG addition to the original 25.
COG_CATEGORY_NAMES = {
    # Information storage and processing
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    # Cellular processes and signaling
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    # Metabolism
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    # Poorly characterized
    "R": "General function prediction only",
    "S": "Function unknown",
    # EggNOG addition
    "X": "Mobilome: prophages, transposons",
}

_COG_UNANNOTATED = "Unannotated"


def cog_category_names(code) -> list:
    """Letter codes in one ``cog_category`` value -> category name(s).

    EggNOG packs every category a protein matches into one string (``"EGP"``);
    each letter becomes a separate name so the protein counts once per category.
    Empty / ``-`` / missing -> ``["Unannotated"]``. Unknown letters map to
    ``"Function unknown"`` (COG S) rather than being silently dropped.
    """
    if is_missing(code):
        return [_COG_UNANNOTATED]
    names = []
    for ch in str(code).strip():
        if ch in ("-", " "):
            continue
        names.append(COG_CATEGORY_NAMES.get(ch.upper(), COG_CATEGORY_NAMES["S"]))
    return names or [_COG_UNANNOTATED]


def kegg_kos(value) -> list:
    """One ``kegg_ko`` value -> KEGG Orthology IDs.

    EggNOG emits ``"ko:K10953;ko:K12516"``; returns ``["K10953", "K12516"]``
    (the ``ko:`` prefix dropped). Empty / ``-`` -> ``[]``.
    """
    if is_missing(value):
        return []
    out = []
    for tok in str(value).strip().split(";"):
        ko = tok.strip().removeprefix("ko:").strip()
        if not is_missing(ko):
            out.append(ko)
    return out


# Bundled KEGG KO -> definition table (rest.kegg.jp/list/ko, gzipped). Lazy-loaded
# once (lru_cache) so figures can show "RTX toxin RtxA" instead of "K10953".
_KEGG_KO_TABLE = os.path.join(os.path.dirname(__file__), "data", "kegg_ko_names.tsv.gz")


@functools.lru_cache(maxsize=1)
def _kegg_ko_names() -> dict:
    names: dict = {}
    try:
        with gzip.open(_KEGG_KO_TABLE, "rt", encoding="utf-8") as fh:
            for line in fh:
                ko, _, defn = line.partition("\t")
                defn = defn.strip()
                if ko and defn:
                    names[ko.strip().upper()] = defn
    except OSError:
        pass  # table not bundled -> fall back to raw KO IDs
    return names


def kegg_ko_name(ko_id: str) -> str:
    """KO id ('K10953') -> short human function name, or the id if unmapped.

    KEGG definitions read "gene1, gene2; function name [EC:...]"; drop the leading
    gene-symbol list and the trailing EC bracket so the label is the function.
    """
    defn = _kegg_ko_names().get(ko_id.strip().upper())
    if not defn:
        return ko_id
    name = defn.split(";", 1)[1].strip() if ";" in defn else defn
    return name.split(" [", 1)[0].strip() or defn


def kegg_descriptions(value) -> list:
    """One ``kegg_ko`` cell -> human KO function names (deduped, order-preserving)."""
    out, seen = [], set()
    for ko in kegg_kos(value):
        name = kegg_ko_name(ko)
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


# The consensus voter's "Secretion system" / "Flagellar" categories are the
# secretion/uptake machinery itself, not cargo: in a secreted-protein functional
# chart they are routed to one labelled bucket rather than read as a function.
_APPARATUS_CATEGORIES = frozenset({"Secretion system", "Flagellar"})
_KNOWN_CATEGORIES = frozenset(name for name, _ in CATEGORY_PATTERNS)
APPARATUS_BUCKET = "Apparatus-associated"
_OTHER_BUCKET = "Other"


def consensus_bucket(value) -> str:
    """Collapse a ``broad_annotation`` / ``broad_consensus_annotation`` value to a
    curated functional bucket.

    Known consensus categories pass through; machinery (``Secretion system`` /
    ``Flagellar``) -> ``"Apparatus-associated"``; blanks and the voter's
    first-3-words fallback noise (anything not a known category) -> ``"Other"``.
    """
    if is_missing(value):
        return _OTHER_BUCKET
    cat = str(value).strip().split(" (", 1)[0].strip()  # drop the "(Tool1, Tool2)" consensus suffix
    if cat in _APPARATUS_CATEGORIES:
        return APPARATUS_BUCKET
    return cat if cat in _KNOWN_CATEGORIES else _OTHER_BUCKET
