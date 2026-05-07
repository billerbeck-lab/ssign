"""Tests for ssign_lib/go_utils.py — GO-slim functional categorization.

Pure-Python surfaces here:

1. `_keyword_fallback` — annotation-text matching when GO terms are
   absent. Categories are joined pipe-delimited; no match → "Unknown".
2. `merge_go_terms_true_path` — True Path Rule: when two GO terms are
   in an ancestor-descendant relationship, keep the descendant
   (more specific) and drop the ancestor.
3. `categorize_protein` — orchestrator; falls through to keyword
   matching when no GO terms or no slim mapping found.
4. `BROAD_CATEGORY_MAP` / `FALLBACK_KEYWORDS` — pinning the shape of
   the canonical lookup tables (no broken IDs / empty values).

OBO download + `goatools.mapslim` paths require network and large data
files; covered by `tests/integration/`.
"""

import os
import sys

import pytest

# go_utils raises RuntimeError on import when obonet or goatools are missing.
# These ship in [extended] extras; in the basic [dev] env they're absent.
pytest.importorskip("obonet")
pytest.importorskip("goatools")


import networkx  # noqa: E402
from ssign_lib.go_utils import (  # noqa: E402
    BROAD_CATEGORY_MAP,
    FALLBACK_KEYWORDS,
    _keyword_fallback,
    categorize_protein,
    merge_go_terms_true_path,
)

# ---------------------------------------------------------------------------
# _keyword_fallback
# ---------------------------------------------------------------------------


class TestKeywordFallback:
    def test_empty_returns_unknown(self):
        result = _keyword_fallback("")
        assert result["func_category_broad"] == "Unknown"
        assert result["categorization_source"] == "keyword_fallback"

    def test_whitespace_only_returns_unknown(self):
        assert _keyword_fallback("   ")["func_category_broad"] == "Unknown"

    def test_no_match_returns_unknown(self):
        assert _keyword_fallback("zorblax bizango")["func_category_broad"] == "Unknown"

    @pytest.mark.parametrize(
        "text, expected_category",
        [
            ("serine protease", "Catalytic"),
            ("ABC transporter", "Transport"),
            ("DNA-binding protein", "Binding"),
            ("flagellin subunit", "Structural"),
            ("GGDEF domain protein", "Signaling"),
            ("hemolysin A", "Virulence"),
            ("biosynthesis pathway", "Metabolism"),
        ],
    )
    def test_keyword_maps_to_broad(self, text, expected_category):
        assert _keyword_fallback(text)["func_category_broad"] == expected_category

    def test_case_insensitive(self):
        assert _keyword_fallback("PROTEASE")["func_category_broad"] == "Catalytic"
        assert _keyword_fallback("ProTeAsE")["func_category_broad"] == "Catalytic"

    def test_multiple_categories_pipe_joined_sorted(self):
        # An annotation that hits both Virulence (toxin) and Catalytic (protease).
        result = _keyword_fallback("toxin protease complex")
        broad = result["func_category_broad"]
        assert "Catalytic" in broad
        assert "Virulence" in broad
        # Sorted: "C" < "V"
        parts = broad.split("|")
        assert parts == sorted(parts)

    def test_specific_and_detail_empty_for_keyword_path(self):
        # Keyword fallback can't supply slim names or detail GO IDs.
        result = _keyword_fallback("protease")
        assert result["func_category_specific"] == ""
        assert result["func_category_detail"] == ""


# ---------------------------------------------------------------------------
# merge_go_terms_true_path
# ---------------------------------------------------------------------------


def _build_chain_graph(terms):
    """Build a chain graph: terms[0] (most specific) → ... → terms[-1] (root).

    Edges go child → parent, matching obonet's convention. So
    `networkx.descendants(graph, child)` returns all ancestors.
    """
    g = networkx.MultiDiGraph()
    for child, parent in zip(terms, terms[1:]):
        g.add_edge(child, parent)
    return g


class TestMergeGoTermsTruePath:
    def test_redundant_ancestor_dropped(self):
        # Chain: A → B → C (A is most specific). Inputs {A, B, C} should
        # collapse to just {A}.
        graph = _build_chain_graph(["A", "B", "C"])
        result = merge_go_terms_true_path({"A", "C"}, {"B"}, graph)
        assert result == {"A"}

    def test_disconnected_terms_both_kept(self):
        # Two unrelated chains: keep both leaves.
        g = networkx.MultiDiGraph()
        g.add_edge("A", "B")
        g.add_edge("X", "Y")
        result = merge_go_terms_true_path({"A"}, {"X"}, g)
        assert result == {"A", "X"}

    def test_terms_not_in_graph_dropped(self):
        graph = _build_chain_graph(["A", "B"])
        result = merge_go_terms_true_path({"A", "ORPHAN"}, set(), graph)
        assert result == {"A"}

    def test_set_union_semantics(self):
        # interpro_terms ∪ other_terms should be honoured (no duplicates).
        graph = _build_chain_graph(["A", "B"])
        result = merge_go_terms_true_path({"A"}, {"A", "B"}, graph)
        assert result == {"A"}

    def test_empty_inputs_yield_empty(self):
        g = networkx.MultiDiGraph()
        assert merge_go_terms_true_path(set(), set(), g) == set()

    def test_single_term_unchanged(self):
        graph = _build_chain_graph(["A", "B"])
        assert merge_go_terms_true_path({"A"}, set(), graph) == {"A"}


# ---------------------------------------------------------------------------
# categorize_protein — keyword-fallback branch
# ---------------------------------------------------------------------------
# (The GO-slim branch requires real GODag instances and is exercised
#  by tests/integration/.)


class TestCategorizeProteinFallback:
    def test_no_go_terms_falls_to_keyword(self):
        # go_terms=[] forces the fallback path; the dag args are unused.
        result = categorize_protein(
            go_terms=[],
            go_dag=None,  # unused on this branch
            slim_dag=None,
            annotation_text="ABC transporter",
        )
        assert result["categorization_source"] == "keyword_fallback"
        assert result["func_category_broad"] == "Transport"

    def test_blank_go_terms_treated_as_empty(self):
        # Strings of whitespace are stripped → no valid_terms → fallback path.
        result = categorize_protein(
            go_terms=["", "   ", "\t"],
            go_dag=None,
            slim_dag=None,
            annotation_text="hemolysin",
        )
        assert result["categorization_source"] == "keyword_fallback"
        assert result["func_category_broad"] == "Virulence"

    def test_no_terms_no_annotation_unknown(self):
        result = categorize_protein(
            go_terms=[],
            go_dag=None,
            slim_dag=None,
            annotation_text="",
        )
        assert result["func_category_broad"] == "Unknown"


# ---------------------------------------------------------------------------
# BROAD_CATEGORY_MAP and FALLBACK_KEYWORDS pinning
# ---------------------------------------------------------------------------
# The constants drive every downstream broad-category assignment. If
# someone refactors a value, every downstream report shifts. Pin the shape.


class TestBroadCategoryMap:
    # Documented broad categories (cells in the report). New additions are
    # fine — but every value must be in this allow-list, otherwise the
    # downstream per-category reports will silently drop the unknown bucket.
    KNOWN_CATEGORIES = {
        "Catalytic",
        "Regulation",
        "Transport",
        "Binding",
        "Signaling",
        "Structural",
        "Stress Response",
        "Virulence",
        "Metabolism",
        "Extracellular",
        "Membrane-associated",
        "Periplasmic",
    }

    def test_every_value_is_known_category(self):
        unknown = set(BROAD_CATEGORY_MAP.values()) - self.KNOWN_CATEGORIES
        assert unknown == set(), f"Unexpected categories: {unknown}"

    def test_all_keys_are_go_ids(self):
        # Every key must look like "GO:NNNNNNN".
        for k in BROAD_CATEGORY_MAP:
            assert k.startswith("GO:"), f"{k!r} does not look like a GO ID"
            assert k[3:].isdigit(), f"{k!r} body is not all digits"

    def test_no_empty_values(self):
        assert all(v.strip() for v in BROAD_CATEGORY_MAP.values())


class TestFallbackKeywords:
    def test_categories_overlap_broad_map(self):
        # Every fallback category should also be a recognised broad category.
        # Otherwise the keyword fallback emits a category that disappears
        # from downstream per-category reports.
        unknown = set(FALLBACK_KEYWORDS.keys()) - TestBroadCategoryMap.KNOWN_CATEGORIES
        assert unknown == set(), f"Fallback uses unknown categories: {unknown}"

    def test_no_empty_keyword_lists(self):
        assert all(v for v in FALLBACK_KEYWORDS.values())

    def test_no_duplicate_keywords_within_category(self):
        for cat, kws in FALLBACK_KEYWORDS.items():
            assert len(kws) == len(set(kws)), f"Duplicate keyword in {cat}"
