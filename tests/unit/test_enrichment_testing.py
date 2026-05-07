"""Tests for enrichment_testing.py.

Three pure-Python surfaces here:

1. `get_functional_category` — priority-ordered fallback through five
   annotation fields. Must skip empty / `-` / "hypothetical" /
   "uncharacterized" sentinels and land on `"Unknown"` if nothing real.
2. `fishers_exact_enrichment` — Fisher's exact (greater) per
   SS-type x category cell, plus Benjamini-Hochberg FDR correction.
3. `circular_permutation_test` — permutation null where category
   labels are shuffled across substrates; deterministic with a seed.

scipy is required for #2; `pytest.importorskip` would mask a missing
dep — instead we assert it's present in the dev env (it is, via the
`[dev]` extras) and let the import error fail loudly if that ever
breaks.
"""

import os
import sys

import pytest


from _helpers import write_tsv  # noqa: E402
from enrichment_testing import (  # noqa: E402
    circular_permutation_test,
    fishers_exact_enrichment,
    get_functional_category,
    load_integrated_csv,
)


def _substrate(
    locus_tag="GENE_001",
    nearby="T1SS",
    broad_annotation="",
    blastp_hit_description="",
    pdb_top1_description="",
    pfam_top1_description="",
    interpro_descriptions="",
):
    """Build a minimal substrate row with the priority-ordered annotation
    fields explicit. Defaults: only nearby_ss_types is populated."""
    return {
        "locus_tag": locus_tag,
        "nearby_ss_types": nearby,
        "broad_annotation": broad_annotation,
        "blastp_hit_description": blastp_hit_description,
        "pdb_top1_description": pdb_top1_description,
        "pfam_top1_description": pfam_top1_description,
        "interpro_descriptions": interpro_descriptions,
    }


# ---------------------------------------------------------------------------
# get_functional_category
# ---------------------------------------------------------------------------


class TestFunctionalCategoryPriority:
    """Priority order: broad_annotation → blastp_hit_description →
    pdb_top1_description → pfam_top1_description → interpro_descriptions →
    'Unknown'. Sentinels filtered: '', '-', 'hypothetical protein',
    'uncharacterized protein' (case-insensitive)."""

    def test_broad_annotation_wins(self):
        assert (
            get_functional_category(
                _substrate(
                    broad_annotation="adhesin",
                    blastp_hit_description="hemolysin",  # would lose
                )
            )
            == "adhesin"
        )

    def test_falls_through_to_blastp(self):
        assert (
            get_functional_category(
                _substrate(
                    broad_annotation="",
                    blastp_hit_description="hemolysin",
                )
            )
            == "hemolysin"
        )

    def test_falls_through_to_pdb(self):
        row = _substrate(pdb_top1_description="autotransporter beta-domain")
        assert get_functional_category(row) == "autotransporter beta-domain"

    def test_falls_through_to_pfam(self):
        row = _substrate(pfam_top1_description="Ig-like fold")
        assert get_functional_category(row) == "Ig-like fold"

    def test_falls_through_to_interpro(self):
        row = _substrate(interpro_descriptions="Outer membrane domain")
        assert get_functional_category(row) == "Outer membrane domain"

    def test_all_empty_returns_unknown(self):
        assert get_functional_category(_substrate()) == "Unknown"

    @pytest.mark.parametrize(
        "sentinel",
        ["", "-", "hypothetical protein", "Hypothetical Protein", "uncharacterized protein", "UNCHARACTERIZED PROTEIN"],
    )
    def test_sentinels_skipped(self, sentinel):
        # All fields are sentinels except the lowest-priority one that has a
        # real value — make sure the parser falls all the way through.
        row = _substrate(
            broad_annotation=sentinel,
            blastp_hit_description=sentinel,
            pdb_top1_description=sentinel,
            pfam_top1_description=sentinel,
            interpro_descriptions="hemolysin",
        )
        assert get_functional_category(row) == "hemolysin"

    def test_whitespace_only_treated_as_empty(self):
        # `.strip()` drops "    " → empty → sentinel-like behaviour
        row = _substrate(broad_annotation="   ", blastp_hit_description="real value")
        assert get_functional_category(row) == "real value"


# ---------------------------------------------------------------------------
# load_integrated_csv
# ---------------------------------------------------------------------------


class TestLoadIntegratedCsv:
    def test_round_trips_dict_rows(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "integrated.csv"),
            ["locus_tag", "nearby_ss_types", "broad_annotation"],
            [
                {"locus_tag": "GENE_001", "nearby_ss_types": "T1SS", "broad_annotation": "adhesin"},
                {"locus_tag": "GENE_002", "nearby_ss_types": "T2SS", "broad_annotation": ""},
            ],
            delimiter=",",
        )
        rows = load_integrated_csv(path)
        assert len(rows) == 2
        assert rows[0]["locus_tag"] == "GENE_001"

    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.csv")
        open(path, "w").close()
        assert load_integrated_csv(path) == []


# ---------------------------------------------------------------------------
# fishers_exact_enrichment
# ---------------------------------------------------------------------------


class TestFishersExactEnrichment:
    def test_strong_enrichment_yields_low_pvalue(self):
        # 10 T1SS substrates all annotated "hemolysin"; 10 T2SS substrates all
        # annotated "starch-binding". Cells should be perfectly diagonal —
        # Fisher's exact (greater) should give a tiny p-value for both.
        subs = [_substrate(locus_tag=f"T1_{i}", nearby="T1SS", broad_annotation="hemolysin") for i in range(10)]
        subs += [_substrate(locus_tag=f"T2_{i}", nearby="T2SS", broad_annotation="starch") for i in range(10)]
        results = fishers_exact_enrichment(subs)

        # Two SS x category cells with non-zero observation; cross-cells (T1SS
        # x starch, T2SS x hemolysin) have a=0 and are dropped.
        assert len(results) == 2
        for r in results:
            assert r["pvalue"] < 0.001
            assert r["significant"]

    def test_no_enrichment_no_signal(self):
        # Categories are uniformly distributed across SS types; OR ≈ 1.
        subs = [_substrate(locus_tag=f"T1_h_{i}", nearby="T1SS", broad_annotation="hemolysin") for i in range(5)]
        subs += [_substrate(locus_tag=f"T1_s_{i}", nearby="T1SS", broad_annotation="starch") for i in range(5)]
        subs += [_substrate(locus_tag=f"T2_h_{i}", nearby="T2SS", broad_annotation="hemolysin") for i in range(5)]
        subs += [_substrate(locus_tag=f"T2_s_{i}", nearby="T2SS", broad_annotation="starch") for i in range(5)]
        results = fishers_exact_enrichment(subs)
        # All 4 cells have a > 0 and OR ≈ 1.0; none should be significant.
        for r in results:
            assert not r["significant"]

    def test_zero_observation_cells_dropped(self):
        # Only one diagonal cell has a > 0; the off-diagonal cell is dropped.
        subs = [
            _substrate(locus_tag="g1", nearby="T1SS", broad_annotation="hemolysin"),
            _substrate(locus_tag="g2", nearby="T2SS", broad_annotation="starch"),
        ]
        results = fishers_exact_enrichment(subs)
        # 2 SS x 2 cat = 4 cells; 2 of those have a=0 and must be skipped.
        assert len(results) == 2

    def test_bh_rank_assignment(self):
        # Sorted by p-value, rank should be 1, 2, 3, ...
        subs = []
        for i in range(5):
            subs.append(_substrate(locus_tag=f"T1_{i}", nearby="T1SS", broad_annotation="hemolysin"))
            subs.append(_substrate(locus_tag=f"T2_{i}", nearby="T2SS", broad_annotation="starch"))
        results = fishers_exact_enrichment(subs)
        ranks = [r["bh_rank"] for r in results]
        assert ranks == sorted(ranks)
        assert ranks[0] == 1

    def test_multiple_ss_types_per_substrate(self):
        # nearby_ss_types is comma-separated — each SS type counts independently
        subs = [
            _substrate(locus_tag="g1", nearby="T1SS,T2SS", broad_annotation="adhesin"),
        ]
        results = fishers_exact_enrichment(subs)
        ss_types = {r["ss_type"] for r in results}
        assert ss_types == {"T1SS", "T2SS"}

    def test_empty_input_returns_empty(self):
        assert fishers_exact_enrichment([]) == []

    def test_blank_nearby_skipped(self):
        # Substrate with no nearby_ss_types contributes nothing
        subs = [_substrate(locus_tag="g1", nearby="", broad_annotation="hemolysin")]
        assert fishers_exact_enrichment(subs) == []


# ---------------------------------------------------------------------------
# circular_permutation_test
# ---------------------------------------------------------------------------


class TestCircularPermutationTest:
    def test_deterministic_with_fixed_seed(self):
        subs = [_substrate(locus_tag=f"g{i}", nearby="T1SS", broad_annotation="hemolysin") for i in range(5)] + [
            _substrate(locus_tag=f"g{i}", nearby="T2SS", broad_annotation="starch") for i in range(5, 10)
        ]
        # Same seed → same p-values
        r1 = circular_permutation_test(subs, gene_orders={}, n_perms=100, seed=42)
        r2 = circular_permutation_test(subs, gene_orders={}, n_perms=100, seed=42)
        assert r1 == r2

    def test_different_seeds_differ(self):
        # Need enough substrates that the shuffle produces different counts
        subs = []
        for i in range(20):
            ss = "T1SS" if i < 10 else "T2SS"
            cat = "hemolysin" if i % 2 == 0 else "starch"
            subs.append(_substrate(locus_tag=f"g{i}", nearby=ss, broad_annotation=cat))
        r1 = circular_permutation_test(subs, gene_orders={}, n_perms=50, seed=1)
        r2 = circular_permutation_test(subs, gene_orders={}, n_perms=50, seed=999)
        # At least one cell's p-value should differ across seeds
        assert any(a["perm_pvalue"] != b["perm_pvalue"] for a, b in zip(r1, r2))

    def test_no_observations_returns_empty(self):
        # Every substrate has a blank nearby_ss_types → observed{} is empty
        subs = [_substrate(locus_tag="g1", nearby="", broad_annotation="x")]
        assert circular_permutation_test(subs, gene_orders={}, n_perms=10) == []

    def test_pvalue_smoothed_minimum(self):
        # +1 smoothing means p-value floor is 1 / (n_perms + 1), never zero.
        subs = [_substrate(locus_tag=f"g{i}", nearby="T1SS", broad_annotation="hemolysin") for i in range(5)]
        results = circular_permutation_test(subs, gene_orders={}, n_perms=10, seed=42)
        assert all(r["perm_pvalue"] >= 1 / 11 for r in results)

    def test_n_permutations_recorded_per_result(self):
        subs = [_substrate(locus_tag="g1", nearby="T1SS", broad_annotation="hemolysin")]
        results = circular_permutation_test(subs, gene_orders={}, n_perms=42, seed=42)
        assert all(r["n_permutations"] == 42 for r in results)

    def test_observed_count_correctly_reported(self):
        # 3 substrates of T1SS-hemolysin → observed count = 3
        subs = [_substrate(locus_tag=f"g{i}", nearby="T1SS", broad_annotation="hemolysin") for i in range(3)]
        results = circular_permutation_test(subs, gene_orders={}, n_perms=10, seed=42)
        cell = next(r for r in results if r["ss_type"] == "T1SS" and r["category"] == "hemolysin")
        assert cell["observed"] == 3
