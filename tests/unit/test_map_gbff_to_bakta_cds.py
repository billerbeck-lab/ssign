"""Unit tests for map_gbff_to_bakta_cds.py.

Exercises the pure-Python reciprocal-overlap matcher that carries user
GenBank `product` annotations across to Bakta's fresh CDS coordinates
for annotation-consensus voting.
"""

import os
import sys

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from map_gbff_to_bakta_cds import (  # noqa: E402, F401
    _read_gene_info,
    _reciprocal_overlap,
    best_gbff_match,
    map_gbff_to_bakta,
)


def _cds(locus_tag, contig, start, end, strand, product=""):
    """Mint a minimal CDS dict matching the output shape of _read_gene_info."""
    return {
        "locus_tag": locus_tag,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "product": product,
    }


class TestReciprocalOverlap:
    """The core geometry: min(overlap/len_a, overlap/len_b) for two CDS."""

    def test_identical_ranges_return_1_0(self):
        assert _reciprocal_overlap(100, 200, 100, 200) == 1.0

    def test_no_overlap_returns_0(self):
        assert _reciprocal_overlap(100, 200, 300, 400) == 0.0

    def test_abutting_not_overlapping(self):
        """Ranges that touch but don't overlap — [100,200) and [200,300)."""
        assert _reciprocal_overlap(100, 200, 200, 300) == 0.0

    def test_a_fully_inside_b_reciprocal_is_shorter_over_longer(self):
        # a = [100,150) len 50; b = [100,200) len 100. Overlap = 50.
        # reciprocal = min(50/50, 50/100) = 0.5
        assert _reciprocal_overlap(100, 150, 100, 200) == 0.5

    def test_partial_overlap_midway(self):
        # a = [100,200); b = [150,250); overlap = 50; each len 100.
        # reciprocal = min(50/100, 50/100) = 0.5
        assert _reciprocal_overlap(100, 200, 150, 250) == 0.5

    def test_near_identical_small_shift(self):
        # a = [100,200); b = [110,210); overlap = 90; each len 100.
        # reciprocal = min(90/100, 90/100) = 0.9
        assert _reciprocal_overlap(100, 200, 110, 210) == 0.9

    def test_degenerate_zero_length_returns_0(self):
        assert _reciprocal_overlap(100, 100, 100, 200) == 0.0


class TestBestGbffMatch:
    """Pick the best-overlapping GenBank CDS on the same contig + strand."""

    def test_single_good_match(self):
        bakta = _cds("B1", "c1", 100, 200, "+")
        candidates = [_cds("G1", "c1", 95, 205, "+", "hemolysin")]
        match = best_gbff_match(bakta, candidates, min_overlap=0.8)
        assert match is not None
        assert match["product"] == "hemolysin"

    def test_best_of_many_wins(self):
        """Three candidates overlap; the one with highest reciprocal wins."""
        bakta = _cds("B1", "c1", 100, 200, "+")
        candidates = [
            _cds("G_low", "c1", 50, 150, "+", "low_overlap"),
            _cds("G_best", "c1", 100, 200, "+", "perfect"),
            _cds("G_mid", "c1", 120, 220, "+", "mid_overlap"),
        ]
        # candidates must be sorted by start for best_gbff_match to work
        candidates.sort(key=lambda r: r["start"])
        match = best_gbff_match(bakta, candidates, min_overlap=0.5)
        assert match["product"] == "perfect"

    def test_below_threshold_returns_none(self):
        bakta = _cds("B1", "c1", 100, 200, "+")
        # overlap = 10, lens 100 and 1000; reciprocal = 0.01
        candidates = [_cds("G1", "c1", 190, 1190, "+", "barely_touches")]
        match = best_gbff_match(bakta, candidates, min_overlap=0.8)
        assert match is None

    def test_strand_mismatch_rejected_by_default(self):
        """Same bp range but opposite strand is a different gene."""
        bakta = _cds("B1", "c1", 100, 200, "+")
        candidates = [_cds("G1", "c1", 100, 200, "-", "antisense_gene")]
        match = best_gbff_match(bakta, candidates, min_overlap=0.8)
        assert match is None

    def test_strand_mismatch_accepted_when_strand_check_disabled(self):
        bakta = _cds("B1", "c1", 100, 200, "+")
        candidates = [_cds("G1", "c1", 100, 200, "-", "whatever")]
        match = best_gbff_match(
            bakta, candidates, min_overlap=0.8, strand_must_match=False
        )
        assert match is not None

    def test_empty_candidates_returns_none(self):
        assert best_gbff_match(_cds("B1", "c1", 100, 200, "+"), []) is None

    def test_candidates_past_bakta_end_terminate_early(self):
        """When sorted candidates pass bakta's end, the search short-circuits."""
        bakta = _cds("B1", "c1", 100, 200, "+")
        candidates = [
            _cds("G1", "c1", 120, 180, "+", "inside"),
            _cds("G2", "c1", 500, 600, "+", "far_away"),  # past bakta end
        ]
        match = best_gbff_match(bakta, candidates, min_overlap=0.5)
        assert match["product"] == "inside"


class TestMapGbffToBakta:
    """Full pipeline: build the annotated Bakta rows from both gene-info dicts."""

    def test_perfect_match_carries_annotation(self, tmp_dir):
        bakta = {
            "c1": [_cds("B1", "c1", 100, 200, "+")],
        }
        gbff = {
            "c1": [_cds("G1", "c1", 100, 200, "+", "alpha-hemolysin")],
        }
        rows = list(map_gbff_to_bakta(bakta, gbff))
        assert rows[0]["gbff_annotation"] == "alpha-hemolysin"

    def test_no_match_leaves_empty_annotation(self, tmp_dir):
        bakta = {"c1": [_cds("B1", "c1", 100, 200, "+")]}
        gbff = {"c1": [_cds("G1", "c1", 500, 600, "+", "far_gene")]}
        rows = list(map_gbff_to_bakta(bakta, gbff))
        assert rows[0]["gbff_annotation"] == ""

    def test_different_contig_no_match(self, tmp_dir):
        """Same coordinates but different contig — must not match."""
        bakta = {"c1": [_cds("B1", "c1", 100, 200, "+")]}
        gbff = {"c2": [_cds("G1", "c2", 100, 200, "+", "unrelated")]}
        rows = list(map_gbff_to_bakta(bakta, gbff))
        assert rows[0]["gbff_annotation"] == ""

    def test_preserves_bakta_order_and_locus_tags(self, tmp_dir):
        bakta = {
            "c1": [
                _cds("B1", "c1", 100, 200, "+"),
                _cds("B2", "c1", 300, 400, "+"),
            ]
        }
        gbff = {
            "c1": [
                _cds("G1", "c1", 100, 200, "+", "first"),
                _cds("G2", "c1", 300, 400, "+", "second"),
            ]
        }
        rows = list(map_gbff_to_bakta(bakta, gbff))
        assert [r["locus_tag"] for r in rows] == ["B1", "B2"]
        assert rows[0]["gbff_annotation"] == "first"
        assert rows[1]["gbff_annotation"] == "second"


class TestReadGeneInfo:
    def test_reads_valid_rows(self, tmp_dir):
        path = os.path.join(tmp_dir, "gene_info.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tcontig\tstart\tend\tstrand\tproduct\n")
            f.write("G1\tcontig_1\t100\t200\t+\themolysin\n")
            f.write("G2\tcontig_1\t300\t400\t-\tautotransporter\n")
            f.write("G3\tcontig_2\t50\t150\t+\ttoxin\n")

        result = _read_gene_info(path)
        assert set(result.keys()) == {"contig_1", "contig_2"}
        assert len(result["contig_1"]) == 2
        assert len(result["contig_2"]) == 1

    def test_sorts_by_start_within_each_contig(self, tmp_dir):
        """Output must be sorted — best_gbff_match relies on early-termination."""
        path = os.path.join(tmp_dir, "gene_info.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tcontig\tstart\tend\tstrand\tproduct\n")
            f.write("late\tc1\t1000\t1100\t+\tfoo\n")
            f.write("early\tc1\t100\t200\t+\tbar\n")
            f.write("mid\tc1\t500\t600\t+\tbaz\n")

        result = _read_gene_info(path)
        starts = [r["start"] for r in result["c1"]]
        assert starts == sorted(starts)

    def test_row_with_bad_coords_skipped(self, tmp_dir):
        path = os.path.join(tmp_dir, "gene_info.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tcontig\tstart\tend\tstrand\tproduct\n")
            f.write("G1\tc1\t100\t200\t+\tgood\n")
            f.write("G2\tc1\tnotanumber\t200\t+\tbad\n")
            f.write("G3\tc1\t300\t400\t+\talso_good\n")

        result = _read_gene_info(path)
        assert len(result["c1"]) == 2
