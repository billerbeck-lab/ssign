"""Tests for run_blastp.py.

Two pure-Python surfaces here:

1. `parse_blast_tabular` — outfmt-6 parser. Per-query bitscore re-sort
   (BLAST's row order is DB-traversal order, not score order, after
   widening `-max_target_seqs`). NCBI concatenates redundant subject
   titles with " >"; the split must happen *before* EXCLUDE_TERMS,
   otherwise a real hit reading "hemolysin >X hypothetical protein"
   gets wrongly dropped — Critical Bug Fix #3.
2. `filter_hits` — pident, qcov, and EXCLUDE_TERMS gating.

The `run_local_blastp` subprocess path requires NCBI BLAST+ on PATH and
is exercised by tests/integration/.
"""

import os
import sys

import pytest


from _helpers import BLAST_OUTFMT_COLS, make_blast_outfmt_row  # noqa: E402
from run_blastp import (  # noqa: E402
    BLAST_OUTFMT,
    EXCLUDE_TERMS,
    filter_hits,
    parse_blast_tabular,
)


@pytest.mark.parametrize(
    "stitle, expected_desc",
    [
        ("hemolysin BL21-like protein >WP_X.1 hypothetical protein", "hemolysin BL21-like protein"),
        ("hemolysin >ID_2 something >ID_3 hypothetical protein", "hemolysin"),
        ("hemolysin", "hemolysin"),
    ],
)
def test_ncbi_title_split_takes_first_segment(stitle, expected_desc):
    hits = parse_blast_tabular(make_blast_outfmt_row(stitle=stitle))
    assert hits["GENE_0001"]["blastp_hit_description"] == expected_desc


def test_description_truncated_to_200_chars():
    long_desc = "x" * 500
    hits = parse_blast_tabular(
        make_blast_outfmt_row(stitle=f"{long_desc} >ID hypothetical protein"),
    )
    assert len(hits["GENE_0001"]["blastp_hit_description"]) == 200


class TestFilterRespectsSplitDescription:
    """Without the split, a real hit whose stitle ends in
    "...>X hypothetical protein" would be wrongly excluded — the unsplit
    string contains "hypothetical protein", an EXCLUDE_TERM."""

    def test_real_hit_kept_when_redundant_tail_says_hypothetical(self):
        hits = parse_blast_tabular(
            make_blast_outfmt_row(stitle="hemolysin >WP_X.1 hypothetical protein"),
        )
        assert "GENE_0001" in filter_hits(hits, min_pident=80, min_qcov=80)

    def test_truly_hypothetical_hit_dropped(self):
        hits = parse_blast_tabular(
            make_blast_outfmt_row(stitle="hypothetical protein >WP_X.1 hemolysin"),
        )
        assert filter_hits(hits, min_pident=80, min_qcov=80) == {}


class TestBitscoreSort:
    """parse_blast_tabular re-sorts hits per query by bitscore — BLAST's row
    order is DB-traversal order, not bitscore order, so without the sort the
    wrong best hit gets reported."""

    def test_best_bitscore_wins(self):
        rows = "\n".join(
            [
                make_blast_outfmt_row(qseqid="Q", sseqid="LOWER", bitscore=100.0),
                make_blast_outfmt_row(qseqid="Q", sseqid="HIGHER", bitscore=900.0),
                make_blast_outfmt_row(qseqid="Q", sseqid="MID", bitscore=500.0),
            ]
        )
        hits = parse_blast_tabular(rows)
        assert hits["Q"]["blastp_hit_accession"] == "HIGHER"


class TestParserResilience:
    def test_empty_output_returns_empty_dict(self):
        assert parse_blast_tabular("") == {}

    def test_blank_lines_ignored(self):
        rows = "\n".join(
            [
                "",
                make_blast_outfmt_row(qseqid="Q1"),
                "",
                make_blast_outfmt_row(qseqid="Q2"),
                "",
            ]
        )
        assert set(parse_blast_tabular(rows).keys()) == {"Q1", "Q2"}

    def test_short_rows_skipped(self):
        # Fewer than the 15 required fields → silently dropped (corrupt output)
        assert parse_blast_tabular("Q1\tWP_X.1\t95.0\t100") == {}


@pytest.mark.parametrize("term", EXCLUDE_TERMS)
def test_every_exclude_term_filters_out(term):
    hits = parse_blast_tabular(make_blast_outfmt_row(stitle=term))
    assert filter_hits(hits, min_pident=0, min_qcov=0) == {}


class TestQueryCoverage:
    """qcov = aln_len / qlen * 100. Defensive on qlen=0."""

    @pytest.mark.parametrize(
        "aln_len, qlen, expected_qcov",
        [
            (200, 200, 100.0),  # full coverage
            (100, 200, 50.0),
            (160, 200, 80.0),
            (50, 200, 25.0),
        ],
    )
    def test_qcov_calculated_from_aln_len_and_qlen(self, aln_len, qlen, expected_qcov):
        hits = parse_blast_tabular(make_blast_outfmt_row(aln_len=aln_len, qlen=qlen))
        assert hits["GENE_0001"]["blastp_qcov"] == expected_qcov

    def test_qlen_zero_yields_zero_qcov_not_div_by_zero(self):
        hits = parse_blast_tabular(make_blast_outfmt_row(aln_len=100, qlen=0))
        assert hits["GENE_0001"]["blastp_qcov"] == 0


class TestFilterThresholds:
    def test_pident_below_threshold_dropped(self):
        hits = parse_blast_tabular(make_blast_outfmt_row(pident=70.0))
        assert filter_hits(hits, min_pident=80, min_qcov=0) == {}

    def test_pident_at_threshold_kept(self):
        # `>=` comparison — equality passes (pin against accidental flip to `>`).
        hits = parse_blast_tabular(make_blast_outfmt_row(pident=80.0))
        assert "GENE_0001" in filter_hits(hits, min_pident=80, min_qcov=0)

    def test_qcov_below_threshold_dropped(self):
        hits = parse_blast_tabular(make_blast_outfmt_row(aln_len=100, qlen=200))
        assert filter_hits(hits, min_pident=0, min_qcov=80) == {}

    def test_qcov_at_threshold_kept(self):
        hits = parse_blast_tabular(make_blast_outfmt_row(aln_len=160, qlen=200))
        assert "GENE_0001" in filter_hits(hits, min_pident=0, min_qcov=80)


class TestMultiQueryIndependence:
    """Each query keeps its own best hit; queries don't pool."""

    def test_two_queries_get_independent_best_hits(self):
        rows = "\n".join(
            [
                make_blast_outfmt_row(qseqid="Q1", sseqid="Q1_BEST", bitscore=900.0),
                make_blast_outfmt_row(qseqid="Q1", sseqid="Q1_OTHER", bitscore=200.0),
                make_blast_outfmt_row(qseqid="Q2", sseqid="Q2_BEST", bitscore=500.0),
                make_blast_outfmt_row(qseqid="Q2", sseqid="Q2_OTHER", bitscore=100.0),
            ]
        )
        hits = parse_blast_tabular(rows)
        assert hits["Q1"]["blastp_hit_accession"] == "Q1_BEST"
        assert hits["Q2"]["blastp_hit_accession"] == "Q2_BEST"


class TestEntryFieldsPopulated:
    """All seven output fields present, with the expected types and rounding."""

    def test_all_fields_present_and_typed(self):
        hits = parse_blast_tabular(
            make_blast_outfmt_row(
                qseqid="Q",
                sseqid="WP_001",
                pident=92.347,
                aln_len=180,
                qlen=200,
                evalue=1e-50,
                bitscore=400.0,
                stitle="hemolysin",
            )
        )
        e = hits["Q"]
        assert e["locus_tag"] == "Q"
        assert e["blastp_hit_accession"] == "WP_001"
        assert e["blastp_hit_description"] == "hemolysin"
        # pident rounded to 1 dp
        assert e["blastp_pident"] == 92.3
        # qcov = 180/200 * 100 = 90.0 (rounded to 1 dp)
        assert e["blastp_qcov"] == 90.0
        assert e["blastp_evalue"] == 1e-50
        assert e["blastp_bitscore"] == 400.0


class TestOutfmtPinning:
    """The BLAST_OUTFMT string is the contract between the subprocess command
    and the parser's column indices. If the format ever drifts out of sync
    with `_helpers.BLAST_OUTFMT_COLS`, every column-index lookup silently
    returns the wrong field. Pin the alignment."""

    def test_outfmt_column_order_matches_helpers(self):
        # BLAST_OUTFMT is "6 col1 col2 ... colN"; first token is the format
        # specifier, rest are columns.
        cols_in_outfmt = BLAST_OUTFMT.split()
        assert cols_in_outfmt[0] == "6"
        assert cols_in_outfmt[1:] == BLAST_OUTFMT_COLS

    def test_outfmt_has_15_columns(self):
        # _BLAST_MIN_FIELDS = 15 in the parser; column-count drift below 15
        # silently drops every row.
        assert len(BLAST_OUTFMT.split()) - 1 == 15
        assert len(BLAST_OUTFMT_COLS) == 15
