"""Tests for extract_gene_order.py.

The output gene_order.tsv is consumed by proximity_analysis to compute
±N-gene windows around each SS component. Two invariants are
load-bearing for the rest of the pipeline:

1. Per-contig 0-based gene_index ascending by `start` — without this,
   proximity arithmetic returns the wrong neighbours.
2. Contigs are output in alphabetical order — keeps the file
   deterministic across runs (helps the golden-output E2E test in 4.i).
"""

import os
import sys
import tempfile

from hypothesis import given, settings
from hypothesis import strategies as st


from _helpers import (  # noqa: E402
    GENE_INFO_FIELDS,
    GENE_ORDER_FIELDS,
    read_tsv_rows,
    run_script_main,
    write_tsv,
)
from extract_gene_order import main as extract_gene_order_main  # noqa: E402


def _gene_row(locus_tag, contig, start, end=None, strand="+"):
    return {
        "locus_tag": locus_tag,
        "protein_id": locus_tag,
        "gene": "",
        "product": "hypothetical protein",
        "contig": contig,
        "start": str(start),
        "end": str(end if end is not None else start + 999),
        "strand": strand,
    }


def _run(monkeypatch, tmp_dir, gene_info_rows):
    gene_info = write_tsv(
        os.path.join(tmp_dir, "gene_info.tsv"),
        GENE_INFO_FIELDS,
        gene_info_rows,
    )
    out = os.path.join(tmp_dir, "gene_order.tsv")
    run_script_main(
        monkeypatch,
        extract_gene_order_main,
        [
            "extract_gene_order",
            "--gene-info",
            gene_info,
            "--output",
            out,
        ],
    )
    return read_tsv_rows(out)


# ---------------------------------------------------------------------------
# Core invariants
# ---------------------------------------------------------------------------


def test_per_contig_gene_index_starts_at_zero(monkeypatch, tmp_dir):
    rows = _run(
        monkeypatch,
        tmp_dir,
        [
            _gene_row("A_0", "contig_A", 0),
            _gene_row("A_1", "contig_A", 1000),
            _gene_row("B_0", "contig_B", 500),
        ],
    )
    by_contig = {}
    for r in rows:
        by_contig.setdefault(r["contig"], []).append(r)
    assert by_contig["contig_A"][0]["gene_index"] == "0"
    assert by_contig["contig_B"][0]["gene_index"] == "0"


def test_unsorted_input_produces_start_sorted_output(monkeypatch, tmp_dir):
    rows = _run(
        monkeypatch,
        tmp_dir,
        [
            _gene_row("A_2", "contig_A", 2000),
            _gene_row("A_0", "contig_A", 0),
            _gene_row("A_1", "contig_A", 1000),
        ],
    )
    starts = [int(r["start"]) for r in rows]
    assert starts == sorted(starts)


def test_gene_index_assigned_in_start_sorted_order(monkeypatch, tmp_dir):
    rows = _run(
        monkeypatch,
        tmp_dir,
        [
            _gene_row("LATE", "contig_A", 2000),
            _gene_row("EARLY", "contig_A", 0),
            _gene_row("MID", "contig_A", 1000),
        ],
    )
    by_index = {int(r["gene_index"]): r["locus_tag"] for r in rows}
    assert by_index == {0: "EARLY", 1: "MID", 2: "LATE"}


def test_contigs_output_alphabetically(monkeypatch, tmp_dir):
    rows = _run(
        monkeypatch,
        tmp_dir,
        [
            _gene_row("Z_0", "contig_Z", 0),
            _gene_row("A_0", "contig_A", 0),
            _gene_row("M_0", "contig_M", 0),
        ],
    )
    contig_order = [r["contig"] for r in rows]
    assert contig_order == ["contig_A", "contig_M", "contig_Z"]


def test_multi_contig_gene_index_resets_per_contig(monkeypatch, tmp_dir):
    rows = _run(
        monkeypatch,
        tmp_dir,
        [
            _gene_row("A_0", "contig_A", 0),
            _gene_row("A_1", "contig_A", 1000),
            _gene_row("A_2", "contig_A", 2000),
            _gene_row("B_0", "contig_B", 500),
            _gene_row("B_1", "contig_B", 1500),
        ],
    )
    a_indexes = [int(r["gene_index"]) for r in rows if r["contig"] == "contig_A"]
    b_indexes = [int(r["gene_index"]) for r in rows if r["contig"] == "contig_B"]
    assert a_indexes == [0, 1, 2]
    assert b_indexes == [0, 1]


def test_output_schema_matches_canonical_fields(monkeypatch, tmp_dir):
    rows = _run(monkeypatch, tmp_dir, [_gene_row("A_0", "contig_A", 0)])
    assert set(rows[0].keys()) == set(GENE_ORDER_FIELDS)


def test_single_gene_per_contig(monkeypatch, tmp_dir):
    rows = _run(
        monkeypatch,
        tmp_dir,
        [
            _gene_row("A_0", "contig_A", 1234),
            _gene_row("B_0", "contig_B", 5678),
        ],
    )
    assert len(rows) == 2
    assert all(r["gene_index"] == "0" for r in rows)


def test_empty_input_writes_header_only(monkeypatch, tmp_dir):
    rows = _run(monkeypatch, tmp_dir, [])
    assert rows == []
    with open(os.path.join(tmp_dir, "gene_order.tsv")) as f:
        assert f.readline().strip().split("\t") == GENE_ORDER_FIELDS


# ---------------------------------------------------------------------------
# Property: gene_index is always the rank of `start` within its contig
# ---------------------------------------------------------------------------


@settings(max_examples=25, deadline=None)
@given(
    starts=st.lists(
        st.integers(min_value=0, max_value=10_000_000),
        min_size=1,
        max_size=20,
        unique=True,
    )
)
def test_gene_index_equals_start_rank(starts):
    """For any set of distinct start positions on one contig, the emitted
    gene_index for each gene must equal that gene's rank (ascending) by start.

    Uses an inline TemporaryDirectory rather than a pytest fixture because
    function-scoped fixtures aren't reset between hypothesis-generated inputs.
    """
    with tempfile.TemporaryDirectory() as td:
        gene_info = write_tsv(
            os.path.join(td, "gene_info.tsv"),
            GENE_INFO_FIELDS,
            [_gene_row(f"GENE_{i}", "contig_A", s) for i, s in enumerate(starts)],
        )
        out = os.path.join(td, "gene_order.tsv")
        saved_argv = sys.argv
        try:
            sys.argv = [
                "extract_gene_order",
                "--gene-info",
                gene_info,
                "--output",
                out,
            ]
            extract_gene_order_main()
        finally:
            sys.argv = saved_argv
        rows = read_tsv_rows(out)

    by_locus = {r["locus_tag"]: int(r["gene_index"]) for r in rows}
    expected = {f"GENE_{i}": rank for rank, (i, _) in enumerate(sorted(enumerate(starts), key=lambda pair: pair[1]))}
    assert by_locus == expected
