"""Tests for extract_neighborhood.py.

Three pure functions plus a CLI driver. The window-arithmetic invariant
(``range(max(0, i-w), min(len, i+w+1))``) is shared with
``proximity_analysis.py`` — the canonical regression for that fix lives
in ``test_proximity_analysis.py``; here we exercise the standalone
``get_neighborhood_proteins`` helper plus the loader's column-alias
tolerance for legacy upstream column names.
"""

import os
import sys

import pytest
from hypothesis import given
from hypothesis import strategies as st


from _helpers import (  # noqa: E402
    GENE_ORDER_FIELDS,
    SS_COMPONENT_FIELDS,
    make_ss_component_row,
    run_script_main,
    write_tsv,
)
from extract_neighborhood import (  # noqa: E402
    get_neighborhood_proteins,
    load_gene_order,
    load_ss_components,
)
from extract_neighborhood import (
    main as neighborhood_main,
)

# ---------------------------------------------------------------------------
# get_neighborhood_proteins — window arithmetic
# ---------------------------------------------------------------------------


def _gene_order(contig, n_genes):
    """Build [(position, locus_tag), ...] for a contig with n_genes."""
    return [(i, f"{contig}_{i:04d}") for i in range(n_genes)]


class TestNeighborhoodWindow:
    def test_window_3_around_index_5_yields_seven_genes(self):
        genes = {"contig_A": _gene_order("contig_A", 10)}
        components = {"contig_A_0005"}
        result = get_neighborhood_proteins(genes, components, window=3)
        # ±3 around index 5 = indices 2..8 inclusive (7 genes)
        assert result == {f"contig_A_{i:04d}" for i in range(2, 9)}

    def test_window_zero_returns_only_the_component(self):
        genes = {"contig_A": _gene_order("contig_A", 10)}
        components = {"contig_A_0005"}
        result = get_neighborhood_proteins(genes, components, window=0)
        assert result == {"contig_A_0005"}

    def test_window_clips_at_contig_start(self):
        genes = {"contig_A": _gene_order("contig_A", 10)}
        components = {"contig_A_0000"}
        result = get_neighborhood_proteins(genes, components, window=3)
        # ±3 around index 0 = indices 0..3 (no negatives)
        assert result == {f"contig_A_{i:04d}" for i in range(0, 4)}

    def test_window_clips_at_contig_end(self):
        genes = {"contig_A": _gene_order("contig_A", 10)}
        components = {"contig_A_0009"}
        result = get_neighborhood_proteins(genes, components, window=3)
        # ±3 around index 9 (last) = indices 6..9
        assert result == {f"contig_A_{i:04d}" for i in range(6, 10)}

    def test_window_does_not_cross_contigs(self):
        genes = {
            "contig_A": _gene_order("contig_A", 5),
            "contig_B": _gene_order("contig_B", 5),
        }
        # Component is at the end of contig_A → must not pull contig_B genes
        components = {"contig_A_0004"}
        result = get_neighborhood_proteins(genes, components, window=3)
        assert all(locus.startswith("contig_A") for locus in result)

    def test_multiple_components_take_union_of_windows(self):
        genes = {"contig_A": _gene_order("contig_A", 20)}
        components = {"contig_A_0005", "contig_A_0015"}
        result = get_neighborhood_proteins(genes, components, window=2)
        # Two non-overlapping windows: {3-7} ∪ {13-17}
        expected = {f"contig_A_{i:04d}" for i in range(3, 8)} | {f"contig_A_{i:04d}" for i in range(13, 18)}
        assert result == expected

    def test_overlapping_components_collapse_correctly(self):
        # Two components, three genes apart, window=3 → windows overlap, no double-count
        genes = {"contig_A": _gene_order("contig_A", 10)}
        components = {"contig_A_0003", "contig_A_0006"}
        result = get_neighborhood_proteins(genes, components, window=3)
        # Union of {0-6} and {3-9} = {0-9}
        assert result == {f"contig_A_{i:04d}" for i in range(0, 10)}

    def test_empty_components_yields_empty(self):
        genes = {"contig_A": _gene_order("contig_A", 10)}
        assert get_neighborhood_proteins(genes, set(), window=3) == set()

    def test_component_not_on_any_contig_is_ignored(self):
        # Robustness: a stale SS_component pointing to an unknown locus
        # must not crash; it just contributes nothing.
        genes = {"contig_A": _gene_order("contig_A", 10)}
        components = {"orphan_locus"}
        assert get_neighborhood_proteins(genes, components, window=3) == set()


# ---------------------------------------------------------------------------
# Property: every result for a given (component, window) lives within the
# expected index range on the component's contig
# ---------------------------------------------------------------------------


@given(
    n_genes=st.integers(min_value=1, max_value=50),
    component_idx=st.integers(min_value=0, max_value=49),
    window=st.integers(min_value=0, max_value=10),
)
def test_window_bounded_by_clipped_range(n_genes, component_idx, window):
    component_idx = component_idx % n_genes  # keep in range
    genes = {"contig_A": _gene_order("contig_A", n_genes)}
    components = {f"contig_A_{component_idx:04d}"}
    result = get_neighborhood_proteins(genes, components, window=window)

    expected_indices = set(
        range(
            max(0, component_idx - window),
            min(n_genes, component_idx + window + 1),
        )
    )
    expected = {f"contig_A_{i:04d}" for i in expected_indices}
    assert result == expected


# ---------------------------------------------------------------------------
# load_gene_order — column-alias tolerance
# ---------------------------------------------------------------------------


class TestLoadGeneOrder:
    def test_canonical_columns(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "gene_order.tsv"),
            GENE_ORDER_FIELDS,
            [
                {"contig": "c_A", "gene_index": "0", "locus_tag": "g0", "start": "0", "end": "100", "strand": "+"},
                {"contig": "c_A", "gene_index": "1", "locus_tag": "g1", "start": "200", "end": "300", "strand": "+"},
            ],
        )
        result = load_gene_order(path)
        assert result["c_A"] == [(0, "g0"), (1, "g1")]

    def test_legacy_column_aliases_replicon_position(self, tmp_dir):
        # Older upstream emitted `replicon` instead of `contig` and `position`
        # instead of `gene_index`. Loader tolerates both.
        path = write_tsv(
            os.path.join(tmp_dir, "gene_order.tsv"),
            ["replicon", "position", "locus_tag"],
            [
                {"replicon": "c_A", "position": "1", "locus_tag": "g1"},
                {"replicon": "c_A", "position": "0", "locus_tag": "g0"},
            ],
        )
        result = load_gene_order(path)
        # Sorted by position regardless of input order
        assert result["c_A"] == [(0, "g0"), (1, "g1")]

    def test_protein_id_alias_for_locus_tag(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "gene_order.tsv"),
            ["contig", "gene_index", "protein_id"],
            [{"contig": "c_A", "gene_index": "0", "protein_id": "p0"}],
        )
        result = load_gene_order(path)
        assert result["c_A"] == [(0, "p0")]

    def test_rows_with_blank_contig_or_locus_skipped(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "gene_order.tsv"),
            GENE_ORDER_FIELDS,
            [
                {"contig": "", "gene_index": "0", "locus_tag": "orphan", "start": "0", "end": "100", "strand": "+"},
                {"contig": "c_A", "gene_index": "0", "locus_tag": "", "start": "0", "end": "100", "strand": "+"},
                {"contig": "c_A", "gene_index": "1", "locus_tag": "valid", "start": "200", "end": "300", "strand": "+"},
            ],
        )
        result = load_gene_order(path)
        assert result == {"c_A": [(1, "valid")]}


class TestLoadSsComponents:
    def test_canonical_locus_tag_column(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "ss_components.tsv"),
            SS_COMPONENT_FIELDS,
            [
                make_ss_component_row("g_T2SS_1", "T2SS"),
                make_ss_component_row("g_T2SS_2", "T2SS"),
            ],
        )
        assert load_ss_components(path) == {"g_T2SS_1", "g_T2SS_2"}

    def test_protein_id_alias(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "ss_components.tsv"),
            ["protein_id", "ss_type"],
            [{"protein_id": "p1", "ss_type": "T2SS"}],
        )
        assert load_ss_components(path) == {"p1"}


# ---------------------------------------------------------------------------
# main() integration — empty SS components and FASTA filtering
# ---------------------------------------------------------------------------


@pytest.fixture
def _proteins_fasta(tmp_dir):
    path = os.path.join(tmp_dir, "proteins.faa")
    with open(path, "w") as f:
        for i in range(10):
            f.write(f">contig_A_{i:04d}\nMKT{i}\n")
    return path


@pytest.fixture
def _ten_gene_order_tsv(tmp_dir):
    """gene_order.tsv with 10 genes on contig_A, locus_tag contig_A_0000..0009."""
    return write_tsv(
        os.path.join(tmp_dir, "gene_order.tsv"),
        GENE_ORDER_FIELDS,
        [
            {
                "contig": "contig_A",
                "gene_index": str(i),
                "locus_tag": f"contig_A_{i:04d}",
                "start": str(i * 1000),
                "end": str(i * 1000 + 999),
                "strand": "+",
            }
            for i in range(10)
        ],
    )


def _run_main(monkeypatch, tmp_dir, gene_order_path, ss_components_path, proteins_path, *, window=3, output_ids=False):
    out_fasta = os.path.join(tmp_dir, "neighborhood.faa")
    argv = [
        "extract_neighborhood",
        "--gene-order",
        gene_order_path,
        "--ss-components",
        ss_components_path,
        "--proteins",
        proteins_path,
        "--window",
        str(window),
        "--output",
        out_fasta,
    ]
    out_ids_path = None
    if output_ids:
        out_ids_path = os.path.join(tmp_dir, "ids.txt")
        argv += ["--output-ids", out_ids_path]
    run_script_main(monkeypatch, neighborhood_main, argv)
    return out_fasta, out_ids_path


def test_main_writes_neighborhood_proteins(
    monkeypatch,
    tmp_dir,
    _proteins_fasta,
    _ten_gene_order_tsv,
):
    ss = write_tsv(
        os.path.join(tmp_dir, "ss_components.tsv"),
        SS_COMPONENT_FIELDS,
        [make_ss_component_row("contig_A_0005", "T2SS")],
    )
    out_fasta, _ = _run_main(
        monkeypatch,
        tmp_dir,
        _ten_gene_order_tsv,
        ss,
        _proteins_fasta,
        window=2,
    )
    with open(out_fasta) as f:
        text = f.read()
    # ±2 around 0005 = 0003..0007
    for i in range(3, 8):
        assert f"contig_A_{i:04d}" in text
    # 0001 (out of window) must not appear
    assert "contig_A_0001" not in text


def test_main_empty_components_writes_empty_fasta(monkeypatch, tmp_dir, _proteins_fasta):
    gene_order = write_tsv(
        os.path.join(tmp_dir, "gene_order.tsv"),
        GENE_ORDER_FIELDS,
        [{"contig": "contig_A", "gene_index": "0", "locus_tag": "g0", "start": "0", "end": "100", "strand": "+"}],
    )
    ss = write_tsv(
        os.path.join(tmp_dir, "ss_components.tsv"),
        SS_COMPONENT_FIELDS,
        [],
    )
    out_fasta, _ = _run_main(monkeypatch, tmp_dir, gene_order, ss, _proteins_fasta)
    assert os.path.exists(out_fasta)
    assert os.path.getsize(out_fasta) == 0


def test_main_writes_id_list_when_requested(
    monkeypatch,
    tmp_dir,
    _proteins_fasta,
    _ten_gene_order_tsv,
):
    ss = write_tsv(
        os.path.join(tmp_dir, "ss_components.tsv"),
        SS_COMPONENT_FIELDS,
        [make_ss_component_row("contig_A_0005", "T2SS")],
    )
    _, ids_path = _run_main(
        monkeypatch,
        tmp_dir,
        _ten_gene_order_tsv,
        ss,
        _proteins_fasta,
        window=2,
        output_ids=True,
    )
    with open(ids_path) as f:
        ids = {line.strip() for line in f if line.strip()}
    assert ids == {f"contig_A_{i:04d}" for i in range(3, 8)}
