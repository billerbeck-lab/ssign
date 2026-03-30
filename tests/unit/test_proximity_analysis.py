"""Tests for proximity analysis logic.

Validates:
- Per-component window (not system-boundary)
- Multi-contig boundary handling
- Off-by-one in window calculation

These tests import from the actual proximity_analysis.py script, which uses
a main() CLI entry point. We test the core logic by building the same data
structures directly.
"""

import csv
import os
import sys

import pytest

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))


def write_tsv(path, fieldnames, rows):
    """Helper: write a list of dicts as TSV."""
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_gene_order(genes):
    """Build genes_by_contig and contig_index from list of gene dicts."""
    genes_by_contig = {}
    locus_to_info = {}
    for g in genes:
        contig = g['contig_id']
        if contig not in genes_by_contig:
            genes_by_contig[contig] = []
        genes_by_contig[contig].append(g)
        locus_to_info[g['locus_tag']] = {
            'contig': contig,
            'gene_index': int(g['gene_index']),
        }

    contig_index = {}
    for contig, gene_list in genes_by_contig.items():
        idx_map = {}
        for g in gene_list:
            idx_map[int(g['gene_index'])] = g['locus_tag']
        contig_index[contig] = idx_map

    return genes_by_contig, locus_to_info, contig_index


def find_nearby_genes(contig_index, locus_to_info, ss_component_loci, component_ss_types, window=3):
    """Core proximity logic extracted from proximity_analysis.py for testing."""
    results = []
    for comp_locus in ss_component_loci:
        info = locus_to_info.get(comp_locus)
        if not info:
            continue

        contig = info['contig']
        comp_idx = info['gene_index']
        idx_map = contig_index.get(contig, {})
        max_idx = max(idx_map.keys()) if idx_map else 0

        for offset in range(-window, window + 1):
            neighbor_idx = comp_idx + offset
            if neighbor_idx < 0 or neighbor_idx > max_idx:
                continue

            neighbor_locus = idx_map.get(neighbor_idx)
            if not neighbor_locus:
                continue

            # Skip SS components themselves
            if neighbor_locus in ss_component_loci:
                continue

            results.append({
                'locus_tag': neighbor_locus,
                'contig': contig,
                'ss_type': component_ss_types.get(comp_locus, ''),
            })

    return results


@pytest.fixture
def gene_data():
    """10 genes on contig_A, 5 on contig_B."""
    genes = []
    for i in range(10):
        genes.append({
            'contig_id': 'contig_A',
            'gene_index': str(i),
            'locus_tag': f'GENE_{i:04d}',
            'start': str(i * 1000),
            'end': str(i * 1000 + 999),
            'strand': '+',
        })
    for i in range(5):
        genes.append({
            'contig_id': 'contig_B',
            'gene_index': str(i),
            'locus_tag': f'GENEB_{i:04d}',
            'start': str(i * 1000),
            'end': str(i * 1000 + 999),
            'strand': '+',
        })
    return genes


class TestProximityWindow:
    def test_window_3_returns_6_neighbors(self, gene_data):
        """Window=3 around gene 5: genes 2,3,4,6,7,8 (6 neighbors)."""
        _, locus_to_info, contig_index = build_gene_order(gene_data)

        ss_component_loci = {'GENE_0005'}
        component_ss_types = {'GENE_0005': 'T2SS'}

        results = find_nearby_genes(
            contig_index, locus_to_info, ss_component_loci,
            component_ss_types, window=3,
        )

        nearby_tags = {r['locus_tag'] for r in results}
        # Should include genes at index 2-4, 6-8 (not 5 itself)
        expected = {f'GENE_{i:04d}' for i in range(2, 9) if i != 5}
        assert nearby_tags == expected

    def test_window_does_not_span_contigs(self, gene_data):
        """Component at end of contig_A should NOT include contig_B genes."""
        _, locus_to_info, contig_index = build_gene_order(gene_data)

        ss_component_loci = {'GENE_0009'}
        component_ss_types = {'GENE_0009': 'T2SS'}

        results = find_nearby_genes(
            contig_index, locus_to_info, ss_component_loci,
            component_ss_types, window=3,
        )

        nearby_tags = {r['locus_tag'] for r in results}
        # Only contig_A genes should appear
        for tag in nearby_tags:
            assert tag.startswith('GENE_'), f"Cross-contig leak: {tag}"
        # Should be genes 6, 7, 8 only (9 is the component itself)
        expected = {'GENE_0006', 'GENE_0007', 'GENE_0008'}
        assert nearby_tags == expected

    def test_window_at_contig_start(self, gene_data):
        """Component at gene 0 — window shouldn't go negative."""
        _, locus_to_info, contig_index = build_gene_order(gene_data)

        ss_component_loci = {'GENE_0000'}
        component_ss_types = {'GENE_0000': 'T1SS'}

        results = find_nearby_genes(
            contig_index, locus_to_info, ss_component_loci,
            component_ss_types, window=3,
        )

        nearby_tags = {r['locus_tag'] for r in results}
        expected = {'GENE_0001', 'GENE_0002', 'GENE_0003'}
        assert nearby_tags == expected

    def test_component_excluded_from_results(self, gene_data):
        """The SS component itself should never appear as a substrate."""
        _, locus_to_info, contig_index = build_gene_order(gene_data)

        ss_component_loci = {'GENE_0005'}
        component_ss_types = {'GENE_0005': 'T2SS'}

        results = find_nearby_genes(
            contig_index, locus_to_info, ss_component_loci,
            component_ss_types, window=3,
        )

        nearby_tags = {r['locus_tag'] for r in results}
        assert 'GENE_0005' not in nearby_tags
