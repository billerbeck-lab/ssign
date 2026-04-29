#!/usr/bin/env python3
"""Extract proteins in the genomic neighborhood of secretion system components.

Given gene order, SS components, and a proximity window, outputs a FASTA file
containing only proteins within ±window genes of any SS component. This allows
downstream tools (DeepLocPro, DeepSecE, SignalP) to run on a focused subset
instead of the full proteome.

Also includes the SS component proteins themselves (needed for localization
validation in cross_validate_predictions.py).
"""

import argparse
import csv
import logging
from collections import defaultdict

from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def load_gene_order(path):
    """Load gene order TSV → dict of contig → [(position, locus_tag), ...]."""
    contigs = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row.get('contig', row.get('replicon', ''))
            locus = row.get('locus_tag', row.get('protein_id', ''))
            pos = int(row.get('position', row.get('gene_index', 0)))
            if contig and locus:
                contigs[contig].append((pos, locus))
    # Sort by position
    for contig in contigs:
        contigs[contig].sort()
    return contigs


def load_ss_components(path):
    """Load SS components TSV → set of locus_tags."""
    components = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            locus = row.get('locus_tag', row.get('protein_id', ''))
            if locus:
                components.add(locus)
    return components


def get_neighborhood_proteins(gene_order, ss_components, window):
    """Find all proteins within ±window of any SS component on same contig.

    Returns set of locus_tags including the components themselves.
    """
    neighborhood = set()

    for contig, genes in gene_order.items():
        for i, (pos, locus) in enumerate(genes):
            if locus in ss_components:
                # Include the component itself
                neighborhood.add(locus)
                # Include ±window neighbors
                for j in range(max(0, i - window), min(len(genes), i + window + 1)):
                    neighborhood.add(genes[j][1])

    return neighborhood


def main():
    parser = argparse.ArgumentParser(
        description="Extract neighborhood proteins around SS components"
    )
    parser.add_argument("--gene-order", required=True,
                        help="Gene order TSV from extract_gene_order.py")
    parser.add_argument("--ss-components", required=True,
                        help="SS components TSV from validate_macsyfinder_systems.py")
    parser.add_argument("--proteins", required=True,
                        help="Full proteome FASTA")
    parser.add_argument("--window", type=int, default=3,
                        help="Proximity window (default: 3 genes)")
    parser.add_argument("--output", required=True,
                        help="Output FASTA with neighborhood proteins only")
    parser.add_argument("--output-ids", default="",
                        help="Optional: output file listing neighborhood protein IDs")
    args = parser.parse_args()

    # Load data
    gene_order = load_gene_order(args.gene_order)
    ss_components = load_ss_components(args.ss_components)

    if not ss_components:
        logger.warning("No SS components found — outputting empty FASTA")
        with open(args.output, 'w') as f:
            pass
        return

    logger.info(f"Loaded {sum(len(v) for v in gene_order.values())} genes "
                f"across {len(gene_order)} contigs")
    logger.info(f"Found {len(ss_components)} SS component proteins")

    # Get neighborhood
    neighborhood = get_neighborhood_proteins(gene_order, ss_components, args.window)
    logger.info(f"Neighborhood (±{args.window} genes): {len(neighborhood)} proteins")

    # Extract sequences
    n_written = 0
    with open(args.output, 'w') as out:
        for rec in SeqIO.parse(args.proteins, 'fasta'):
            if rec.id in neighborhood:
                SeqIO.write(rec, out, 'fasta')
                n_written += 1

    logger.info(f"Wrote {n_written} neighborhood proteins to {args.output}")

    # Optionally write ID list
    if args.output_ids:
        with open(args.output_ids, 'w') as f:
            for pid in sorted(neighborhood):
                f.write(pid + '\n')

    # Summary
    total_proteins = sum(1 for _ in SeqIO.parse(args.proteins, 'fasta'))
    pct = 100 * n_written / max(total_proteins, 1)
    logger.info(f"Reduction: {total_proteins} → {n_written} proteins ({pct:.1f}% of proteome)")


if __name__ == '__main__':
    main()
