#!/usr/bin/env python3
"""Produce contig-sorted gene order from gene_info.tsv.

The gene order file is used by proximity_analysis.py to find proteins
within +/- N genes of secretion system components on the same contig.

Output columns:
  contig, gene_index, locus_tag, start, end, strand
"""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description="Extract gene order from gene_info TSV")
    parser.add_argument("--gene-info", required=True, help="Input gene_info.tsv")
    parser.add_argument("--output", required=True, help="Output gene_order.tsv")
    args = parser.parse_args()

    # Read gene info
    genes_by_contig = {}
    with open(args.gene_info) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row["contig"]
            if contig not in genes_by_contig:
                genes_by_contig[contig] = []
            genes_by_contig[contig].append(row)

    # Sort genes by start position within each contig
    for contig in genes_by_contig:
        genes_by_contig[contig].sort(key=lambda g: int(g["start"]))

    # Write gene order
    fieldnames = ["contig", "gene_index", "locus_tag", "start", "end", "strand"]
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for contig in sorted(genes_by_contig.keys()):
            for idx, gene in enumerate(genes_by_contig[contig]):
                writer.writerow({
                    "contig": contig,
                    "gene_index": idx,
                    "locus_tag": gene["locus_tag"],
                    "start": gene["start"],
                    "end": gene["end"],
                    "strand": gene["strand"],
                })


if __name__ == '__main__':
    main()
