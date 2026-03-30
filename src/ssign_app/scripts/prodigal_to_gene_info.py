#!/usr/bin/env python3
"""Convert Prodigal output to ssign's standard gene_info format.

Prodigal FASTA headers contain coordinates:
  >contig_1 # start # end # strand # ID=1_1;partial=...

This script parses those headers and produces:
  - Cleaned proteins FASTA with simple locus_tag headers
  - gene_info.tsv matching the format from extract_proteins.py
"""

import argparse
import csv
import re


def main():
    parser = argparse.ArgumentParser(description="Convert Prodigal output to gene_info")
    parser.add_argument("--proteins", required=True, help="Prodigal proteins FASTA")
    parser.add_argument("--gff", required=True, help="Prodigal GFF output")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--out-proteins", required=True, help="Output cleaned FASTA")
    parser.add_argument("--out-gene-info", required=True, help="Output gene_info.tsv")
    args = parser.parse_args()

    entries = []
    sequences = {}

    # Parse Prodigal FASTA — headers have coordinates embedded
    with open(args.proteins) as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                # Parse Prodigal header: >contig_N # start # end # strand # attrs
                parts = line[1:].split(' # ')
                raw_id = parts[0].strip()
                start = int(parts[1]) if len(parts) > 1 else 0
                end = int(parts[2]) if len(parts) > 2 else 0
                strand_code = int(parts[3]) if len(parts) > 3 else 1
                strand = '+' if strand_code == 1 else '-'

                # Extract contig from the ID (Prodigal appends _N to contig name)
                # e.g., "contig_1_5" means contig "contig_1", gene 5
                match = re.match(r'^(.+)_(\d+)$', raw_id)
                if match:
                    contig = match.group(1)
                    gene_num = match.group(2)
                else:
                    contig = raw_id
                    gene_num = "0"

                # Create a proper locus_tag
                locus_tag = f"{args.sample}_{gene_num.zfill(5)}"
                current_id = locus_tag

                entries.append({
                    "locus_tag": locus_tag,
                    "protein_id": "",
                    "gene": "",
                    "product": "hypothetical protein",
                    "contig": contig,
                    "start": start - 1,  # Convert to 0-based
                    "end": end,
                    "strand": strand,
                })
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)

    # Write cleaned FASTA
    with open(args.out_proteins, 'w') as out:
        for e in entries:
            seq = sequences.get(e["locus_tag"], "")
            if not seq:
                continue
            # Remove trailing * from Prodigal translations
            seq = seq.rstrip('*')
            out.write(f">{e['locus_tag']}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

    # Write gene_info TSV
    fieldnames = ["locus_tag", "protein_id", "gene", "product",
                  "contig", "start", "end", "strand"]
    with open(args.out_gene_info, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for e in entries:
            writer.writerow(e)


if __name__ == '__main__':
    main()
