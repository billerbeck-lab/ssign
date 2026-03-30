#!/usr/bin/env python3
"""Extract protein sequences and gene info from GenBank or GFF3 files.

Supports:
  - GenBank/GBFF files (Bakta, NCBI, Prokka output)
  - GFF3 + FASTA pairs

Produces:
  - proteins.faa: FASTA of all CDS translations
  - gene_info.tsv: tab-separated metadata per CDS

Adapted from pipeline/lib/id_mapping.py::extract_accessions_from_gbff()
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def extract_from_genbank(gbff_path: str, sample_id: str):
    """Extract CDS features from a GenBank/GBFF file.

    Yields dicts with keys: locus_tag, protein_id, product, contig,
    start, end, strand, sequence.
    """
    for record in SeqIO.parse(gbff_path, "genbank"):
        contig_id = record.id
        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers

            # Skip pseudogenes
            if "pseudo" in qualifiers or "pseudogene" in qualifiers:
                continue

            # Get translation
            translation = qualifiers.get("translation", [None])[0]
            if not translation:
                continue

            # Get identifiers
            locus_tag = qualifiers.get("locus_tag", [None])[0]
            protein_id = qualifiers.get("protein_id", [None])[0]
            product = qualifiers.get("product", ["hypothetical protein"])[0]
            gene = qualifiers.get("gene", [""])[0]

            # Use locus_tag as primary ID, fallback to protein_id
            primary_id = locus_tag or protein_id
            if not primary_id:
                continue

            # Location
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'

            yield {
                "locus_tag": primary_id,
                "protein_id": protein_id or "",
                "gene": gene,
                "product": product,
                "contig": contig_id,
                "start": start,
                "end": end,
                "strand": strand,
                "sequence": translation,
            }


def extract_from_gff3(gff_path: str, fasta_path: str, sample_id: str):
    """Extract CDS features from GFF3 + FASTA pair.

    Uses BioPython to parse GFF3 and translate CDS from the FASTA.
    """
    from Bio.Seq import Seq

    # Load genome sequences
    genome_seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome_seqs[record.id] = str(record.seq)

    # Parse GFF3
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue

            contig = parts[0]
            start = int(parts[3]) - 1  # GFF is 1-based
            end = int(parts[4])
            strand = parts[6]

            # Parse attributes
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attrs[key] = val

            locus_tag = attrs.get('locus_tag', attrs.get('ID', ''))
            protein_id = attrs.get('protein_id', '')
            product = attrs.get('product', 'hypothetical protein')
            gene = attrs.get('gene', '')

            if not locus_tag:
                continue

            # Get or translate sequence
            translation = attrs.get('translation', '')
            if not translation and contig in genome_seqs:
                nuc_seq = genome_seqs[contig][start:end]
                if strand == '-':
                    nuc_seq = str(Seq(nuc_seq).reverse_complement())
                try:
                    translation = str(Seq(nuc_seq).translate(to_stop=True))
                except Exception:
                    continue

            if not translation:
                continue

            yield {
                "locus_tag": locus_tag,
                "protein_id": protein_id,
                "gene": gene,
                "product": product.replace('%20', ' ').replace('%2C', ','),
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "sequence": translation,
            }


def main():
    parser = argparse.ArgumentParser(description="Extract proteins from genome annotation")
    parser.add_argument("--input", required=True, help="GenBank/GBFF or GFF3 file")
    parser.add_argument("--fasta", help="FASTA file (required for GFF3 input)")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--out-proteins", required=True, help="Output FASTA path")
    parser.add_argument("--out-gene-info", required=True, help="Output gene info TSV path")
    args = parser.parse_args()

    input_path = Path(args.input)
    ext = input_path.suffix.lower()

    # Determine parser
    if ext in ('.gbff', '.gbk', '.gb'):
        entries = list(extract_from_genbank(args.input, args.sample))
    elif ext in ('.gff', '.gff3'):
        if not args.fasta:
            logger.error("GFF3 input requires --fasta argument")
            sys.exit(1)
        entries = list(extract_from_gff3(args.input, args.fasta, args.sample))
    else:
        logger.error(f"Unsupported format: {ext}")
        sys.exit(1)

    # Deduplicate by locus_tag (keep first)
    seen = set()
    unique_entries = []
    for e in entries:
        if e["locus_tag"] not in seen:
            seen.add(e["locus_tag"])
            unique_entries.append(e)

    logger.info(f"Extracted {len(unique_entries)} CDS from {input_path.name} "
                f"({len(entries) - len(unique_entries)} duplicates removed)")

    # Write FASTA
    with open(args.out_proteins, 'w') as fasta_out:
        for e in unique_entries:
            fasta_out.write(f">{e['locus_tag']}\n")
            seq = e["sequence"]
            for i in range(0, len(seq), 80):
                fasta_out.write(seq[i:i+80] + "\n")

    # Write gene info TSV
    fieldnames = ["locus_tag", "protein_id", "gene", "product",
                  "contig", "start", "end", "strand"]
    with open(args.out_gene_info, 'w', newline='') as tsv_out:
        writer = csv.DictWriter(tsv_out, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for e in unique_entries:
            writer.writerow(e)


if __name__ == '__main__':
    main()
