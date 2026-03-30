#!/usr/bin/env python3
"""Run Bakta genome annotation and convert output to ssign standard format.

Bakta (https://github.com/oschwengers/bakta) provides rich annotation
compared to Prodigal — includes functional annotations, gene names, and
database cross-references. Requires a database download (~30GB full,
~2GB light).

Usage:
    bakta --db /path/to/db --output outdir --prefix sample contigs.fasta

Output files used:
    {prefix}.faa  — protein sequences
    {prefix}.tsv  — annotation table with locus tags, coordinates, products

Database download:
    bakta_db download --output /path/to/db --type light   # ~2GB
    bakta_db download --output /path/to/db --type full    # ~30GB
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


# Bakta TSV column names (tab-separated, with header)
# Sequence Id | Type | Start | Stop | Strand | Locus Tag | Gene | Product | DbXrefs
_COL_SEQ_ID = "Sequence Id"
_COL_TYPE = "Type"
_COL_START = "Start"
_COL_STOP = "Stop"
_COL_STRAND = "Strand"
_COL_LOCUS_TAG = "Locus Tag"
_COL_GENE = "Gene"
_COL_PRODUCT = "Product"

# Only these feature types contain proteins
_PROTEIN_TYPES = {"CDS", "sORF"}


def run_bakta(contigs_fasta, db_path, sample_id, output_dir, threads=4):
    """Run Bakta on a genome FASTA file.

    Returns:
        tuple: (proteins_faa_path, tsv_path)
    """
    cmd = [
        "bakta",
        "--db", db_path,
        "--output", output_dir,
        "--prefix", sample_id,
        "--threads", str(threads),
        "--force",     # Overwrite if exists
        "--skip-plot", # Skip plot (not needed in pipeline)
        contigs_fasta,
    ]

    logger.info(f"Running Bakta: bakta --db {db_path} --prefix {sample_id} ...")
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=14400,
        )
        if result.returncode != 0:
            logger.error(f"Bakta failed:\n{result.stderr[:1000]}")
            raise RuntimeError(f"Bakta exit code {result.returncode}")
    except subprocess.TimeoutExpired:
        raise RuntimeError("Bakta timed out after 4 hours")

    proteins_faa = os.path.join(output_dir, f"{sample_id}.faa")
    tsv_path = os.path.join(output_dir, f"{sample_id}.tsv")

    if not os.path.exists(proteins_faa):
        raise FileNotFoundError(f"Bakta proteins not found: {proteins_faa}")
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Bakta TSV not found: {tsv_path}")

    return proteins_faa, tsv_path


def parse_bakta_tsv(tsv_path):
    """Parse Bakta TSV annotation table.

    Returns list of dicts matching gene_info.tsv format:
        locus_tag, protein_id, gene, product, contig, start, end, strand
    """
    entries = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            feat_type = row.get(_COL_TYPE, "").strip()
            if feat_type not in _PROTEIN_TYPES:
                continue

            locus_tag = row.get(_COL_LOCUS_TAG, "").strip()
            if not locus_tag:
                continue

            try:
                start = int(row.get(_COL_START, 0)) - 1  # Convert to 0-based
                end = int(row.get(_COL_STOP, 0))
            except (ValueError, TypeError):
                start, end = 0, 0

            entries.append({
                "locus_tag": locus_tag,
                "protein_id": "",
                "gene": row.get(_COL_GENE, "").strip(),
                "product": row.get(_COL_PRODUCT, "hypothetical protein").strip() or "hypothetical protein",
                "contig": row.get(_COL_SEQ_ID, "").strip(),
                "start": start,
                "end": end,
                "strand": row.get(_COL_STRAND, "+").strip(),
            })

    return entries


def write_proteins_fasta(bakta_faa_path, entries, out_fasta):
    """Write cleaned proteins FASTA using Bakta locus tags.

    Bakta's .faa headers use locus tags directly:
    >{locus_tag} {product} [{contig}]
    """
    # Build set of locus tags we want
    wanted = {e["locus_tag"] for e in entries}

    # Parse Bakta FASTA, keeping only CDS/sORF proteins
    seqs = {}
    current_tag = None
    current_seq = []

    with open(bakta_faa_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_tag and current_tag in wanted:
                    seqs[current_tag] = "".join(current_seq).rstrip('*')
                # Bakta header: >{locus_tag} {description}
                current_tag = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_tag and current_tag in wanted:
            seqs[current_tag] = "".join(current_seq).rstrip('*')

    n_written = 0
    with open(out_fasta, 'w') as f:
        for e in entries:
            tag = e["locus_tag"]
            seq = seqs.get(tag, "")
            if not seq:
                continue
            f.write(f">{tag}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
            n_written += 1

    logger.info(f"Wrote {n_written} protein sequences")
    return n_written


def main():
    parser = argparse.ArgumentParser(description="Run Bakta and convert output to ssign format")
    parser.add_argument("--input", required=True, help="Input genome FASTA (contigs)")
    parser.add_argument("--db", required=True, help="Path to Bakta database")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--threads", type=int, default=4, help="CPU threads")
    parser.add_argument("--out-proteins", required=True, help="Output proteins FASTA")
    parser.add_argument("--out-gene-info", required=True, help="Output gene_info.tsv")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        # Run Bakta
        bakta_faa, bakta_tsv = run_bakta(
            args.input, args.db, args.sample, tmpdir, args.threads
        )

        # Parse annotation table
        entries = parse_bakta_tsv(bakta_tsv)
        logger.info(f"Parsed {len(entries)} CDS/sORF features from Bakta TSV")

        # Write cleaned proteins FASTA
        write_proteins_fasta(bakta_faa, entries, args.out_proteins)

    # Write gene_info TSV
    fieldnames = ["locus_tag", "protein_id", "gene", "product",
                  "contig", "start", "end", "strand"]
    with open(args.out_gene_info, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for e in entries:
            writer.writerow(e)

    logger.info(f"Done: {len(entries)} proteins annotated for {args.sample}")


if __name__ == '__main__':
    main()
