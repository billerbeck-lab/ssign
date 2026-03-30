#!/usr/bin/env python3
"""Run BLASTp in local or remote mode.

Local mode: runs BLAST+ blastp locally against a user-provided database.
Remote mode: submits sequences to NCBI BLASTp web API.

CRITICAL BUG FIX PRESERVED:
  NCBI concatenates multiple hit descriptions with " >".
  We only check the PRIMARY description (before first " >"):
    primary_desc = hit_desc.split(" >")[0]

Adapted from pipeline/scripts/run_blastp.py
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
import time

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

from ssign_lib.fasta_io import read_fasta

# Terms to exclude from BLASTp hits
EXCLUDE_TERMS = [
    "hypothetical protein", "uncharacterized protein",
    "domain of unknown function", "duf",
    "unnamed protein product", "predicted protein",
]


def load_substrate_ids(substrates_path):
    """Load substrate locus_tags from substrate TSV."""
    ids = set()
    with open(substrates_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ids.add(row['locus_tag'])
    return ids


def run_local_blastp(query_fasta, db_path, evalue, exclude_taxid, num_threads=4):
    """Run local BLAST+ blastp and return parsed hits."""
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen"

    cmd = [
        "blastp",
        "-query", query_fasta,
        "-db", db_path,
        "-outfmt", outfmt,
        "-evalue", str(evalue),
        "-max_target_seqs", "10",
        "-num_threads", str(num_threads),
    ]
    if exclude_taxid:
        cmd.extend(["-negative_taxids", str(exclude_taxid)])

    logger.info(f"Running local BLASTp: {' '.join(cmd[:6])}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)

    if result.returncode != 0:
        logger.error(f"BLASTp failed: {result.stderr[:500]}")
        raise RuntimeError(f"BLASTp exit code {result.returncode}")

    return parse_blast_tabular(result.stdout)


def run_remote_blastp(query_fasta, evalue, exclude_taxid):
    """Submit sequences to NCBI BLASTp web API and parse results."""
    from Bio.Blast import NCBIWWW, NCBIXML

    sequences = read_fasta(query_fasta)
    all_hits = {}
    batch_size = 15
    protein_ids = list(sequences.keys())

    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i + batch_size]
        batch_fasta = "\n".join(f">{pid}\n{sequences[pid]}" for pid in batch)

        logger.info(f"Submitting batch {i // batch_size + 1} "
                     f"({len(batch)} proteins) to NCBI BLASTp...")

        try:
            result_handle = NCBIWWW.qblast(
                "blastp", "nr", batch_fasta,
                expect=evalue,
                hitlist_size=10,
                entrez_query=f"NOT txid{exclude_taxid}[ORGN]" if exclude_taxid else "",
            )
            records = NCBIXML.parse(result_handle)

            for record in records:
                query_id = record.query.split()[0]
                query_len = record.query_length

                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        pident = (hsp.identities / hsp.align_length) * 100
                        qcov = ((hsp.query_end - hsp.query_start + 1) / query_len) * 100

                        # CRITICAL: only check PRIMARY description
                        hit_desc = alignment.hit_def
                        primary_desc = hit_desc.split(" >")[0]

                        all_hits[query_id] = {
                            'locus_tag': query_id,
                            'blastp_hit_accession': alignment.accession,
                            'blastp_hit_description': primary_desc[:200],
                            'blastp_pident': round(pident, 1),
                            'blastp_qcov': round(qcov, 1),
                            'blastp_evalue': hsp.expect,
                        }
                        break  # Best HSP only
                    break  # Best alignment only

        except Exception as e:
            logger.warning(f"Batch {i // batch_size + 1} failed: {e}")

        if i + batch_size < len(protein_ids):
            time.sleep(15)  # Rate limiting

    return all_hits


def parse_blast_tabular(output_text):
    """Parse BLAST tabular output (outfmt 6)."""
    hits = {}
    for line in output_text.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 15:
            continue

        query_id = parts[0]
        if query_id in hits:
            continue  # Keep first (best) hit

        pident = float(parts[2])
        qlen = int(parts[13])
        aln_len = int(parts[3])
        qcov = (aln_len / qlen * 100) if qlen > 0 else 0

        # CRITICAL: only check PRIMARY description (before first " >")
        hit_desc = parts[12]
        primary_desc = hit_desc.split(" >")[0]

        hits[query_id] = {
            'locus_tag': query_id,
            'blastp_hit_accession': parts[1],
            'blastp_hit_description': primary_desc[:200],
            'blastp_pident': round(pident, 1),
            'blastp_qcov': round(qcov, 1),
            'blastp_evalue': float(parts[10]),
        }

    return hits


def filter_hits(hits, min_pident, min_qcov):
    """Filter BLASTp hits by identity, coverage, and description."""
    filtered = {}
    for pid, hit in hits.items():
        if hit['blastp_pident'] < min_pident:
            continue
        if hit['blastp_qcov'] < min_qcov:
            continue

        desc_lower = hit['blastp_hit_description'].lower()
        if any(term in desc_lower for term in EXCLUDE_TERMS):
            continue

        filtered[pid] = hit

    return filtered


def main():
    parser = argparse.ArgumentParser(description="Run BLASTp (local or remote)")
    parser.add_argument("--mode", choices=['local', 'remote'], required=True)
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--db", default="")
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min-pident", type=float, default=80)
    parser.add_argument("--min-qcov", type=float, default=80)
    parser.add_argument("--exclude-taxid", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Extract substrate sequences into temp FASTA
    substrate_ids = load_substrate_ids(args.substrates)
    all_seqs = read_fasta(args.proteins)
    sub_seqs = {k: v for k, v in all_seqs.items() if k in substrate_ids}

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        for pid, seq in sub_seqs.items():
            tmp.write(f">{pid}\n{seq}\n")
        tmp_path = tmp.name

    try:
        if args.mode == 'local':
            if not args.db:
                logger.error("Local mode requires --db")
                sys.exit(1)
            hits = run_local_blastp(tmp_path, args.db, args.evalue, args.exclude_taxid)
        else:
            hits = run_remote_blastp(tmp_path, args.evalue, args.exclude_taxid)
    finally:
        os.unlink(tmp_path)

    filtered = filter_hits(hits, args.min_pident, args.min_qcov)
    logger.info(f"{len(filtered)}/{len(hits)} hits pass filters for {args.sample}")

    fieldnames = ['locus_tag', 'blastp_hit_accession', 'blastp_hit_description',
                  'blastp_pident', 'blastp_qcov', 'blastp_evalue']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for hit in filtered.values():
            writer.writerow(hit)


if __name__ == '__main__':
    main()
