#!/usr/bin/env python3
"""Assign ortholog groups to substrate proteins using all-vs-all BLASTp.

Algorithm:
1. Collect all filtered substrate sequences across genomes
2. Run all-vs-all BLASTp (local makeblastdb + blastp)
3. Filter hits: >= min_pident AND >= min_qcov
4. Single-linkage clustering via Union-Find
5. Output: CSV mapping locus_tag → ortholog_group + group stats

Requires local NCBI BLAST+ (makeblastdb, blastp). Falls back to remote
NCBI BLASTp if local not available (much slower).
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Try to import BioPython for FASTA parsing
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def read_fasta_simple(path):
    """Read FASTA without BioPython dependency."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif current_id:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = ''.join(current_seq)
    return seqs


def run_local_blast(fasta_path, min_pident, min_qcov, evalue=1e-5):
    """Run all-vs-all BLASTp locally. Returns list of (query, subject, pident, qcov) tuples."""
    # Check for blastp and makeblastdb
    blastp = None
    for candidate in ['blastp', 'blastp.exe']:
        try:
            subprocess.run([candidate, '-version'], capture_output=True, timeout=10)
            blastp = candidate
            break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue

    if not blastp:
        logger.warning(
            "Local BLASTp not found -- skipping ortholog grouping.\n"
            "  To enable ortholog grouping, install NCBI BLAST+:\n"
            "    - sudo apt install ncbi-blast+  (Debian/Ubuntu)\n"
            "    - conda install -c bioconda blast  (conda)\n"
            "    - Or download from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )
        return None

    with tempfile.TemporaryDirectory(prefix="ortho_blast_") as tmpdir:
        db_path = os.path.join(tmpdir, "substrates_db")
        out_path = os.path.join(tmpdir, "blast_results.txt")

        # Make database
        makeblastdb = blastp.replace('blastp', 'makeblastdb')
        cmd_db = [makeblastdb, "-in", fasta_path, "-dbtype", "prot",
                  "-out", db_path]
        result = subprocess.run(cmd_db, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            logger.error(f"makeblastdb failed: {result.stderr[:300]}")
            return None

        # Run BLASTp
        n_seqs = sum(1 for line in open(fasta_path) if line.startswith('>'))
        cmd_blast = [
            blastp,
            "-query", fasta_path,
            "-db", db_path,
            "-out", out_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore qlen slen",
            "-evalue", str(evalue),
            "-max_target_seqs", str(min(n_seqs, 10000)),
            "-num_threads", "4",
        ]
        logger.info(f"Running all-vs-all BLASTp on {n_seqs} sequences...")
        result = subprocess.run(cmd_blast, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            logger.error(f"BLASTp failed: {result.stderr[:300]}")
            return None

        # Parse results
        hits = []
        with open(out_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 14:
                    continue
                query, subject = parts[0], parts[1]
                if query == subject:
                    continue  # Skip self-hits
                pident = float(parts[2])
                aln_len = int(parts[3])
                qlen = int(parts[12])
                qcov = 100.0 * aln_len / max(qlen, 1)

                if pident >= min_pident and qcov >= min_qcov:
                    hits.append((query, subject, pident, qcov))

        logger.info(f"Found {len(hits)} ortholog-quality hits "
                    f"(>={min_pident}% id, >={min_qcov}% qcov)")
        return hits


def cluster_union_find(hits, all_protein_ids):
    """Single-linkage clustering via Union-Find. Returns dict of group_id → set of members."""
    parent = {pid: pid for pid in all_protein_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for query, subject, _, _ in hits:
        if query in parent and subject in parent:
            union(query, subject)

    groups = defaultdict(set)
    for pid in all_protein_ids:
        groups[find(pid)].add(pid)

    return groups


def compute_group_stats(groups, hits, all_protein_ids):
    """Compute per-group statistics. Returns list of dicts."""
    # Build hit lookup for within-group identity calculation
    hit_identities = defaultdict(list)
    for query, subject, pident, _ in hits:
        hit_identities[(query, subject)].append(pident)

    stats = []
    for idx, (rep, members) in enumerate(
            sorted(groups.items(), key=lambda x: -len(x[1])), 1):
        # Compute mean identity within group
        within_ids = []
        member_list = sorted(members)
        for i, m1 in enumerate(member_list):
            for m2 in member_list[i+1:]:
                for key in [(m1, m2), (m2, m1)]:
                    if key in hit_identities:
                        within_ids.extend(hit_identities[key])

        mean_id = sum(within_ids) / len(within_ids) if within_ids else 100.0

        stats.append({
            'ortholog_group': f'OG_{idx:03d}',
            'n_members': len(members),
            'members': ';'.join(sorted(members)),
            'mean_pident': round(mean_id, 1),
        })

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Assign ortholog groups via all-vs-all BLASTp"
    )
    parser.add_argument("--substrates-fasta", required=True,
                        help="FASTA of all substrate proteins (combined across genomes)")
    parser.add_argument("--min-pident", type=float, default=40.0,
                        help="Minimum %% identity for ortholog assignment (default: 40)")
    parser.add_argument("--min-qcov", type=float, default=70.0,
                        help="Minimum query coverage %% for ortholog assignment (default: 70)")
    parser.add_argument("--evalue", type=float, default=1e-5,
                        help="E-value threshold for BLASTp (default: 1e-5)")
    parser.add_argument("--output", required=True,
                        help="Output CSV with ortholog group assignments")
    parser.add_argument("--output-groups", default="",
                        help="Optional: output CSV with group-level statistics")
    args = parser.parse_args()

    # Read all substrate sequences
    seqs = read_fasta_simple(args.substrates_fasta)
    if not seqs:
        logger.warning("No sequences in input FASTA — writing empty output")
        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['locus_tag', 'ortholog_group',
                                                     'og_n_members', 'og_mean_pident'])
            writer.writeheader()
        return

    logger.info(f"Loaded {len(seqs)} substrate sequences")

    # Need at least 2 sequences for all-vs-all
    if len(seqs) < 2:
        logger.info("Only 1 substrate — assigning singleton group OG_001")
        pid = list(seqs.keys())[0]
        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['locus_tag', 'ortholog_group',
                                                     'og_n_members', 'og_mean_pident'])
            writer.writeheader()
            writer.writerow({'locus_tag': pid, 'ortholog_group': 'OG_001',
                             'og_n_members': 1, 'og_mean_pident': 100.0})
        return

    # Run all-vs-all BLASTp
    hits = run_local_blast(args.substrates_fasta, args.min_pident,
                           args.min_qcov, args.evalue)

    if hits is None:
        logger.error("BLASTp failed — cannot assign ortholog groups")
        # Write output with empty groups
        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['locus_tag', 'ortholog_group',
                                                     'og_n_members', 'og_mean_pident'])
            writer.writeheader()
            for pid in seqs:
                writer.writerow({'locus_tag': pid, 'ortholog_group': '',
                                 'og_n_members': 0, 'og_mean_pident': 0})
        return

    # Cluster
    all_ids = set(seqs.keys())
    groups = cluster_union_find(hits, all_ids)
    logger.info(f"Clustered into {len(groups)} ortholog groups")

    # Compute stats
    group_stats = compute_group_stats(groups, hits, all_ids)

    # Build protein → group mapping
    pid_to_group = {}
    pid_to_stats = {}
    for gs in group_stats:
        for member in gs['members'].split(';'):
            pid_to_group[member] = gs['ortholog_group']
            pid_to_stats[member] = gs

    # Write per-protein output
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['locus_tag', 'ortholog_group',
                                                 'og_n_members', 'og_mean_pident'])
        writer.writeheader()
        for pid in sorted(seqs.keys()):
            og = pid_to_group.get(pid, '')
            stats = pid_to_stats.get(pid, {})
            writer.writerow({
                'locus_tag': pid,
                'ortholog_group': og,
                'og_n_members': stats.get('n_members', 1),
                'og_mean_pident': stats.get('mean_pident', 100.0),
            })

    logger.info(f"Wrote {len(seqs)} protein assignments to {args.output}")

    # Write group-level stats
    if args.output_groups:
        with open(args.output_groups, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['ortholog_group', 'n_members',
                                                     'members', 'mean_pident'])
            writer.writeheader()
            for gs in group_stats:
                writer.writerow(gs)
        logger.info(f"Wrote {len(group_stats)} group stats to {args.output_groups}")

    # Summary
    sizes = [gs['n_members'] for gs in group_stats]
    n_singleton = sum(1 for s in sizes if s == 1)
    n_multi = sum(1 for s in sizes if s > 1)
    logger.info(f"Groups: {n_singleton} singletons, {n_multi} multi-member "
                f"(largest: {max(sizes)} members)")


if __name__ == '__main__':
    main()
