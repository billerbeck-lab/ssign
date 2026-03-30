#!/usr/bin/env python3
"""Deduplicate protein sequences before sending to remote tools.

Groups proteins by identical sequence, sends only unique representatives,
then maps annotations back to all copies. This saves API calls and time.

Usage:
    # Deduplicate
    unique_fasta, seq_groups = deduplicate_fasta(input_fasta)

    # Run tool on unique_fasta...

    # Map results back
    full_results = expand_results(tool_results, seq_groups)
"""

import csv
import hashlib
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

import os as _os, sys as _sys
_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)
from ssign_lib.fasta_io import read_fasta


def deduplicate_dict(seqs):
    """Deduplicate a dict of {protein_id: sequence}.

    Returns:
        unique_seqs: dict of {representative_id: sequence} (one per unique seq)
        seq_groups: dict of {representative_id: [all_ids_with_same_seq]}
    """
    hash_to_ids = {}
    hash_to_rep = {}
    hash_to_seq = {}

    for pid, seq in seqs.items():
        seq_hash = hashlib.md5(seq.encode()).hexdigest()
        if seq_hash not in hash_to_ids:
            hash_to_ids[seq_hash] = []
            hash_to_rep[seq_hash] = pid
            hash_to_seq[seq_hash] = seq
        hash_to_ids[seq_hash].append(pid)

    unique_seqs = {hash_to_rep[h]: hash_to_seq[h] for h in hash_to_ids}
    seq_groups = {hash_to_rep[h]: hash_to_ids[h] for h in hash_to_ids}

    n_total = len(seqs)
    n_unique = len(unique_seqs)
    n_dups = n_total - n_unique
    if n_dups > 0:
        logger.info(f"Deduplicated: {n_total} -> {n_unique} unique sequences ({n_dups} duplicates removed)")

    return unique_seqs, seq_groups


def expand_results_dict(results, seq_groups, id_key='locus_tag'):
    """Expand results dict from unique representatives to all proteins.

    Args:
        results: dict of {representative_id: {id_key: ..., ...}}
        seq_groups: dict from deduplicate_dict
        id_key: field name for the protein ID

    Returns:
        expanded dict with entries for all duplicate members
    """
    expanded = {}
    for rep_id, entry in results.items():
        if rep_id in seq_groups:
            for member_id in seq_groups[rep_id]:
                new_entry = dict(entry)
                new_entry[id_key] = member_id
                expanded[member_id] = new_entry
        else:
            expanded[rep_id] = entry
    return expanded


def deduplicate_fasta(input_fasta, output_fasta):
    """Read a FASTA, write only unique sequences, return grouping.

    Returns:
        seq_groups: dict mapping representative_id -> [all_ids_with_same_seq]
                    The representative is the first ID encountered.
    """
    sequences = read_fasta(input_fasta)

    # Group by sequence hash
    hash_to_ids = {}  # seq_hash -> [list of protein IDs]
    hash_to_rep = {}  # seq_hash -> representative ID (first seen)
    hash_to_seq = {}  # seq_hash -> actual sequence

    for pid, seq in sequences.items():
        seq_hash = hashlib.md5(seq.encode()).hexdigest()
        if seq_hash not in hash_to_ids:
            hash_to_ids[seq_hash] = []
            hash_to_rep[seq_hash] = pid
            hash_to_seq[seq_hash] = seq
        hash_to_ids[seq_hash].append(pid)

    # Build groups: representative -> all IDs
    seq_groups = {}
    for seq_hash, ids in hash_to_ids.items():
        rep = hash_to_rep[seq_hash]
        seq_groups[rep] = ids

    # Write unique sequences
    n_unique = len(seq_groups)
    n_total = len(sequences)
    n_dups = n_total - n_unique

    with open(output_fasta, 'w') as f:
        for seq_hash in hash_to_ids:
            rep = hash_to_rep[seq_hash]
            seq = hash_to_seq[seq_hash]
            f.write(f">{rep}\n{seq}\n")

    if n_dups > 0:
        logger.info(f"Deduplicated: {n_total} -> {n_unique} unique sequences ({n_dups} duplicates)")
    else:
        logger.info(f"No duplicates found in {n_total} sequences")

    return seq_groups


def expand_results_tsv(input_tsv, output_tsv, seq_groups, id_column='locus_tag'):
    """Expand a TSV of results from unique representatives to all proteins.

    For each row where the ID matches a representative, create copies
    for all other members of that group.
    """
    if not os.path.exists(input_tsv):
        return

    with open(input_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        rows = list(reader)

    expanded = []
    for row in rows:
        rep_id = row.get(id_column, '')
        if rep_id in seq_groups:
            # Add a row for each member of the group
            for member_id in seq_groups[rep_id]:
                new_row = dict(row)
                new_row[id_column] = member_id
                expanded.append(new_row)
        else:
            expanded.append(row)

    n_added = len(expanded) - len(rows)
    if n_added > 0:
        logger.info(f"Expanded {len(rows)} -> {len(expanded)} rows ({n_added} from dedup)")

    with open(output_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in expanded:
            writer.writerow(row)


def expand_results_csv(input_csv, output_csv, seq_groups, id_column='locus_tag'):
    """Same as expand_results_tsv but for comma-separated files."""
    if not os.path.exists(input_csv):
        return

    with open(input_csv) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    expanded = []
    for row in rows:
        rep_id = row.get(id_column, '')
        if rep_id in seq_groups:
            for member_id in seq_groups[rep_id]:
                new_row = dict(row)
                new_row[id_column] = member_id
                expanded.append(new_row)
        else:
            expanded.append(row)

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in expanded:
            writer.writerow(row)
