#!/usr/bin/env python3
"""Detect input genome file format for ssign pipeline.

Prints one of: genbank, gff3, fasta_contigs
"""

import sys
from pathlib import Path


def detect_format(filepath: str) -> str:
    """Detect genome file format by extension and content inspection."""
    p = Path(filepath)
    ext = p.suffix.lower()

    # Extension-based detection
    if ext in ('.gbff', '.gbk', '.gb'):
        return 'genbank'
    if ext in ('.gff', '.gff3', '.gtf'):
        return 'gff3'
    if ext == '.faa':
        return 'protein_fasta'
    if ext in ('.fasta', '.fna', '.fa'):
        # Could be annotated proteins or raw contigs — check content
        return _inspect_fasta(filepath)

    # Fallback: inspect first lines
    return _inspect_content(filepath)


def _inspect_fasta(filepath: str) -> str:
    """Check if FASTA contains protein sequences or nucleotide contigs."""
    with open(filepath) as f:
        seq_chars = set()
        lines_read = 0
        for line in f:
            if line.startswith('>'):
                continue
            seq_chars.update(line.strip().upper())
            lines_read += 1
            if lines_read > 50:
                break

    # If mostly ATGCN, it's nucleotide contigs
    nuc_chars = {'A', 'T', 'G', 'C', 'N'}
    if seq_chars and (seq_chars - nuc_chars) == set():
        return 'fasta_contigs'

    # If it has amino acid characters, it might be proteins
    # but for pipeline purposes we treat bare FASTA as contigs
    return 'fasta_contigs'


def _inspect_content(filepath: str) -> str:
    """Inspect first few lines to guess format."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('LOCUS'):
                return 'genbank'
            if line.startswith('##gff-version'):
                return 'gff3'
            if line.startswith('>'):
                return 'fasta_contigs'
            break

    raise ValueError(f"Cannot determine format of {filepath}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: detect_input_format.py <file>", file=sys.stderr)
        sys.exit(1)
    print(detect_format(sys.argv[1]))
