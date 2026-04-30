#!/usr/bin/env python3
"""Compute physicochemical properties using BioPython ProtParam.

Computes: molecular weight, isoelectric point, GRAVY (hydrophobicity),
instability index, aromaticity, and charge at pH 7.

Adapted from analyze_substrate_annotations.py::compute_physicochemical_properties()
"""

import argparse
import csv
import logging

from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

import os as _os
import sys as _sys

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)
from ssign_lib.fasta_io import read_fasta
from ssign_lib.substrates import load_substrate_ids


def main():
    parser = argparse.ArgumentParser(description="Compute physicochemical properties")
    parser.add_argument("--substrates", required=True, help="Substrate TSV")
    parser.add_argument("--proteins", required=True, help="Proteins FASTA")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    substrate_ids = load_substrate_ids(args.substrates)

    # Load sequences
    sequences = read_fasta(args.proteins)

    # Compute properties
    fieldnames = ["locus_tag", "mw_da", "isoelectric_point", "gravy", "instability_index", "aromaticity", "charge_ph7"]

    n_ok = 0
    n_fail = 0
    results = []

    for locus_tag in sorted(substrate_ids):
        seq = sequences.get(locus_tag, "")
        if not seq:
            n_fail += 1
            continue

        # Remove non-standard amino acids
        clean_seq = "".join(c for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        if len(clean_seq) < 10:
            n_fail += 1
            continue

        try:
            pa = ProteinAnalysis(clean_seq)
            results.append(
                {
                    "locus_tag": locus_tag,
                    "mw_da": round(pa.molecular_weight(), 2),
                    "isoelectric_point": round(pa.isoelectric_point(), 4),
                    "gravy": round(pa.gravy(), 4),
                    "instability_index": round(pa.instability_index(), 4),
                    "aromaticity": round(pa.aromaticity(), 4),
                    "charge_ph7": round(pa.charge_at_pH(7.0), 4),
                }
            )
            n_ok += 1
        except Exception as e:
            logger.warning(f"ProtParam failed for {locus_tag}: {e}")
            n_fail += 1

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    logger.info(f"Computed properties for {n_ok}/{n_ok + n_fail} proteins")


if __name__ == "__main__":
    main()
