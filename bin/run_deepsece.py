#!/usr/bin/env python3
"""Run DeepSecE and parse output to standard format.

DeepSecE predicts which secretion system type a protein is secreted by.
MIT licensed — can be pip-installed.

Column mapping from pipeline/scripts/parse_deepsece.py.

Output columns: locus_tag, dse_ss_type, dse_max_prob, plus per-type probabilities.
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

# Column name mapping from raw DeepSecE output
_COLUMN_MAP = {
    "protein_id": "locus_tag",
    "deepsece_prediction": "deepsece_prediction",
    "deepsece_ss_type": "dse_ss_type",
    "max_prob": "dse_max_prob",
    "nonsec_prob": "dse_nonsec_prob",
    "T1_prob": "dse_T1_prob",
    "T2_prob": "dse_T2_prob",
    "T3_prob": "dse_T3_prob",
    "T4_prob": "dse_T4_prob",
    "T6_prob": "dse_T6_prob",
}


def run_deepsece(input_fasta, output_dir):
    """Run DeepSecE prediction."""
    # Try deepsece CLI first
    cmd = ["deepsece", "predict", "-i", input_fasta, "-o", output_dir]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode == 0:
            return find_output(output_dir)
    except FileNotFoundError:
        pass

    # Fallback: try python -m deepsece
    cmd = [sys.executable, "-m", "deepsece", "predict", "-i", input_fasta, "-o", output_dir]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        logger.error(f"DeepSecE failed: {result.stderr}")
        raise RuntimeError(f"DeepSecE exit code {result.returncode}")

    return find_output(output_dir)


def find_output(output_dir):
    """Find DeepSecE output CSV."""
    for fname in os.listdir(output_dir):
        if fname.endswith('.csv') or fname.endswith('.tsv'):
            return os.path.join(output_dir, fname)
    raise FileNotFoundError(f"No output file in {output_dir}")


def parse_deepsece_output(results_path):
    """Parse DeepSecE output into standardized format."""
    entries = []

    for sep in [',', '\t']:
        try:
            with open(results_path) as f:
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    entry = {}
                    for raw_col, std_col in _COLUMN_MAP.items():
                        entry[std_col] = row.get(raw_col, '')
                    entries.append(entry)
                if entries:
                    return entries
        except Exception:
            continue

    return entries


def main():
    parser = argparse.ArgumentParser(description="Run DeepSecE")
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        results_path = run_deepsece(args.input, tmpdir)
        entries = parse_deepsece_output(results_path)

    logger.info(f"Parsed {len(entries)} DeepSecE predictions for {args.sample}")

    fieldnames = list(_COLUMN_MAP.values())
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for e in entries:
            writer.writerow(e)


if __name__ == '__main__':
    main()
