#!/usr/bin/env python3
"""Run SignalP 6.0 and parse output to standard format.

Two modes:
- local: requires DTU academic license download.
- remote: uses BioLib wrapper (pybiolib package) — no license needed.

SignalP predicts signal peptides: Sec/SPI, Sec/SPII, Tat/SPI, Tat/SPII, Sec/SPIII.

Output columns: locus_tag, signalp_prediction, signalp_probability, signalp_cs_position
"""

import argparse
import csv
import logging
import os
import subprocess
import tempfile

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


# ── Local mode ──

def run_local_signalp(input_fasta, signalp_path, output_dir):
    """Run SignalP 6.0 CLI locally."""
    signalp_bin = os.path.join(signalp_path, "signalp6") if signalp_path else "signalp6"

    cmd = [
        signalp_bin,
        "--fastafile", input_fasta,
        "--output_dir", output_dir,
        "--organism", "gram-",
        "--format", "txt",
    ]

    logger.info(f"Running local SignalP: {' '.join(cmd[:4])}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        logger.error(f"SignalP failed: {result.stderr[:500]}")
        raise RuntimeError(f"SignalP exit code {result.returncode}")

    return find_output_file(output_dir)


# ── Remote mode (BioLib) ──

def run_remote_signalp(input_fasta, output_dir):
    """Run SignalP via BioLib (no DTU license needed).

    Uses pybiolib package: pip install pybiolib
    BioLib app: DTU/SignalP-6
    """
    try:
        import biolib
    except ImportError:
        logger.error(
            "pybiolib not installed. Install with: pip install pybiolib\n"
            "Or use --mode local with a DTU-licensed SignalP install."
        )
        raise RuntimeError("pybiolib not installed")

    logger.info("Submitting to BioLib SignalP 6.0...")
    app = biolib.load("DTU/SignalP-6")

    result = app.cli(
        args=f"--fastafile {input_fasta} --organism gram- --format txt",
    )

    result.save_files(output_dir)
    return find_output_file(output_dir)


def find_output_file(output_dir):
    """Find SignalP output in directory."""
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            if fname.endswith('.txt') or fname.endswith('.signalp5') or fname.endswith('.gff3'):
                return os.path.join(root, fname)

    raise FileNotFoundError(f"No SignalP output in {output_dir}")


def parse_signalp_output(results_path):
    """Parse SignalP output file.

    SignalP 6.0 output (tab-separated, with header starting with #):
    # ID  Prediction  SP(Sec/SPI)  TAT(Tat/SPI)  LIPO(Sec/SPII)  ...  CS Position
    """
    entries = []

    with open(results_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            protein_id = parts[0]
            prediction = parts[1]

            # Find CS position (column that contains "CS pos:")
            cs_position = ""
            for part in parts:
                if "CS pos:" in part:
                    cs_position = part

            # Find max probability among signal peptide types
            max_prob = 0.0
            for part in parts[2:]:
                try:
                    prob = float(part)
                    if prob > max_prob:
                        max_prob = prob
                except ValueError:
                    continue

            entries.append({
                'locus_tag': protein_id,
                'signalp_prediction': prediction,
                'signalp_probability': round(max_prob, 4),
                'signalp_cs_position': cs_position,
            })

    return entries


def main():
    parser = argparse.ArgumentParser(description="Run SignalP 6.0 (local or remote)")
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", choices=['local', 'remote'], default='remote',
                        help="local: DTU license needed. remote: uses BioLib (free).")
    parser.add_argument("--signalp-path", default="", help="Path to SignalP install (local mode)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        if args.mode == 'local':
            results_path = run_local_signalp(args.input, args.signalp_path, tmpdir)
        else:
            results_path = run_remote_signalp(args.input, tmpdir)

        entries = parse_signalp_output(results_path)

    logger.info(f"Parsed {len(entries)} SignalP predictions for {args.sample}")
    n_sp = sum(1 for e in entries if e['signalp_prediction'] != 'OTHER')
    logger.info(f"{n_sp}/{len(entries)} proteins have signal peptides")

    fieldnames = ['locus_tag', 'signalp_prediction', 'signalp_probability',
                  'signalp_cs_position']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for e in entries:
            writer.writerow(e)


if __name__ == '__main__':
    main()
