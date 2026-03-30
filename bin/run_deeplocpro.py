#!/usr/bin/env python3
"""Run DeepLocPro and parse output to standard format.

Two modes:
- local: requires DTU academic license, user provides install path.
- remote: uses BioLib wrapper (pybiolib package) — no license needed.

Output columns: locus_tag, predicted_localization, extracellular_prob,
periplasmic_prob, outer_membrane_prob, cytoplasmic_prob, product
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


# ── Local mode ──

def run_local_deeplocpro(input_fasta, deeplocpro_path, output_dir, organism="gram-"):
    """Run DeepLocPro CLI locally and return path to results file."""
    # Try common entry points
    candidates = ["deeplocpro", "predict.py", "deeploc"]
    cmd_base = None

    if deeplocpro_path and os.path.isdir(deeplocpro_path):
        for entry in candidates:
            candidate = os.path.join(deeplocpro_path, entry)
            if os.path.exists(candidate):
                cmd_base = candidate
                break

    if not cmd_base:
        # Try as a command on PATH
        cmd_base = "deeplocpro"

    cmd = [cmd_base, "-f", input_fasta, "-o", output_dir, "-g", "negative"]

    logger.info(f"Running local DeepLocPro: {' '.join(cmd[:4])}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

    if result.returncode != 0:
        logger.error(f"DeepLocPro failed: {result.stderr[:500]}")
        raise RuntimeError(f"DeepLocPro exit code {result.returncode}")

    return find_output_file(output_dir)


# ── Remote mode (BioLib) ──

def run_remote_deeplocpro(input_fasta, output_dir):
    """Run DeepLocPro via BioLib (no DTU license needed).

    Uses pybiolib package: pip install pybiolib
    BioLib app: KU/DeepLocPro
    """
    try:
        import biolib
    except ImportError:
        logger.error(
            "pybiolib not installed. Install with: pip install pybiolib\n"
            "Or use --mode local with a DTU-licensed DeepLocPro install."
        )
        raise RuntimeError("pybiolib not installed")

    logger.info("Submitting to BioLib DeepLocPro...")
    app = biolib.load("KU/DeepLocPro")

    result = app.cli(
        args=f"--fasta {input_fasta} --organism gram-",
    )

    # Save output files
    result.save_files(output_dir)

    return find_output_file(output_dir)


def find_output_file(output_dir):
    """Find DeepLocPro output CSV/TSV in output directory."""
    for fname in os.listdir(output_dir):
        if fname.endswith('.csv') or fname.endswith('.tsv'):
            return os.path.join(output_dir, fname)

    # Check subdirectories
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            if fname.endswith('.csv') or fname.endswith('.tsv'):
                return os.path.join(root, fname)

    raise FileNotFoundError(f"No output file found in {output_dir}")


def parse_deeplocpro_output(results_path):
    """Parse DeepLocPro output CSV into standardized format."""
    entries = []

    for sep in [',', '\t']:
        try:
            with open(results_path) as f:
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    # Handle various column naming conventions
                    protein_id = (row.get('protein_id') or row.get('Protein_ID')
                                  or row.get('ID') or row.get('id', ''))

                    ext_prob = float(row.get('Extracellular', row.get('extracellular', 0)))
                    peri_prob = float(row.get('Periplasm', row.get('periplasmic',
                                     row.get('Periplasmic', 0))))
                    om_prob = float(row.get('OuterMembrane', row.get('outer_membrane',
                                   row.get('Outer Membrane', 0))))
                    cyto_prob = float(row.get('Cytoplasm', row.get('cytoplasmic',
                                     row.get('Cytoplasmic', 0))))

                    probs = {
                        'Extracellular': ext_prob,
                        'Periplasmic': peri_prob,
                        'Outer Membrane': om_prob,
                        'Cytoplasmic': cyto_prob,
                    }
                    predicted = max(probs, key=probs.get)

                    entries.append({
                        'locus_tag': protein_id,
                        'predicted_localization': predicted,
                        'extracellular_prob': round(ext_prob, 4),
                        'periplasmic_prob': round(peri_prob, 4),
                        'outer_membrane_prob': round(om_prob, 4),
                        'cytoplasmic_prob': round(cyto_prob, 4),
                        'product': row.get('annotation', row.get('product', '')),
                    })

                if entries:
                    return entries
        except Exception:
            continue

    return entries


def main():
    parser = argparse.ArgumentParser(description="Run DeepLocPro (local or remote)")
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", choices=['local', 'remote'], default='remote',
                        help="local: DTU license needed. remote: uses BioLib (free).")
    parser.add_argument("--deeplocpro-path", default="", help="Path to DeepLocPro install (local mode)")
    parser.add_argument("--conf-threshold", type=float, default=0.8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        if args.mode == 'local':
            results_path = run_local_deeplocpro(args.input, args.deeplocpro_path, tmpdir)
        else:
            results_path = run_remote_deeplocpro(args.input, tmpdir)

        entries = parse_deeplocpro_output(results_path)

    logger.info(f"Parsed {len(entries)} proteins from DeepLocPro")

    fieldnames = ['locus_tag', 'predicted_localization', 'extracellular_prob',
                  'periplasmic_prob', 'outer_membrane_prob', 'cytoplasmic_prob', 'product']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for e in entries:
            writer.writerow(e)

    n_ext = sum(1 for e in entries if e['extracellular_prob'] >= args.conf_threshold)
    logger.info(f"{n_ext}/{len(entries)} proteins are extracellular (>= {args.conf_threshold})")


if __name__ == '__main__':
    main()
