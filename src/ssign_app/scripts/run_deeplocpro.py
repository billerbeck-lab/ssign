#!/usr/bin/env python3
"""Run DeepLocPro and parse output to standard format.

Two modes:
- local: requires DTU academic license, user provides install path.
- remote: submits directly to DTU web server (free, no license needed).

Output columns: locus_tag, predicted_localization, extracellular_prob,
periplasmic_prob, outer_membrane_prob, cytoplasmic_prob, product
"""

import argparse
import csv
import logging
import os
import re
import subprocess
import tempfile
import time

import requests

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


# ── Local mode ──

def run_local_deeplocpro(input_fasta, deeplocpro_path, output_dir, organism="gram-"):
    """Run DeepLocPro CLI locally and return path to results file."""
    # Handle 0-sequence input gracefully: write empty output and return
    with open(input_fasta) as _f:
        n_seqs = sum(1 for line in _f if line.startswith('>'))
    if n_seqs == 0:
        logger.info("0 sequences in input FASTA — writing empty output")
        empty_out = os.path.join(output_dir, "deeplocpro_results.csv")
        with open(empty_out, 'w') as f:
            f.write("locus_tag\tpredicted_localization\textracellular_prob\t"
                    "periplasmic_prob\touter_membrane_prob\tcytoplasmic_prob\t"
                    "cytoplasmic_membrane_prob\tmax_localization\tmax_probability\n")
        return empty_out

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
    # FRAGILE: subprocess call requires DeepLocPro binary on PATH or at deeplocpro_path
    # If this breaks: install DeepLocPro locally or switch to --mode remote
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"DeepLocPro binary not found: {e}\n"
            f"  Common causes:\n"
            f"    - DeepLocPro is not installed or not on PATH\n"
            f"    - Wrong --deeplocpro-path specified\n"
            f"  How to fix:\n"
            f"    - Install DeepLocPro (requires DTU academic license)\n"
            f"    - Or use --mode remote (free, no license needed)"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"DeepLocPro timed out after 2 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input sequences\n"
            f"    - Or use --mode remote to offload computation to DTU servers"
        ) from e

    if result.returncode != 0:
        logger.error(f"DeepLocPro failed: {result.stderr[:500]}")
        raise RuntimeError(f"DeepLocPro exit code {result.returncode}")

    return find_output_file(output_dir)


# ── Remote mode (DTU web server) ──

DTU_SUBMIT_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
DTU_RESULTS_BASE = "https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp"
DTU_MAX_POLL = 360  # 30 minutes at 5s intervals (enough for large batches)
DTU_POLL_INTERVAL = 5


DTU_BATCH_SIZE = 500  # DTU accepts max 500 sequences per submission


def _split_fasta_bytes(fasta_content, batch_size):
    """Split FASTA content into batches of batch_size sequences."""
    lines = fasta_content.decode().split('\n')
    batches = []
    current_batch = []
    current_count = 0

    for line in lines:
        if line.startswith('>'):
            if current_count >= batch_size:
                batches.append('\n'.join(current_batch).encode())
                current_batch = []
                current_count = 0
            current_count += 1
        current_batch.append(line)

    if current_batch:
        batches.append('\n'.join(current_batch).encode())

    return batches


def _submit_and_poll_dtu(fasta_bytes, batch_num, total_batches, max_retries=3):
    """Submit one batch to DTU and poll for results, with retry on failure."""
    for attempt in range(1, max_retries + 1):
        try:
            return _submit_and_poll_dtu_once(fasta_bytes, batch_num, total_batches)
        except RuntimeError as e:
            if attempt < max_retries:
                wait = 30 * attempt  # 30s, 60s, 90s
                logger.warning(
                    f"Batch {batch_num}/{total_batches} attempt {attempt}/{max_retries} "
                    f"failed: {e}. Retrying in {wait}s..."
                )
                time.sleep(wait)
            else:
                raise


def _submit_and_poll_dtu_once(fasta_bytes, batch_num, total_batches):
    """Submit one batch to DTU and poll for results. Returns results dict."""
    files = {
        "uploadfile": ("input.fasta", fasta_bytes, "text/plain"),
    }
    data = {
        "configfile": "/var/www/services/services/DeepLocPro-1.0/webface.cf",
        "fasta": "",
        "group": "negative",
        "format": "long",
    }

    # FRAGILE: DTU web server submission can fail due to network issues or server maintenance
    # If this breaks: check https://services.healthtech.dtu.dk status, or use --mode local
    try:
        resp = requests.post(DTU_SUBMIT_URL, data=data, files=files, timeout=60)
    except requests.ConnectionError as e:
        raise RuntimeError(
            f"Cannot connect to DTU DeepLocPro server: {e}\n"
            f"  Common causes:\n"
            f"    - DTU server is down for maintenance\n"
            f"    - Network/firewall blocking outbound HTTPS\n"
            f"  How to fix:\n"
            f"    - Check https://services.healthtech.dtu.dk manually\n"
            f"    - Or use --mode local with a DTU academic license"
        ) from e
    except requests.Timeout as e:
        raise RuntimeError(
            f"DTU DeepLocPro server timed out during submission: {e}\n"
            f"  Common causes:\n"
            f"    - Server overloaded or slow\n"
            f"  How to fix:\n"
            f"    - Retry later, or use --mode local"
        ) from e
    if resp.status_code != 200:
        raise RuntimeError(f"DTU server returned HTTP {resp.status_code}")

    # FRAGILE: job ID regex parsing depends on DTU response format
    # If this breaks: DTU may have changed their web interface HTML/redirect format
    job_match = re.search(r"jobid=([A-F0-9]+)", resp.url)
    if not job_match:
        job_match = re.search(r"jobid=([A-F0-9]+)", resp.text)
    if not job_match:
        raise RuntimeError(
            "Could not parse job ID from DTU response.\n"
            "  Common causes:\n"
            "    - DTU changed their web interface or redirect format\n"
            "    - The response HTML no longer contains 'jobid=<HEX>'\n"
            "  How to fix:\n"
            "    - Check DTU website manually and report to ssign maintainers\n"
            "    - Or use --mode local with a DTU academic license"
        )

    job_id = job_match.group(1)
    logger.info(f"Batch {batch_num}/{total_batches}: DTU job {job_id}")

    for poll_num in range(DTU_MAX_POLL):
        time.sleep(DTU_POLL_INTERVAL)
        try:
            ajax_resp = requests.get(
                f"{DTU_SUBMIT_URL}?ajax=1&jobid={job_id}", timeout=15
            )
            status_data = ajax_resp.json()
            status = status_data.get("status", "unknown")
            runtime = status_data.get("runtime", 0)

            if status == "finished":
                logger.info(f"Batch {batch_num}/{total_batches}: completed in {runtime}s")
                break
            elif status in ("failed", "error"):
                raise RuntimeError(
                    f"DTU job failed after {runtime}s (batch {batch_num})"
                )
            else:
                if poll_num % 12 == 0 and poll_num > 0:
                    logger.info(f"Batch {batch_num}/{total_batches}: {status} ({runtime}s)")
        except requests.RequestException as e:
            logger.warning(f"Poll error: {e}")
    else:
        raise RuntimeError(f"DTU job {job_id} timed out")

    # Fetch results JSON
    results_url = f"{DTU_RESULTS_BASE}/{job_id}/results.json"
    results_resp = requests.get(results_url, timeout=30)
    if results_resp.status_code != 200:
        raise RuntimeError(f"Could not fetch results: HTTP {results_resp.status_code}")

    return results_resp.json()


def run_remote_deeplocpro(input_fasta, output_dir):
    """Submit to DTU DeepLocPro web server directly.

    CRITICAL: FASTA must be sent as file upload (uploadfile field), NOT as
    textarea text (fasta field). Sequences are batched in groups of 500.
    """
    logger.info("Submitting to DTU DeepLocPro web server...")

    with open(input_fasta, 'rb') as f:
        fasta_content = f.read()

    n_seqs = sum(1 for line in fasta_content.decode().split('\n') if line.startswith('>'))

    # Handle 0-sequence input gracefully: write empty output and return
    if n_seqs == 0:
        logger.info("0 sequences in input FASTA — writing empty output")
        empty_out = os.path.join(output_dir, "deeplocpro_results.csv")
        with open(empty_out, 'w') as f:
            f.write("locus_tag\tpredicted_localization\textracellular_prob\t"
                    "periplasmic_prob\touter_membrane_prob\tcytoplasmic_prob\t"
                    "cytoplasmic_membrane_prob\tmax_localization\tmax_probability\n")
        return empty_out

    n_batches = (n_seqs + DTU_BATCH_SIZE - 1) // DTU_BATCH_SIZE
    logger.info(f"Submitting {n_seqs} sequences in {n_batches} batch(es)")

    # Split into batches of 500
    batches = _split_fasta_bytes(fasta_content, DTU_BATCH_SIZE)
    n_batches = len(batches)

    all_sequences = {}
    localizations = []

    for batch_num, batch_bytes in enumerate(batches, 1):
        batch_n = sum(1 for line in batch_bytes.decode().split('\n') if line.startswith('>'))
        logger.info(f"Batch {batch_num}/{n_batches}: {batch_n} sequences")

        try:
            results = _submit_and_poll_dtu(batch_bytes, batch_num, n_batches)
            batch_seqs = results.get("sequences", {})
            all_sequences.update(batch_seqs)
            if not localizations:
                localizations = results.get("Localization", [])
            logger.info(f"Batch {batch_num}/{n_batches}: got {len(batch_seqs)} predictions")
        except Exception as e:
            logger.warning(f"Batch {batch_num}/{n_batches} failed: {e}")
            # Continue with remaining batches

    logger.info(f"Total: {len(all_sequences)} predictions from {n_batches} batches")

    if not all_sequences:
        raise RuntimeError("All DTU batches failed — no predictions returned")

    # Write combined CSV from all batch results
    csv_out = os.path.join(output_dir, "deeplocpro_results.csv")
    with open(csv_out, 'w') as f:
        # Header
        header = ["Protein_ID", "Prediction"] + localizations
        f.write(",".join(header) + "\n")
        for name, seq_data in all_sequences.items():
            probs = seq_data.get("Probability", [])
            row = [name, seq_data.get("Prediction", "")] + [str(p) for p in probs]
            f.write(",".join(row) + "\n")

    return csv_out


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
                    # DTU web server uses 'ACC', local uses 'protein_id' or 'Protein_ID'
                    protein_id = (row.get('ACC') or row.get('protein_id')
                                  or row.get('Protein_ID') or row.get('ID')
                                  or row.get('id') or row.get('Name', ''))

                    if not protein_id:
                        continue

                    ext_prob = float(row.get('Extracellular', row.get('extracellular', 0)))
                    peri_prob = float(row.get('Periplasmic', row.get('Periplasm',
                                     row.get('periplasmic', 0))))
                    om_prob = float(row.get('Outer Membrane', row.get('OuterMembrane',
                                   row.get('outer_membrane', 0))))
                    cyto_prob = float(row.get('Cytoplasmic', row.get('Cytoplasm',
                                     row.get('cytoplasmic', 0))))

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

    try:
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
    except RuntimeError:
        raise  # Already has a user-friendly message
    except Exception as e:
        alt_mode = "remote" if args.mode == "local" else "local"
        raise RuntimeError(
            f"DeepLocPro pipeline failed ({args.mode} mode): {e}\n"
            f"  How to fix:\n"
            f"    - Check the error details above\n"
            f"    - Try --mode {alt_mode} as an alternative"
        ) from e


if __name__ == '__main__':
    main()
