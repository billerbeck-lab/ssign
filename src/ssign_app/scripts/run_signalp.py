#!/usr/bin/env python3
"""Run SignalP 6.0 and parse output to standard format.

Two modes:
- local: requires DTU academic license download.
- remote: submits directly to DTU web server (free, no license needed).

SignalP predicts signal peptides: Sec/SPI, Sec/SPII, Tat/SPI, Tat/SPII, Sec/SPIII.

Output columns: locus_tag, signalp_prediction, signalp_probability, signalp_cs_position
"""

import argparse
import csv
import logging
import os
import re
import subprocess
import tempfile
import time

import requests as http_requests

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

DTU_SUBMIT_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
DTU_RESULTS_BASE = "https://services.healthtech.dtu.dk/services/SignalP-6.0/tmp"
DTU_MAX_POLL = 360  # 30 minutes at 5s intervals (enough for large submissions)
DTU_POLL_INTERVAL = 5


# ── Local mode ──

def run_local_signalp(input_fasta, signalp_path, output_dir):
    """Run SignalP 6.0 CLI locally."""
    # Handle 0-sequence input gracefully: write empty output and return
    with open(input_fasta) as _f:
        n_seqs = sum(1 for line in _f if line.startswith('>'))
    if n_seqs == 0:
        logger.info("0 sequences in input FASTA — writing empty output")
        empty_out = os.path.join(output_dir, "signalp_results.tsv")
        with open(empty_out, 'w') as f:
            f.write("locus_tag\tsignalp_prediction\tsignalp_probability\tsignalp_cs_position\n")
        return empty_out

    signalp_bin = os.path.join(signalp_path, "signalp6") if signalp_path else "signalp6"

    cmd = [
        signalp_bin,
        "--fastafile", input_fasta,
        "--output_dir", output_dir,
        # SignalP v6 collapsed v5's gram-/gram+/arch into a single
        # "other" group (the protein language model no longer needs
        # organism information); only `eukarya`/`euk`/`other` are valid.
        "--organism", "other",
        "--format", "txt",
        "--mode", "fast",
        # Cap threads — upstream default of 8 oversubscribes on small
        # machines (and on multi-genome runs where each PipelineRunner
        # would spawn its own 8-thread pool).
        "--torch_num_threads", str(min(4, os.cpu_count() or 1)),
    ]

    logger.info(f"Running local SignalP: {' '.join(cmd[:4])}...")
    # FRAGILE: subprocess call requires signalp6 binary on PATH or at signalp_path
    # If this breaks: install SignalP 6.0 locally or switch to --mode remote
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"SignalP binary (signalp6) not found: {e}\n"
            f"  Common causes:\n"
            f"    - SignalP 6.0 is not installed or not on PATH\n"
            f"    - Wrong --signalp-path specified\n"
            f"  How to fix:\n"
            f"    - Install SignalP 6.0 (requires DTU academic license)\n"
            f"    - Or use --mode remote (free, no license needed)"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"SignalP timed out after 2 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input sequences\n"
            f"    - Or use --mode remote to offload computation to DTU servers"
        ) from e

    if result.returncode != 0:
        logger.error(f"SignalP failed: {result.stderr[:500]}")
        raise RuntimeError(f"SignalP exit code {result.returncode}")

    return find_output_file(output_dir)


# ── Remote mode (DTU web server) ──

def run_remote_signalp(input_fasta, output_dir, max_retries=3):
    """Submit to DTU SignalP 6.0 web server with retry on failure."""
    for attempt in range(1, max_retries + 1):
        try:
            return _run_remote_signalp_once(input_fasta, output_dir)
        except RuntimeError as e:
            if attempt < max_retries:
                wait = 30 * attempt
                logger.warning(
                    f"DTU SignalP attempt {attempt}/{max_retries} failed: {e}. "
                    f"Retrying in {wait}s..."
                )
                time.sleep(wait)
            else:
                raise


def _run_remote_signalp_once(input_fasta, output_dir):
    """Submit to DTU SignalP 6.0 web server directly.

    Same approach as DeepLocPro: FASTA sent via file upload field.
    """
    logger.info("Submitting to DTU SignalP 6.0 web server...")

    with open(input_fasta, 'rb') as f:
        fasta_content = f.read()

    n_seqs = sum(1 for line in fasta_content.decode().split('\n') if line.startswith('>'))

    # Handle 0-sequence input gracefully: write empty output and return
    if n_seqs == 0:
        logger.info("0 sequences in input FASTA — writing empty output")
        empty_out = os.path.join(output_dir, "signalp_results.tsv")
        with open(empty_out, 'w') as f:
            f.write("locus_tag\tsignalp_prediction\tsignalp_probability\tsignalp_cs_position\n")
        return empty_out

    logger.info(f"Submitting {n_seqs} sequences to DTU SignalP (max 5000)")

    if n_seqs > 5000:
        raise RuntimeError(
            f"DTU SignalP accepts max 5000 sequences, got {n_seqs}. "
            "Use --mode local for larger datasets."
        )

    # Submit via file upload
    files = {
        "uploadfile": ("input.fasta", fasta_content, "text/plain"),
    }
    data = {
        "configfile": "/var/www/services/services/SignalP-6.0/webface.cf",
        "fasta": "",  # textarea left empty — use file upload instead
        "organism": "Other",  # "Other" = non-eukaryote (gram-neg/pos/archaea)
        "format": "short",
        "mode": "fast",
    }

    # FRAGILE: DTU web server submission can fail due to network issues or server maintenance
    # If this breaks: check https://services.healthtech.dtu.dk status, or use --mode local
    try:
        resp = http_requests.post(DTU_SUBMIT_URL, data=data, files=files, timeout=60)
    except http_requests.ConnectionError as e:
        raise RuntimeError(
            f"Cannot connect to DTU SignalP server: {e}\n"
            f"  Common causes:\n"
            f"    - DTU server is down for maintenance\n"
            f"    - Network/firewall blocking outbound HTTPS\n"
            f"  How to fix:\n"
            f"    - Check https://services.healthtech.dtu.dk manually\n"
            f"    - Or use --mode local with a DTU academic license"
        ) from e
    except http_requests.Timeout as e:
        raise RuntimeError(
            f"DTU SignalP server timed out during submission: {e}\n"
            f"  Common causes:\n"
            f"    - Server overloaded or slow\n"
            f"  How to fix:\n"
            f"    - Retry later, or use --mode local"
        ) from e

    if resp.status_code != 200:
        raise RuntimeError(f"DTU server returned HTTP {resp.status_code}")

    # FRAGILE: job ID regex parsing depends on DTU response redirect format
    # If this breaks: DTU may have changed their web interface HTML/redirect format
    job_match = re.search(r"jobid=([A-F0-9]+)", resp.url)
    if not job_match:
        job_match = re.search(r"jobid=([A-F0-9]+)", resp.text)
    if not job_match:
        raise RuntimeError(
            "Could not parse job ID from DTU SignalP response.\n"
            "  Common causes:\n"
            "    - DTU changed their web interface or redirect format\n"
            "    - The response HTML no longer contains 'jobid=<HEX>'\n"
            "  How to fix:\n"
            "    - Check DTU website manually and report to ssign maintainers\n"
            "    - Or use --mode local with a DTU academic license"
        )

    job_id = job_match.group(1)
    logger.info(f"DTU SignalP job submitted: {job_id}")

    # Poll for completion
    for poll_num in range(DTU_MAX_POLL):
        time.sleep(DTU_POLL_INTERVAL)
        try:
            ajax_resp = http_requests.get(
                f"{DTU_SUBMIT_URL}?ajax=1&jobid={job_id}", timeout=15
            )
            status_data = ajax_resp.json()
            status = status_data.get("status", "unknown")
            runtime = status_data.get("runtime", 0)

            if status == "finished":
                logger.info(f"DTU SignalP job completed in {runtime}s")
                break
            elif status in ("failed", "error"):
                raise RuntimeError(
                    f"DTU SignalP job failed after {runtime}s. "
                    "Try again later or use local mode."
                )
            else:
                if poll_num % 6 == 0:
                    logger.info(f"DTU SignalP job {status} (runtime={runtime}s)")
        except http_requests.RequestException as e:
            logger.warning(f"Poll error: {e}")
    else:
        raise RuntimeError(f"DTU SignalP job {job_id} timed out")

    # Try to find output files in the job directory
    # SignalP output is typically a TSV with predictions
    job_dir_url = f"{DTU_RESULTS_BASE}/{job_id}/"
    dir_resp = http_requests.get(job_dir_url, timeout=15)

    if dir_resp.status_code == 200:
        # Parse directory listing for output files
        file_links = re.findall(r'href="([^"]+\.(?:txt|tsv|json|signalp))"', dir_resp.text)
        for fname in file_links:
            file_url = f"{DTU_RESULTS_BASE}/{job_id}/{fname}"
            file_resp = http_requests.get(file_url, timeout=30)
            if file_resp.status_code == 200:
                local_path = os.path.join(output_dir, fname)
                with open(local_path, 'w') as f:
                    f.write(file_resp.text)
                logger.info(f"Downloaded: {fname} ({len(file_resp.text)} bytes)")

    # Also try the prediction_results.txt standard name
    for candidate in ["prediction_results.txt", "output.txt", f"{job_id}_summary.signalp5"]:
        candidate_url = f"{DTU_RESULTS_BASE}/{job_id}/{candidate}"
        resp = http_requests.get(candidate_url, timeout=15)
        if resp.status_code == 200 and len(resp.text) > 10:
            local_path = os.path.join(output_dir, candidate)
            with open(local_path, 'w') as f:
                f.write(resp.text)
            logger.info(f"Downloaded: {candidate}")

    return find_output_file(output_dir)


def find_output_file(output_dir):
    """Find SignalP prediction output in directory.

    Prefers prediction_results.txt over plot files.
    """
    # Priority 1: prediction_results.txt (DTU web server output)
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            if fname == "prediction_results.txt":
                return os.path.join(root, fname)

    # Priority 2: any file with "prediction" or "summary" in name
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            if ("prediction" in fname or "summary" in fname) and fname.endswith('.txt'):
                fpath = os.path.join(root, fname)
                if os.path.getsize(fpath) > 10:
                    return fpath

    # Priority 3: .signalp5 files (local install output)
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            if fname.endswith('.signalp5'):
                return os.path.join(root, fname)

    raise FileNotFoundError(f"No SignalP output in {output_dir}")


_SIGNALP6_COLS = (
    "ID", "Prediction", "OTHER",
    "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS Position",
)
_COL_ID, _COL_PRED, _COL_OTHER = 0, 1, 2
_COL_SP, _COL_LIPO, _COL_TAT, _COL_TATLIPO, _COL_PILIN = 3, 4, 5, 6, 7
_COL_CS = 8

# CS position cell looks like "CS pos: 20-21. Pr: 0.9756" — extract just
# the position range so downstream consumers don't have to re-parse.
_CS_POS_RE = re.compile(r"CS pos:\s*(\d+-\d+)")


def parse_signalp_output(results_path):
    """Parse SignalP 6.0 prediction_results.txt.

    Format (header line begins with `#`):
        # ID	Prediction	OTHER	SP(Sec/SPI)	LIPO(Sec/SPII)	TAT(Tat/SPI)
            TATLIPO(Tat/SPII)	PILIN(Sec/SPIII)	CS Position

    The ID column carries the full FASTA header (locus_tag + description);
    we keep only the first whitespace-delimited token so the locus_tag
    matches the rest of the pipeline.

    `signalp_probability` is the max across the five signal-peptide-type
    probabilities (SP/LIPO/TAT/TATLIPO/PILIN) — the OTHER column is
    excluded so non-secreted proteins report ~0 instead of ~1.
    """
    entries = []
    with open(results_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#'):
                # Two '#' lines precede the data: a metadata banner
                # (`# SignalP-6.0  Organism: Other  Timestamp: ...`) and
                # the column-header line (`# ID  Prediction  OTHER  ...`).
                # Validate the latter against the expected v6 schema so
                # format drift fails loudly rather than silently zeroing
                # all probabilities.
                header = tuple(
                    c.strip() for c in line.lstrip('#').strip().split('\t')
                )
                if header and header[0] == 'ID' and header != _SIGNALP6_COLS:
                    logger.warning(
                        "SignalP header drift — expected %s, got %s",
                        _SIGNALP6_COLS, header,
                    )
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) <= _COL_PILIN:
                continue

            locus_tag = parts[_COL_ID].split()[0] if parts[_COL_ID] else ""
            if not locus_tag:
                continue

            sp_probs = []
            for col in (_COL_SP, _COL_LIPO, _COL_TAT, _COL_TATLIPO, _COL_PILIN):
                try:
                    sp_probs.append(float(parts[col]))
                except ValueError:
                    logger.warning(
                        "Non-float in SignalP col %d for %s: %r",
                        col, locus_tag, parts[col],
                    )
                    sp_probs.append(0.0)
            sp_max = max(sp_probs) if sp_probs else 0.0

            cs_raw = parts[_COL_CS] if len(parts) > _COL_CS else ""
            cs_match = _CS_POS_RE.search(cs_raw)
            cs_position = cs_match.group(1) if cs_match else ""

            entries.append({
                'locus_tag': locus_tag,
                'signalp_prediction': parts[_COL_PRED],
                'signalp_probability': round(sp_max, 4),
                'signalp_cs_position': cs_position,
            })

    return entries


def main():
    parser = argparse.ArgumentParser(description="Run SignalP 6.0 (local or remote)")
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", choices=['local', 'remote'], default='remote',
                        help="local: DTU license needed. remote: uses DTU web server (free).")
    parser.add_argument("--signalp-path", default="", help="Path to SignalP install (local mode)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    try:
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
    except RuntimeError:
        raise  # Already has a user-friendly message
    except Exception as e:
        alt_mode = "remote" if args.mode == "local" else "local"
        raise RuntimeError(
            f"SignalP pipeline failed ({args.mode} mode): {e}\n"
            f"  How to fix:\n"
            f"    - Check the error details above\n"
            f"    - Try --mode {alt_mode} as an alternative"
        ) from e


if __name__ == '__main__':
    main()
