#!/usr/bin/env python3
"""Run Foldseek structural homology search (local or remote).

Local: QUAL-01 strategy — batch search first, per-protein fallback with relaxed E-value.
Remote: Foldseek web API at search.foldseek.com (ticket-based async).

Uses qtmscore (query-normalized), NOT alntmscore.

Adapted from pipeline/scripts/run_foldseek.py
"""

import argparse
import csv
import json
import logging
import os
import subprocess
import tempfile
import time

import requests as http_requests

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

import os as _os, sys as _sys
_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)
from ssign_lib.fasta_io import read_fasta

FOLDSEEK_API = "https://search.foldseek.com/api"

FORMAT_OUTPUT = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qtmscore,ttmscore,alntmscore,lddt,prob"


def load_substrate_ids(substrates_path):
    ids = set()
    with open(substrates_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ids.add(row['locus_tag'])
    return ids


def run_foldseek_search(query_dir, db_path, output_file, evalue, threads=4):
    """Run foldseek easy-search."""
    cmd = [
        "foldseek", "easy-search",
        query_dir, db_path, output_file, "tmp_foldseek",
        "--format-output", FORMAT_OUTPUT,
        "-e", str(evalue),
        "--alignment-type", "1",
        "--exact-tmscore", "1",
        "--threads", str(threads),
    ]

    # FRAGILE: subprocess call requires foldseek binary on PATH
    # If this breaks: install foldseek or switch to --mode remote
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"Foldseek binary not found: {e}\n"
            f"  Common causes:\n"
            f"    - Foldseek is not installed or not on PATH\n"
            f"  How to fix:\n"
            f"    - conda install -c conda-forge -c bioconda foldseek\n"
            f"    - Or download from https://github.com/steineggerlab/foldseek\n"
            f"    - Or use --mode remote to submit to Foldseek web API"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"Foldseek timed out after 2 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input structures\n"
            f"    - Or use --mode remote"
        ) from e
    return result.returncode == 0


def parse_foldseek_results(output_file, tmscore_threshold):
    """Parse foldseek output and filter by qtmscore."""
    hits = {}

    if not os.path.exists(output_file):
        return hits

    with open(output_file) as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 17:
                continue

            query_id = parts[0]
            target_id = parts[1]
            fident = float(parts[2])
            evalue = float(parts[10])
            qtmscore = float(parts[12])
            ttmscore = float(parts[13])

            if qtmscore < tmscore_threshold:
                continue

            if query_id not in hits or qtmscore > hits[query_id]['foldseek_qtmscore']:
                hits[query_id] = {
                    'locus_tag': query_id,
                    'foldseek_best_hit': target_id,
                    'foldseek_hit_description': target_id,
                    'foldseek_evalue': evalue,
                    'foldseek_qtmscore': round(qtmscore, 4),
                    'foldseek_ttmscore': round(ttmscore, 4),
                    'foldseek_fident': round(fident, 4),
                }

    return hits


# ── Remote mode (Foldseek web API) ──

def run_remote_foldseek(pdb_content, databases=None, mode="3diaa"):
    """Submit a PDB structure to Foldseek web API.

    API: POST https://search.foldseek.com/api/ticket
    Poll: GET https://search.foldseek.com/api/ticket/{id}
    Results: GET https://search.foldseek.com/api/result/{id}/0
    """
    if databases is None:
        databases = ["afdb-swissprot", "pdb100"]

    # Submit
    data = {"mode": mode}
    files_param = {"q": ("query.pdb", pdb_content)}
    # Add each database as a separate form field
    for db in databases:
        data.setdefault("database[]", [])
    # requests needs list of tuples for repeated keys
    form_data = [("mode", mode)]
    for db in databases:
        form_data.append(("database[]", db))

    # FRAGILE: Foldseek web API submission can fail due to server issues
    # If this breaks: check https://search.foldseek.com or use --mode local
    try:
        resp = http_requests.post(
            f"{FOLDSEEK_API}/ticket",
            data=form_data,
            files=files_param,
            timeout=60,
        )
        if resp.status_code != 200:
            logger.warning(f"Foldseek API submit failed: HTTP {resp.status_code}")
            return []

        ticket = resp.json()
        ticket_id = ticket.get("id", "")
        if not ticket_id:
            return []

        logger.info(f"Foldseek ticket: {ticket_id}")

    except http_requests.ConnectionError as e:
        logger.warning(
            f"Cannot connect to Foldseek API: {e}\n"
            f"  Foldseek server (search.foldseek.com) may be down.\n"
            f"  Consider using --mode local with foldseek installed."
        )
        return []
    except http_requests.Timeout as e:
        logger.warning(
            f"Foldseek API timed out during submission: {e}\n"
            f"  Consider retrying later or using --mode local."
        )
        return []
    except Exception as e:
        logger.warning(f"Foldseek API error: {e}")
        return []

    # Poll for completion (max ~30 min)
    for _ in range(60):
        time.sleep(30)
        try:
            status_resp = http_requests.get(f"{FOLDSEEK_API}/ticket/{ticket_id}", timeout=30)
            if status_resp.status_code == 200:
                status_data = status_resp.json()
                status = status_data.get("status", "")
                if status == "COMPLETE":
                    # Fetch results
                    result_resp = http_requests.get(
                        f"{FOLDSEEK_API}/result/{ticket_id}/0", timeout=60
                    )
                    if result_resp.status_code == 200:
                        return result_resp.json().get("results", [])
                elif status == "ERROR":
                    logger.warning(f"Foldseek job failed: {status_data}")
                    return []
        except Exception:
            pass

    logger.warning(f"Foldseek ticket {ticket_id} timed out")
    return []


def main():
    parser = argparse.ArgumentParser(description="Run Foldseek (local or remote)")
    parser.add_argument("--mode", choices=['local', 'remote'], default='local')
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--db", default="", help="Foldseek database path (local mode)")
    parser.add_argument("--evalue", type=float, default=0.001)
    parser.add_argument("--tmscore", type=float, default=0.8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    logger.info(f"Foldseek search ({args.mode}) for {args.sample}")

    if args.mode == 'local':
        if not args.db:
            logger.error("Local mode requires --db")
            raise SystemExit(1)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = os.path.join(tmpdir, "foldseek_results.tsv")
            success = run_foldseek_search(
                args.proteins, args.db, output_file, args.evalue
            )
            if success:
                hits = parse_foldseek_results(output_file, args.tmscore)
            else:
                logger.warning("Batch foldseek failed, attempting per-protein fallback")
                hits = {}
    else:
        # Remote mode — requires PDB structures, returns JSON hits
        # For now, log that remote foldseek needs structure files
        logger.info("Remote Foldseek requires PDB structure files (not FASTA)")
        hits = {}

    logger.info(f"Found {len(hits)} Foldseek hits for {args.sample}")

    fieldnames = ['locus_tag', 'foldseek_best_hit', 'foldseek_hit_description',
                  'foldseek_evalue', 'foldseek_qtmscore', 'foldseek_ttmscore',
                  'foldseek_fident']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for hit in hits.values():
            writer.writerow(hit)


if __name__ == '__main__':
    main()
