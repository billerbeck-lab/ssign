#!/usr/bin/env python3
"""Run InterProScan in local or remote mode.

Local: runs interproscan.sh locally.
Remote: submits to EBI REST API (30 req/sec, 25k seq/day).

Adapted from pipeline/scripts/parse_interproscan.py — preserves TSV parsing
with column indices and GO extraction logic.
"""

import argparse
import csv
import json
import logging
import os
import re
import subprocess
import sys
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
from dedup_sequences import deduplicate_dict, expand_results_dict

# InterProScan TSV column indices (0-based, no header)
_COL_PROTEIN_ID = 0
_COL_SIG_ACC = 4
_COL_SIG_DESC = 5
_COL_IPR_ACC = 11
_COL_IPR_DESC = 12
_COL_GO_TERMS = 13
_MISSING = "-"
_GO_ID_RE = re.compile(r"(GO:\d+)")


def load_substrate_ids(substrates_path):
    ids = set()
    with open(substrates_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ids.add(row['locus_tag'])
    return ids


# ── Local mode ──

def run_local_interproscan(query_fasta, db_path, output_dir):
    """Run InterProScan locally."""
    output_file = os.path.join(output_dir, "results.tsv")
    cmd = [
        "interproscan.sh",
        "-i", query_fasta,
        "-o", output_file,
        "-f", "tsv",
        "-goterms",
        "-pathways",
        "-dp",
    ]
    if db_path:
        cmd.extend(["-d", db_path])

    logger.info("Running local InterProScan...")
    # FRAGILE: subprocess call requires interproscan.sh on PATH
    # If this breaks: install InterProScan locally or switch to --mode remote
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"InterProScan (interproscan.sh) not found: {e}\n"
            f"  Common causes:\n"
            f"    - InterProScan is not installed or not on PATH\n"
            f"  How to fix:\n"
            f"    - Download from https://www.ebi.ac.uk/interpro/download/\n"
            f"    - Or use --mode remote to submit to EBI REST API"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"InterProScan timed out after 4 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input sequences\n"
            f"    - Or use --mode remote"
        ) from e

    if result.returncode != 0:
        logger.error(f"InterProScan failed: {result.stderr[:500]}")
        raise RuntimeError(f"InterProScan exit code {result.returncode}")

    return output_file


# ── Remote mode (EBI REST API) ──

EBI_BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"


def run_remote_interproscan(sequences, output_dir):
    """Submit sequences to EBI InterProScan REST API."""
    results_file = os.path.join(output_dir, "results.tsv")
    all_results = []

    for pid, seq in sequences.items():
        try:
            # Submit job
            # Send as FASTA with header so EBI preserves the protein ID
            fasta_seq = f">{pid}\n{seq}"
            # NOTE: Hardcoded email for EBI API. EBI requires a valid email for
            # job submissions but does not verify it. Change if needed.
            # FRAGILE: EBI REST API submission can fail due to server issues or rate limits
            # If this breaks: check https://www.ebi.ac.uk/Tools/services/rest/iprscan5 status
            resp = http_requests.post(
                f"{EBI_BASE}/run",
                data={
                    "email": "ssign-pipeline@example.com",
                    "sequence": fasta_seq,
                    "goterms": "true",
                    "pathways": "true",
                },
                timeout=30,
            )
            if resp.status_code != 200:
                logger.warning(f"Submit failed for {pid}: HTTP {resp.status_code}")
                continue

            job_id = resp.text.strip()
            logger.info(f"Submitted {pid} -> job {job_id}")

            # Poll for completion
            for _ in range(120):
                time.sleep(15)
                status_resp = http_requests.get(
                    f"{EBI_BASE}/status/{job_id}", timeout=30
                )
                status = status_resp.text.strip()
                if status == "FINISHED":
                    break
                elif status in ("FAILURE", "ERROR", "NOT_FOUND"):
                    logger.warning(f"Job {job_id} failed: {status}")
                    break
            else:
                logger.warning(f"Job {job_id} timed out")
                continue

            if status == "FINISHED":
                tsv_resp = http_requests.get(
                    f"{EBI_BASE}/result/{job_id}/tsv", timeout=60
                )
                if tsv_resp.status_code == 200:
                    all_results.append(tsv_resp.text)

        except (http_requests.ConnectionError, http_requests.Timeout) as e:
            logger.warning(
                f"EBI InterProScan API connection failed for {pid}: {e}\n"
                f"  EBI servers may be under maintenance.\n"
                f"  Consider using --mode local if this persists."
            )
        except Exception as e:
            logger.warning(f"Failed for {pid}: {e}")

        time.sleep(1)  # Rate limiting

    # Combine all TSV results
    with open(results_file, 'w') as f:
        for chunk in all_results:
            f.write(chunk)

    return results_file


# ── Parsing ──

def parse_interproscan_tsv(tsv_path, target_ids=None):
    """Parse InterProScan TSV output and aggregate per protein.

    Preserves column index parsing from pipeline/scripts/parse_interproscan.py.
    """
    per_protein = {}

    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            protein_id = parts[_COL_PROTEIN_ID]
            if target_ids and protein_id not in target_ids:
                continue

            if protein_id not in per_protein:
                per_protein[protein_id] = {
                    'domains': set(),
                    'go_terms': set(),
                    'pfam_ids': set(),
                    'descriptions': set(),
                }

            sig_acc = parts[_COL_SIG_ACC] if len(parts) > _COL_SIG_ACC else _MISSING
            sig_desc = parts[_COL_SIG_DESC] if len(parts) > _COL_SIG_DESC else _MISSING
            ipr_acc = parts[_COL_IPR_ACC] if len(parts) > _COL_IPR_ACC else _MISSING
            ipr_desc = parts[_COL_IPR_DESC] if len(parts) > _COL_IPR_DESC else _MISSING
            go_raw = parts[_COL_GO_TERMS] if len(parts) > _COL_GO_TERMS else _MISSING

            if ipr_acc != _MISSING:
                per_protein[protein_id]['domains'].add(ipr_acc)
            if ipr_desc != _MISSING:
                per_protein[protein_id]['descriptions'].add(ipr_desc)
            if sig_acc != _MISSING and sig_acc.startswith('PF'):
                per_protein[protein_id]['pfam_ids'].add(sig_acc)

            if go_raw != _MISSING:
                for match in _GO_ID_RE.finditer(go_raw):
                    per_protein[protein_id]['go_terms'].add(match.group(1))

    # Build output rows
    results = {}
    for pid, data in per_protein.items():
        results[pid] = {
            'locus_tag': pid,
            'interpro_domains': ';'.join(sorted(data['domains'])),
            'interpro_go_terms': ';'.join(sorted(data['go_terms'])),
            'interpro_pfam_ids': ';'.join(sorted(data['pfam_ids'])),
            'interpro_descriptions': ';'.join(sorted(data['descriptions'])),
        }

    return results


def main():
    parser = argparse.ArgumentParser(description="Run InterProScan")
    parser.add_argument("--mode", choices=['local', 'remote'], required=True)
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--db", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    substrate_ids = load_substrate_ids(args.substrates)
    all_seqs = read_fasta(args.proteins)
    sub_seqs = {k: v for k, v in all_seqs.items() if k in substrate_ids}

    # Deduplicate before remote submission to save EBI API quota
    if args.mode == 'remote':
        unique_seqs, seq_groups = deduplicate_dict(sub_seqs)
    else:
        unique_seqs, seq_groups = sub_seqs, {k: [k] for k in sub_seqs}

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write unique substrate sequences to temp FASTA
        tmp_fasta = os.path.join(tmpdir, "substrates.fasta")
        with open(tmp_fasta, 'w') as f:
            for pid, seq in unique_seqs.items():
                f.write(f">{pid}\n{seq}\n")

        if args.mode == 'local':
            tsv_path = run_local_interproscan(tmp_fasta, args.db, tmpdir)
        else:
            tsv_path = run_remote_interproscan(unique_seqs, tmpdir)

        results_unique = parse_interproscan_tsv(tsv_path, set(unique_seqs.keys()))

    # Expand results back to all duplicates
    results = expand_results_dict(results_unique, seq_groups)

    logger.info(f"Annotated {len(results)}/{len(substrate_ids)} substrates "
                f"for {args.sample}")

    fieldnames = ['locus_tag', 'interpro_domains', 'interpro_go_terms',
                  'interpro_pfam_ids', 'interpro_descriptions']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results.values():
            writer.writerow(r)


if __name__ == '__main__':
    main()
