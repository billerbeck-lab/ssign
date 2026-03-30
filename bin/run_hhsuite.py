#!/usr/bin/env python3
"""Run HH-suite (local hhblits/hhsearch) or MPI Toolkit API (remote).

Local mode: runs hhblits for MSA generation, then hhsearch against Pfam-A and PDB70.
Remote mode: submits to MPI Toolkit HHpred API.

CRITICAL: Remote mode uses "alignment" parameter, NOT "sequence"!
Rate limit: 200 jobs/hour, 2000/day. Use >=20s delay.

Adapted from submit_passenger_hhpred.py
"""

import argparse
import csv
import json
import logging
import os
import subprocess
import sys
import tempfile
import time

import requests

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

from ssign_lib.fasta_io import read_fasta

MPI_BASE_URL = "https://toolkit.tuebingen.mpg.de"
MPI_DELAY = 20.0  # Seconds between submissions
MPI_POLL_INTERVAL = 60

# MPI HHpred defaults
HHPRED_DEFAULTS = {
    "msa_gen_method": "UniRef30",
    "msa_gen_max_iter": "3",
    "ss_scoring": "2",
    "alignmacmode": "loc",
    "macthreshold": "0.3",
    "desc": "250",
    "min_cov": "20",
    "pmin": "20",
}

STATUS_NAMES = {
    1: "submitted", 2: "queued", 3: "running", 4: "error",
    5: "done", 6: "warning", 7: "MSA_building", 8: "searching",
}


def load_substrate_ids(substrates_path):
    ids = set()
    with open(substrates_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ids.add(row['locus_tag'])
    return ids


# ── Local mode ──

def run_local_hhsuite(sequences, pfam_db, pdb70_db, uniclust_db, output_dir):
    """Run hhblits + hhsearch locally."""
    results = {}

    for pid, seq in sequences.items():
        query_file = os.path.join(output_dir, f"{pid}.fasta")
        with open(query_file, 'w') as f:
            f.write(f">{pid}\n{seq}\n")

        # Step 1: Generate MSA with hhblits
        a3m_file = os.path.join(output_dir, f"{pid}.a3m")
        cmd_hhblits = [
            "hhblits",
            "-i", query_file,
            "-d", uniclust_db,
            "-n", "2",
            "-o", "/dev/null",
            "-oa3m", a3m_file,
        ]

        try:
            subprocess.run(cmd_hhblits, capture_output=True, text=True,
                           timeout=600, check=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.warning(f"hhblits failed for {pid}: {e}")
            continue

        # Step 2: Search against Pfam and PDB70
        for db_name, db_path in [("pfam", pfam_db), ("pdb", pdb70_db)]:
            if not db_path:
                continue
            hhr_file = os.path.join(output_dir, f"{pid}_{db_name}.hhr")
            cmd_hhsearch = [
                "hhsearch",
                "-i", a3m_file,
                "-d", db_path,
                "-o", hhr_file,
            ]
            try:
                subprocess.run(cmd_hhsearch, capture_output=True, text=True,
                               timeout=600, check=True)
                hits = parse_hhr(hhr_file, db_name)
                if hits:
                    if pid not in results:
                        results[pid] = {'locus_tag': pid}
                    results[pid].update(hits)
            except Exception as e:
                logger.warning(f"hhsearch {db_name} failed for {pid}: {e}")

    return results


def parse_hhr(hhr_path, db_prefix):
    """Parse HHR result file and extract top hit."""
    if not os.path.exists(hhr_path):
        return {}

    with open(hhr_path) as f:
        in_hits = False
        for line in f:
            if line.startswith(' No Hit'):
                in_hits = True
                continue
            if in_hits and line.strip() and not line.startswith(' No'):
                # Parse hit line
                parts = line.split()
                if len(parts) >= 8:
                    hit_id = parts[1]
                    # Find probability and e-value
                    try:
                        prob = float(parts[-5]) if len(parts) > 5 else 0
                        evalue = float(parts[-4]) if len(parts) > 4 else 999
                        score = float(parts[-3]) if len(parts) > 3 else 0
                    except (ValueError, IndexError):
                        continue

                    desc = ' '.join(parts[1:-5]) if len(parts) > 6 else parts[1]

                    return {
                        f'{db_prefix}_top1_id': hit_id,
                        f'{db_prefix}_top1_description': desc[:200],
                        f'{db_prefix}_top1_probability': prob,
                        f'{db_prefix}_top1_evalue': evalue,
                        f'{db_prefix}_top1_score': score,
                    }
                break
    return {}


# ── Remote mode (MPI Toolkit API) ──

def run_remote_hhpred(sequences, db_name="pdb"):
    """Submit sequences to MPI Toolkit HHpred API."""
    db_configs = {
        "pdb": {"hhsuitedb": "mmcif70/pdb70"},
        "pfam": {"hhsuitedb": "pfama/pfama"},
    }

    results = {}
    job_ids = {}

    # Submit jobs
    for pid, seq in sequences.items():
        payload = {
            **HHPRED_DEFAULTS,
            "alignment": seq,  # CRITICAL: "alignment" NOT "sequence"!
            **db_configs.get(db_name, db_configs["pdb"]),
        }

        try:
            resp = requests.post(
                f"{MPI_BASE_URL}/api/jobs/?toolName=hhpred",
                data=payload,
                timeout=30,
            )
            if resp.status_code == 200:
                data = resp.json()
                job_id = data.get("id", data.get("jobID", ""))
                if job_id:
                    job_ids[pid] = job_id
                    logger.info(f"Submitted {pid} -> job {job_id}")
            else:
                logger.warning(f"Submit failed for {pid}: HTTP {resp.status_code}")
        except Exception as e:
            logger.warning(f"Submit failed for {pid}: {e}")

        time.sleep(MPI_DELAY)

    # Poll and collect results
    for pid, job_id in job_ids.items():
        hits = poll_and_collect(job_id, pid, db_name)
        if hits:
            results[pid] = {'locus_tag': pid, **hits}

    return results


def poll_and_collect(job_id, protein_id, db_prefix, max_polls=120):
    """Poll MPI Toolkit job and collect results."""
    for _ in range(max_polls):
        try:
            resp = requests.get(f"{MPI_BASE_URL}/api/jobs/{job_id}", timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                status = data.get("status", 0)

                if status == 5:  # Done
                    return collect_hhpred_result(job_id, db_prefix)
                elif status == 4:  # Error
                    logger.warning(f"Job {job_id} ({protein_id}) failed")
                    return {}
        except Exception:
            pass

        time.sleep(MPI_POLL_INTERVAL)

    logger.warning(f"Job {job_id} ({protein_id}) timed out")
    return {}


def collect_hhpred_result(job_id, db_prefix):
    """Download and parse HHpred result from MPI Toolkit."""
    try:
        resp = requests.get(
            f"{MPI_BASE_URL}/results/{job_id}/hhpred.hhr",
            timeout=60,
        )
        if resp.status_code != 200:
            return {}

        # Parse HHR content
        in_hits = False
        for line in resp.text.split('\n'):
            if line.startswith(' No Hit'):
                in_hits = True
                continue
            if in_hits and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        return {
                            f'{db_prefix}_top1_id': parts[1],
                            f'{db_prefix}_top1_description': ' '.join(parts[1:-5])[:200],
                            f'{db_prefix}_top1_probability': float(parts[-5]),
                            f'{db_prefix}_top1_evalue': float(parts[-4]),
                            f'{db_prefix}_top1_score': float(parts[-3]),
                        }
                    except (ValueError, IndexError):
                        pass
                break
    except Exception as e:
        logger.warning(f"Failed to collect {job_id}: {e}")

    return {}


def main():
    parser = argparse.ArgumentParser(description="Run HH-suite (local or remote)")
    parser.add_argument("--mode", choices=['local', 'remote'], required=True)
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pfam-db", default="")
    parser.add_argument("--pdb70-db", default="")
    parser.add_argument("--uniclust-db", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    substrate_ids = load_substrate_ids(args.substrates)
    all_seqs = read_fasta(args.proteins)
    sub_seqs = {k: v for k, v in all_seqs.items() if k in substrate_ids}

    logger.info(f"Processing {len(sub_seqs)} substrate proteins for {args.sample}")

    if args.mode == 'local':
        with tempfile.TemporaryDirectory() as tmpdir:
            results = run_local_hhsuite(
                sub_seqs, args.pfam_db, args.pdb70_db, args.uniclust_db, tmpdir
            )
    else:
        # Run against both PDB and Pfam remotely
        results_pdb = run_remote_hhpred(sub_seqs, "pdb")
        results_pfam = run_remote_hhpred(sub_seqs, "pfam")
        # Merge
        results = {}
        for pid in sub_seqs:
            entry = {'locus_tag': pid}
            if pid in results_pdb:
                entry.update(results_pdb[pid])
            if pid in results_pfam:
                entry.update(results_pfam[pid])
            if len(entry) > 1:
                results[pid] = entry

    # Determine output columns from results
    all_cols = set()
    for r in results.values():
        all_cols.update(r.keys())
    fieldnames = ['locus_tag'] + sorted(c for c in all_cols if c != 'locus_tag')

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results.values():
            writer.writerow(r)

    logger.info(f"Wrote {len(results)} HH-suite annotations for {args.sample}")


if __name__ == '__main__':
    main()
