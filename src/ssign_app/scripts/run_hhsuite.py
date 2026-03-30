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
import re
import subprocess
import sys
import tempfile
import time

import requests

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

import os as _os, sys as _sys
_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)
from ssign_lib.fasta_io import read_fasta
from dedup_sequences import deduplicate_dict, expand_results_dict

MPI_BASE_URL = "https://toolkit.tuebingen.mpg.de"
MPI_DELAY = 20.0  # Seconds between submissions
MPI_POLL_INTERVAL = 10  # 10 seconds between polls (jobs take ~90s)

# MPI HHpred defaults — ALL fields required (missing fields cause silent failure)
# Discovered via browser curl capture: the session cookie + full payload are critical
HHPRED_DEFAULTS = {
    "msa_gen_method": "UniRef30",
    "msa_gen_max_iter": "3",
    "hhpred_incl_eval": "1e-3",
    "min_seqid_query": "0",
    "min_cov": "20",
    "ss_scoring": "2",
    "alignmacmode": "loc",
    "macthreshold": "0.3",
    "desc": "250",
    "pmin": "20",
    "proteomes": "",
}

# Regex for HHR hit table lines — works for both PDB (desc ends with ;) and Pfam (multiple ;)
# Fields (from right): (TLen) TRange QRange Cols SS Score Pval Eval Prob
_HHR_HIT_RE = re.compile(
    r"^\s*\d+\s+(\S+)\s+(.*?)\s+"
    r"(\d[\d.E+-]+)\s+"    # Prob
    r"([\d.E+-]+)\s+"      # E-value
    r"[\d.E+-]+\s+"        # P-value
    r"([\d.E+-]+)\s+"      # Score
    r"[\d.E+-]+\s+"        # SS
    r"\d+\s+"              # Cols
    r"\d+-\d+\s+"          # Query HMM range
    r"\d+-\d+\s+"          # Template HMM range
    r"\(\d+\)\s*$"         # Template length
)

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

        # FRAGILE: subprocess call requires hhblits binary on PATH
        # If this breaks: install HH-suite or switch to --mode remote
        try:
            subprocess.run(cmd_hhblits, capture_output=True, text=True,
                           timeout=600, check=True)
        except FileNotFoundError as e:
            raise RuntimeError(
                f"hhblits binary not found: {e}\n"
                f"  Common causes:\n"
                f"    - HH-suite is not installed or not on PATH\n"
                f"  How to fix:\n"
                f"    - conda install -c bioconda hhsuite\n"
                f"    - Or use --mode remote to submit to MPI Toolkit API"
            ) from e
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
            # FRAGILE: subprocess call requires hhsearch binary on PATH
            # If this breaks: install HH-suite or switch to --mode remote
            try:
                subprocess.run(cmd_hhsearch, capture_output=True, text=True,
                               timeout=600, check=True)
                hits = parse_hhr(hhr_file, db_name)
                if hits:
                    if pid not in results:
                        results[pid] = {'locus_tag': pid}
                    results[pid].update(hits)
            except FileNotFoundError as e:
                raise RuntimeError(
                    f"hhsearch binary not found: {e}\n"
                    f"  Common causes:\n"
                    f"    - HH-suite is not installed or not on PATH\n"
                    f"  How to fix:\n"
                    f"    - conda install -c bioconda hhsuite\n"
                    f"    - Or use --mode remote to submit to MPI Toolkit API"
                ) from e
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
            if in_hits and line.strip():
                m = _HHR_HIT_RE.match(line)
                if m:
                    hit_id, desc, prob, evalue, score = m.group(1, 2, 3, 4, 5)
                    return {
                        f'{db_prefix}_top1_id': hit_id,
                        f'{db_prefix}_top1_description': desc.strip('; ')[:200],
                        f'{db_prefix}_top1_probability': float(prob),
                        f'{db_prefix}_top1_evalue': float(evalue),
                        f'{db_prefix}_top1_score': float(score),
                    }
                break
    return {}


# ── Remote mode (MPI Toolkit API) ──

JOB_TIMEOUT = 600  # 10 min max per job (status 7 can hang indefinitely)
MAX_RETRIES = 1     # Retry once if stuck at MSA building


def _submit_one(session, pid, seq, db_config):
    """Submit a single HHpred job. Returns job_id or None."""
    if not seq.startswith(">"):
        fasta_seq = f">{pid}\n{seq}"
    else:
        fasta_seq = seq

    payload = {
        **HHPRED_DEFAULTS,
        "alignment": fasta_seq,  # CRITICAL: "alignment" NOT "sequence"!
        **db_config,
    }

    try:
        resp = session.post(
            f"{MPI_BASE_URL}/api/jobs/?toolName=hhpred",
            json=payload,  # MUST be json= not data= (415 error otherwise)
            timeout=30,
        )
        if resp.status_code == 200:
            data = resp.json()
            job_id = data.get("id", data.get("jobID", ""))
            if job_id:
                return job_id
        else:
            logger.warning(f"Submit failed for {pid}: HTTP {resp.status_code}")
    except Exception as e:
        logger.warning(f"Submit failed for {pid}: {e}")
    return None


def run_remote_hhpred(sequences, db_name="pdb"):
    """Submit sequences to MPI Toolkit HHpred API.

    Uses parallel polling: submit all jobs first, then poll all concurrently.
    Includes per-job timeout (10 min) and retry logic for stuck MSA building.

    CRITICAL: Must use requests.Session() for cookie handling.
    Without session cookies, jobs submit but silently fail during execution.
    """
    db_configs = {
        "pdb": {"hhsuitedb": "mmcif70/pdb70"},
        "pfam": {"hhsuitedb": "pfama/pfama"},
    }

    db_config = db_configs.get(db_name, db_configs["pdb"])
    results = {}

    # CRITICAL: Use Session for cookie handling — without this, jobs fail silently
    session = requests.Session()
    # FRAGILE: Session cookie acquisition is CRITICAL for MPI Toolkit.
    # Without the session cookie obtained from this GET request, job submissions
    # will appear to succeed (HTTP 200) but the jobs will silently fail during
    # execution on the server side. This is an undocumented MPI Toolkit requirement.
    # If this breaks: MPI Toolkit may be down or blocking requests
    try:
        session.get(f"{MPI_BASE_URL}/tools/hhpred", timeout=30)  # Get session cookie
    except requests.ConnectionError as e:
        raise RuntimeError(
            f"Cannot connect to MPI Toolkit to obtain session cookie: {e}\n"
            f"  Common causes:\n"
            f"    - MPI Toolkit (toolkit.tuebingen.mpg.de) is down\n"
            f"    - Network/firewall blocking outbound HTTPS\n"
            f"  How to fix:\n"
            f"    - Check https://toolkit.tuebingen.mpg.de manually\n"
            f"    - Or use --mode local with HH-suite installed:\n"
            f"      conda install -c bioconda hhsuite"
        ) from e
    except requests.Timeout as e:
        raise RuntimeError(
            f"MPI Toolkit timed out during session initialization: {e}\n"
            f"  How to fix:\n"
            f"    - Retry later, or use --mode local with HH-suite installed"
        ) from e

    # Phase 1: Submit all jobs
    job_ids = {}  # pid → job_id
    job_start = {}  # pid → start_time
    for pid, seq in sequences.items():
        job_id = _submit_one(session, pid, seq, db_config)
        if job_id:
            job_ids[pid] = job_id
            job_start[pid] = time.time()
            logger.info(f"Submitted {pid} -> job {job_id}")
        time.sleep(MPI_DELAY)

    if not job_ids:
        return results

    # Phase 2: Poll all jobs in parallel (round-robin)
    pending = dict(job_ids)  # pid → job_id (jobs still waiting)
    retried = set()  # pids that have been retried

    while pending:
        to_remove = []
        for pid, job_id in pending.items():
            elapsed = time.time() - job_start[pid]

            # Per-job timeout
            if elapsed > JOB_TIMEOUT:
                logger.warning(f"Job {job_id} ({pid}) timed out after {elapsed:.0f}s")

                # Retry once if we haven't already
                if pid not in retried:
                    retried.add(pid)
                    logger.info(f"Retrying {pid} (attempt 2)...")
                    # Get fresh session cookie before retry
                    try:
                        session.get(f"{MPI_BASE_URL}/tools/hhpred", timeout=30)
                    except Exception:
                        pass
                    new_job_id = _submit_one(session, pid, sequences[pid], db_config)
                    if new_job_id:
                        pending[pid] = new_job_id
                        job_start[pid] = time.time()
                        logger.info(f"Resubmitted {pid} -> job {new_job_id}")
                        time.sleep(MPI_DELAY)
                        continue
                to_remove.append(pid)
                continue

            # Poll
            try:
                resp = session.get(f"{MPI_BASE_URL}/api/jobs/{job_id}", timeout=15)
                if resp.status_code == 200:
                    data = resp.json()
                    status = data.get("status", 0)

                    if status == 5:  # Done
                        hits = collect_hhpred_result(session, job_id, db_name)
                        if hits:
                            results[pid] = {'locus_tag': pid, **hits}
                        to_remove.append(pid)
                        logger.info(f"Collected {pid} ({elapsed:.0f}s)")
                    elif status == 4:  # Error
                        logger.warning(f"Job {job_id} ({pid}) errored after {elapsed:.0f}s")
                        to_remove.append(pid)
            except Exception:
                pass  # Network blip, retry next round

        for pid in to_remove:
            pending.pop(pid, None)

        if pending:
            time.sleep(MPI_POLL_INTERVAL)

    return results


def poll_and_collect(session, job_id, protein_id, db_prefix, max_polls=60):
    """Poll a single MPI Toolkit job. Used for legacy compatibility.

    Timeout: 60 polls × 10s = 10 minutes max.
    """
    http = session if session else requests
    for _ in range(max_polls):
        try:
            resp = http.get(f"{MPI_BASE_URL}/api/jobs/{job_id}", timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                status = data.get("status", 0)

                if status == 5:  # Done
                    return collect_hhpred_result(http, job_id, db_prefix)
                elif status == 4:  # Error
                    logger.warning(f"Job {job_id} ({protein_id}) failed")
                    return {}
        except Exception:
            pass

        time.sleep(MPI_POLL_INTERVAL)

    logger.warning(f"Job {job_id} ({protein_id}) timed out after {max_polls * MPI_POLL_INTERVAL}s")
    return {}


def collect_hhpred_result(http, job_id, db_prefix):
    """Download and parse HHpred result from MPI Toolkit.

    URL discovered from toolkit JS source:
    /api/jobs/{jobID}/results/files/{jobID}.hhr
    (NOT /results/{jobID}/hhpred.hhr — that's the SPA HTML)
    """
    try:
        resp = http.get(
            f"{MPI_BASE_URL}/api/jobs/{job_id}/results/files/{job_id}.hhr",
            timeout=60,
        )
        if resp.status_code != 200:
            return {}

        # Parse HHR hit table using regex — handles both PDB and Pfam formats
        in_hits = False
        for line in resp.text.split('\n'):
            if line.startswith(' No Hit'):
                in_hits = True
                continue
            if in_hits and line.strip():
                m = _HHR_HIT_RE.match(line)
                if m:
                    hit_id, desc, prob, evalue, score = m.group(1, 2, 3, 4, 5)
                    return {
                        f'{db_prefix}_top1_id': hit_id,
                        f'{db_prefix}_top1_description': desc.strip('; ')[:200],
                        f'{db_prefix}_top1_probability': float(prob),
                        f'{db_prefix}_top1_evalue': float(evalue),
                        f'{db_prefix}_top1_score': float(score),
                    }
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
        # Deduplicate before remote submission to save MPI API quota
        unique_seqs, seq_groups = deduplicate_dict(sub_seqs)
        # Run against both PDB and Pfam remotely
        results_pdb = run_remote_hhpred(unique_seqs, "pdb")
        results_pfam = run_remote_hhpred(unique_seqs, "pfam")
        # Merge unique results
        results_unique = {}
        for pid in unique_seqs:
            entry = {'locus_tag': pid}
            if pid in results_pdb:
                entry.update(results_pdb[pid])
            if pid in results_pfam:
                entry.update(results_pfam[pid])
            if len(entry) > 1:
                results_unique[pid] = entry
        # Expand back to all duplicates
        results = expand_results_dict(results_unique, seq_groups)

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
