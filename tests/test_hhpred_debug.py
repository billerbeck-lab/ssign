#!/usr/bin/env python3
"""Debug HHpred: test each protein × each database individually.

Identifies which combinations hang at status 7 (MSA building).
"""
import sys, os, time, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'ssign_app', 'scripts'))

import requests
import csv
from run_hhsuite import (
    MPI_BASE_URL, MPI_DELAY, MPI_POLL_INTERVAL,
    HHPRED_DEFAULTS, STATUS_NAMES, _HHR_HIT_RE,
)
from ssign_lib.fasta_io import read_fasta

# Load the 4 substrate proteins
SUBSTRATES_TSV = "/tmp/ssign_multi_test/xanthobacter/substrates_filtered.tsv"
PROTEINS_FAA = "/tmp/ssign_multi_test/xanthobacter/proteins.faa"

sub_ids = set()
with open(SUBSTRATES_TSV) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        sub_ids.add(row['locus_tag'])

all_seqs = read_fasta(PROTEINS_FAA)
sub_seqs = {k: v for k, v in all_seqs.items() if k in sub_ids}

print(f"Testing {len(sub_seqs)} proteins × 2 databases")
print(f"Proteins: {list(sub_seqs.keys())}")
print()

DB_CONFIGS = {
    "pdb": {"hhsuitedb": "mmcif70/pdb70"},
    "pfam": {"hhsuitedb": "pfama/pfama"},
}

# Test each protein × each database individually
session = requests.Session()
session.get(f"{MPI_BASE_URL}/tools/hhpred", timeout=30)
print("Got session cookie")

results = {}
TIMEOUT_PER_JOB = 600  # 10 min max per job

for db_name, db_config in DB_CONFIGS.items():
    print(f"\n{'='*60}")
    print(f"DATABASE: {db_name}")
    print(f"{'='*60}")

    for pid, seq in sub_seqs.items():
        print(f"\n--- {pid} vs {db_name} ---")
        fasta_seq = f">{pid}\n{seq}"

        payload = {
            **HHPRED_DEFAULTS,
            "alignment": fasta_seq,
            **db_config,
        }

        # Submit
        t0 = time.time()
        try:
            resp = session.post(
                f"{MPI_BASE_URL}/api/jobs/?toolName=hhpred",
                json=payload, timeout=30,
            )
            data = resp.json()
            job_id = data.get("id", data.get("jobID", ""))
            print(f"  Submitted -> {job_id} (HTTP {resp.status_code})")
        except Exception as e:
            print(f"  SUBMIT FAILED: {e}")
            continue

        # Poll with detailed status tracking
        last_status = None
        status_history = []
        done = False
        for poll in range(int(TIMEOUT_PER_JOB / 10)):
            time.sleep(10)
            elapsed = time.time() - t0

            try:
                resp = session.get(f"{MPI_BASE_URL}/api/jobs/{job_id}", timeout=30)
                data = resp.json()
                status = data.get("status", 0)
                status_name = STATUS_NAMES.get(status, f"unknown({status})")

                if status != last_status:
                    print(f"  [{elapsed:.0f}s] Status: {status} ({status_name})")
                    status_history.append((elapsed, status, status_name))
                    last_status = status

                if status == 5:  # Done
                    # Collect result
                    resp2 = session.get(
                        f"{MPI_BASE_URL}/api/jobs/{job_id}/results/files/{job_id}.hhr",
                        timeout=60,
                    )
                    if resp2.status_code == 200:
                        # Parse top hit
                        in_hits = False
                        for line in resp2.text.split('\n'):
                            if line.startswith(' No Hit'):
                                in_hits = True
                                continue
                            if in_hits and line.strip():
                                m = _HHR_HIT_RE.match(line)
                                if m:
                                    hit_id, desc, prob = m.group(1), m.group(2), m.group(3)
                                    print(f"  TOP HIT: {hit_id} | {desc[:60]} | prob={prob}")
                                    break
                    print(f"  DONE in {elapsed:.0f}s")
                    results[(pid, db_name)] = {"status": "done", "time": elapsed}
                    done = True
                    break
                elif status == 4:  # Error
                    print(f"  ERROR after {elapsed:.0f}s")
                    results[(pid, db_name)] = {"status": "error", "time": elapsed}
                    done = True
                    break

            except Exception as e:
                print(f"  [{elapsed:.0f}s] Poll error: {e}")

            # Hard timeout
            if elapsed > TIMEOUT_PER_JOB:
                print(f"  TIMEOUT after {elapsed:.0f}s (last status: {last_status})")
                results[(pid, db_name)] = {
                    "status": "timeout",
                    "time": elapsed,
                    "last_status": last_status,
                    "history": status_history,
                }
                done = True
                break

        if not done:
            print(f"  LOOP EXHAUSTED")
            results[(pid, db_name)] = {"status": "loop_exhausted", "time": time.time() - t0}

        # Rate limit delay before next submission
        time.sleep(MPI_DELAY)

# Summary
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
for (pid, db), info in sorted(results.items()):
    status = info['status']
    t = info['time']
    extra = ""
    if status == "timeout":
        extra = f" (stuck at status {info.get('last_status', '?')})"
    print(f"  {pid:20s} x {db:5s} -> {status:12s} ({t:.0f}s){extra}")
