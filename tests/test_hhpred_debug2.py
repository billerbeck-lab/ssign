#!/usr/bin/env python3
"""Debug HHpred: test proteins individually against PDB and Pfam.

Tests with short and long proteins to identify the hanging pattern.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'ssign_app', 'scripts'))

import requests
from run_hhsuite import (
    MPI_BASE_URL, MPI_DELAY, HHPRED_DEFAULTS, STATUS_NAMES, _HHR_HIT_RE,
)

# Test proteins: 2 short, 2 long (representing typical substrates)
TEST_PROTEINS = {
    "SHORT_90aa": "MKKTAIAIAVALAGFATVAQAATENGHTGASATTLNNETATSFKVGDTSTRSSSTLGRFKSNQDNGNTTPVSSQTDGTSTQTTGAGKYDYENT",
    "MEDIUM_200aa": (
        "MRRFLAFVLGALALGQALPAVAAHTAAETDTQAANGKDAADTKAAAQ"
        "KASDEATKAATQAAEEAKTAAAQKAAETAAAQKAAEQKAAETAAAQKA"
        "AEQKAAETAAAQKAAEQKAAETAAAQKAAEQKAAETAAAQKAAEQKAA"
        "ETAAAQKAAEQKAAETAAAQKAAEQKAAETAAAQKAAEQKAAETKAAT"
        "KAATQAAEEAKTAAAQ"
    ),
    "LONG_600aa": (
        "MKYLLPTAAAGLLLLAAQPAMAQVQLVESGGGLVQAGGSLRLSCAASG"
        "FTFSNYAMSWVRQAPGKGLEWVSAISSNGGSTYYADSVKGRFTISRDN"
        "AKNTVYLQMNSLRAEDTAVYYCAARGRGSTFSGYYRGQVTVSSASTKG"
        "PSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTF"
        "PAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEP"
        "KSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVD"
        "VSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDW"
        "LNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQ"
        "VSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLT"
        "VDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGKAATGGPSVFPLAP"
        "SSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGL"
        "YSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKAAETAAAQK"
    ),
}

DB_CONFIGS = {
    "pdb": {"hhsuitedb": "mmcif70/pdb70"},
    "pfam": {"hhsuitedb": "pfama/pfama"},
}

print(f"Testing {len(TEST_PROTEINS)} proteins × {len(DB_CONFIGS)} databases")
print(f"Max timeout per job: 600s (10 min)")
print()

# Get session cookie
session = requests.Session()
session.get(f"{MPI_BASE_URL}/tools/hhpred", timeout=30)
print("Got MPI session cookie\n")

results = {}

for db_name, db_config in DB_CONFIGS.items():
    print(f"\n{'='*60}")
    print(f"DATABASE: {db_name}")
    print(f"{'='*60}")

    for pid, seq in TEST_PROTEINS.items():
        print(f"\n--- {pid} ({len(seq)}aa) vs {db_name} ---")
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
            if resp.status_code != 200:
                print(f"  SUBMIT FAILED: HTTP {resp.status_code}: {resp.text[:200]}")
                results[(pid, db_name)] = {"status": "submit_failed", "time": 0}
                time.sleep(MPI_DELAY)
                continue
            data = resp.json()
            job_id = data.get("id", data.get("jobID", ""))
            print(f"  Submitted -> {job_id}")
        except Exception as e:
            print(f"  SUBMIT ERROR: {e}")
            results[(pid, db_name)] = {"status": "submit_error", "time": 0}
            time.sleep(MPI_DELAY)
            continue

        # Poll with status tracking
        last_status = None
        status_times = {}
        for _ in range(60):  # 60 × 10s = 600s max
            time.sleep(10)
            elapsed = time.time() - t0

            try:
                resp = session.get(f"{MPI_BASE_URL}/api/jobs/{job_id}", timeout=30)
                data = resp.json()
                status = data.get("status", 0)
                status_name = STATUS_NAMES.get(status, f"unknown({status})")

                if status != last_status:
                    print(f"  [{elapsed:5.0f}s] Status {status} ({status_name})")
                    if status not in status_times:
                        status_times[status] = elapsed
                    last_status = status

                if status == 5:  # Done
                    # Collect top hit
                    try:
                        resp2 = session.get(
                            f"{MPI_BASE_URL}/api/jobs/{job_id}/results/files/{job_id}.hhr",
                            timeout=60,
                        )
                        in_hits = False
                        for line in resp2.text.split('\n'):
                            if line.startswith(' No Hit'):
                                in_hits = True
                                continue
                            if in_hits and line.strip():
                                m = _HHR_HIT_RE.match(line)
                                if m:
                                    print(f"  HIT: {m.group(1)} | prob={m.group(3)}")
                                    break
                    except Exception:
                        pass
                    print(f"  DONE in {elapsed:.0f}s")
                    results[(pid, db_name)] = {
                        "status": "done", "time": elapsed,
                        "status_times": status_times,
                    }
                    break

                elif status == 4:  # Error
                    print(f"  ERROR after {elapsed:.0f}s")
                    results[(pid, db_name)] = {
                        "status": "error", "time": elapsed,
                        "status_times": status_times,
                    }
                    break

            except Exception as e:
                print(f"  [{elapsed:5.0f}s] Poll error: {e}")

            if elapsed > 600:
                print(f"  TIMEOUT at {elapsed:.0f}s (last: {status_name})")
                results[(pid, db_name)] = {
                    "status": "timeout", "time": elapsed,
                    "last_status": last_status,
                    "status_times": status_times,
                }
                break
        else:
            elapsed = time.time() - t0
            print(f"  LOOP END at {elapsed:.0f}s (last: {STATUS_NAMES.get(last_status, '?')})")
            results[(pid, db_name)] = {
                "status": "loop_end", "time": elapsed,
                "last_status": last_status,
                "status_times": status_times,
            }

        time.sleep(MPI_DELAY)

# Summary
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
print(f"{'Protein':20s} {'DB':6s} {'Status':12s} {'Time':>6s}  Status flow")
print("-" * 70)
for (pid, db), info in sorted(results.items()):
    flow = " → ".join(
        f"{STATUS_NAMES.get(s, s)}"
        for s in sorted(info.get('status_times', {}).keys())
    )
    print(f"{pid:20s} {db:6s} {info['status']:12s} {info['time']:5.0f}s  {flow}")
