#!/usr/bin/env python3
"""Test full DTU DeepLocPro web form submission flow."""
import requests
import re
import time
import json
import sys

DTU_SUBMIT = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
DTU_RESULTS = "https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp"

fasta = (
    ">NJFHAN_00001\n"
    "MPDLLLELFSEEIPARMQAQAAEALRKLVTDKLVERGLIYEGAKAFVTPRRLALSVHGLP\n"
    "GRQADQKEEKKGPRVGAPEALLLIALLREEGIDPELREDLLPHKASFLLDAGFDPIQRHL\n"
    "DTLSAQEAELFRQMMQRLGYSPEQLEAMLLQFNRHFPDVLADTRSLFAEMTKEVFWQLVG\n"
    "EAAKAGQTVTISGDITDDNHDFKRTGYRYGFCTDAWSFDARLRRTFDEACSAGCADMVFS\n"
)

data = {
    "configfile": "/var/www/services/services/DeepLocPro-1.0/webface.cf",
    "fasta": fasta,
    "group": "negative",
    "format": "short",
}

# Add browser-like headers
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
    "Origin": "https://services.healthtech.dtu.dk",
    "Referer": "https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/",
}

print("Step 1: Submitting...")
resp = requests.post(DTU_SUBMIT, data=data, headers=headers, timeout=30)
print(f"  Status: {resp.status_code}, URL: {resp.url}")

job_match = re.search(r"jobid=([A-F0-9]+)", resp.url)
if not job_match:
    print("ERROR: No job ID")
    print(resp.text[:500])
    sys.exit(1)

job_id = job_match.group(1)
print(f"  Job ID: {job_id}")

print("Step 2: Polling with ajax=1...")
for i in range(60):
    time.sleep(5)
    ajax_url = f"{DTU_SUBMIT}?ajax=1&jobid={job_id}"
    r = requests.get(ajax_url, headers=headers, timeout=15)

    try:
        status_data = r.json()
        status = status_data.get("status", "unknown")
        print(f"  Poll {i+1}: status={status}")

        if status == "finished":
            print("Step 3: Fetching results.json...")
            results_url = f"{DTU_RESULTS}/{job_id}/results.json"
            rr = requests.get(results_url, timeout=15)
            if rr.status_code == 200:
                results = rr.json()
                print(f"  Sequences: {results.get('info', {}).get('size', 0)}")
                for name, seq_data in results.get("sequences", {}).items():
                    print(f"  {name}: {seq_data['Prediction']}")
                    locs = results.get("Localization", [])
                    probs = seq_data.get("Probability", [])
                    for loc, prob in zip(locs, probs):
                        print(f"    {loc}: {prob}")

                # Also fetch CSV
                csv_path = results.get("csv_file", "")
                if csv_path:
                    csv_url = f"https://services.healthtech.dtu.dk{csv_path}"
                    cr = requests.get(csv_url, timeout=15)
                    print(f"\n  CSV ({cr.status_code}):")
                    print(f"  {cr.text[:500]}")
            else:
                print(f"  Results fetch failed: {rr.status_code}")
            break

        elif status in ("error", "failed"):
            print(f"  JOB FAILED: {status_data}")
            break

    except ValueError:
        # Not JSON — still loading
        print(f"  Poll {i+1}: not JSON yet (len={len(r.text)})")
else:
    print("TIMEOUT")

print("\nDone.")
