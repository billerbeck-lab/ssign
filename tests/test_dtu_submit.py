#!/usr/bin/env python3
"""Test DTU DeepLocPro web form submission."""
import requests
import re
import time
import sys

# Real protein from test genome
fasta = (
    ">NJFHAN_00001\n"
    "MPDLLLELFSEEIPARMQAQAAEALRKLVTDKLVERGLIYEGAKAFVTPRRLALSVHGLP\n"
    "GRQADQKEEKKGPRVGAPEALLLIALLREEGIDPELREDLLPHKASFLLDAGFDPIQRHL\n"
    "DTLSAQEAELFRQMMQRLGYSPEQLEAMLLQFNRHFPDVLADTRSLFAEMTKEVFWQLVG\n"
    "EAAKAGQTVTISGDITDDNHDFKRTGYRYGFCTDAWSFDARLRRTFDEACSAGCADMVFS\n"
)

url = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
data = {
    "configfile": "/var/www/services/services/DeepLocPro-1.0/webface.cf",
    "fasta": fasta,
    "group": "negative",
    "format": "short",
}

print("Submitting to DTU DeepLocPro...")
resp = requests.post(url, data=data, timeout=30)
print(f"Response URL: {resp.url}")

job_match = re.search(r"jobid=([A-F0-9]+)", resp.url)
if not job_match:
    job_match = re.search(r"([A-F0-9]{24,})", resp.text)

if not job_match:
    print("ERROR: No job ID found")
    print(resp.text[:500])
    sys.exit(1)

job_id = job_match.group(1)
print(f"Job ID: {job_id}")

# Poll for completion
for i in range(40):
    time.sleep(10)
    poll_url = f"https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi?jobid={job_id}&wait=20"
    r = requests.get(poll_url, timeout=30)
    title = re.search(r"<title>(.*?)</title>", r.text)
    title_text = title.group(1) if title else "unknown"
    print(f"Poll {i+1}: title={title_text}, len={len(r.text)}")

    if "failed" in title_text.lower():
        print("JOB FAILED!")
        # Look for error message
        err_idx = r.text.lower().find("error")
        if err_idx > 0:
            print(r.text[max(0, err_idx-50):err_idx+200])
        break

    if "active" in r.text and "launchcheck" in r.text:
        continue

    # Job might be done - try to find output files
    # DTU stores results at a predictable path
    for path_pattern in [
        f"https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp/{job_id}/output.csv",
        f"https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp/{job_id}/results.csv",
        f"https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp/{job_id}/",
    ]:
        try:
            fr = requests.get(path_pattern, timeout=10)
            if fr.status_code == 200 and len(fr.text) > 10:
                print(f"FOUND: {path_pattern}")
                print(fr.text[:500])
                print("---")
        except Exception:
            pass

    # Check if Angular JSON endpoint exists
    json_url = f"https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp/{job_id}/output.json"
    try:
        jr = requests.get(json_url, timeout=10)
        if jr.status_code == 200:
            print(f"JSON FOUND: {json_url}")
            print(jr.text[:500])
    except Exception:
        pass

    # If page is long, results might be embedded
    if len(r.text) > 5000 and "csv_file" in r.text:
        print("Results page detected - looking for data URL...")
        csv_ref = re.search(r'JSON\.csv_file\s*=\s*["\']([^"\']+)', r.text)
        if csv_ref:
            print(f"CSV ref: {csv_ref.group(1)}")

        # Find any URL with the job ID
        all_urls = re.findall(f'["\']([^"\']*{job_id}[^"\']*)["\']', r.text)
        for u in all_urls[:5]:
            print(f"Job URL found: {u}")
        break

print("Done.")
