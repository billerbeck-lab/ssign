#!/usr/bin/env python3
"""Test DTU DeepLocPro via file upload."""
import requests
import re
import time

fasta_content = b""">NJFHAN_00001
MPDLLLELFSEEIPARMQAQAAEALRKLVTDKLVERGLIYEGAKAFVTPRRLALSVHGLP
GRQADQKEEKKGPRVGAPEALLLIALLREEGIDPELREDLLPHKASFLLDAGFDPIQRHL
DTLSAQEAELFRQMMQRLGYSPEQLEAMLLQFNRHFPDVLADTRSLFAEMTKEVFWQLVG
EAAKAGQTVTISGDITDDNHDFKRTGYRYGFCTDAWSFDARLRRTFDEACSAGCADMVFS
"""

url = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"

# Submit via FILE UPLOAD instead of textarea
files = {
    "uploadfile": ("input.fasta", fasta_content, "text/plain"),
}
data = {
    "configfile": "/var/www/services/services/DeepLocPro-1.0/webface.cf",
    "fasta": "",
    "group": "negative",
    "format": "long",
}

print("Submitting via file upload...")
resp = requests.post(url, data=data, files=files, timeout=30)
print(f"URL: {resp.url}")

job_match = re.search(r"jobid=([A-F0-9]+)", resp.url)
if not job_match:
    print("No job")
    print(resp.text[:300])
else:
    job_id = job_match.group(1)
    print(f"Job: {job_id}")
    for i in range(30):
        time.sleep(5)
        r = requests.get(f"{url}?ajax=1&jobid={job_id}", timeout=15)
        try:
            d = r.json()
            st = d.get("status")
            rt = d.get("runtime", 0)
            print(f"Poll {i+1}: {st} (rt={rt}s)")
            if st == "finished":
                rr = requests.get(
                    f"https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/tmp/{job_id}/results.json"
                )
                print(rr.text[:500])
                print("SUCCESS!")
                break
            elif st in ("failed", "error"):
                break
        except Exception:
            print(f"Poll {i+1}: wait")
