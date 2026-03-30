#!/usr/bin/env python3
"""Quick HHpred server health check - single protein, PDB only."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'ssign_app', 'scripts'))

import requests
from run_hhsuite import MPI_BASE_URL, HHPRED_DEFAULTS, STATUS_NAMES

session = requests.Session()
session.get(f'{MPI_BASE_URL}/tools/hhpred', timeout=30)

seq = '>TEST\nMKKTAIAIAVALAGFATVAQAATENGHTGASATTLNNETATSFKVGDTSTRSSSTLGRFKSNQDNGNTTPVSSQTDGTSTQTTGAGKYDYENT'
payload = {**HHPRED_DEFAULTS, 'alignment': seq, 'hhsuitedb': 'mmcif70/pdb70'}
resp = session.post(f'{MPI_BASE_URL}/api/jobs/?toolName=hhpred', json=payload, timeout=30)
data = resp.json()
job_id = data.get('id', '')
print(f'Submitted: {job_id}')

for i in range(30):
    time.sleep(10)
    elapsed = (i + 1) * 10
    resp = session.get(f'{MPI_BASE_URL}/api/jobs/{job_id}', timeout=15)
    data = resp.json()
    # Handle both dict and list response
    if isinstance(data, list):
        data = data[0] if data else {}
    status = data.get('status', 0)
    sname = STATUS_NAMES.get(status, 'unknown')
    print(f'  [{elapsed}s] Status {status} ({sname})')
    if status == 5:
        print('DONE! Server is healthy.')
        break
    if status == 4:
        print('ERROR - job failed')
        break
else:
    print(f'Still stuck at status {status} after 300s - server may be overloaded')
