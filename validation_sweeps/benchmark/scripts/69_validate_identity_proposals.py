#!/usr/bin/env python3
"""Phase A second-pass review: deterministically re-check every accession the identity agents proposed.

The agents propose UniProt accessions; this script does NOT trust their reported length/organism. It reads
the reconciliation, pulls each distinct proposed accession straight from the UniProt /accessions endpoint,
and re-derives length-consistency (vs the gold CDS span) and organism-match itself. That makes the
adjudication evidence machine-checked, not agent-asserted -- the same safety net scripts/66 gave the inputs.

Pass/fail per proposed accession:
  - exists      : the accession is live in UniProtKB
  - len_ok      : UniProt length within +/-15% of the gold row's expected_aa (CDS span/3 - 1)
  - organism_ok : the gold organism's genus+species tokens appear in the UniProt organism string

Inputs : gold_review2/sweep1_identity/reconciliation.tsv  (proposed_uniprot column)
         gold_list_final.tsv                              (expected_aa + organism per row)
Outputs: gold_review2/sweep1_identity/proposal_validation.tsv  (+ printed summary)
Run    : .venv/bin/python scripts/69_validate_identity_proposals.py
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
RECON = VD / "gold_review2" / "sweep1_identity" / "reconciliation.tsv"
OUT = VD / "gold_review2" / "sweep1_identity" / "proposal_validation.tsv"
CACHE = VD / "gold_review2" / "sweep1_identity" / "_proposal_cache.json"
API = "https://rest.uniprot.org/uniprotkb/accessions"
FIELDS = "accession,reviewed,length,organism_name"
LEN_LO, LEN_HI = 0.85, 1.15


def fetch(accs: list[str]) -> dict[str, dict]:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    out: dict[str, dict] = {}
    for i in range(0, len(accs), 80):
        chunk = accs[i : i + 80]
        url = f"{API}?{urllib.parse.urlencode({'accessions': ','.join(chunk), 'fields': FIELDS, 'format': 'json'})}"
        with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 (trusted host)
            for e in json.load(r).get("results", []):
                out[e["primaryAccession"]] = {
                    "reviewed": e.get("entryType", ""),
                    "length": e.get("sequence", {}).get("length"),
                    "organism": e.get("organism", {}).get("scientificName", ""),
                }
    CACHE.write_text(json.dumps(out, indent=0))
    return out


def main() -> int:
    if not RECON.exists():
        print(f"no reconciliation yet: {RECON} (run scripts/68 sweep1_identity first)", file=sys.stderr)
        return 1
    gold = {r["instance_id"]: r for r in read_tsv(GOLD)}
    recon = read_tsv(RECON)
    proposals = [(r["row_id"], r["proposed_uniprot"].strip()) for r in recon if r.get("proposed_uniprot", "").strip()]
    accs = sorted({a for _, a in proposals})
    if not accs:
        print("no proposed accessions in the reconciliation -- nothing to validate")
        return 0
    meta = fetch(accs)

    rows, n_fail = [], 0
    for row_id, acc in proposals:
        g = gold.get(row_id, {})
        try:
            exp_aa = round((int(g["stop"]) - int(g["start"]) + 1) / 3) - 1
        except (ValueError, KeyError):
            exp_aa = None
        m = meta.get(acc)
        up_len = m["length"] if m else None
        ratio = round(up_len / exp_aa, 3) if (up_len and exp_aa) else ""
        gtoks = " ".join((g.get("organism", "").split() + ["", ""])[:2]).lower()
        up_org = (m["organism"] if m else "").lower()
        exists = m is not None
        len_ok = isinstance(ratio, float) and LEN_LO <= ratio <= LEN_HI
        organism_ok = bool(gtoks) and all(t in up_org for t in gtoks.split())
        ok = exists and len_ok and organism_ok
        n_fail += not ok
        rows.append(
            {
                "row_id": row_id,
                "proposed_uniprot": acc,
                "exists": "yes" if exists else "no",
                "reviewed": "yes" if m and m["reviewed"] == "UniProtKB reviewed (Swiss-Prot)" else "no",
                "uniprot_len": up_len or "",
                "expected_aa": exp_aa if exp_aa is not None else "",
                "len_ratio": ratio,
                "len_ok": "yes" if len_ok else "no",
                "uniprot_organism": m["organism"] if m else "",
                "organism_ok": "yes" if organism_ok else "no",
                "verdict": "PASS" if ok else "CHECK",
            }
        )

    rows.sort(key=lambda r: (r["verdict"] == "PASS", r["row_id"]))
    write_tsv(OUT, list(rows[0].keys()), rows)
    print(f"proposal_validation.tsv: {len(rows)} proposed accessions, {n_fail} need a manual look (verdict=CHECK)")
    for r in rows:
        print(
            f"  {r['verdict']:5} {r['row_id']:10} {r['proposed_uniprot']:12} "
            f"len {r['uniprot_len']}/{r['expected_aa']} ratio={r['len_ratio']} "
            f"rev={r['reviewed']} org_ok={r['organism_ok']}  {r['uniprot_organism'][:40]}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
