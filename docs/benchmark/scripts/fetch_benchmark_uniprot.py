#!/usr/bin/env python3
"""Fetch sequences and accessions for every protein in the ssign benchmarking list.

Writes:
  work/benchmark_uniprot.json   per-protein record (sequence, secondary ACs, xrefs)
  work/benchmark_proteins.fasta sequences keyed by instance_id

84 of the 85 resolve from UniProt. YezP (T6SS_13) has no UniProtKB entry at all,
so it is taken from the RefSeq CDS at the coordinates recorded in the
benchmarking list; that fallback is handled here rather than as a manual patch,
so the whole audit reproduces from this one script.
"""

import csv
import json
import os
import re
import time
import urllib.error
import urllib.request

from _paths import BENCHMARK_CSV, WORK

UA = {"User-Agent": "ssign-training-overlap-audit"}
# 400 bad accession, 404 unknown, 410 deleted: all mean "no entry", not "retry".
NOT_FOUND_CODES = {400, 404, 410}


def get(url, tries=4, parse_json=True):
    """GET with retries on transient failures only."""
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                body = r.read().decode()
                return (json.loads(body) if parse_json else body), r.geturl()
        except urllib.error.HTTPError as e:
            if e.code in NOT_FOUND_CODES:
                return None, url
            last = e
            if e.code == 429:
                time.sleep(float(e.headers.get("Retry-After", 5)))
                continue
        except Exception as e:  # network hiccup, timeout
            last = e
        time.sleep(2 * (i + 1))
    raise RuntimeError(f"giving up on {url}: {last}")


def uniprot_record(acc):
    # The full JSON entry is fetched rather than a ?fields= projection because
    # secondaryAccessions is not exposed by the field-selected endpoints.
    j, final_url = get(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
    if not j:
        return None
    rec = {
        "found": True,
        "primaryAccession": j.get("primaryAccession"),
        "secondaryAccessions": j.get("secondaryAccessions", []),
        "uniProtkbId": j.get("uniProtkbId"),
        "entryType": j.get("entryType"),
        "sequence": j.get("sequence", {}).get("value"),
        "seq_len": j.get("sequence", {}).get("length"),
        "seq_md5": j.get("sequence", {}).get("md5"),
        "organism": j.get("organism", {}).get("scientificName"),
        "taxonId": j.get("organism", {}).get("taxonId"),
        "lineage": j.get("organism", {}).get("lineage", []),
        "proteinExistence": j.get("proteinExistence"),
        "entryAudit": j.get("entryAudit", {}),
    }
    # A merged accession silently redirects to the surviving entry, which would
    # swap in a different protein's sequence without any signal. Record it.
    if rec["primaryAccession"] != acc:
        rec["merged_from"] = acc
        print(f"    ! {acc} redirects to {rec['primaryAccession']} ({final_url})")
    xrefs = {}
    embl_prot = []
    for x in j.get("uniProtKBCrossReferences", []):
        xrefs.setdefault(x["database"], []).append(x["id"])
        if x["database"] == "EMBL":
            for p in x.get("properties", []):
                if p.get("key") == "ProteinId" and p.get("value") not in ("-", None):
                    embl_prot.append(p["value"])
    rec["xref_refseq"] = xrefs.get("RefSeq", [])
    rec["xref_embl"] = xrefs.get("EMBL", [])
    rec["xref_pdb"] = xrefs.get("PDB", [])
    rec["embl_protein_ids"] = embl_prot
    return rec


def refseq_cds(contig, start, stop, strand):
    """Translated CDS from a RefSeq region, for proteins with no UniProt entry."""
    gb, _ = get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
        f"&id={contig}&seq_start={start}&seq_stop={stop}"
        f"&strand={2 if str(strand).startswith('-') else 1}&rettype=gb&retmode=text",
        parse_json=False,
    )
    if not gb:
        return None
    m = re.search(r'/translation="([^"]+)"', gb, re.S)
    if not m:
        return None
    pid = re.search(r'/protein_id="([^"]+)"', gb)
    seq = "".join(m.group(1).split())
    return {
        "found": True,
        "primaryAccession": pid.group(1) if pid else f"{contig}:{start}-{stop}",
        "secondaryAccessions": [],
        "uniProtkbId": None,
        "entryType": "RefSeq CDS (no UniProtKB entry)",
        "sequence": seq,
        "seq_len": len(seq),
        "xref_refseq": [pid.group(1)] if pid else [],
        "xref_embl": [],
        "embl_protein_ids": [],
        "source_note": f"translated from {contig}:{start}-{stop}({strand}); not in UniProt",
    }


rows = list(csv.DictReader(open(BENCHMARK_CSV)))
out = []
for n, r in enumerate(rows, 1):
    acc = r["uniprot"].strip()
    rec = {
        "instance_id": r["instance_id"],
        "ss_type": r["ss_type"],
        "subtype": r["subtype"],
        "gene": r["gene"],
        "organism_csv": r["organism"],
        "uniprot_csv": acc,
        "found": False,
    }
    got = uniprot_record(acc) if acc and acc != "-" else None
    if got is None and r.get("contig") and r.get("start"):
        got = refseq_cds(r["contig"], r["start"], r["stop"], r["strand"])
        if got:
            got.setdefault("organism", r["organism"])
    if got:
        rec.update(got)
    out.append(rec)
    print(
        f"[{n:>3}/{len(rows)}] {rec['instance_id']:<10} {acc:<14} "
        f"{'OK len=%s' % rec.get('seq_len') if rec['found'] else 'NOT FOUND'}",
        flush=True,
    )
    time.sleep(0.15)

json.dump(out, open(os.path.join(WORK, "benchmark_uniprot.json"), "w"), indent=1)
with open(os.path.join(WORK, "benchmark_proteins.fasta"), "w") as fh:
    for r in out:
        if r.get("sequence"):
            fh.write(f">{r['instance_id']}|{r['primaryAccession']}|{r['gene']}|{r['ss_type']}\n")
            for i in range(0, len(r["sequence"]), 60):
                fh.write(r["sequence"][i : i + 60] + "\n")

n_ok = sum(1 for r in out if r["found"])
missing = [r["instance_id"] for r in out if not r["found"]]
print(f"\nresolved {n_ok}/{len(out)}" + (f"; missing: {missing}" if missing else ""))
