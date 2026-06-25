#!/usr/bin/env python3
"""Phase A (benchmark-final-validation 1.4 + 1.3 backbone): deterministic mechanical verification.

For every ledger effector this resolves the two falsifiable, content-filter-immune axes with no agents:

  1. primary_ref DOI -> does it resolve? (Crossref) + the paper title/journal/year, so a human or a
     blind agent can eyeball topic match.
  2. uniprot accession -> is it live (UniProt REST) and does its gene name match the ledger gene?

This is the spine of the re-verification: it catches dead DOIs and stale/obsolete accessions
mechanically, leaving only the secretion-claim JUDGEMENT to the blind agents (1.3). Both lookups are
cached to JSON so reruns are offline.

Inputs : data/phase2/verification_phase_a/phasea_ledger.tsv
Outputs: data/phase2/verification_phase_a/mechanical_check.tsv
         data/phase2/verification_phase_a/{crossref_cache.json,uniprot_cache.json}
Run    : .venv/bin/python scripts/61_phasea_mechanical_check.py
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VDIR = BENCH / "data" / "phase2" / "verification_phase_a"
LEDGER = VDIR / "phasea_ledger.tsv"
OUT = VDIR / "mechanical_check.tsv"
CR_CACHE = VDIR / "crossref_cache.json"
UP_CACHE = VDIR / "uniprot_cache.json"
UA = "ssign-benchmark/1.0 (mailto:teoreid@gmail.com)"
FIELDS = [
    "row_id",
    "ss_type",
    "gene",
    "uniprot",
    "doi",
    "doi_resolves",
    "paper_title",
    "paper_journal",
    "paper_year",
    "uniprot_live",
    "uniprot_gene",
    "uniprot_organism",
    "gene_matches",
    "notes",
]


def _norm_doi(ref: str) -> str:
    ref = (ref or "").strip()
    for p in ("https://doi.org/", "http://doi.org/", "doi:", "DOI:"):
        if ref.startswith(p):
            ref = ref[len(p) :]
    return ref.strip()


def _get_json(url: str) -> dict | None:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:  # noqa: S310 (fixed Crossref/UniProt hosts)
            return json.load(r)
    except Exception:
        return None


def crossref(doi: str, cache: dict) -> dict:
    if doi in cache:
        return cache[doi]
    d = _get_json(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}")
    m = (d or {}).get("message", {})
    rec = {
        "resolves": bool(d),
        "title": (m.get("title") or [""])[0][:200],
        "journal": (m.get("container-title") or [""])[0][:120],
        "year": str((m.get("issued", {}).get("date-parts") or [[""]])[0][0] or ""),
    }
    cache[doi] = rec
    time.sleep(0.05)  # polite pool
    return rec


def uniprot(acc: str, cache: dict) -> dict:
    if acc in cache:
        return cache[acc]
    url = f"https://rest.uniprot.org/uniprotkb/{acc}?fields=accession,gene_names,organism_name,protein_name&format=json"
    d = _get_json(url)
    genes = []
    if d:
        for g in d.get("genes", []):
            for k in ("geneName", "orderedLocusNames", "orfNames"):
                v = g.get(k)
                if isinstance(v, dict):
                    genes.append(v.get("value", ""))
                elif isinstance(v, list):
                    genes += [x.get("value", "") for x in v]
            genes += [s.get("value", "") for s in g.get("synonyms", [])]
    rec = {
        "live": bool(d),
        "genes": [g for g in genes if g],
        "organism": (d or {}).get("organism", {}).get("scientificName", "")[:80] if d else "",
    }
    cache[acc] = rec
    time.sleep(0.05)
    return rec


def main() -> int:
    rows = read_tsv(LEDGER)
    cr = json.loads(CR_CACHE.read_text()) if CR_CACHE.exists() else {}
    up = json.loads(UP_CACHE.read_text()) if UP_CACHE.exists() else {}
    out = []
    for i, r in enumerate(rows, 1):
        doi = _norm_doi(r["primary_ref"])
        c = crossref(doi, cr) if doi else {"resolves": False, "title": "", "journal": "", "year": ""}
        acc = (r.get("uniprot") or "").strip()
        has_acc = acc and acc != "-"
        u = uniprot(acc, up) if has_acc else {"live": None, "genes": [], "organism": ""}
        gene = (r.get("gene") or "").strip().lower()
        gmatch = "na" if not has_acc else ("yes" if any(gene == g.lower() for g in u["genes"]) else "no")
        out.append(
            {
                "row_id": r["row_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "uniprot": acc,
                "doi": doi,
                "doi_resolves": "yes" if c["resolves"] else "no",
                "paper_title": c["title"],
                "paper_journal": c["journal"],
                "paper_year": c["year"],
                "uniprot_live": "" if u["live"] is None else ("yes" if u["live"] else "no"),
                "uniprot_gene": ";".join(u["genes"][:4]),
                "uniprot_organism": u["organism"],
                "gene_matches": gmatch,
                "notes": "",
            }
        )
        if i % 50 == 0:
            CR_CACHE.write_text(json.dumps(cr))
            UP_CACHE.write_text(json.dumps(up))
            print(f"  ...{i}/{len(rows)}")
    CR_CACHE.write_text(json.dumps(cr))
    UP_CACHE.write_text(json.dumps(up))
    write_tsv(OUT, FIELDS, out)

    dead = [r for r in out if r["doi_resolves"] == "no"]
    stale = [r for r in out if r["uniprot_live"] == "no"]
    gmis = [r for r in out if r["gene_matches"] == "no"]
    print(f"\nwrote {OUT.relative_to(BENCH)} ({len(out)} rows)")
    print(f"  DOI does NOT resolve: {len(dead)}  -> {[r['row_id'] + '/' + r['gene'] for r in dead]}")
    print(f"  UniProt accession dead: {len(stale)} -> {[r['row_id'] + '/' + r['uniprot'] for r in stale]}")
    print(
        f"  gene != UniProt gene: {len(gmis)} (review; may be naming) -> {[r['row_id'] + '/' + r['gene'] for r in gmis][:25]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
