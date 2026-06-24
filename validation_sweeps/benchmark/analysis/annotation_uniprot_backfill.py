#!/usr/bin/env python3
"""Backfill UniProt ground-truth for the accession-less effectors in the
annotation-accuracy sheet, so they can be graded against ssign's calls.

Of the 51 emitted_secreted effectors, ~17 carry no UniProt accession in the
corpus (uniprot = '-' or blank), so `annotation_accuracy_sheet.py` never
fetched a UniProt name/family/function for them. Many of these DO have a
UniProt entry that simply was not linked. For each, we:

  1. Inherit: if another panel row with the SAME locus_tag already has an
     accession (e.g. the PAO1 INSDC/RefSeq duplicate PA1510), reuse it.
  2. Search UniProt by gene token(s) + organism, reviewed (Swiss-Prot) first,
     then by locus tag, then gene-only.

We write a CANDIDATES TSV for manual verification. We do NOT overwrite the
graded sheet: gene+organism matching can be wrong, and the corpus is hand-
audited, so a human confirms each accession before it is trusted.

    .venv/bin/python annotation_uniprot_backfill.py

Network: queries https://rest.uniprot.org/uniprotkb/search once per distinct
query, cached to disk so reruns are offline.
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SHEET = os.path.join(HERE, "annotation_accuracy_sheet.tsv")
OUT = os.path.join(HERE, "annotation_uniprot_backfill_candidates.tsv")
CACHE = os.path.join(HERE, "annotation_uniprot_search_cache.json")

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,reviewed,protein_name,gene_names,protein_families,cc_function,organism_name"
OUT_FIELDS = [
    "gene",
    "organism",
    "ss_type",
    "ssign_locus",
    "method",
    "cand_accession",
    "reviewed",
    "n_hits",
    "cand_name",
    "cand_families",
    "cand_function",
    "other_accessions",
    "query",
]


def simplify_organism(org: str) -> str:
    """Genus + species, dropping strain/serovar noise UniProt won't match on."""
    org = re.sub(r"\(.*?\)", "", org).strip()
    parts = org.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else org


def gene_tokens(gene: str) -> list[str]:
    """Candidate real-gene tokens from a possibly-synthetic corpus label.

    'TplE_alias_Tle4' -> ['TplE', 'Tle4']; 'Tae_Atu4347' -> ['Tae', 'Atu4347'].
    Drops the 'alias' joiner and de-dups while preserving order.
    """
    raw = [gene] + [p for p in gene.split("_") if p and p.lower() != "alias"]
    seen, out = set(), []
    for t in raw:
        if t and t.lower() not in seen:
            seen.add(t.lower())
            out.append(t)
    return out


def _load_cache() -> dict:
    return json.load(open(CACHE)) if os.path.exists(CACHE) else {}


def uniprot_search(query: str, cache: dict) -> list[dict]:
    if query in cache:
        return cache[query]
    params = urllib.parse.urlencode({"query": query, "fields": FIELDS, "format": "tsv", "size": "15"})
    url = f"{SEARCH}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310 (fixed UniProt host)
            lines = resp.read().decode().splitlines()
    except Exception as exc:  # transient failure: return empty but do NOT cache, so a rerun retries
        print(f"  ! query failed ({query!r}): {exc}", file=sys.stderr)
        return []
    hits = list(csv.DictReader(lines, delimiter="\t"))  # a genuine 0-hit result IS cached
    cache[query] = hits
    json.dump(cache, open(CACHE, "w"), indent=0)
    return hits


def rank(hits: list[dict], org_simple: str, gene_toks: list[str]) -> list[dict]:
    """Reviewed first, then organism match, then exact gene match, then has-function."""
    tok_l = {t.lower() for t in gene_toks}

    def key(h):
        reviewed = h.get("Reviewed", "") == "reviewed"
        org_ok = org_simple.lower() in (h.get("Organism", "") or "").lower()
        genes = (h.get("Gene Names", "") or "").lower().split()
        gene_ok = bool(tok_l & set(genes))
        has_fn = bool((h.get("Function [CC]", "") or "").strip())
        return (reviewed, org_ok, gene_ok, has_fn)

    return sorted(hits, key=key, reverse=True)


def main():
    rows = list(csv.DictReader(open(SHEET), delimiter="\t"))
    have_acc = {r["ssign_locus"]: r["uniprot"] for r in rows if r["uniprot"].strip() and r["uniprot"] != "-"}
    missing = [r for r in rows if not (r["uniprot"].strip() and r["uniprot"] != "-")]
    print(f"{len(missing)}/{len(rows)} effectors are accession-less\n")

    cache = _load_cache()
    out_rows, n_recovered, n_reviewed, n_none, n_lowconf = [], 0, 0, 0, 0
    for r in missing:
        gene, org = r.get("gene", ""), r.get("organism", "")
        locus = r.get("ssign_locus", "")
        org_simple = simplify_organism(org)
        toks = gene_tokens(gene)
        rec = {
            "gene": gene,
            "organism": org,
            "ss_type": r.get("ss_type", ""),
            "ssign_locus": locus,
            "method": "none",
            "cand_accession": "",
            "reviewed": "",
            "n_hits": 0,
            "cand_name": "",
            "cand_families": "",
            "cand_function": "",
            "other_accessions": "",
            "query": "",
        }

        # 1. inherit from a same-locus row that already has an accession
        if locus in have_acc:
            rec.update(method="internal-locus", cand_accession=have_acc[locus])
            out_rows.append(rec)
            n_recovered += 1
            print(f"  {gene or locus:16s} -> {have_acc[locus]} (same-locus inherit)")
            continue

        # 2. locus+organism (most specific), then gene+organism. Always constrained
        # to Bacteria + the source organism: an org-less gene search matches a
        # eukaryote sharing the symbol (Tle1->human TLE1, Tde1->human SERINC3), so
        # we never query without the organism.
        queries = []
        if org_simple:
            tax = "(taxonomy_name:Bacteria)"
            if locus:
                queries.append(f'({locus}) AND (organism_name:"{org_simple}") AND {tax}')
            queries += [f'(gene:{t}) AND (organism_name:"{org_simple}") AND {tax}' for t in toks]

        hits, used_q = [], ""
        for q in queries:
            hits = uniprot_search(q, cache)
            if hits:
                used_q = q
                break

        if not hits:
            n_none += 1
            out_rows.append(rec)
            print(f"  {gene or locus:16s} -> NONE")
            continue

        ranked = rank(hits, org_simple, toks)
        best = ranked[0]
        reviewed = best.get("Reviewed", "")
        fn = (best.get("Function [CC]", "") or "").replace("FUNCTION: ", "").strip()
        # genus agreement: the best hit's organism must share the source genus,
        # else the gene symbol matched an unrelated bacterium -> flag for review.
        expect_genus = org_simple.split()[0].lower() if org_simple else ""
        genus_ok = expect_genus and expect_genus in (best.get("Organism", "") or "").lower()
        rec.update(
            method="uniprot-search" if genus_ok else "low-confidence-genus-mismatch",
            cand_accession=best.get("Entry", ""),
            reviewed=reviewed,
            n_hits=len(hits),
            cand_name=best.get("Protein names", ""),
            cand_families=best.get("Protein families", ""),
            cand_function=fn[:300] + ("..." if len(fn) > 300 else ""),
            other_accessions=",".join(h.get("Entry", "") for h in ranked[1:6]),
            query=used_q,
        )
        out_rows.append(rec)
        if genus_ok:
            n_recovered += 1
            if reviewed == "reviewed":
                n_reviewed += 1
        else:
            n_lowconf += 1
        flag = ("REVIEWED" if reviewed == "reviewed" else "unreviewed") if genus_ok else "GENUS-MISMATCH"
        print(f"  {gene or locus:16s} -> {best.get('Entry', '')} [{flag}] {best.get('Protein names', '')[:50]}")

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"\nwrote {OUT}")
    print(f"  recovered (genus-confirmed): {n_recovered}/{len(missing)}  (reviewed/Swiss-Prot: {n_reviewed})")
    print(f"  low-confidence (verify):     {n_lowconf}/{len(missing)}")
    print(f"  no hit (manual/lit):         {n_none}/{len(missing)}")


if __name__ == "__main__":
    main()
