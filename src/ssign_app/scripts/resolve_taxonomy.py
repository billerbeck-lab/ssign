#!/usr/bin/env python3
"""Resolve organism names to NCBI taxonomy IDs.

Uses NCBI Entrez esearch API (no API key needed for low-volume queries).
Returns both species-level and genus-level taxids so the user can choose.

Rate-limited to stay under NCBI's 3 req/sec limit.
Results cached in-memory to avoid re-queries on Streamlit reruns.

Examples:
    resolve_organism("Xanthomonas campestris")
    -> {"species": {"name": "Xanthomonas campestris", "taxid": "339"},
        "genus": {"name": "Xanthomonas", "taxid": "338"}}
"""

import json
import logging
import sys
import time
import urllib.parse
import urllib.request

logger = logging.getLogger(__name__)

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Rate limiting — NCBI allows 3 req/sec without API key
_last_ncbi_call = 0.0


def _ncbi_delay():
    """Ensure at least 0.4s between NCBI API calls (< 3 req/sec)."""
    global _last_ncbi_call
    elapsed = time.time() - _last_ncbi_call
    if elapsed < 0.4:
        time.sleep(0.4 - elapsed)
    _last_ncbi_call = time.time()


# In-memory cache to prevent re-queries on Streamlit reruns
_cache = {}


def _esearch_taxonomy(query):
    """Search NCBI Taxonomy for a name. Returns list of taxids."""
    _ncbi_delay()
    params = urllib.parse.urlencode({
        "db": "taxonomy",
        "term": query,
        "retmode": "json",
        "retmax": "5",
    })
    url = f"{ESEARCH_URL}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
            return data.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        logger.warning(
            f"NCBI taxonomy search failed for {query!r}: {e}\n"
            f"  Note: NCBI Entrez has scheduled maintenance windows (usually weekends).\n"
            f"  If this persists, check https://www.ncbi.nlm.nih.gov/Status/ for outages."
        )
        return []


def _esummary_taxonomy(taxids):
    """Get summary info for taxonomy IDs. Returns dict of {taxid: {name, rank, ...}}."""
    if not taxids:
        return {}
    _ncbi_delay()
    params = urllib.parse.urlencode({
        "db": "taxonomy",
        "id": ",".join(taxids),
        "retmode": "json",
    })
    url = f"{ESUMMARY_URL}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
            result = data.get("result", {})
            out = {}
            for tid in taxids:
                if tid in result:
                    info = result[tid]
                    out[tid] = {
                        "taxid": tid,
                        "name": info.get("scientificname", ""),
                        "rank": info.get("rank", ""),
                        "division": info.get("division", ""),
                    }
            return out
    except Exception as e:
        logger.warning(
            f"NCBI taxonomy summary failed: {e}\n"
            f"  Note: NCBI Entrez has scheduled maintenance windows (usually weekends).\n"
            f"  If this persists, check https://www.ncbi.nlm.nih.gov/Status/ for outages."
        )
        return {}


def _get_lineage(taxid):
    """Get lineage taxids for a given taxid using efetch XML."""
    _ncbi_delay()
    params = urllib.parse.urlencode({
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml",
    })
    url = f"{EFETCH_URL}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            import xml.etree.ElementTree as ET
            tree = ET.parse(resp)
            root = tree.getroot()
            lineage_items = []
            for taxon in root.iter("LineageEx"):
                for item in taxon.findall("Taxon"):
                    tid = item.findtext("TaxId", "")
                    name = item.findtext("ScientificName", "")
                    rank = item.findtext("Rank", "")
                    lineage_items.append({"taxid": tid, "name": name, "rank": rank})
            return lineage_items
    except Exception as e:
        logger.warning(
            f"NCBI lineage fetch failed for {taxid}: {e}\n"
            f"  Note: NCBI Entrez has scheduled maintenance windows (usually weekends).\n"
            f"  If this persists, check https://www.ncbi.nlm.nih.gov/Status/ for outages."
        )
        return []


def resolve_organism(organism_name):
    """Resolve an organism name to species and genus taxids.

    Args:
        organism_name: e.g. "Xanthomonas campestris pv. campestris"

    Returns:
        dict with keys:
            species: {"name": ..., "taxid": ...} or None
            genus: {"name": ..., "taxid": ...} or None
    """
    if not organism_name or not organism_name.strip():
        return {"species": None, "genus": None}

    organism_name = organism_name.strip()

    # Check cache first
    if organism_name in _cache:
        return _cache[organism_name]

    result = {"species": None, "genus": None}

    # Search for the full organism name
    taxids = _esearch_taxonomy(f'"{organism_name}"[Scientific Name]')
    if not taxids:
        taxids = _esearch_taxonomy(f'{organism_name}[Scientific Name]')

    if taxids:
        summaries = _esummary_taxonomy(taxids)
        # Find the best species-level match — only accept actual species ranks
        # (NOT "no rank" which can match genus-level entries)
        for tid, info in summaries.items():
            if info["rank"] in ("species", "subspecies", "strain"):
                result["species"] = {"name": info["name"], "taxid": tid}
                break

        # Fallback: if no species-rank match, accept if name contains 2+ words
        # (genus names are single words, species names have 2+)
        if not result["species"] and summaries:
            for tid, info in summaries.items():
                if " " in info["name"]:
                    result["species"] = {"name": info["name"], "taxid": tid}
                    break

    # Find genus from lineage (most reliable method)
    if result["species"]:
        lineage = _get_lineage(result["species"]["taxid"])
        for item in lineage:
            if item["rank"] == "genus":
                result["genus"] = {"name": item["name"], "taxid": item["taxid"]}
                break

    # Fallback genus: extract from first word of organism name
    if not result["genus"]:
        genus_name = organism_name.split()[0]
        genus_taxids = _esearch_taxonomy(
            f'"{genus_name}"[Scientific Name] AND "genus"[Rank]'
        )
        if genus_taxids:
            genus_summaries = _esummary_taxonomy(genus_taxids[:1])
            if genus_summaries:
                info = list(genus_summaries.values())[0]
                result["genus"] = {"name": info["name"], "taxid": info["taxid"]}

    # Cache result
    _cache[organism_name] = result
    return result


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: resolve_taxonomy.py <organism_name>", file=sys.stderr)
        sys.exit(1)
    name = " ".join(sys.argv[1:])
    result = resolve_organism(name)
    print(json.dumps(result, indent=2))
