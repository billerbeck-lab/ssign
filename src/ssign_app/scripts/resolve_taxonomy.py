#!/usr/bin/env python3
"""Resolve organism names to NCBI taxonomy IDs using a local taxdump.

Reads NCBI's `nodes.dmp` and `names.dmp` via the `taxopy` library. The
taxdump is bundled with ssign's base-tier Zenodo deposit and fetched by
`scripts/fetch_databases.sh`. Default location: ~/.ssign/taxdump/.
Override with the SSIGN_TAXDUMP_DIR environment variable.

Examples:
    resolve_organism("Xanthomonas campestris")
    -> {"species": {"name": "Xanthomonas campestris", "taxid": "339"},
        "genus": {"name": "Xanthomonas", "taxid": "338"}}
"""

import json
import logging
import os
import sys
from typing import Any

logger = logging.getLogger(__name__)

DEFAULT_TAXDUMP_DIR = os.path.join(os.path.expanduser("~"), ".ssign", "taxdump")
TAXDUMP_DIR = os.environ.get("SSIGN_TAXDUMP_DIR", DEFAULT_TAXDUMP_DIR)

# A {"name": ..., "taxid": ...} entry for a single rank, or None when no
# match was found. Used for both species-level and genus-level results.
TaxonEntry = dict[str, str] | None
TaxonomyResult = dict[str, TaxonEntry]

# Module-level taxdump handle and resolution cache. _taxdb is set on first
# call to _load_taxdb(); annotated as Any because taxopy is ignored by mypy.
_taxdb: Any = None
_cache: dict[str, TaxonomyResult] = {}


def _load_taxdb():
    """Load the NCBI taxdump on first use (parses ~1.5 GB uncompressed).

    Raises RuntimeError if taxopy isn't installed or the dump files are
    missing — callers should catch and degrade gracefully.
    """
    global _taxdb
    if _taxdb is not None:
        return _taxdb

    try:
        import taxopy
    except ImportError as e:
        raise RuntimeError(
            "taxopy not installed (required for local taxonomy resolution).\n  How to fix:\n    - pip install taxopy"
        ) from e

    nodes_dmp = os.path.join(TAXDUMP_DIR, "nodes.dmp")
    names_dmp = os.path.join(TAXDUMP_DIR, "names.dmp")
    missing = [p for p in (nodes_dmp, names_dmp) if not os.path.exists(p)]
    if missing:
        raise RuntimeError(
            f"NCBI taxdump files not found under {TAXDUMP_DIR}.\n"
            f"  Missing: {missing}\n"
            f"  How to fix:\n"
            f"    - Run: bash scripts/fetch_databases.sh --tier base\n"
            f"    - Or set SSIGN_TAXDUMP_DIR to an existing taxdump directory\n"
            f"    - Manual download: "
            f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz "
            f"(extract nodes.dmp + names.dmp into the target dir)"
        )

    logger.info(f"Loading NCBI taxdump from {TAXDUMP_DIR} (one-time, ~15-30 s, ~2 GB resident RAM)")
    _taxdb = taxopy.TaxDb(nodes_dmp=nodes_dmp, names_dmp=names_dmp, keep_files=True)
    return _taxdb


def _find_species_and_genus(taxids, taxdb):
    """Pick the best species-level match and its genus from a list of taxids."""
    import taxopy

    result: TaxonomyResult = {"species": None, "genus": None}
    if not taxids:
        return result

    # Constructing Taxon objects walks the parent chain — do it once per id.
    # Prefer an exact species/subspecies/strain rank match; fall back to any
    # taxon whose name has 2+ words (genera are single-word).
    taxa = [taxopy.Taxon(int(tid), taxdb) for tid in taxids]
    species_taxon = next(
        (t for t in taxa if t.rank in ("species", "subspecies", "strain")),
        None,
    )
    if species_taxon is None:
        species_taxon = next((t for t in taxa if " " in t.name), None)

    if species_taxon is not None:
        result["species"] = {
            "name": species_taxon.name,
            "taxid": str(species_taxon.taxid),
        }
        genus_taxid = species_taxon.rank_taxid_dictionary.get("genus")
        if genus_taxid is not None:
            genus_taxon = taxopy.Taxon(int(genus_taxid), taxdb)
            result["genus"] = {
                "name": genus_taxon.name,
                "taxid": str(genus_taxid),
            }

    return result


def resolve_organism(organism_name):
    """Resolve an organism name to species and genus taxids.

    Args:
        organism_name: e.g. "Xanthomonas campestris pv. campestris"

    Returns:
        dict with keys:
            species: {"name": ..., "taxid": ...} or None
            genus:   {"name": ..., "taxid": ...} or None

    On missing taxdump or taxopy, returns {"species": None, "genus": None}
    and logs a warning — callers should treat taxonomy as best-effort.
    """
    if not organism_name or not organism_name.strip():
        return {"species": None, "genus": None}

    organism_name = organism_name.strip()
    if organism_name in _cache:
        return _cache[organism_name]

    try:
        taxdb = _load_taxdb()
    except RuntimeError as e:
        logger.warning(f"Taxonomy resolution unavailable: {e}")
        _cache[organism_name] = {"species": None, "genus": None}
        return _cache[organism_name]

    import taxopy

    taxids = taxopy.taxid_from_name(organism_name, taxdb) or []
    result = _find_species_and_genus(taxids, taxdb)

    # Fallback: if we still don't have a genus, try resolving the first
    # word of the organism name directly as a genus.
    if result["genus"] is None:
        genus_name = organism_name.split()[0]
        genus_taxids = taxopy.taxid_from_name(genus_name, taxdb) or []
        for gtid in genus_taxids:
            taxon = taxopy.Taxon(int(gtid), taxdb)
            if taxon.rank == "genus":
                result["genus"] = {"name": taxon.name, "taxid": str(gtid)}
                break

    _cache[organism_name] = result
    return result


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: resolve_taxonomy.py <organism_name>", file=sys.stderr)
        sys.exit(1)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    name = " ".join(sys.argv[1:])
    print(json.dumps(resolve_organism(name), indent=2))
