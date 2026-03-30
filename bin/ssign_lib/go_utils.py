"""Gene Ontology utility library for GO-slim-based hierarchical categorization.

Provides GO DAG loading with local caching, GO-to-slim mapping via goatools
``mapslim``, True Path Rule enforcement when merging GO terms from multiple
annotation tools, and a three-level hierarchical categorization engine
(broad -> specific -> detail).  Includes a keyword-based fallback for
proteins with no GO terms.

Adapted from the original pipeline/lib/go_utils.py for use in the ssign
Nextflow pipeline.  OBO downloads use HTTPS.
"""

import logging
import os
import urllib.request
from typing import Any

import networkx
import obonet
from goatools.base import download_go_basic_obo
from goatools.mapslim import mapslim
from goatools.obo_parser import GODag

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Broad category map: GO slim IDs -> human-readable broad categories
# Covers Molecular Function (MF), Biological Process (BP), and Cellular
# Component (CC) namespaces from goslim_metagenomics.obo.
# ---------------------------------------------------------------------------

BROAD_CATEGORY_MAP: dict[str, str] = {
    # -- Molecular Function slim terms --
    "GO:0003824": "Catalytic",           # catalytic activity
    "GO:0016787": "Catalytic",           # hydrolase activity
    "GO:0016740": "Catalytic",           # transferase activity
    "GO:0016491": "Catalytic",           # oxidoreductase activity
    "GO:0016829": "Catalytic",           # lyase activity
    "GO:0016853": "Catalytic",           # isomerase activity
    "GO:0016874": "Catalytic",           # ligase activity
    "GO:0008233": "Catalytic",           # peptidase activity
    "GO:0004857": "Regulation",          # enzyme inhibitor activity
    "GO:0005215": "Transport",           # transporter activity
    "GO:0022857": "Transport",           # transmembrane transporter activity
    "GO:0046872": "Binding",             # metal ion binding
    "GO:0043167": "Binding",             # ion binding
    "GO:0005515": "Binding",             # protein binding
    "GO:0030246": "Binding",             # carbohydrate binding
    "GO:0003676": "Binding",             # nucleic acid binding
    "GO:0003677": "Binding",             # DNA binding
    "GO:0003723": "Binding",             # RNA binding
    "GO:0000166": "Binding",             # nucleotide binding
    "GO:0005488": "Binding",             # binding
    "GO:0038023": "Signaling",           # signaling receptor activity
    "GO:0004871": "Signaling",           # signal transducer activity
    "GO:0060089": "Signaling",           # molecular transducer activity
    "GO:0003700": "Regulation",          # DNA-binding transcription factor activity
    "GO:0003735": "Structural",          # structural constituent of ribosome
    "GO:0005198": "Structural",          # structural molecule activity
    "GO:0016209": "Stress Response",     # antioxidant activity
    "GO:0030234": "Regulation",          # enzyme regulator activity
    "GO:0008289": "Binding",             # lipid binding
    "GO:0019825": "Binding",             # oxygen binding
    "GO:0140110": "Regulation",          # transcription regulator activity
    # -- Biological Process slim terms --
    "GO:0006810": "Transport",           # transport
    "GO:0055085": "Transport",           # transmembrane transport
    "GO:0007165": "Signaling",           # signal transduction
    "GO:0006508": "Catalytic",           # proteolysis
    "GO:0009403": "Virulence",           # toxin biosynthetic process
    "GO:0009404": "Virulence",           # toxin metabolic process
    "GO:0005975": "Metabolism",          # carbohydrate metabolic process
    "GO:0006629": "Metabolism",          # lipid metabolic process
    "GO:0009058": "Metabolism",          # biosynthetic process
    "GO:0006950": "Stress Response",     # response to stress
    "GO:0009372": "Signaling",           # quorum sensing
    "GO:0006259": "Metabolism",          # DNA metabolic process
    "GO:0006412": "Metabolism",          # translation
    "GO:0006396": "Metabolism",          # RNA processing
    "GO:0006351": "Metabolism",          # DNA-templated transcription
    "GO:0006355": "Regulation",          # regulation of DNA-templated transcription
    "GO:0050896": "Signaling",           # response to stimulus
    "GO:0007154": "Signaling",           # cell communication
    "GO:0044238": "Metabolism",          # primary metabolic process
    "GO:0043170": "Metabolism",          # macromolecule metabolic process
    "GO:0009056": "Metabolism",          # catabolic process
    "GO:0008152": "Metabolism",          # metabolic process
    "GO:0006807": "Metabolism",          # nitrogen compound metabolic process
    "GO:0009607": "Stress Response",     # response to biotic stimulus
    "GO:0006928": "Structural",          # movement of cell or subcellular component
    "GO:0071554": "Structural",          # cell wall organization or biogenesis
    "GO:0030001": "Transport",           # metal ion transport
    "GO:0006811": "Transport",           # ion transport
    "GO:0044419": "Virulence",           # biological process involved in interspecies interaction
    "GO:0009405": "Virulence",           # pathogenesis
    "GO:0051704": "Virulence",           # multi-organism process
    "GO:0007155": "Structural",          # cell adhesion
    "GO:0030031": "Structural",          # cell projection assembly
    "GO:0044764": "Metabolism",          # multi-organism cellular process
    "GO:0019835": "Metabolism",          # cytolysis
    # -- Cellular Component slim terms --
    "GO:0005576": "Extracellular",       # extracellular region
    "GO:0019867": "Membrane-associated", # outer membrane
    "GO:0042597": "Periplasmic",         # periplasmic space
    "GO:0005886": "Membrane-associated", # plasma membrane
    "GO:0016020": "Membrane-associated", # membrane
    "GO:0005737": "Structural",          # cytoplasm
    "GO:0005622": "Structural",          # intracellular anatomical structure
    "GO:0005694": "Structural",          # chromosome
    "GO:0005840": "Structural",          # ribosome
    "GO:0009279": "Membrane-associated", # cell outer membrane
    "GO:0030312": "Structural",          # external encapsulating structure
    "GO:0009274": "Structural",          # peptidoglycan-based cell wall
    "GO:0009288": "Structural",          # bacterial-type flagellum
    "GO:0110165": "Structural",          # cellular anatomical entity
    "GO:0043226": "Structural",          # organelle
    "GO:0032991": "Structural",          # protein-containing complex
    "GO:0005615": "Extracellular",       # extracellular space
}

# ---------------------------------------------------------------------------
# Keyword fallback categories for proteins with no GO terms
# ---------------------------------------------------------------------------

FALLBACK_KEYWORDS: dict[str, list[str]] = {
    "Catalytic": [
        "protease", "peptidase", "hydrolase", "kinase", "transferase",
        "oxidoreductase", "lyase", "ligase", "isomerase", "dehydrogenase",
        "reductase", "synthase", "synthetase",
    ],
    "Transport": [
        "transporter", "permease", "porin", "channel", "pump", "efflux",
        "ABC", "TonB",
    ],
    "Binding": [
        "binding protein", "receptor",
    ],
    "Structural": [
        "flagellin", "pilin", "fimbrial", "pilus",
    ],
    "Signaling": [
        "sensor", "response regulator", "histidine kinase", "GGDEF", "EAL",
    ],
    "Virulence": [
        "effector", "toxin", "hemolysin", "virulence",
    ],
    "Metabolism": [
        "metabolic", "biosynthesis", "catabolic", "degradation",
    ],
}


# ===================================================================
# Public API
# ===================================================================


def load_go_dags(data_dir: str) -> tuple[GODag, GODag]:
    """Load the full GO DAG and the metagenomics slim DAG.

    Downloads OBO files into *data_dir* if they are not already cached.

    Args:
        data_dir: Directory to store/look for OBO files.

    Returns:
        ``(go_dag, slim_dag)`` tuple of :class:`GODag` objects.
    """
    os.makedirs(data_dir, exist_ok=True)

    obo_path = os.path.join(data_dir, "go-basic.obo")
    slim_path = os.path.join(data_dir, "goslim_metagenomics.obo")

    # Download go-basic.obo if missing
    if not os.path.exists(obo_path):
        logger.info("Downloading go-basic.obo to %s ...", obo_path)
        download_go_basic_obo(obo_path)
    else:
        logger.debug("go-basic.obo already cached at %s", obo_path)

    # Download goslim_metagenomics.obo if missing (HTTPS)
    if not os.path.exists(slim_path):
        url = (
            "https://current.geneontology.org/ontology/subsets/"
            "goslim_metagenomics.obo"
        )
        logger.info("Downloading goslim_metagenomics.obo from %s ...", url)
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as resp, open(slim_path, "wb") as out:
            out.write(resp.read())
    else:
        logger.debug("goslim_metagenomics.obo already cached at %s", slim_path)

    go_dag = GODag(obo_path)
    slim_dag = GODag(slim_path)

    logger.info("Loaded GO DAG: %d terms", len(go_dag))
    logger.info("Loaded GO slim: %d terms", len(slim_dag))

    return go_dag, slim_dag


def load_go_graph(obo_path: str) -> networkx.MultiDiGraph:
    """Load an OBO file as a networkx graph via obonet.

    Needed for True Path Rule operations.  Edge direction: child -> parent.

    Args:
        obo_path: Path to the ``.obo`` file.

    Returns:
        :class:`networkx.MultiDiGraph` with is_a edges.
    """
    graph = obonet.read_obo(obo_path)
    logger.info("Loaded GO graph via obonet: %d nodes", graph.number_of_nodes())
    return graph


def merge_go_terms_true_path(
    interpro_terms: set[str],
    other_terms: set[str],
    go_graph: networkx.MultiDiGraph,
) -> set[str]:
    """Merge GO terms from multiple sources respecting True Path Rule.

    Takes the union of term sets, filters to terms present in the GO graph,
    then removes redundant ancestor terms — keeping only the most specific.

    In obonet, ``networkx.descendants(graph, term)`` returns ANCESTORS
    (more general terms) because edges go child -> parent.

    Args:
        interpro_terms: GO IDs from InterProScan.
        other_terms: GO IDs from any other source.
        go_graph: networkx graph from :func:`load_go_graph`.

    Returns:
        Set of most-specific GO IDs after redundancy removal.
    """
    all_terms = interpro_terms | other_terms

    # Filter to terms actually present in the DAG
    valid_terms = {t for t in all_terms if t in go_graph}
    skipped = len(all_terms) - len(valid_terms)
    if skipped:
        logger.debug(
            "Skipped %d GO terms not found in graph (out of %d)",
            skipped,
            len(all_terms),
        )

    # Remove redundant ancestors: keep only the most specific terms.
    most_specific = set(valid_terms)
    for term in valid_terms:
        try:
            parents = networkx.descendants(go_graph, term)
        except networkx.NetworkXError:
            continue
        most_specific -= parents

    removed = len(valid_terms) - len(most_specific)
    if removed:
        logger.debug(
            "Removed %d redundant ancestor terms (kept %d most-specific)",
            removed,
            len(most_specific),
        )

    return most_specific


def map_go_to_slim(
    go_terms: list[str],
    go_dag: GODag,
    slim_dag: GODag,
) -> dict[str, set[str]]:
    """Map a list of GO terms to their metagenomics slim categories.

    Args:
        go_terms: List of GO IDs.
        go_dag: Full GO DAG.
        slim_dag: Slim GO DAG.

    Returns:
        Dict with ``"direct_slim_ids"`` and ``"all_slim_ids"`` sets.
    """
    direct_slim_ids: set[str] = set()
    all_slim_ids: set[str] = set()

    for term in go_terms:
        if term not in go_dag:
            logger.warning("GO term %s not found in GO DAG, skipping", term)
            continue
        try:
            direct, all_anc = mapslim(term, go_dag, slim_dag)
            direct_slim_ids |= direct
            all_slim_ids |= all_anc
        except (KeyError, ValueError) as exc:
            logger.warning("mapslim failed for %s: %s", term, exc)

    return {
        "direct_slim_ids": direct_slim_ids,
        "all_slim_ids": all_slim_ids,
    }


def categorize_protein(
    go_terms: list[str],
    go_dag: GODag,
    slim_dag: GODag,
    annotation_text: str = "",
) -> dict[str, str]:
    """Assign a three-level hierarchical category to a protein.

    If *go_terms* are available, uses GO-slim mapping.  Otherwise falls
    back to keyword matching on *annotation_text*.

    Three-level hierarchy:
        * ``func_category_broad``: pipe-delimited sorted broad categories
        * ``func_category_specific``: pipe-delimited sorted slim term names
        * ``func_category_detail``: pipe-delimited sorted original GO IDs

    Args:
        go_terms: List of GO IDs for the protein.
        go_dag: Full GO DAG.
        slim_dag: Slim GO DAG.
        annotation_text: Free-text annotation used as fallback.

    Returns:
        Dict with keys ``func_category_broad``, ``func_category_specific``,
        ``func_category_detail``, and ``categorization_source``.
    """
    valid_terms = [t.strip() for t in go_terms if t and t.strip()]

    if not valid_terms:
        return _keyword_fallback(annotation_text)

    slim_result = map_go_to_slim(valid_terms, go_dag, slim_dag)
    direct_slim_ids = slim_result["direct_slim_ids"]

    if not direct_slim_ids:
        return _keyword_fallback(annotation_text)

    broad_categories: set[str] = set()
    specific_names: set[str] = set()

    for slim_id in direct_slim_ids:
        broad_cat = BROAD_CATEGORY_MAP.get(slim_id, "Other")
        broad_categories.add(broad_cat)

        if slim_id in go_dag:
            specific_names.add(go_dag[slim_id].name)
        elif slim_id in slim_dag:
            specific_names.add(slim_dag[slim_id].name)

    return {
        "func_category_broad": "|".join(sorted(broad_categories)),
        "func_category_specific": "|".join(sorted(specific_names)),
        "func_category_detail": "|".join(sorted(valid_terms)),
        "categorization_source": "go_slim",
    }


def _keyword_fallback(annotation_text: str) -> dict[str, str]:
    """Categorize a protein by keyword matching on annotation text."""
    if not annotation_text or not annotation_text.strip():
        return {
            "func_category_broad": "Unknown",
            "func_category_specific": "",
            "func_category_detail": "",
            "categorization_source": "keyword_fallback",
        }

    text_lower = annotation_text.lower()
    matched_categories: set[str] = set()

    for category, keywords in FALLBACK_KEYWORDS.items():
        for kw in keywords:
            if kw.lower() in text_lower:
                matched_categories.add(category)
                break

    if not matched_categories:
        return {
            "func_category_broad": "Unknown",
            "func_category_specific": "",
            "func_category_detail": "",
            "categorization_source": "keyword_fallback",
        }

    return {
        "func_category_broad": "|".join(sorted(matched_categories)),
        "func_category_specific": "",
        "func_category_detail": "",
        "categorization_source": "keyword_fallback",
    }
