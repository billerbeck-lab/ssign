#!/usr/bin/env python3
"""Compute annotation consensus across tools (v2, tool-weighted).

Each substrate gets ONE broad functional category by a weighted vote over its
per-tool annotation descriptions:

  - Tool credibility sets the vote weight. Functional-name tools (BLASTp/Swiss-
    Prot, EggNOG, Bakta, GBFF) = 3, domain tools (HHpred_Pfam, InterProScan) = 2,
    structure-only tools (pLM-BLAST/ECOD, HHpred_PDB) = 1.
  - A description that IS a structural fold / superfamily name is capped to
    weight 1 regardless of source tool, so a fold ("Haem peroxidase domain
    superfamily") can never outvote a real functional call.
  - Each description classifies to a single category by a most-specific-first
    ordered rule set (first match wins); no title-case fallback is minted.
  - Machinery / translocator components are detected by component identity and
    map to ``Apparatus-associated``, which competes in the weighted vote.
  - ``Hypothetical`` is a floor (never beats a functional call). When no tool
    yields a functional, apparatus, or hypothetical call the result is
    ``Unclassified`` (honest, not an invented label).

Output dict keys are unchanged so ``integrate_annotations`` and the CSV schema
stay the same:
  - broad_annotation: winning broad category
  - broad_consensus_annotation: "Category (Tool1, Tool2, ...)"
  - detailed_annotation: distinct specific terms from the tool descriptions
  - detailed_consensus_annotation: same as broad_consensus_annotation
  - n_tools_agreeing: tools whose description classified to the winner
  - n_tools_with_hits: tools with a usable description
  - concordance_ratio: n_agreeing / n_with_hits
  - confidence_tier: from the winner's weighted support (High/Medium/Low/no_hits)
  - evidence_keywords: "Category[Tool1,Tool2]; Category2[Tool3]"

See docs/explanation/design_decisions.md § 4.1 for the rationale.
"""

import re
from collections import Counter

APPARATUS = "Apparatus-associated"
UNCLASSIFIED = "Unclassified"
HYPOTHETICAL = "Hypothetical"  # the sole scoring floor (see _FLOORS)

# Ordered MOST-SPECIFIC -> generic; first match wins per description. The
# "__APPARATUS__" sentinel routes translocator/secretion context to
# `Apparatus-associated` (resolved in classify_description). Taxonomy grounded in
# the T1/T3/T5/T6 effector reviews cited in the audit; RTX is deliberately NOT a
# category (it is a T1SS secretion motif spanning proteases/cytolysins/cyclases):
# serralysins fall through to Protease, cytolysins/CyaA to Pore-forming.
CATEGORY_PATTERNS = [
    (
        "Phage/mobile element",
        r"\bphage\b|prophage|\bintegrase\b|transposase|insertion sequence"
        r"|capsid|tail fib|terminase|mobile element|conjugativ",
    ),
    # translocator/apparatus context BEFORE the toxin rules: "haemolysin secretion protein" is machinery, not a toxin
    (
        "__APPARATUS__",
        r"(ha?emolysin|leukotoxin|filamentous.?ha?emagglutinin).{0,25}(secretion|activation|transport|export)"
        r"|\bshlb\b|\bfhac\b|\bhecb\b|\btpsb\b|two.?partner.*(transport|secretion protein)"
        r"|patatin.?like.*secretion",
    ),
    (
        "Peptidoglycan hydrolase",
        r"amidase|muramidase|glucosaminidase|lytic transglycosylase|murein.*hydrolase"
        r"|peptidoglycan.*(hydrolase|amidase)|n-acetylmuram|\btae\d|\btge\d",
    ),
    ("ADP-ribosyltransferase", r"adp.?ribosyl|mono.?adp|\bmartx\b|nad.*glycohydrolase|\bnadase\b|tox-art"),
    ("Glycosyltransferase", r"glycosyltransferase|glcnac|glucosyltransferase|arginine.?glcnac|\bnleb\b"),
    ("Phosphothreonine lyase", r"phosphothreonine.?lyase|phospho.?lyase|\bospf\b|\bsple\b"),
    ("Ubiquitin-pathway", r"ubiquitin.*ligase|e3 ligase|deubiquitinase|ubiquitin.?specific|\bnel\b"),
    ("GTPase modulator", r"gtpase.?activating|\bgap domain\b|guanine.?nucleotide|\bgef\b|rho.*gtpase"),
    ("Kinase/Phosphatase", r"\bkinase\b|phosphatase|phosphotransferase"),
    ("Acyl/Acetyltransferase", r"acetyltransferase|acyltransferase|n-acyl"),
    ("Beta-lactamase", r"beta.?lactamase|\bpenicillin.?binding|carbapenemase"),
    (
        "Protease/Peptidase",
        r"protease|peptidase|proteinase|endopeptidase|metalloprotease|serralysin|subtilisin|\bspate\b|collagenase|transpeptidase|caspase|ulp1",
    ),
    ("Lipase/Phospholipase", r"lipase|esterase|phospholipase|\btle\d|\bpla[12]\b|patatin|alpha/beta hydrolase"),
    ("Nuclease", r"nuclease|\bdnase\b|\brnase\b|endonuclease|exonuclease|his-me finger|\bntox\d|pd-\(d/e\)xk"),
    (
        "Glycoside hydrolase/CAZy",
        r"glycos[iy][dl]e?\s*hydrolase|cellulase|chitinase|amylase|lysozyme|pectate.?lyase|pectin.?lyase|xylanase|glucanase",
    ),
    ("Oxidoreductase", r"oxidoreductase|dehydrogenase|\boxidase\b|reductase|peroxidase|catalase"),
    (
        "Pore-forming toxin",
        r"\bha?emolysin\b|cytolysin|pore.?forming|leukocidin|aerolysin|\bapx[iv]|leukotoxin|adenylate.?cyclase.?toxin",
    ),
    ("Toxin (other)", r"\btoxin\b|colicin|bacteriocin|\bvasx\b"),
    (
        "Adhesin",
        r"adhesin|ha?emagglutinin|fimb|pilin|\bpili\b|invasin|intimin|trimeric.?autotransporter|two.?partner|filamentous.?ha?em|\bhia\b|\byada\b",
    ),
    (
        "Hemophore/metal uptake",
        r"hemophore|haemophore|heme.?binding|haem.?binding|\bhasa\b|siderophore|tonb.?dependent|ferric.*receptor",
    ),
    ("S-layer", r"s.?layer|surface.?layer|paracrystalline|\brsaa\b"),
    ("Lectin/CBM", r"\blectin\b|carbohydrate.?binding module|\bcbm\d"),
    ("Autotransporter passenger", r"autotransporter|passenger.?domain"),
    ("Chaperone", r"chaperone|foldase|\beag\b|duf4123"),
    ("Regulatory", r"regulat|transcription|repressor|activator|two.?component|response.?regulator|sensor.?kinase"),
    (
        "Transporter/channel",
        r"transporter|permease|efflux|\bchannel\b|porin|abc.?transporter|\bmfs\b|solute.?binding|substrate.?binding",
    ),
    (HYPOTHETICAL, r"hypothetical|uncharacter|domain.*unknown|\bduf\d"),  # FLOOR
]

_COMPILED = [(cat, re.compile(pat, re.IGNORECASE)) for cat, pat in CATEGORY_PATTERNS]

# Machinery detected by COMPONENT IDENTITY (word boundaries), not by naming a
# secretion-system type: an effector annotated "T6SS amidase (Tae)" must NOT be
# hijacked to machinery. VgrG/Hcp/PAAR/Tss*/Gsp*/Vir[BD]/Sec/Tat/flagellar etc.
_MACH = re.compile(
    r"\b(vgrg|hcp|paar|tss[a-m]|clpv|icmf|dotu|gsp[c-o]|pul[c-o]|sec[a-fy]|tat[abc]|flg[a-z]|fli[a-z]|flh[a-z]|mot[ab])\b|vir[bd]\d|flagell|translocon",
    re.IGNORECASE,
)

# Structural fold / superfamily names (not functional annotations). A description
# that IS one of these gets its tool weight capped to 1 so a fold name can't
# outvote a real functional call.
_IS_FOLD = re.compile(
    r"superfamily|-like domain|domain[- ]containing protein|\bfold\b|"
    r"beta.?(helix|barrel|propeller|sandwich)|jelly.?roll|tim.?barrel",
    re.IGNORECASE,
)

_FLOORS = {HYPOTHETICAL}

# Tool credibility weight: functional-name tools > domain tools > structure-only.
# Keys match integrate_annotations.TOOL_HIT_COLUMNS; default 2 for any new tool.
TOOL_WEIGHT = {
    "BLASTp": 3,
    "EggNOG": 3,
    "Bakta": 3,
    "GBFF": 3,
    "HHpred_Pfam": 2,
    "InterProScan": 2,
    "pLM-BLAST": 1,
    "HHpred_PDB": 1,
}
_DEFAULT_WEIGHT = 2

# The labels classify_description can emit (for the figure-side vocabulary). The
# __APPARATUS__ sentinel resolves to APPARATUS; machinery (_MACH) also maps there.
CATEGORY_NAMES = list(
    dict.fromkeys((APPARATUS if cat == "__APPARATUS__" else cat) for cat, _ in CATEGORY_PATTERNS).keys()
)
if APPARATUS not in CATEGORY_NAMES:
    CATEGORY_NAMES.append(APPARATUS)


def _norm(desc: str) -> str:
    """Normalise a description for dedup: lowercase, collapse non-word runs."""
    return re.sub(r"\W+", " ", (desc or "").lower()).strip()


def classify_description(description: str) -> str | None:
    """Classify one tool description into a single broad category.

    Most-specific-first ordered rules, first match wins. Machinery/translocator
    context returns ``Apparatus-associated``. Returns ``None`` when nothing
    matches (no invented title-case fallback).
    """
    if not description or not description.strip():
        return None
    for cat, pattern in _COMPILED:
        if pattern.search(description):
            return APPARATUS if cat == "__APPARATUS__" else cat
    if _MACH.search(description):
        return APPARATUS
    return None


def _tool_weight(tool: str, description: str) -> int:
    """Vote weight for a tool's call: credibility tier, capped to 1 if the
    description is itself a structural fold/superfamily name."""
    if _IS_FOLD.search(description or ""):
        return 1
    return TOOL_WEIGHT.get(tool, _DEFAULT_WEIGHT)


def _confidence_tier(weighted_support: int) -> str:
    """Confidence in the winning functional call from its weighted support.

    High >= 5 (e.g. two Tier-1 tools, or Tier-1 + Tier-2), Medium >= 3 (one
    Tier-1, or two Tier-2), Low otherwise. Floors (Hypothetical/Unclassified)
    call this with 0 -> reported as Low (no confident functional call).
    """
    if weighted_support >= 5:
        return "High"
    if weighted_support >= 3:
        return "Medium"
    return "Low"


def compute_consensus(tool_descriptions: dict[str, str]) -> dict:
    """Compute the tool-weighted annotation consensus for one protein.

    Args:
        tool_descriptions: {tool_name: description} for tools with hits.

    Returns:
        dict with the consensus fields (keys documented in the module header).
    """
    if not tool_descriptions:
        # Sentinel for "no annotation tools produced hits for this protein".
        # "no_hits" (not the string "None") avoids pandas.read_csv coercing the
        # cell to NaN on the master-CSV round-trip.
        return {
            "broad_annotation": "",
            "broad_consensus_annotation": "",
            "detailed_annotation": "",
            "detailed_consensus_annotation": "",
            "evidence_keywords": "",
            "n_tools_agreeing": 0,
            "n_tools_with_hits": 0,
            "concordance_ratio": 0.0,
            "confidence_tier": "no_hits",
        }

    n_tools = len(tool_descriptions)

    # Per-tool category + weight.
    tool_category: dict[str, str | None] = {}
    tool_weight: dict[str, int] = {}
    for tool, desc in tool_descriptions.items():
        tool_category[tool] = classify_description(desc)
        tool_weight[tool] = _tool_weight(tool, desc)

    # Weighted vote, deduping identical description strings (keep the highest
    # weight) so two tools echoing the same call count once. Functional
    # categories AND Apparatus compete by weight; Hypothetical is a floor.
    best: dict[str, tuple[int, str | None]] = {}
    for tool, desc in tool_descriptions.items():
        key = _norm(desc)
        if not key:
            continue
        w = tool_weight[tool]
        if key not in best or w > best[key][0]:
            best[key] = (w, tool_category[tool])

    scores: Counter[str] = Counter()
    floors: Counter[str] = Counter()
    for w, cat in best.values():
        if not cat:
            continue
        (floors if cat in _FLOORS else scores)[cat] += w

    if scores:
        broad = scores.most_common(1)[0][0]
        weighted_support = scores[broad]
    elif floors.get(HYPOTHETICAL):
        broad = HYPOTHETICAL
        weighted_support = 0  # a hypothetical winner is no confident functional call -> Low tier
    else:
        broad = UNCLASSIFIED
        weighted_support = 0

    # Tools whose description classified to each category (raw, all tools).
    category_tools: dict[str, list[str]] = {}
    for tool in tool_descriptions:
        cat = tool_category[tool]
        if cat:
            category_tools.setdefault(cat, []).append(tool)

    supporting = sorted(category_tools.get(broad, []))
    n_agreeing = len(supporting)

    # Evidence: every category with its supporting tools, most-supported first.
    def _score_of(cat: str) -> int:
        return scores.get(cat, floors.get(cat, 0))

    ordered = sorted(category_tools, key=lambda c: (-_score_of(c), c))
    evidence = "; ".join(f"{c}[{','.join(sorted(category_tools[c]))}]" for c in ordered)

    # Detailed: distinct specific terms lifted from each tool description
    # (independent of the category vote).
    specific_terms = set()
    for desc in tool_descriptions.values():
        for part in re.split(r"[;|]", desc):
            term = part.strip()
            if term and len(term) > 3 and term.lower() not in ("hypothetical protein", "uncharacterized protein", ""):
                if len(term) > 60:
                    term = term[:60].rsplit(" ", 1)[0] + "..."
                specific_terms.add(term)
    detailed = " | ".join(sorted(specific_terms)[:15])

    concordance = n_agreeing / n_tools if n_tools > 0 else 0.0
    tier = _confidence_tier(weighted_support)

    # broad is always non-empty here (a category, HYPOTHETICAL, or UNCLASSIFIED);
    # append the supporting tools when there are any, else emit the bare label.
    consensus = f"{broad} ({', '.join(supporting)})" if supporting else broad

    return {
        "broad_annotation": broad,
        "broad_consensus_annotation": consensus,
        "detailed_annotation": detailed,
        "detailed_consensus_annotation": consensus,
        "evidence_keywords": evidence,
        "n_tools_agreeing": n_agreeing,
        "n_tools_with_hits": n_tools,
        "concordance_ratio": round(concordance, 3),
        "confidence_tier": tier,
    }
