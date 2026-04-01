#!/usr/bin/env python3
"""Compute annotation consensus across tools.

Classifies tool descriptions into broad functional categories using
keyword matching, then computes:
  - broad_annotation: most supported broad category
  - broad_consensus_annotation: "Category (Tool1, Tool2, ...)"
  - detailed_annotation: all unique categories found, pipe-separated
  - n_tools_agreeing: tools supporting the broad_annotation
  - concordance_ratio: n_agreeing / n_with_hits
  - confidence_tier: High (≥3 tools agree), Medium (2), Low (1), None (0)
  - evidence_keywords: "Category[Tool1,Tool2]; Category2[Tool3]"
"""

import re
from collections import Counter

# Broad functional categories with keyword patterns.
# Order matters — first match wins for ambiguous descriptions.
CATEGORY_PATTERNS = [
    ('Adhesin', r'adhesin|hemagg?lutinin|fimb|pili|pilin|attachment'),
    ('Autotransporter', r'autotransporter|passenger.*domain|barr?el.*domain'),
    ('Protease', r'protease|peptidase|proteinase|endopeptidase|metalloprotease'),
    ('Lipase/Esterase', r'lipase|esterase|phospholipase|acyltransferase'),
    ('Nuclease', r'nuclease|dnase|rnase|endonuclease|exonuclease'),
    ('Glycoside hydrolase', r'glycos[iy][dl]e?\s*hydrolase|cellulase|chitinase|amylase|lysozyme'),
    ('Toxin', r'toxin|hemolysin|cytolysin|leukocidin|colicin'),
    ('Transporter', r'transporter|permease|efflux|channel|porin|substrate.binding|abc.transporter|mfs'),
    ('Secretion system', r'secretion.*system|type.*secretion|t[1-9]ss|vir[bd]'),
    ('Flagellar', r'flagell|flg[a-z]|fli[a-z]|mot[ab]|hook|basal.body'),
    ('Oxidoreductase', r'oxidoreductase|dehydrogenase|oxidase|reductase|peroxidase|catalase'),
    ('Transferase', r'transferase|kinase|methyltransferase|acetyltransferase'),
    ('Chaperone', r'chaperone|foldase|isomerase|usher'),
    ('Binding protein', r'binding.*protein|receptor|lectin'),
    ('Structural', r'structural|outer.*membrane|lipo.*protein|murein|peptidoglycan'),
    ('Regulatory', r'regulat|transcription|repressor|activator|sensor|response.*regulator'),
    ('Hypothetical', r'hypothetical|uncharacterized|domain.*unknown|duf\d'),
]

_COMPILED = [(cat, re.compile(pat, re.IGNORECASE)) for cat, pat in CATEGORY_PATTERNS]


def classify_description(description: str) -> list[str]:
    """Classify a tool description into broad categories.

    Returns list of matching category names (may be multiple).
    """
    if not description or not description.strip():
        return []

    desc = description.strip()
    categories = []
    for cat, pattern in _COMPILED:
        if pattern.search(desc):
            categories.append(cat)

    # If nothing matched, use the first few words as a fallback category
    if not categories:
        words = desc.split()[:3]
        fallback = ' '.join(w for w in words if len(w) > 2).title()
        if fallback and fallback.lower() not in ('nan', 'none', ''):
            categories.append(fallback)

    return categories


def compute_consensus(tool_descriptions: dict[str, str]) -> dict:
    """Compute annotation consensus from multiple tool descriptions.

    Args:
        tool_descriptions: {tool_name: description} for tools with hits.

    Returns:
        dict with consensus fields.
    """
    if not tool_descriptions:
        return {
            'broad_annotation': '',
            'broad_consensus_annotation': '',
            'detailed_annotation': '',
            'detailed_consensus_annotation': '',
            'evidence_keywords': '',
            'n_tools_agreeing': 0,
            'n_tools_with_hits': 0,
            'concordance_ratio': 0.0,
            'confidence_tier': 'None',
        }

    n_tools = len(tool_descriptions)

    # Classify each tool's description
    tool_categories = {}  # tool → list of categories
    all_categories = Counter()  # category → count of tools supporting it
    category_tools = {}  # category → list of tools

    for tool, desc in tool_descriptions.items():
        cats = classify_description(desc)
        tool_categories[tool] = cats
        for cat in cats:
            all_categories[cat] += 1
            category_tools.setdefault(cat, []).append(tool)

    # Find the most supported broad category
    if all_categories:
        broad, n_agreeing = all_categories.most_common(1)[0]
    else:
        broad, n_agreeing = '', 0

    # Build evidence keywords: "Category[Tool1,Tool2]; Category2[Tool3]"
    evidence_parts = []
    for cat, count in all_categories.most_common():
        tools = sorted(category_tools[cat])
        evidence_parts.append(f"{cat}[{','.join(tools)}]")
    evidence = '; '.join(evidence_parts)

    # Detailed: all unique categories found
    all_cats = sorted(set(cat for cats in tool_categories.values() for cat in cats))
    detailed = ' | '.join(all_cats)

    # Concordance
    concordance = n_agreeing / n_tools if n_tools > 0 else 0.0

    # Confidence tier
    if n_agreeing >= 3:
        tier = 'High'
    elif n_agreeing == 2:
        tier = 'Medium'
    elif n_agreeing == 1:
        tier = 'Low'
    else:
        tier = 'None'

    # Consensus annotation with supporting tools
    if broad:
        supporting = sorted(category_tools.get(broad, []))
        consensus = f"{broad} ({', '.join(supporting)})"
    else:
        consensus = ''

    return {
        'broad_annotation': broad,
        'broad_consensus_annotation': consensus,
        'detailed_annotation': detailed,
        'detailed_consensus_annotation': consensus,
        'evidence_keywords': evidence,
        'n_tools_agreeing': n_agreeing,
        'n_tools_with_hits': n_tools,
        'concordance_ratio': round(concordance, 3),
        'confidence_tier': tier,
    }
