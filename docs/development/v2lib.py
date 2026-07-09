"""v2 matcher (prototype, tool-weighted) — importable so analysis scripts reuse it."""

import re
from collections import Counter

# Ordered MOST-SPECIFIC -> generic; first match wins per description.
V2_PATTERNS = [
    (
        "Phage/mobile element",
        r"\bphage\b|prophage|\bintegrase\b|transposase|insertion sequence|capsid|tail fib|terminase|mobile element|conjugativ",
    ),
    # translocator/apparatus context BEFORE the toxin rules: "haemolysin secretion protein" is machinery, not a toxin
    (
        "__APPARATUS__",
        r"(ha?emolysin|leukotoxin|filamentous.?ha?emagglutinin).{0,25}(secretion|activation|transport|export)|\bshlb\b|\bfhac\b|\bhecb\b|\btpsb\b|two.?partner.*(transport|secretion protein)|patatin.?like.*secretion",
    ),
    (
        "Peptidoglycan hydrolase",
        r"amidase|muramidase|glucosaminidase|lytic transglycosylase|murein.*hydrolase|peptidoglycan.*(hydrolase|amidase)|n-acetylmuram|\btae\d|\btge\d",
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
    # RTX dropped as a category (it's a T1SS secretion motif spanning proteases/cytolysins/cyclases,
    # not a function): serralysins fall through to Protease, cytolysins/CyaA to Pore-forming.
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
    ("Hypothetical", r"hypothetical|uncharacter|domain.*unknown|\bduf\d"),  # FLOOR
]
_V2 = [(c, re.compile(p, re.IGNORECASE)) for c, p in V2_PATTERNS]
_MACH = re.compile(
    r"\b(vgrg|hcp|paar|tss[a-m]|clpv|icmf|dotu|gsp[c-o]|pul[c-o]|sec[a-fy]|tat[abc]|flg[a-z]|fli[a-z]|flh[a-z]|mot[ab])\b|vir[bd]\d|flagell|translocon",
    re.IGNORECASE,
)
_FLOORS = {"Hypothetical"}
APPARATUS = "Apparatus-associated"
# Structural fold / superfamily names (not functional annotations). A description that
# IS one of these gets its tool weight capped to 1 so a fold name ("Haem peroxidase
# domain superfamily", "beta-helix") can't outvote a real functional call.
_IS_FOLD = re.compile(
    r"superfamily|-like domain|domain[- ]containing protein|\bfold\b|"
    r"beta.?(helix|barrel|propeller|sandwich)|jelly.?roll|tim.?barrel",
    re.IGNORECASE,
)

# tool credibility weight: functional-name tools > structural-homology tools
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


def classify_v2(desc):
    if not desc or not desc.strip():
        return None
    for cat, pat in _V2:
        if pat.search(desc):
            return APPARATUS if cat == "__APPARATUS__" else cat
    if _MACH.search(desc):
        return APPARATUS
    return None


def consensus_v2(tool_descs, weighted=True):
    # dedup identical strings, keeping the highest-weight tool that produced each;
    # then weighted vote — functional categories beat floors/apparatus.
    best = {}
    for tool, desc in tool_descs.items():
        key = re.sub(r"\W+", " ", (desc or "").lower()).strip()
        if not key:
            continue
        w = TOOL_WEIGHT.get(tool, 2) if weighted else 1
        if _IS_FOLD.search(desc or ""):
            w = 1  # structural fold name, not a functional annotation
        if key not in best or w > best[key][0]:
            best[key] = (w, desc)
    # Functional categories AND Apparatus compete by weight; only Hypothetical is a floor.
    scores, floors = Counter(), Counter()
    for w, desc in best.values():
        c = classify_v2(desc)
        if not c:
            continue
        (floors if c in _FLOORS else scores)[c] += w
    if scores:
        return scores.most_common(1)[0][0]
    if floors.get("Hypothetical"):
        return "Hypothetical"
    return "Unclassified"
