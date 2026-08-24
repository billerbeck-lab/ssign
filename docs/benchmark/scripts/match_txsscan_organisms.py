#!/usr/bin/env python3
"""Match benchmark organisms against the strains that seeded / validated TXSScan.

MacSyFinder has no training set. The nearest equivalent exposure is Abby et al.
2016 Table S1 (reference systems whose components seeded the in-house HMM
profiles) and Table S2 (strains the finished models were validated on). If a
benchmark protein's genome AND its system type both appear there, ssign's
machinery call for that pair is not independent of how TXSScan was built.

Three wrinkles the naive match gets wrong, all handled here:
  * the supplement is a multi-column PDF, so extraction interleaves the organism
    column with neighbouring text and splits "Genus species" across lines.
    Matching is genus + epithet within a character window, not adjacent tokens.
  * several benchmark organisms carry post-2016 names (Caulobacter vibrioides =
    C. crescentus, Dickeya chrysanthemi = Erwinia chrysanthemi, Agrobacterium
    fabrum = A. tumefaciens, Edwardsiella piscicida = split from E. tarda).
  * an organism can appear under a system type that is irrelevant to the protein
    in hand (Neisseria meningitidis is listed only under T4P, an appendage ssign
    excludes by default). Matching is therefore done per system-type section,
    and same-type hits are reported separately from any-type hits.
"""

import json
import re

from _paths import WORK

W = WORK
WINDOW = 260

SYNONYMS = {
    "Caulobacter vibrioides": ["Caulobacter crescentus"],
    "Dickeya chrysanthemi": ["Erwinia chrysanthemi", "Dickeya dadantii"],
    "Dickeya dadantii": ["Erwinia chrysanthemi", "Dickeya dadantii"],
    "Agrobacterium fabrum": ["Agrobacterium tumefaciens"],
    "Salmonella typhimurium": ["Salmonella enterica"],
    "Photorhabdus laumondii": ["Photorhabdus luminescens"],
    "Edwardsiella piscicida": ["Edwardsiella tarda"],
}
SYS_RE = r"\b(T1SS|T2SS|T3SS|T4SS|T4P|Tad|T5aSS|T5bSS|T5cSS|T6SS|T9SS|Flagellum)\b"


def sections(path):
    """Split a table into (system_type, text) chunks.

    Consecutive markers naming the SAME system are merged. A row's comment text
    often mentions its own system again ("T2SS Legionella T2SS scattered in 5
    loci ... pneumophila"), and cutting there would strand the genus in one chunk
    and the epithet in the next, so the organism is never matched. Merging is
    safe because both chunks carry the same system tag, which is all this uses.
    """
    flat = re.sub(r"\s+", " ", open(path).read())
    marks = [(m.start(), m.group(1)) for m in re.finditer(SYS_RE, flat)]
    merged = [m for i, m in enumerate(marks) if i == 0 or m[1] != marks[i - 1][1]]
    out = []
    for i, (pos, sysname) in enumerate(merged):
        end = merged[i + 1][0] if i + 1 < len(merged) else len(flat)
        out.append((sysname, flat[pos:end]))
    return out


TABLES = {
    "S1_reference": sections(f"{W}/txsscan_tableS1_reference_systems_full.txt"),
    "S2_validation": sections(f"{W}/txsscan_tableS2_validation_strains_full.txt"),
}


def species_of(org):
    org = re.sub(r"\(.*?\)", " ", org)
    toks = [t for t in re.split(r"[\s,]+", org) if t]
    if len(toks) < 2:
        return None
    g, s = toks[0], toks[1]
    if not re.match(r"^[A-Z][a-z]+$", g) or not re.match(r"^[a-z]+$", s):
        return None
    return f"{g} {s}"


def found_in(text, sp):
    genus, epi = sp.split()
    gpat = rf"(?:{re.escape(genus)}|{genus[0]}\.)"
    for gm in re.finditer(gpat, text):
        seg = text[gm.start() : gm.start() + WINDOW]
        if re.search(rf"\b{re.escape(epi)}\b", seg, re.I):
            return text[max(0, gm.start() - 25) : gm.start() + WINDOW].strip()
    return None


def wanted_sections(ss_type, subtype):
    if ss_type == "T5SS":
        # TXSScan v1.1.4 models only T5aSS, T5bSS and T5cSS. T5dSS and T5eSS have
        # no model at all (Abby 2016: "not matched by the T5SS profiles"), so a
        # T5d/T5e protein has no same-system reference entry to overlap with.
        return [subtype] if subtype in ("T5aSS", "T5bSS", "T5cSS") else []
    return [ss_type]


recs = json.load(open(f"{W}/benchmark_uniprot.json"))
rows = []
for r in recs:
    org = r.get("organism") or r.get("organism_csv") or ""
    sp = species_of(org)
    names = ([sp] + SYNONYMS.get(sp, [])) if sp else []
    want = wanted_sections(r["ss_type"], r["subtype"])
    hit = {
        "instance_id": r["instance_id"],
        "ss_type": r["ss_type"],
        "subtype": r["subtype"],
        "gene": r["gene"],
        "organism": org[:60],
        "species": sp or "-",
        "matched_as": "",
    }
    for tag in ("S1_reference", "S2_validation"):
        same, any_, ctx, othersys = "no", "no", "", set()
        for sysname, text in TABLES[tag]:
            for nm in names:
                c = found_in(text, nm)
                if not c:
                    continue
                any_ = "yes"
                if nm != sp:
                    hit["matched_as"] = nm
                if sysname in want:
                    same = "yes"
                    if not ctx:
                        ctx = c[:190]
                else:
                    othersys.add(sysname)
                break
        hit[f"{tag}_same_system"] = same
        hit[f"{tag}_any_system"] = any_
        hit[f"{tag}_other_systems"] = ",".join(sorted(othersys))
        hit[f"{tag}_context"] = ctx
    rows.append(hit)

cols = [
    "instance_id",
    "ss_type",
    "subtype",
    "gene",
    "organism",
    "species",
    "matched_as",
    "S1_reference_same_system",
    "S1_reference_any_system",
    "S1_reference_other_systems",
    "S2_validation_same_system",
    "S2_validation_any_system",
    "S2_validation_other_systems",
    "S1_reference_context",
    "S2_validation_context",
]
with open(f"{W}/txsscan_organism_overlap.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")

n = len(rows)  # denominators follow the list, which loses rows as types are retired
s1s = sum(r["S1_reference_same_system"] == "yes" for r in rows)
s2s = sum(r["S2_validation_same_system"] == "yes" for r in rows)
both = sum(r["S1_reference_same_system"] == "yes" or r["S2_validation_same_system"] == "yes" for r in rows)
anyt = sum(r["S1_reference_any_system"] == "yes" or r["S2_validation_any_system"] == "yes" for r in rows)
print(f"SAME system type + same species, Table S1 (profile seeds):     {s1s}/{n}")
print(f"SAME system type + same species, Table S2 (validation):        {s2s}/{n}")
print(f"SAME system type in either table:                              {both}/{n}")
print(f"species present anywhere in either table (any system type):    {anyt}/{n}\n")
print(f"{'id':<10} {'type':<7} {'species':<38} S1 S2")
for r in rows:
    if r["S1_reference_same_system"] == "yes" or r["S2_validation_same_system"] == "yes":
        print(
            f"{r['instance_id']:<10} {r['subtype'] or r['ss_type']:<7} {r['species']:<38} "
            f"{r['S1_reference_same_system']:<3}{r['S2_validation_same_system']}"
        )
