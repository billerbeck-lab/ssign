#!/usr/bin/env python3
"""Grade the consensus annotation logic + word-matching against real data.

Sources:
  1. Xanthobacter 74-genome panel combined_results.csv (593 substrates) -> in-the-wild behaviour
  2. benchmark annotation_accuracy_sheet.tsv (55 known secreted proteins, tool descs + known function)
  3. benchmark annotation_groundtruth_fill.tsv (17 effectors, known families)
  4. benchmark gold_verified_readable.tsv (338 known secreted proteins)

Findings-only; no code changes. Prints a structured report.
"""

import csv
import os
import sys
from collections import Counter

ROOT = "/home/teo/Desktop/billerbeck_lab/ssign_package"
sys.path.insert(0, os.path.join(ROOT, "src/ssign_app/scripts"))
from annotation_consensus import _COMPILED, CATEGORY_PATTERNS, classify_description, compute_consensus  # noqa: E402
from ssign_lib.functional_vocab import consensus_bucket  # noqa: E402

KNOWN = {name for name, _ in CATEGORY_PATTERNS}
XPANEL = os.path.expanduser("~/Desktop/cx3_runs/batched_RTX6000_20260706_195242_3232160/combined_results.csv")
BENCH = os.path.join(ROOT, "validation_sweeps/benchmark/analysis")
GOLD = os.path.join(ROOT, "validation_sweeps/benchmark/data/dataset/gold_verified_readable.tsv")

# Xantho column -> tool name (mirrors integrate_annotations TOOL_COLUMNS; blastp absent at extended tier)
XCOLS = {
    "Bakta": "product",
    "InterProScan": "interpro_descriptions",
    "EggNOG": "eggnog_description",
    "pLM-BLAST": "ecod_top1_description",
}


def _load(path, delim):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter=delim))


def _missing(v):
    return not v or str(v).strip().lower() in ("", "nan", "none", "-")


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


# ─────────────────────────────────────────────────────────────────────────
# PART 1 — Xantho in-the-wild behaviour (593 substrates)
# ─────────────────────────────────────────────────────────────────────────
hr("PART 1 — Xanthobacter 74-genome panel (593 substrates): behaviour")
xrows = _load(XPANEL, ",")
print(f"rows: {len(xrows)}")

broad_dist = Counter()
tier_dist = Counter()
bucket_dist = Counter()
n_fallback = 0  # broad_annotation not a known category (title-case fallback noise)
n_hypothetical = 0
n_blank = 0
n_multi_tool_multi_cat = 0  # >=1 tool matched >1 category
suppressed = []  # Hypothetical/fallback won but a real functional cat was present in some tool
reproduced = 0
mismatch_repro = 0

for r in xrows:
    tool_descs = {t: str(r.get(c, "")).strip() for t, c in XCOLS.items() if not _missing(r.get(c, ""))}
    res = compute_consensus(tool_descs)
    broad = res["broad_annotation"]
    stored = (r.get("broad_annotation") or "").strip()
    if broad == stored:
        reproduced += 1
    else:
        mismatch_repro += 1
    broad_dist[broad or "(blank)"] += 1
    tier_dist[res["confidence_tier"]] += 1
    bucket_dist[consensus_bucket(broad)] += 1
    if not broad:
        n_blank += 1
    elif broad not in KNOWN and broad not in ("Secretion system", "Flagellar"):
        n_fallback += 1
    if broad == "Hypothetical":
        n_hypothetical += 1
    # multi-cat per tool
    real_cats_present = set()
    for t, d in tool_descs.items():
        cats = classify_description(d)
        if len([c for c in cats if c in KNOWN]) > 1:
            n_multi_tool_multi_cat += 1
        for c in cats:
            if c in KNOWN and c not in ("Hypothetical",):
                real_cats_present.add(c)
    # suppression: a real functional category existed but Hypothetical/fallback/blank won
    if (broad == "Hypothetical" or broad not in KNOWN or not broad) and real_cats_present:
        suppressed.append(
            (
                r.get("source_genome", "")[:22],
                broad,
                sorted(real_cats_present),
                tool_descs.get("Bakta", "")[:48],
                tool_descs.get("EggNOG", "")[:48],
            )
        )

print(f"reproduced stored broad_annotation: {reproduced}/{len(xrows)}  (mismatch {mismatch_repro})")
print("\nbroad_annotation distribution (top 20):")
for cat, n in broad_dist.most_common(20):
    flag = "" if cat in KNOWN or cat in ("(blank)", "Secretion system", "Flagellar") else "  <- FALLBACK noise"
    print(f"  {n:4d}  {cat}{flag}")
print(f"\nconfidence_tier: {dict(tier_dist)}")
print(f"consensus_bucket (figure vocab): {dict(bucket_dist)}")
print(
    f"\nfallback-category (title-case noise) broad calls: {n_fallback}/{len(xrows)} ({100 * n_fallback / len(xrows):.0f}%)"
)
print(f"'Hypothetical' won: {n_hypothetical}/{len(xrows)} ({100 * n_hypothetical / len(xrows):.0f}%)")
print(f"blank broad: {n_blank}/{len(xrows)}")
print(f"rows where >=1 tool matched >1 known category (vote split): {n_multi_tool_multi_cat}")
print(f"\nSUPPRESSION (real functional cat present but Hypothetical/fallback/blank won): {len(suppressed)}")
for g, b, cats, bak, egg in suppressed[:15]:
    print(f"  won='{b}'  had={cats}\n      bakta='{bak}'  eggnog='{egg}'")


# ─────────────────────────────────────────────────────────────────────────
# PART 2 — Accuracy vs known truth (55-row accuracy sheet)
# ─────────────────────────────────────────────────────────────────────────
hr("PART 2 — accuracy vs known function (annotation_accuracy_sheet.tsv, 55 known proteins)")
arows = _load(os.path.join(BENCH, "annotation_accuracy_sheet.tsv"), "\t")
acols = {
    "InterProScan": "ssign_interpro",
    "EggNOG": "ssign_eggnog",
    "pLM-BLAST": "ssign_plmblast_ecod",
    "Bakta": "ssign_pfam",
}
print(f"rows: {len(arows)}")
print(f"\n{'gene':<10}{'ss_type':<8}{'consensus_call':<22}{'known_family / function':<44}")
print("-" * 84)
acc_other = 0
for r in arows:
    tool_descs = {t: str(r.get(c, "")).strip() for t, c in acols.items() if not _missing(r.get(c, ""))}
    res = compute_consensus(tool_descs)
    broad = res["broad_annotation"] or "(blank)"
    known = (r.get("known_family") or r.get("known_uniprot_name") or r.get("known_uniprot_function") or "")[:42]
    if broad not in KNOWN:
        acc_other += 1
    print(f"{(r.get('gene') or '')[:9]:<10}{(r.get('ss_type') or '')[:7]:<8}{broad[:21]:<22}{known:<44}")
print(f"\nconsensus call landed OUTSIDE the known vocabulary (Other/fallback/blank): {acc_other}/{len(arows)}")


# ─────────────────────────────────────────────────────────────────────────
# PART 3 — Vocabulary coverage of known effector families (lit-curated)
# ─────────────────────────────────────────────────────────────────────────
hr("PART 3 — vocabulary coverage of KNOWN effector families (gold + groundtruth_fill)")
fam_rows = _load(os.path.join(BENCH, "annotation_groundtruth_fill.tsv"), "\t")
print("Known effector family/function -> what CATEGORY_PATTERNS assigns (coverage gap = -> Other):")
print(f"\n{'gene':<10}{'ss_type':<8}{'classify(known_family)':<26}{'known_family':<40}")
print("-" * 84)
gap = 0
for r in fam_rows:
    fam = (r.get("known_family") or "") + " " + (r.get("known_function") or "")
    cats = classify_description(fam)
    hit = next((c for c in cats if c in KNOWN), None) or "(no known-cat match -> Other)"
    if hit.startswith("(no"):
        gap += 1
    print(
        f"{(r.get('gene') or '')[:9]:<10}{(r.get('ss_type') or '')[:7]:<8}{hit[:25]:<26}{(r.get('known_family') or '')[:39]:<40}"
    )
print(f"\nknown effector families with NO matching category: {gap}/{len(fam_rows)}")

# gold list families
grows = _load(GOLD, "\t")
by_type_gap = Counter()
by_type_tot = Counter()
for r in grows:
    st = r.get("ss_type", "")
    by_type_tot[st] += 1
    # gold list has gene + subtype; classify the gene name + subtype as a weak proxy
    probe = f"{r.get('gene', '')} {r.get('subtype', '')}"
    cats = [c for c in classify_description(probe) if c in KNOWN]
    if not cats:
        by_type_gap[st] += 1
print("\ngold list (338): gene+subtype probes with no category (expected high — gene names aren't descriptions):")
for st in sorted(by_type_tot):
    print(f"  {st:<8} {by_type_gap[st]:>3}/{by_type_tot[st]:<3} no-match")


# ─────────────────────────────────────────────────────────────────────────
# PART 4 — per-rule usage + regex false-match probes
# ─────────────────────────────────────────────────────────────────────────
hr("PART 4 — rule usage across Xantho descriptions + false-match probes")
rule_hits = Counter()
all_descs = []
for r in xrows:
    for c in XCOLS.values():
        v = r.get(c, "")
        if not _missing(v):
            all_descs.append(str(v))
for d in all_descs:
    for cat, pat in _COMPILED:
        if pat.search(d):
            rule_hits[cat] += 1
print(f"descriptions scanned: {len(all_descs)}")
print("rule hit counts (0 = never fired on this corpus):")
for cat, _ in CATEGORY_PATTERNS:
    print(f"  {rule_hits[cat]:5d}  {cat}")

hr("false-match probes (does the regex catch things it shouldn't?)")
probes = [
    "lipid II flippase MurJ",
    "flippase family protein",
    "response regulator receiver",
    "ATP synthase subunit",
    "substrate-binding protein of ABC transporter",
    "outer membrane protease OmpT",
    "cell division protein FtsI",
    "conflict resolution",
    "amino acid transferase",
    "nucleotide sugar epimerase",
]
for p in probes:
    cats = classify_description(p)
    print(f"  '{p}' -> {cats}")
