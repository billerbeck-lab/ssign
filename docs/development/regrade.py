#!/usr/bin/env python3
"""v2-weighted re-grade on Xantho 593 + 55 gold + full-tier 24."""

import csv
import os
import sys
from collections import Counter

sys.path.insert(0, "/home/teo/Desktop/billerbeck_lab/ssign_package/src/ssign_app/scripts")
sys.path.insert(0, os.path.dirname(__file__))
from annotation_consensus import compute_consensus as v1c  # noqa
from v2lib import consensus_v2, classify_v2  # noqa

CX = os.path.expanduser("~/Desktop/cx3_runs")
XPANEL = f"{CX}/batched_RTX6000_20260706_195242_3232160/combined_results.csv"
FT3 = f"{CX}/batched_RTX6000_20260709_124244_3256040/combined_results.csv"
BENCH = "/home/teo/Desktop/billerbeck_lab/ssign_package/validation_sweeps/benchmark/analysis"
V1_KNOWN = {
    "Adhesin",
    "Autotransporter",
    "Protease",
    "Lipase/Esterase",
    "Nuclease",
    "Glycoside hydrolase",
    "Toxin",
    "Transporter",
    "Secretion system",
    "Flagellar",
    "Oxidoreductase",
    "Transferase",
    "Chaperone",
    "Binding protein",
    "Structural",
    "Regulatory",
    "Hypothetical",
}


def miss(v):
    return not v or str(v).strip().lower() in ("", "nan", "none", "-")


def load(p, d=","):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter=d))


def td_from(r, m):
    return {t: str(r.get(c, "")).strip() for t, c in m.items() if not miss(r.get(c, ""))}


# ---- Xantho 593 (extended tools + GBFF) ----
XMAP = {
    "Bakta": "product",
    "InterProScan": "interpro_descriptions",
    "EggNOG": "eggnog_description",
    "pLM-BLAST": "ecod_top1_description",
    "GBFF": "gbff_annotation",
}
xr = load(XPANEL)
d1, d2 = Counter(), Counter()
fb = ap = 0
for r in xr:
    td = td_from(r, XMAP)
    b1 = v1c(td)["broad_annotation"] or "(blank)"
    b2 = consensus_v2(td)
    d1[b1] += 1
    d2[b2] += 1
    if b1 not in V1_KNOWN and b1 != "(blank)":
        fb += 1
n = len(xr)
print("=" * 74, "\nXANTHO 593 — v1 vs v2-WEIGHTED\n" + "=" * 74, sep="")
print(
    f"max single-cat share: v1 {d1.most_common(1)[0][0]} {100 * d1.most_common(1)[0][1] / n:.0f}%  ->  v2 {d2.most_common(1)[0][0]} {100 * d2.most_common(1)[0][1] / n:.0f}%"
)
print(f"v1 fallback-noise: {fb}  ->  v2 0")
print("v2 distribution (top 16):")
for c, k in d2.most_common(16):
    print(f"  {k:4d}  {c}")

# ---- full-tier 24: does weighting fix Thermolabile-hemolysin / ShlB? ----
FMAP = {
    "EggNOG": "eggnog_description",
    "HHpred_Pfam": "pfam_top1_description",
    "HHpred_PDB": "pdb_top1_description",
    "InterProScan": "interpro_descriptions",
    "pLM-BLAST": "ecod_top1_description",
    "Bakta": "product",
    "BLASTp": "blastp_hit_description",
}
fr = load(FT3)
print("\n" + "=" * 74, "\nFULL-TIER 24 — v1 -> v2-weighted (the PDB false-name cases)\n" + "=" * 74, sep="")
for r in fr:
    td = td_from(r, FMAP)
    b1 = v1c(td)["broad_annotation"] or "(blank)"
    b2 = consensus_v2(td)
    if b1 != b2:
        pf = r.get("pfam_top1_description", "")
        pdb = r.get("pdb_top1_description", "")
        print(f"  v1 '{b1}' -> v2 '{b2}'   [Pfam='{pf[:26]}' PDB='{pdb[:26]}']")

# ---- 55 gold accuracy (correct tool names: ssign_pfam = HHpred_Pfam) ----
AMAP = {
    "InterProScan": "ssign_interpro",
    "EggNOG": "ssign_eggnog",
    "pLM-BLAST": "ssign_plmblast_ecod",
    "HHpred_Pfam": "ssign_pfam",
}
ar = load(os.path.join(BENCH, "annotation_accuracy_sheet.tsv"), "\t")
print("\n" + "=" * 74, "\n55 GOLD — v1 -> v2-weighted (changed rows only)\n" + "=" * 74, sep="")
v1o = v2o = 0
for r in ar:
    td = td_from(r, AMAP)
    b1 = v1c(td)["broad_annotation"] or "(blank)"
    b2 = consensus_v2(td)
    if b1 not in V1_KNOWN:
        v1o += 1
    if b2 == "Unclassified":
        v2o += 1
    if b1 != b2:
        print(
            f"  {(r.get('gene') or '')[:9]:<10}{(r.get('ss_type') or '')[:5]:<6} v1 '{b1[:20]}' -> v2 '{b2[:24]}'  [{(r.get('known_family') or r.get('known_uniprot_name') or '')[:26]}]"
        )
print(f"\noutside-vocab/unclassified: v1 {v1o}/55 -> v2 {v2o}/55")
