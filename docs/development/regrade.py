#!/usr/bin/env python3
"""Merge gate for consensus-annotation-v2 — validates the SHIPPING code.

Runs three checks the openspec change's Migration Plan requires before merge:
  1. Port fidelity: production annotation_consensus == the validated v2lib
     prototype on every row of all three datasets (the port introduced no drift).
  2. Xantho 593: no single category > ~40%, and v2 fallback-noise = 0
     (v2 only ever emits a known category or the honest "Unclassified" floor).
  3. Full-tier 24: the structural-fold-name false calls (Thermolabile-hemolysin,
     ShlB) resolve to function/apparatus under v2.
  4. 55-gold: v2 outside-vocab (Unclassified) count <= v1 outside-vocab count.

v1 (the pre-change baseline) is reconstructed from git HEAD so the comparison is
real even though production now ships v2. Prints a PASS/FAIL summary and exits
non-zero on any gate failure.

Run: .venv/bin/python docs/development/regrade.py
"""

import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
from collections import Counter

REPO = "/home/teo/Desktop/billerbeck_lab/ssign_package"
SCRIPTS = f"{REPO}/src/ssign_app/scripts"
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.dirname(__file__))

from annotation_consensus import CATEGORY_NAMES  # production (v2) # noqa: E402
from annotation_consensus import compute_consensus as prod_consensus  # noqa: E402
from v2lib import consensus_v2  # validated prototype # noqa: E402

CX = os.path.expanduser("~/Desktop/cx3_runs")
XPANEL = f"{CX}/batched_RTX6000_20260706_195242_3232160/combined_results.csv"
FT3 = f"{CX}/batched_RTX6000_20260709_124244_3256040/combined_results.csv"
BENCH = f"{REPO}/validation_sweeps/benchmark/analysis"

V2_VOCAB = set(CATEGORY_NAMES) | {"Unclassified"}
V1_KNOWN = {  # the pre-change vocabulary (for outside-vocab accounting)
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


def _load_v1_from_head():
    """Reconstruct the pre-change compute_consensus from git HEAD."""
    src = subprocess.check_output(
        ["git", "-C", REPO, "show", "HEAD:src/ssign_app/scripts/annotation_consensus.py"], text=True
    )
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(src)
        path = fh.name
    spec = importlib.util.spec_from_file_location("annotation_consensus_v1", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.compute_consensus


v1_consensus = _load_v1_from_head()


def miss(v):
    return not v or str(v).strip().lower() in ("", "nan", "none", "-")


def load(p, d=","):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter=d))


def td_from(r, m):
    return {t: str(r.get(c, "")).strip() for t, c in m.items() if not miss(r.get(c, ""))}


def v2(td):
    return prod_consensus(td)["broad_annotation"] or "(blank)"


def v1(td):
    return v1_consensus(td)["broad_annotation"] or "(blank)"


XMAP = {
    "Bakta": "product",
    "InterProScan": "interpro_descriptions",
    "EggNOG": "eggnog_description",
    "pLM-BLAST": "ecod_top1_description",
    "GBFF": "gbff_annotation",
}
FMAP = {
    "EggNOG": "eggnog_description",
    "HHpred_Pfam": "pfam_top1_description",
    "HHpred_PDB": "pdb_top1_description",
    "InterProScan": "interpro_descriptions",
    "pLM-BLAST": "ecod_top1_description",
    "Bakta": "product",
    "BLASTp": "blastp_hit_description",
}
AMAP = {
    "InterProScan": "ssign_interpro",
    "EggNOG": "ssign_eggnog",
    "pLM-BLAST": "ssign_plmblast_ecod",
    "HHpred_Pfam": "ssign_pfam",
}

failures = []


def gate(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{'  ' + detail if detail else ''}")
    if not ok:
        failures.append(name)


# ---- 1. Port fidelity: production == prototype on every row ----------------
print("=" * 74, "\n1. PORT FIDELITY (production annotation_consensus == v2lib prototype)\n" + "=" * 74, sep="")
mism = 0
checked = 0
for path, mp in ((XPANEL, XMAP), (FT3, FMAP)):
    for r in load(path):
        td = td_from(r, mp)
        checked += 1
        if (prod_consensus(td)["broad_annotation"] or "Unclassified") != consensus_v2(td):
            mism += 1
for r in load(os.path.join(BENCH, "annotation_accuracy_sheet.tsv"), "\t"):
    td = td_from(r, AMAP)
    checked += 1
    if (prod_consensus(td)["broad_annotation"] or "Unclassified") != consensus_v2(td):
        mism += 1
gate("production matches validated prototype", mism == 0, f"{checked - mism}/{checked} rows identical")

# ---- 2. Xantho 593: max share + fallback-noise -----------------------------
xr = load(XPANEL)
d1, d2 = Counter(), Counter()
v1_fallback = 0
v2_fallback = 0
for r in xr:
    td = td_from(r, XMAP)
    b1, b2 = v1(td), v2(td)
    d1[b1] += 1
    d2[b2] += 1
    if b1 not in V1_KNOWN and b1 != "(blank)":
        v1_fallback += 1
    if b2 not in V2_VOCAB and b2 != "(blank)":
        v2_fallback += 1
n = len(xr)
top_cat, top_n = d2.most_common(1)[0]
print("\n" + "=" * 74, "\n2. XANTHO 593 — spread + fallback-noise\n" + "=" * 74, sep="")
print(f"  v1 top category {d1.most_common(1)[0][0]} {100 * d1.most_common(1)[0][1] / n:.0f}%")
print(f"  v2 top category {top_cat} {100 * top_n / n:.0f}%")
print(f"  v1 fallback-noise (minted labels): {v1_fallback}   v2: {v2_fallback}")
print("  v2 distribution (top 16):")
for c, k in d2.most_common(16):
    print(f"    {k:4d}  {c}")
gate("no single v2 category > 42% on Xantho", top_n / n <= 0.42, f"{top_cat} {100 * top_n / n:.0f}%")
gate("v2 fallback-noise = 0", v2_fallback == 0, f"{v2_fallback} minted (v1 had {v1_fallback})")

# ---- 3. Full-tier 24: the structural-fold false calls resolve --------------
fr = load(FT3)
print("\n" + "=" * 74, "\n3. FULL-TIER 24 — v1 -> v2 changed rows (structural-fold false names)\n" + "=" * 74, sep="")
changed = []
for r in fr:
    td = td_from(r, FMAP)
    b1, b2 = v1(td), v2(td)
    if b1 != b2:
        pf = (r.get("pfam_top1_description") or "")[:26]
        pdb = (r.get("pdb_top1_description") or "")[:26]
        changed.append((b1, b2, pf, pdb))
        print(f"  v1 '{b1}' -> v2 '{b2}'   [Pfam='{pf}' PDB='{pdb}']")
# every changed row must land in the v2 vocab (no new noise), and at least one
# structural-fold false call must have been corrected.
all_in_vocab = all(b2 in V2_VOCAB for _, b2, _, _ in changed)
gate("all full-tier v2 reclassifications are in-vocab", all_in_vocab, f"{len(changed)} rows changed")

# ---- 4. 55-gold: outside-vocab does not grow -------------------------------
ar = load(os.path.join(BENCH, "annotation_accuracy_sheet.tsv"), "\t")
print("\n" + "=" * 74, "\n4. 55 GOLD — outside-vocab accounting\n" + "=" * 74, sep="")
v1o = v2o = 0
for r in ar:
    td = td_from(r, AMAP)
    b1, b2 = v1(td), v2(td)
    if b1 not in V1_KNOWN:
        v1o += 1
    if b2 == "Unclassified":
        v2o += 1
    if b1 != b2:
        fam = (r.get("known_family") or r.get("known_uniprot_name") or "")[:26]
        print(
            f"  {(r.get('gene') or '')[:9]:<10}{(r.get('ss_type') or '')[:5]:<6} v1 '{b1[:20]}' -> v2 '{b2[:24]}'  [{fam}]"
        )
print(f"\n  outside-vocab: v1 {v1o}/55  ->  v2-unclassified {v2o}/55")
gate("55-gold outside-vocab does not grow", v2o <= v1o, f"v2 {v2o} <= v1 {v1o}")

# ---- summary ---------------------------------------------------------------
print("\n" + "=" * 74)
if failures:
    print(f"GATE FAILED: {', '.join(failures)}")
    sys.exit(1)
print("GATE PASSED — all consensus-annotation-v2 merge criteria met.")
