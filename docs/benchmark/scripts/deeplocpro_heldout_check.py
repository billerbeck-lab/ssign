#!/usr/bin/env python3
"""Would DeepLocPro still call these proteins secreted if it had NOT been trained on them?

DeepLocPro publishes held-out predictions from its nested cross-validation: each
of the 5 folds is a test fold 4 times, so every training protein has 4 predictions
made by models that never saw it. That is an unbiased read of what DeepLocPro
would say about a protein it had not memorised.

ssign flags a protein as DLP-secreted when extracellular_prob >= 0.8
(constants.CONF_THRESHOLD, applied in cross_validate_predictions._dlp_flag).
The wider T5 rule in T5SS_COMPONENT_RULES is keyed on (subtype, *TXSScan
component name*), not on the benchmark subtype, so it applies only to proteins
TXSScan actually models as a component: the T5aSS/T5cSS autotransporters, where
passenger and translocator are one chain. The T5bSS entries in this list are
TpsA passengers, which TXSScan does not model (its T5bSS profile is the TpsB
pore), so ssign gates them with the ordinary extracellular rule and so does this
script. Getting that wrong flips cdrA, hxuA and shlA.
"""

import csv
import glob
import json
import os

from _paths import RAW, WORK

CV = os.path.join(RAW, "deeplocpro/github_DeepLocPro_model/data/cross_validation_held_out_predictions")
# "CYtoplasmicMembrane" is DeepLocPro's own capitalisation in its released files.
# Do not "fix" it; the CSV headers and FASTA labels both spell it this way.
CLASSES = ["Cellwall", "Extracellular", "Cytoplasmic", "CYtoplasmicMembrane", "OuterMembrane", "Periplasmic"]
THRESH = 0.8
# TXSScan profiles whose hit means the benchmark protein IS the modelled
# component, so ssign's wider extracellular-or-outer-membrane rule applies.
EXT_OR_OM_PROFILES = {"T5aSS_PF03797", "T5cSS_PF03895"}

preds = {}
for f in sorted(glob.glob(os.path.join(CV, "*.csv"))):
    with open(f) as fh:
        for row in csv.DictReader(fh):
            acc = row[""] if "" in row else next(iter(row.values()))
            preds.setdefault(acc, []).append({c: float(row[c]) for c in CLASSES})
print(
    f"held-out predictions loaded for {len(preds)} training proteins "
    f"({sum(len(v) for v in preds.values())} rows across 20 CV models)"
)

label = {}
for line in open(os.path.join(RAW, "deeplocpro/full_dataset.fasta")):
    if line.startswith(">"):
        p = line[1:].strip().split("|")
        label[p[0]] = p[1]

self_component = set()
for r in csv.DictReader(open(os.path.join(WORK, "txsscan_profile_hits.tsv")), delimiter="\t"):
    if r["txsscan_profile"] in EXT_OR_OM_PROFILES:
        self_component.add(r["instance_id"])

recs = json.load(open(os.path.join(WORK, "benchmark_uniprot.json")))
# Restrict the overlap table to proteins still on the list: match_training_set.py's
# output is kept as-is when rows are retired from the benchmark, so counting over
# all of it would report a training-set size for a list that no longer exists.
listed = {r["instance_id"] for r in recs}
over = {
    r["instance_id"]: r
    for r in csv.DictReader(open(os.path.join(WORK, "overlap_deeplocpro.tsv")), delimiter="\t")
    if r["instance_id"] in listed
}

rows = []
for r in recs:
    iid = r["instance_id"]
    o = over.get(iid, {})
    if o.get("verdict") != "IN_TRAINING":
        continue
    # A protein can be in the training set under a different accession, matched
    # by identical sequence. Prefer the accession of the record actually matched;
    # without this, 4 of the 38 are silently skipped.
    candidates = [o.get("matched_record", "").split("|")[0], r.get("primaryAccession")]
    acc = next((a for a in candidates if a in preds), None)
    if acc is None:
        print(f"  ! {iid}: in training but no held-out prediction (dropped by DeepLocPro's homology partitioning)")
        continue
    ps = preds[acc]
    mean = {c: sum(p[c] for p in ps) / len(ps) for c in CLASSES}
    ext, om = mean["Extracellular"], mean["OuterMembrane"]
    if iid in self_component:
        rule, val = "ext_or_OM", max(ext, om)
    else:
        rule, val = "ext_only", ext
    rows.append(
        {
            "instance_id": iid,
            "ss_type": r["ss_type"],
            "subtype": r["subtype"],
            "gene": r["gene"],
            "benchmark_uniprot": r.get("primaryAccession"),
            "deeplocpro_accession": acc,
            "n_heldout_models": len(ps),
            "training_label": label.get(acc, "?"),
            "heldout_top_class": max(mean, key=mean.get),
            "heldout_extracellular_prob": round(ext, 3),
            "heldout_outer_membrane_prob": round(om, 3),
            "ssign_rule": rule,
            "heldout_rule_value": round(val, 3),
            "heldout_passes_ssign_threshold": "yes" if val >= THRESH else "no",
        }
    )

if not rows:
    raise SystemExit("no benchmark proteins found in DeepLocPro's training set")

cols = list(rows[0].keys())
with open(os.path.join(WORK, "deeplocpro_heldout.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

n = len(rows)
npass = sum(r["heldout_passes_ssign_threshold"] == "yes" for r in rows)
n_in_training = sum(1 for o in over.values() if o["verdict"] == "IN_TRAINING")
print(f"\nbenchmark proteins in DeepLocPro's training set: {n_in_training}; with held-out predictions: {n}")
print(f"still pass ssign's >= {THRESH} rule when held out: {npass}/{n} ({100 * npass / n:.0f}%)")
print(f"would FLIP to negative without having been trained on: {n - npass}/{n}\n")
print(f"{'id':<10} {'type':<7} {'gene':<20} {'label':<16} {'heldout top':<20} {'rule':<10} {'val':>6}  pass")
for r in sorted(rows, key=lambda x: (x["ss_type"], x["instance_id"])):
    print(
        f"{r['instance_id']:<10} {r['subtype'] or r['ss_type']:<7} {r['gene'][:20]:<20} "
        f"{r['training_label']:<16} {r['heldout_top_class']:<20} {r['ssign_rule']:<10} "
        f"{r['heldout_rule_value']:>6}  {r['heldout_passes_ssign_threshold']}"
    )
