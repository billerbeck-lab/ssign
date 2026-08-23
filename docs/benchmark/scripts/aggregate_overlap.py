#!/usr/bin/env python3
"""Combine the per-tool overlap results into one table + a summary.

Emits:
  docs/benchmark/training_data_overlap.csv   one row per benchmark protein
  work/summary.txt                           counts used to write the report
"""

import csv
import os
from collections import Counter, defaultdict

from _paths import (
    BENCHMARK_CSV,
    CLOSE_HOMOLOG,
    DISTANT_HOMOLOG,
    IN_TRAINING,
    NEAR_IDENTICAL_HOMOLOG,
    NO_MEANINGFUL_MATCH,
    OUT_CSV,
    WORK,
)

W = WORK
bench = {r["instance_id"]: r for r in csv.DictReader(open(BENCHMARK_CSV))}

def load(path):
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}")
    return {r["instance_id"]: r for r in csv.DictReader(open(path), delimiter="\t")}

tools = {
    "deeplocpro":       load(f"{W}/overlap_deeplocpro.tsv"),
    "signalp6":         load(f"{W}/overlap_signalp6_train.tsv"),
    "signalp6_complete":load(f"{W}/overlap_signalp6_complete.tsv"),
    "deepsece_train":   load(f"{W}/overlap_deepsece_train.tsv"),
    "deepsece_holdout": load(f"{W}/overlap_deepsece_holdout.tsv"),
}

# MacSyFinder: profile hits on the protein itself
prof = defaultdict(list)
for r in csv.DictReader(open(f"{W}/txsscan_profile_hits.tsv"), delimiter="\t"):
    prof[r["instance_id"]].append((r["txsscan_profile"], r["i_evalue"], r["profile_coverage"]))
org = {r["instance_id"]: r for r in csv.DictReader(open(f"{W}/txsscan_organism_overlap.tsv"), delimiter="\t")}

SEEN = {IN_TRAINING}
HOMOLOG = {NEAR_IDENTICAL_HOMOLOG, CLOSE_HOMOLOG}

# signalp6_complete is the pre-partitioning candidate set. It is emitted to the
# CSV for reference but deliberately excluded from the leakage verdict and the
# summary tables, so SignalP is not counted twice.
TOOL_LABELS = [
    ("deeplocpro", "DeepLocPro"),
    ("signalp6", "SignalP 6.0"),
    ("deepsece_train", "DeepSecE (train)"),
    ("deepsece_holdout", "DeepSecE (hold-out test)"),
]

rows = []
missing = set()
for iid, b in bench.items():
    o = {
        "instance_id": iid, "ss_type": b["ss_type"], "subtype": b["subtype"],
        "gene": b["gene"], "uniprot": b["uniprot"], "organism": b["organism"],
        "reachable_within_3": b["reachable_within_3"], "found_by_ssign": b["found_by_ssign"],
    }
    for tag, d in tools.items():
        # A protein that failed to resolve upstream is absent from the tool TSVs.
        # Report it as a gap rather than crashing the whole aggregation.
        r = d.get(iid)
        if r is None:
            missing.add(iid)
            r = {"verdict": "UNRESOLVED", "best_pct_identity": "",
                 "best_identical_residue_fraction": "", "matched_record": "",
                 "best_hit_header": ""}
        o[f"{tag}_status"] = r["verdict"]
        o[f"{tag}_pct_identity"] = r["best_pct_identity"]
        o[f"{tag}_identical_residue_pct"] = r["best_identical_residue_fraction"]
        o[f"{tag}_matched_record"] = r["matched_record"] if r["verdict"] == IN_TRAINING else r["best_hit_header"]
    # MacSyFinder / TXSScan
    p = prof.get(iid, [])
    o["txsscan_profile_matches_protein"] = "yes" if p else "no"
    o["txsscan_profile"] = ";".join(x[0] for x in p) if p else "-"
    g = org.get(iid, {})
    o["txsscan_system_in_profile_seed_set"] = g.get("S1_reference_same_system", "no")
    o["txsscan_system_in_validation_set"] = g.get("S2_validation_same_system", "no")

    trained_on = [t for t in ("deeplocpro", "signalp6", "deepsece_train")
                  if o[f"{t}_status"] in SEEN]
    o["n_predictors_trained_on_it"] = len(trained_on)
    o["predictors_trained_on_it"] = ";".join(trained_on) or "-"
    o["in_deepsece_holdout_test"] = "yes" if o["deepsece_holdout_status"] in SEEN else "no"

    if trained_on:
        o["substrate_leakage"] = "direct"
    elif any(o[f"{t}_status"] in HOMOLOG for t in ("deeplocpro", "signalp6", "deepsece_train")):
        o["substrate_leakage"] = "homolog_only"
    else:
        o["substrate_leakage"] = "none"

    if p:
        o["machinery_leakage"] = "protein_is_a_modelled_component"
    elif o["txsscan_system_in_profile_seed_set"] == "yes":
        o["machinery_leakage"] = "its_system_seeded_the_profiles"
    elif o["txsscan_system_in_validation_set"] == "yes":
        o["machinery_leakage"] = "its_system_was_a_validation_case"
    else:
        o["machinery_leakage"] = "none"
    rows.append(o)

rows.sort(key=lambda r: r["instance_id"])
cols = list(rows[0].keys())
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
with open(OUT_CSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

# ---------------- summary ----------------
L = []
def P(s=""):
    L.append(s)
    print(s)

P(f"benchmark proteins: {len(rows)}")
P()
P("=== per-tool: benchmark proteins present in the tool's TRAINING data ===")
for tag, label in TOOL_LABELS:
    c = Counter(r[f"{tag}_status"] for r in rows)
    P(f"{label:<26} IN_TRAINING={c[IN_TRAINING]:<3} near-identical={c[NEAR_IDENTICAL_HOMOLOG]:<3} "
      f"close={c[CLOSE_HOMOLOG]:<3} distant={c[DISTANT_HOMOLOG]:<3} none={c[NO_MEANINGFUL_MATCH]}")
P()
P("=== per-tool IN_TRAINING, broken down by SS type ===")
types = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS"]
P(f"{'tool':<26}" + "".join(f"{t:>8}" for t in types) + f"{'total':>8}")
for tag, label in TOOL_LABELS:
    cells = []
    for t in types:
        n = sum(1 for r in rows if r["ss_type"] == t and r[f"{tag}_status"] == IN_TRAINING)
        d = sum(1 for r in rows if r["ss_type"] == t)
        cells.append(f"{n}/{d}")
    tot = sum(1 for r in rows if r[f"{tag}_status"] == IN_TRAINING)
    P(f"{label:<26}" + "".join(f"{c:>8}" for c in cells) + f"{str(tot)+'/85':>8}")
P()
P("=== overall substrate-level leakage (DLP / SignalP / DeepSecE training sets) ===")
c = Counter(r["substrate_leakage"] for r in rows)
for k in ("direct", "homolog_only", "none"):
    P(f"  {k:<14} {c[k]:>3}/85  ({100*c[k]/85:.0f}%)")
P(f"  seen by >=2 predictors: {sum(1 for r in rows if r['n_predictors_trained_on_it'] >= 2)}/85")
P(f"  seen by all 3:          {sum(1 for r in rows if r['n_predictors_trained_on_it'] == 3)}/85")
P()
P("=== MacSyFinder / TXSScan (rule-based; no training set) ===")
c = Counter(r["machinery_leakage"] for r in rows)
for k, lbl in [("protein_is_a_modelled_component", "protein IS a modelled component"),
               ("its_system_seeded_the_profiles", "its system seeded the HMM profiles"),
               ("its_system_was_a_validation_case", "its system was a validation case"),
               ("none", "no documented exposure")]:
    P(f"  {lbl:<38} {c[k]:>3}/85")
P()
P("=== the number that matters: leakage among the proteins ssign actually predicts ===")
found = [r for r in rows if r["found_by_ssign"] == "yes"]
notfound = [r for r in rows if r["found_by_ssign"] != "yes"]
for label, grp in [("predicted by ssign", found), ("not predicted", notfound)]:
    c = Counter(r["substrate_leakage"] for r in grp)
    P(f"  {label:<20} n={len(grp):<3} direct={c['direct']:<3} homolog_only={c['homolog_only']:<3} none={c['none']}")
P()
P("  per SS type, among ssign's predictions:")
P(f"  {'type':<6} {'predicted':>9} {'direct-leak':>12} {'clean':>7}")
for t in types:
    g = [r for r in found if r["ss_type"] == t]
    d = sum(1 for r in g if r["substrate_leakage"] == "direct")
    P(f"  {t:<6} {len(g):>9} {d:>12} {len(g)-d:>7}")

open(f"{W}/summary.txt", "w").write("\n".join(L) + "\n")
print(f"\nwrote {OUT_CSV} and {W}/summary.txt")
