#!/usr/bin/env python3
"""Per benchmark protein: did each tool CALL it secreted, and was it TRAINED on it?

This is the question the aggregate leakage counts cannot answer. Knowing a protein
sits in DeepLocPro's training set only matters if DeepLocPro is what actually
carried ssign's call for it. A call made by a tool that never saw the protein is
independent evidence no matter what the other tools were trained on.

Proteins are bridged to the run by genome COORDINATES, not by sequence or locus
tag: Bakta renames every locus on re-annotation, and the assemblies used are often
a different strain from the UniProt reference (prtA differs by 5.7%), so exact
sequence matching silently loses real hits. The repo's own validated join is
reused (benchmark/scripts/rerun_coords.RerunIndex), which also reconciles the
three RefSeq/INSDC contig aliases and prefers the `rerun_fullasm/` unit for the
four genomes re-run on full assemblies.

The "called" flags re-derive ssign's gates from the emitted evidence columns
(cross_validate_predictions._dlp_flag and friends); the run does not record which
tool fired, so this is a reconstruction, not a log.
"""

import csv
import io
import os
import sys

from _paths import BENCHMARK_CSV, ROOT, WORK

sys.path.insert(0, os.path.join(ROOT, "validation_sweeps", "benchmark", "scripts"))
from rerun_coords import RerunIndex, _emitted_path, _unit_dirs  # noqa: E402

CONF = 0.8               # constants.CONF_THRESHOLD
VALID_SEC = ("SP", "LIPO")   # constants.VALID_SEC_SIGNAL_TYPES: Sec/SPI + Sec/SPII
# TXSScan profiles whose hit means the protein IS the modelled component, so the
# wider extracellular-or-outer-membrane rule applies (T5SS_COMPONENT_RULES).
EXT_OR_OM = {"T5aSS_PF03797", "T5cSS_PF03895"}


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def secreted_rows(path):
    lines = open(path, newline="").read().split("\n")
    start = next((i for i, ln in enumerate(lines) if ln.strip().lower().startswith("# secreted")), None)
    if start is None:
        return []
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == ""), len(lines))
    block = lines[start + 1:end]
    return list(csv.DictReader(io.StringIO("\n".join(block)))) if block else []


evidence = {}
for unit in _unit_dirs():
    p = _emitted_path(unit)
    if p:
        for r in secreted_rows(p):
            evidence[(unit.name, r["locus_tag"])] = r

OVERLAP_CSV = os.path.join(ROOT, "docs", "benchmark", "training_data_overlap.csv")
overlap = {r["instance_id"]: r for r in csv.DictReader(open(OVERLAP_CSV))}
self_component = set()
with open(os.path.join(WORK, "txsscan_profile_hits.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["txsscan_profile"] in EXT_OR_OM:
            self_component.add(r["instance_id"])

idx = RerunIndex()
out = []
for b in csv.DictReader(open(BENCHMARK_CSV)):
    iid = b["instance_id"]
    o = overlap[iid]
    j = idx.emitted_overlap(b["contig"], int(b["start"]), int(b["stop"]), b["strand"])
    row = {
        "instance_id": iid, "ss_type": b["ss_type"], "subtype": b["subtype"],
        "gene": b["gene"], "uniprot": b["uniprot"], "organism": b["organism"],
        "found_by_ssign": b["found_by_ssign"],
        "run_unit": j["unit"] if j else "", "emitted_locus": j["emitted_locus"] if j else "",
        "join_reason": j["reason"] if j else "contig_absent_from_run",
    }
    ev = evidence.get((row["run_unit"], row["emitted_locus"])) if row["emitted_locus"] else None
    if ev:
        wide = iid in self_component
        dlp_val = max(f(ev.get("dlp_extracellular_prob")), f(ev.get("outer_membrane_prob")) if wide else 0.0)
        dlp = dlp_val >= CONF
        dse = (bool((ev.get("dse_ss_type") or "").strip())
               and f(ev.get("dse_max_prob")) >= CONF
               and str(ev.get("dse_type_match", "")).strip().lower() in ("true", "yes", "1"))
        sp = any((ev.get("signalp_prediction") or "").upper().startswith(k) for k in VALID_SEC)
        row.update({
            "deeplocpro_called": "yes" if dlp else "no",
            "deeplocpro_value": round(dlp_val, 3),
            "deepsece_called": "yes" if dse else "no",
            "deepsece_value": round(f(ev.get("dse_max_prob")), 3),
            "signalp_called": "yes" if sp else "no",
            "signalp_prediction": (ev.get("signalp_prediction") or "").strip(),
        })
    else:
        row.update({"deeplocpro_called": "", "deeplocpro_value": "", "deepsece_called": "",
                    "deepsece_value": "", "signalp_called": "", "signalp_prediction": ""})
    row["deeplocpro_trained_on"] = "yes" if o["deeplocpro_status"] == "IN_TRAINING" else "no"
    row["deepsece_trained_on"] = "yes" if o["deepsece_train_status"] == "IN_TRAINING" else "no"
    row["signalp_trained_on"] = "yes" if o["signalp6_status"] == "IN_TRAINING" else "no"
    row["txsscan_models_this_protein"] = o["txsscan_profile_matches_protein"]

    fired = [t for t in ("deeplocpro", "deepsece", "signalp") if row[f"{t}_called"] == "yes"]
    clean = [t for t in fired if row[f"{t}_trained_on"] == "no"]
    row["tools_that_called_it"] = ";".join(fired) or "-"
    row["tools_that_called_it_and_were_clean"] = ";".join(clean) or "-"
    row["has_independent_support"] = "yes" if clean else ("no" if fired else "")
    out.append(row)

cols = list(out[0].keys())
dest = os.path.join(ROOT, "docs", "benchmark", "per_tool_call_vs_training.csv")
with open(dest, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(out)

pred = [r for r in out if r["found_by_ssign"] == "yes"]
joined = [r for r in pred if r["tools_that_called_it"] != "-"]
indep = [r for r in joined if r["has_independent_support"] == "yes"]
print(f"wrote {dest}")
print(f"\npredicted by ssign: {len(pred)}; evidence recovered for {len(joined)}")
print(f"  at least one calling tool never saw it: {len(indep)}")
print(f"  every calling tool was trained on it:   {len(joined) - len(indep)}")
missing = [r["instance_id"] for r in pred if r["tools_that_called_it"] == "-"]
if missing:
    print(f"  no evidence recovered for: {missing}")
