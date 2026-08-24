#!/usr/bin/env python3
"""Which benchmark system types can MacSyFinder + TXSScan detect at all?

The benchmarking list contains secretion-system types that TXSScan v1.1.4 has no
model for, so ssign cannot report them however good the rest of the pipeline is.
Scoring those rows as misses understates recall against a limit that is not
ssign's to fix. This script separates three different reasons a row can be
undetectable, which need different treatment:

  no_model     TXSScan ships no model for that subtype. Nothing can be found.
  not_detected A model exists, but no system of a compatible type was detected
               anywhere in that genome. The system is real but too diverged for
               TXSScan's profiles.
  low_wholeness The system WAS detected, but below ssign's wholeness_threshold
               (0.8), so validate_macsyfinder_systems.py drops it. This is an
               ssign parameter, not a TXSScan limit, and is tunable.

The first two are grounds for removing a row from the list. The third is a real
result about ssign's defaults and should stay.

Model inventory is read from the installed TXSScan package rather than the paper,
so it tracks whatever version is actually pinned. Detection is read from the
benchmark run's own results CSVs (both "Secretion Systems" chunks), joined to the
list by genome coordinates via the same RerunIndex the recall table uses.

`--diagnose` re-runs MacSyFinder over the genomes with no compatible detection to
tell not_detected from low_wholeness. It needs macsyfinder + the hmmsearch shim
on PATH and takes a few minutes per genome, so it is opt-in; without it those
rows are reported as `undetected` without a cause. The re-run reads the input
GenBank's own CDS translations, where the benchmark run re-annotated with Bakta,
so a wholeness figure can differ by a component either way; it is a diagnosis of
which side of the threshold a system falls on, not a reproduction of the run.
"""

import argparse
import collections
import csv
import glob
import io
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

from _paths import BENCHMARK_CSV, ROOT, WORK

sys.path.insert(0, os.path.join(ROOT, "validation_sweeps", "benchmark", "scripts"))
from rerun_coords import RerunIndex, _emitted_path, _unit_dirs  # noqa: E402

TXSSCAN_DIR = os.path.join(os.path.expanduser("~"), ".macsyfinder", "models", "TXSScan")
INPUTS_GB = os.path.join(ROOT, "validation_sweeps", "benchmark", "inputs_gb")
WHOLENESS_THRESHOLD = 0.8  # constants / validate_macsyfinder_systems.py --wholeness-threshold

# Benchmark type (subtype where the list gives one) -> TXSScan models that could
# carry it. The T4SS split is the one that is easy to get backwards: TXSScan
# names its two protein-secretion T4SS models after the MPF class of the
# conjugative ancestor, so pT4SSt is VirB/VirD4 (Ti-plasmid-like, Agrobacterium,
# Brucella, Bartonella) and pT4SSi is Dot/Icm (IncI-like, Legionella, Coxiella).
# Both are checked for every T4SS row; the list does not subtype them.
COMPATIBLE_MODELS = {
    "T1SS": {"T1SS"},
    "T2SS": {"T2SS"},
    "T3SS": {"T3SS"},
    "T4SS": {"pT4SSi", "pT4SSt"},
    "T5aSS": {"T5aSS"},
    "T5bSS": {"T5bSS"},
    "T5cSS": {"T5cSS"},
    "T5dSS": set(),
    "T5eSS": set(),
    "T6SS": {"T6SSi", "T6SSii", "T6SSiii"},
}


def installed_models():
    """Every system TXSScan defines, with the version it came from."""
    if not os.path.isdir(TXSSCAN_DIR):
        sys.exit(f"TXSScan not installed at {TXSSCAN_DIR}\n  Install: macsydata install --user TXSScan==1.1.4")
    vers = "?"
    meta = os.path.join(TXSSCAN_DIR, "metadata.yml")
    if os.path.exists(meta):
        for line in open(meta):
            if line.startswith("vers:"):
                vers = line.split(":", 1)[1].strip()
    models = {}
    for path in glob.glob(os.path.join(TXSSCAN_DIR, "definitions", "**", "*.xml"), recursive=True):
        name = os.path.basename(path)[:-4]
        root = ET.parse(path).getroot()
        models[name] = {
            "kingdom": os.path.basename(os.path.dirname(path)),
            "min_genes_required": root.get("min_genes_required"),
            "min_mandatory_genes_required": root.get("min_mandatory_genes_required"),
            "genes": tuple(g.get("name") for g in root.iter("gene") if g.get("name")),
        }
    return vers, models


def excluded_systems():
    """ssign's default excluded_systems, read from the package it ships."""
    sys.path.insert(0, os.path.join(ROOT, "src", "ssign_app", "scripts"))
    from ssign_lib.constants import DEFAULT_EXCLUDED_SYSTEMS

    return set(DEFAULT_EXCLUDED_SYSTEMS)


def systems_in_results(path):
    """Every system row from both "# Secretion Systems" chunks of a results CSV.

    Chunk 2 holds systems that have a secreted protein nearby and chunk 3 the
    rest; reading only the first understates detection (it hides every system
    that was found but produced no substrate, which is exactly the case here).
    """
    lines = open(path).read().splitlines()
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith("# Secretion Systems"):
            j, buf = i + 1, []
            while j < len(lines) and not lines[j].startswith("#"):
                if lines[j].strip():
                    buf.append(lines[j])
                j += 1
            if buf:
                out += [r for r in csv.DictReader(io.StringIO("\n".join(buf))) if r.get("record_type") == "system"]
            i = j
        else:
            i += 1
    return out


def detected_per_unit():
    types = {}
    for unit in _unit_dirs():
        p = _emitted_path(unit)
        found = systems_in_results(p) if p else []
        types[unit.name] = collections.Counter(s["ss_type"] for s in found)
    return types


def _profile_hits(msf_dir):
    """{profile: [gene positions]} for every profile with a passing HMM hit.

    MacSyFinder writes one `<profile>.res_hmm_extract` per profile, already
    filtered to its i-evalue and coverage thresholds, with the gene's ordinal
    position on the replicon in column 3. This is the evidence needed to tell
    "the profiles never matched" from "they matched but were too scattered to
    form a system", which are different limits with different answers.
    """
    hits = {}
    for f in glob.glob(os.path.join(msf_dir, "hmmer_results", "*.res_hmm_extract")):
        name = os.path.basename(f)[: -len(".res_hmm_extract")]
        pos = [int(ln.split("\t")[2]) for ln in open(f) if not ln.startswith("#") and "\t" in ln]
        if pos:
            hits[name] = sorted(pos)
    return hits


def _rejections(msf_dir, model):
    """Distinct reasons MacSyFinder rejected candidate clusters of one model."""
    p = os.path.join(msf_dir, "rejected_candidates.tsv")
    if not os.path.exists(p):
        return []
    out = []
    for ln in open(p):
        if f"/{model}\t" not in ln:
            continue
        why = ln.rstrip("\n").split("\t")[-1]
        for part in why.split("/"):
            part = part.strip()
            if part and part not in out:
                out.append(part)
    return out


def rerun_macsyfinder(contig, outdir):
    """Re-run MacSyFinder on one genome; return (wholeness, profile hits, dir).

    Reads all_systems.tsv, not best_solution.tsv: the point is to see systems
    ssign never got the chance to reject, including ones the best-solution
    scoring dropped in favour of an overlapping model.
    """
    gb = os.path.join(INPUTS_GB, contig + ".gbff")
    if not os.path.exists(gb):
        return None, {}, None
    faa = os.path.join(outdir, contig + ".faa")
    if not os.path.exists(faa):
        from Bio import SeqIO

        n = 0
        with open(faa, "w") as fh:
            for rec in SeqIO.parse(gb, "genbank"):
                for ft in rec.features:
                    if ft.type == "CDS" and ft.qualifiers.get("translation"):
                        lt = ft.qualifiers.get("locus_tag", [f"cds{n}"])[0]
                        fh.write(f">{lt}\n{ft.qualifiers['translation'][0]}\n")
                        n += 1
        if not n:
            return None, {}, None
    msf = os.path.join(outdir, "msf_" + contig)
    tsv = os.path.join(msf, "all_systems.tsv")
    if not os.path.exists(tsv):
        cmd = [
            "macsyfinder",
            "--sequence-db",
            faa,
            "--db-type",
            "ordered_replicon",
            "--models",
            "TXSScan",
            "all",
            "--out-dir",
            msf,
            "-w",
            "2",
            "--mute",
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(tsv):
            print(f"    ! macsyfinder failed on {contig}: {r.stderr.strip().splitlines()[-1:]}", file=sys.stderr)
            return None, {}, None
    rows = [ln.rstrip("\n") for ln in open(tsv) if ln.strip() and not ln.startswith("#")]
    if len(rows) < 2:
        return {}, _profile_hits(msf), msf
    hdr = rows[0].split("\t")
    best = {}
    for ln in rows[1:]:
        r = dict(zip(hdr, ln.split("\t")))
        if not r.get("sys_id"):
            continue
        model = r.get("model_fqn", "").split("/")[-1]
        try:
            w = float(r.get("sys_wholeness", 0))
        except ValueError:
            w = 0.0
        best[model] = max(best.get(model, 0.0), w)
    return best, _profile_hits(msf), msf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagnose", action="store_true", help="re-run MacSyFinder on undetected genomes")
    args = ap.parse_args()

    vers, models = installed_models()
    excluded = excluded_systems()
    in_scope = sorted(set(models) - excluded)
    print(f"TXSScan {vers}: {len(models)} models, {len(excluded)} excluded by ssign default, {len(in_scope)} in scope")
    print("  in scope : " + ", ".join(in_scope))
    print("  excluded : " + ", ".join(sorted(excluded)))

    model_genes = {m: v["genes"] for m, v in models.items()}
    detected = detected_per_unit()
    idx = RerunIndex()
    rows = list(csv.DictReader(open(BENCHMARK_CSV)))

    out, need_diag = [], {}
    for b in rows:
        key = b["subtype"] or b["ss_type"]
        compat = COMPATIBLE_MODELS[key]
        missing = sorted(m for m in compat if m not in models)
        j = idx.join(b["contig"], int(b["start"]), int(b["stop"]))
        unit = os.path.basename(str(j["unit"])) if j else ""
        hits = sorted(m for m in compat if detected.get(unit, {}).get(m))
        if not compat:
            status = "no_model"
        elif hits:
            status = "detected"
        else:
            status = "undetected"
            need_diag[b["contig"]] = None
        out.append(
            {
                "instance_id": b["instance_id"],
                "ss_type": b["ss_type"],
                "subtype": b["subtype"],
                "gene": b["gene"],
                "organism": b["organism"][:60],
                "contig": b["contig"],
                "run_unit": unit,
                "compatible_models": ",".join(sorted(compat)) or "-",
                "models_missing_from_txsscan": ",".join(missing) or "-",
                "detected_in_genome": ",".join(hits) or "-",
                "status": status,
                "cause": "no TXSScan model" if status == "no_model" else "",
                "component_profiles_hit": "",
                "best_wholeness": "",
                "reachable_within_3": b["reachable_within_3"],
                "found_by_ssign": b["found_by_ssign"],
            }
        )

    if args.diagnose and need_diag:
        print(f"\ndiagnosing {len(need_diag)} genomes with no compatible detection (re-running MacSyFinder)")
        for contig in sorted(need_diag):
            need_diag[contig] = rerun_macsyfinder(contig, WORK)
            best, _, msf = need_diag[contig]
            print(f"  {contig}: {'unavailable' if msf is None else best}")
        for r in out:
            if r["status"] != "undetected":
                continue
            best, hits, msf = need_diag.get(r["contig"], (None, {}, None))
            if msf is None:
                continue
            compat = COMPATIBLE_MODELS[r["subtype"] or r["ss_type"]]
            near = {m: best[m] for m in compat if m in best}
            if near:
                # Report every compatible model that was found, not just the
                # highest: for Legionella the biologically right call (pT4SSi,
                # Dot/Icm) scores lower than the incidental pT4SSt, and naming
                # only the winner would misstate what TXSScan actually saw.
                ranked = sorted(near.items(), key=lambda kv: -kv[1])
                r["best_wholeness"] = f"{ranked[0][1]:.3f}"
                r["status"] = "low_wholeness"
                r["cause"] = (
                    ", ".join(f"{m} at {w:.3f}" for m, w in ranked) + f"; all below ssign's {WHOLENESS_THRESHOLD}"
                )
                continue
            # No system of a compatible type. Separate the two ways that
            # happens, because only the first is a profile gap: either the
            # model's profiles never matched, or they matched but MacSyFinder
            # could not assemble a cluster that met the model's quorum.
            matched = sorted({g for m in compat for g in model_genes.get(m, ()) if g in hits})
            if not matched:
                r["status"] = "not_detected"
                r["cause"] = "no component profile of a compatible model matched anywhere"
            else:
                r["status"] = "not_colocalised"
                whys = [w for m in compat for w in _rejections(msf, m)]
                seen, reasons = set(), []
                for w in whys:
                    if w not in seen:
                        seen.add(w)
                        reasons.append(w)
                r["component_profiles_hit"] = ",".join(matched)
                r["cause"] = f"{len(matched)} component profiles matched but no cluster passed: " + "; ".join(
                    reasons[:2]
                )

    dest = os.path.join(ROOT, "docs", "benchmark", "txsscan_model_coverage.csv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    print(f"\n{'type':<8}{'n':>4}{'detected':>10}{'no_model':>10}{'undetected':>12}{'low_whole':>11}{'found':>7}")
    agg = collections.defaultdict(collections.Counter)
    for r in out:
        k = r["subtype"] or r["ss_type"]
        agg[k]["n"] += 1
        agg[k][r["status"]] += 1
        if r["found_by_ssign"].strip().lower() in ("yes", "true", "1"):
            agg[k]["found"] += 1
    for k in sorted(agg):
        a = agg[k]
        print(
            f"{k:<8}{a['n']:>4}{a['detected']:>10}{a['no_model']:>10}"
            f"{a['not_detected'] + a['not_colocalised'] + a['undetected']:>9}"
            f"{a['low_wholeness']:>11}{a['found']:>7}"
        )

    print(f"\nrows ssign's machinery step cannot reach ({'diagnosed' if args.diagnose else 'undiagnosed'}):")
    for r in out:
        if r["status"] != "detected":
            print(
                f"  {r['instance_id']:<10}{(r['subtype'] or r['ss_type']):<7}{r['gene'][:16]:<17}"
                f"{r['status']:<17}{r['cause'] or '-'}"
            )
    print(f"\nwrote {dest}")

    untested = sorted(set(in_scope) - {m for k in COMPATIBLE_MODELS for m in COMPATIBLE_MODELS[k]})
    if untested:
        print(f"\nin-scope models with no benchmark row at all: {', '.join(untested)}")


if __name__ == "__main__":
    main()
