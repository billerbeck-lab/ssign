#!/usr/bin/env python3
"""Build the corrected genome panel for the tier-2 CX3 rerun (change benchmark-final-validation).

The current fleet staged 67 assembly units (`inputs_gb/`), built for the proximity benchmark only.
Two corrections are needed for the rerun:

  1. Drop the staged units that carry NO testable answer-key system. A unit is kept for proximity if
     it has >=1 testable proximity SYSTEM INSTANCE (T1/T2/T3/T4/T6), computed with the SAME instance
     logic + audit-drop filter the recall figure uses (script 52). That yields 51 keep-units; the
     remaining staged units are drop candidates -- UNLESS they also carry a T5SS effector (then they
     are rescued into the T5SS set, not dropped).

  2. Add the T5SS genomes missing from the fleet. T5SS was excluded from the proximity panel, so most
     of its 20 genomes were never staged. Each T5SS genome (t5ss_effectors.tsv `refseq_genome`, given
     as a RefSeq sequence accession OR a GCF_ assembly accession) is resolved to its sequence-accession
     bases and classified:
        staged          - a base is a replicon of a staged inputs_gb unit (already in the run)
        cache-not-staged- a base is in data/refseq_cache but not staged (stage from cache, no fetch)
        needs-fetch     - neither (resolve + fetch the assembly, then stage)
     RefSeq (NC_) and INSDC (AE/BX/CP/U...) accessions name the same genome, so each accession is
     resolved to its full base-alias set (GCF_ via the Datasets API, NC_ via eutils esummary
     `assemblyacc`), cached to _t5_acc_cache.json so reruns are offline. Without this, a T5SS genome
     listed as NC_002929.2 but staged as BX470248 is wrongly counted as missing.

The final panel = the proximity keep-units UNION the T5SS-carrying units. The manifest records, per
genome, its role (proximity / T5SS / both), staged status, and the action needed before submitting.

Inputs : data/phase2/actual_per_effector.panel_genbank_{default,t3ss}.tsv, data/phase1/ceiling_per_effector.tsv
         data/phase2/panel_manifest.tsv, data/dataset/t5ss_effectors.tsv, data/refseq_cache/*.gb
Outputs: data/phase2/rerun_panel_manifest.tsv   (genome, accession, role, staged_status, action, unit_id)
Run    : <repo>/.venv/bin/python scripts/57_build_rerun_panel.py
"""

from __future__ import annotations

import csv
import json
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import clean_dataset  # noqa: E402
from bench_index import accession_base as norm  # noqa: E402  (version+prefix-stripped base, drift-tolerant join)
from bench_io import read_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
CEILING = BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"
MANIFEST = P2 / "panel_manifest.tsv"
T5SS = BENCH / "data" / "dataset" / "t5ss_effectors.tsv"
CACHE = BENCH / "data" / "refseq_cache"
ACC_CACHE = P2 / "_t5_acc_cache.json"
OUT = P2 / "rerun_panel_manifest.tsv"
PROX = {"T1SS", "T2SS", "T3SS", "T4SS", "T6SS"}


def proximity_keep_units() -> set[str]:
    """Staged units with >=1 testable proximity SYSTEM INSTANCE (exact script-52 logic + audit filter)."""
    dropped = clean_dataset.dropped_id()
    cei = read_tsv(CEILING)
    by_loc = {r["effector_locus"]: r.get("instance_id", "") for r in cei if r.get("effector_locus", "").strip()}
    by_uni = {
        r["uniprot"]: r.get("instance_id", "") for r in cei if r.get("uniprot", "").strip() and r["uniprot"] != "-"
    }

    def genomes(tag: str, types: set[str]) -> set[str]:
        inst: dict[tuple, dict] = defaultdict(lambda: {"testable": False, "genome": None})
        for r in clean_dataset.load_clean_actual(P2 / f"actual_per_effector.{tag}.tsv"):
            if (r["gene"], r["uniprot"]) in dropped or (r["gene"], r.get("effector_locus", "")) in dropped:
                continue
            iid = by_loc.get(r.get("effector_locus", "")) or by_uni.get(r.get("uniprot", ""))
            if not iid or r["ss_type"] not in types:
                continue
            v = inst[(r["ss_type"], iid)]
            v["genome"] = v["genome"] or r.get("unit_id")
            if r["testable"] == "yes":
                v["testable"] = True
        return {v["genome"] for v in inst.values() if v["testable"] and v["genome"]}

    return genomes("panel_genbank_default", PROX - {"T3SS"}) | genomes("panel_genbank_t3ss", {"T3SS"})


def _http_json(url: str):
    """GET + parse JSON, throttled + retried for NCBI's 3-req/s anonymous limit (HTTP 429)."""
    import time

    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 (fixed NCBI hosts)
                d = json.load(r)
            time.sleep(0.4)  # stay under 3 req/s for the next call
            return d
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < 4:
                time.sleep(2 * (attempt + 1))
                continue
            raise


def genome_bases(acc: str, cache: dict) -> list[str]:
    """Any genome accession -> all its sequence-accession bases (RefSeq + INSDC), cached.

    RefSeq (NC_) and INSDC (AE/BX/CP/U/AM/L...) accessions name the SAME genome, and the fleet stages a
    genome under whichever spelling the corpus used. So a naive base match misses a T5SS genome listed
    as NC_002929.2 that is actually staged as BX470248. We resolve both spellings:
      - GCF_ assembly  -> Datasets API sequence_reports (every replicon's refseq + genbank accession)
      - sequence acc   -> its own base PLUS the INSDC alias from eutils esummary `assemblyacc`
    """
    if acc in cache:
        return cache[acc]
    bases: set[str] = set()
    try:
        if acc.upper().startswith("GCF_"):
            d = _http_json(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/sequence_reports")
            for rep in d.get("reports", []):
                for k in ("refseq_accession", "genbank_accession"):
                    if rep.get(k):
                        bases.add(norm(rep[k]))
        else:
            bases.add(norm(acc))
            d = _http_json(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={acc}&retmode=json"
            )
            uid = d["result"]["uids"][0]
            alias = d["result"][uid].get("assemblyacc", "")  # the INSDC sequence accession alias
            if alias:
                bases.add(norm(alias))
    except Exception as e:  # noqa: BLE001
        print(f"  WARN could not resolve {acc}: {e}")
        return [norm(acc)] if not acc.upper().startswith("GCF_") else []  # transient: not cached, retried next run
    cache[acc] = sorted(bases)
    return cache[acc]


def main() -> int:
    # staged replicon base -> unit_id, from the existing panel manifest (authoritative "is staged")
    rep2unit: dict[str, str] = {}
    for r in read_tsv(MANIFEST):
        for rep in r["replicons"].split(","):
            rep2unit[norm(rep)] = r["unit_id"]
    staged_units = {r["unit_id"] for r in read_tsv(MANIFEST)}
    cache_bases = {norm(p.stem) for p in CACHE.glob("*.gb")}

    keep = proximity_keep_units()

    # T5SS genomes -> staged status
    acc_cache = json.loads(ACC_CACHE.read_text()) if ACC_CACHE.exists() else {}
    t5_genomes = sorted({r["refseq_genome"].strip() for r in read_tsv(T5SS) if r.get("refseq_genome", "").strip()})
    t5_status: dict[str, tuple[str, str]] = {}  # genome -> (status, unit_or_blank)
    for g in t5_genomes:
        bases = genome_bases(g, acc_cache)
        unit = next((rep2unit[b] for b in bases if b in rep2unit), "")
        if unit:
            t5_status[g] = ("staged", unit)
        elif any(b in cache_bases for b in bases):
            t5_status[g] = ("cache-not-staged", "")
        else:
            t5_status[g] = ("needs-fetch", "")
    ACC_CACHE.write_text(json.dumps(acc_cache, indent=0))

    t5_staged_units = {u for s, u in t5_status.values() if s == "staged" and u}

    # final panel + drop set
    panel_units = keep | t5_staged_units  # staged units in the rerun
    drop_units = staged_units - panel_units  # staged but neither proximity-testable nor T5SS-bearing
    t5_to_add = {g: s for g, (s, u) in t5_status.items() if s != "staged"}  # cached or needs-fetch

    # write manifest
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "role", "staged_status", "action", "unit_id"])
        for u in sorted(keep):
            role = "proximity+T5SS" if u in t5_staged_units else "proximity"
            w.writerow([u, role, "staged", "run", u])
        for u in sorted(t5_staged_units - keep):
            w.writerow([u, "T5SS", "staged", "run", u])
        for g in sorted(t5_to_add):
            act = "stage-from-cache" if t5_to_add[g] == "cache-not-staged" else "fetch+stage"
            w.writerow([g, "T5SS", t5_to_add[g], act, ""])
        for u in sorted(drop_units):
            w.writerow([u, "none", "staged", "DROP (no testable answer-key system)", u])

    # summary
    n_panel = len(panel_units) + len(t5_to_add)
    print(f"proximity keep-units            : {len(keep)}")
    print(f"T5SS genomes (distinct)         : {len(t5_genomes)}")
    print(
        f"  already in a staged unit      : {len(t5_staged_units)} units ({sum(1 for s, _ in t5_status.values() if s == 'staged')} genomes)"
    )
    print(f"  cache-not-staged (stage only) : {sum(1 for s in t5_to_add.values() if s == 'cache-not-staged')}")
    print(f"  needs-fetch                   : {sum(1 for s in t5_to_add.values() if s == 'needs-fetch')}")
    print(f"staged units DROPPED            : {len(drop_units)}")
    print(
        f"FINAL PANEL                     : {n_panel} genomes "
        f"({len(panel_units)} already staged + {len(t5_to_add)} to add)"
    )
    print(f"\nwrote {OUT.relative_to(BENCH)}")
    print("\nT5SS to ADD:")
    for g in sorted(t5_to_add):
        print(f"  {g:18s} {t5_to_add[g]}")
    print("\nDROPPED units:", ", ".join(sorted(drop_units)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
