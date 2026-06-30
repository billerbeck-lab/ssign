#!/usr/bin/env python3
"""Adversarial identity audit of every gold row via its UniProt DOMAIN, not just its length.

The high-churn re-review caught two rows (T6SS_17 Tae4, T6SS_18 Tle1) that had been anchored onto the WRONG
gene: an earlier swap_up assigned an accession whose COORDINATES matched (so the length check passed) but whose
PROTEIN was a structural/regulatory neighbour (TagF regulator / TagL OmpA accessory), not the effector. Only the
Pfam domain revealed it. That check was run reactively, on two rows. This runs it on all 90.

For every row with a UniProt accession, fetch the entry and flag:
  RED_DOMAIN  - Pfam/name matches a class that is machinery/structural/mobile-element, never a secreted substrate
                (OmpA/TagL, TagF, transposase/integrase/IS, pilin/flagellin, ribosomal, porin-as-machinery, ...).
                This is the systematic version of the T6SS_17/18 catch.
  LEN         - |uniprot_len / (cds_span/3 - 1) - 1| > 0.15 (coords encode a very different-size protein).
  ORG         - UniProt organism genus != the row's organism genus (possible wrong-strain/wrong-species anchor).
  NO_PFAM     - no Pfam at all (can't assess by domain; informational, eyeball the name).
NO_ACCESSION rows (blank uniprot) are listed separately -- they were never given an identity to check.

The RED_DOMAIN denylist is high-precision (only unambiguous non-effector classes) to keep false positives near
zero; it is NOT a guarantee of completeness -- a wrong gene with a plausible effector domain won't trip it, so
the printed full table is meant to be eyeballed too. Read-only; writes domain_audit.tsv.
Run: .venv/bin/python scripts/74_domain_audit.py
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

VD = Path(__file__).resolve().parents[1] / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
OUT = VD / "domain_audit.tsv"
CACHE = VD / "_uniprot_entry_cache.json"
LEN_TOL = 0.15

# Unambiguous non-substrate protein classes. A gold effector whose UniProt domain/name hits one of these is
# almost certainly anchored on the wrong gene (machinery, structural accessory, or a mobile/housekeeping ORF).
RED_FLAGS = (
    "ompa",
    "tagf",
    "tagl",
    "transposase",
    "integrase",
    "insertion element",
    "is3",
    "is4",
    "is5",
    "is630",
    "recombinase",
    "ribosomal protein",
    "trna ",
    "flagellin",
    "flagellar",
    "pilin",
    "type iv pilus",
    "dna polymerase",
    "rna polymerase",
    "elongation factor",
    "cell division protein",
    "membrane protein insertase",
    "outer membrane protein assembly",
    "bama",
    "secretion system protein",
    "atp synthase",
)


def gene_tokens(gene: str) -> list[str]:
    """Meaningful name tokens from a row's gene label, e.g. 'TplE_alias_Tle4' -> ['tple','tle4'],
    'YopJ_YopP' -> ['yopj','yopp'], 'YezP_Zn_binding' -> ['yezp']. Used to ask whether the assigned
    UniProt actually names this effector (the tell that exposed the T6SS_17/18 wrong-gene anchors)."""
    drop = {"alias", "stm", "sci1", "zn", "binding", "related", "like", "protease", "nuclease", "ec"}
    toks = []
    for part in gene.replace("/", "_").split("_"):
        p = part.strip().lower()
        if p and p not in drop and not p.isdigit():
            toks.append(p)
    return toks or [gene.strip().lower()]


def fetch(acc: str, cache: dict) -> dict | None:
    acc = acc.strip()
    if not acc or acc == "-":
        return None
    if acc in cache:
        return cache[acc]
    try:
        with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=25) as fh:
            j = json.load(fh)
        genes = [g.get("geneName", {}).get("value", "") for g in j.get("genes", [])]
        names = [j.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")]
        names += [
            a.get("fullName", {}).get("value", "") for a in j.get("proteinDescription", {}).get("alternativeNames", [])
        ]
        pfam = [
            f"{x['id']}:{x.get('properties', [{}])[0].get('value', '')}"
            for x in j.get("uniProtKBCrossReferences", [])
            if x["database"] == "Pfam"
        ]
        kw = [k.get("name", "") for k in j.get("keywords", [])]
        entry = {
            "len": j.get("sequence", {}).get("length"),
            "names": [n for n in names if n],
            "genes": [g for g in genes if g],
            "pfam": pfam,
            "keywords": kw,
            "organism": j.get("organism", {}).get("scientificName", ""),
        }
    except Exception as e:  # noqa: BLE001
        print(f"  WARN fetch failed {acc}: {e}", file=sys.stderr)
        return None
    cache[acc] = entry
    return entry


def main() -> int:
    gold = read_tsv(GOLD)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    no_acc, out = [], []
    for r in gold:
        acc = r["uniprot"].strip()
        if not acc or acc == "-":
            no_acc.append((r["instance_id"], r["ss_type"], r["gene"], r["effector_locus"]))
            continue
        e = fetch(acc, cache)
        if e is None:
            out.append(
                {
                    "instance_id": r["instance_id"],
                    "ss_type": r["ss_type"],
                    "gene": r["gene"],
                    "uniprot": acc,
                    "uniprot_name": "",
                    "pfam": "",
                    "flags": "FETCH_FAIL",
                    "uniprot_len": "",
                    "expected_aa": "",
                    "uniprot_organism": "",
                    "row_organism": r["organism"],
                }
            )
            continue
        blob = " ".join(e["names"] + e["pfam"]).lower()
        flags = []
        red = [t for t in RED_FLAGS if t in blob]
        if red:
            flags.append("RED_DOMAIN:" + "/".join(red))
        expected_aa = (int(r["stop"]) - int(r["start"]) + 1) // 3 - 1
        if e["len"] and expected_aa and abs(e["len"] / expected_aa - 1) > LEN_TOL:
            flags.append(f"LEN({e['len']}/{expected_aa})")
        rg = r["organism"].split()[0].lower() if r["organism"].split() else ""
        ug = e["organism"].split()[0].lower() if e["organism"].split() else ""
        if rg and ug and rg != ug:
            flags.append(f"ORG({ug}!={rg})")
        if not e["pfam"]:
            flags.append("NO_PFAM")
        # does the UniProt entry actually name this effector? (entries often use a locus tag as the name,
        # so a miss is a candidate-for-eyeball, not proof of a wrong anchor -- but it is exactly the T6SS_17/18 tell)
        name_blob = " ".join(e["genes"] + e["names"]).lower()
        if name_blob and not any(t in name_blob for t in gene_tokens(r["gene"])):
            flags.append("NAME_MISMATCH")
        out.append(
            {
                "instance_id": r["instance_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "uniprot": acc,
                "uniprot_name": "; ".join(e["names"])[:80],
                "pfam": "; ".join(e["pfam"])[:90],
                "uniprot_len": e["len"],
                "expected_aa": expected_aa,
                "uniprot_organism": e["organism"],
                "row_organism": r["organism"],
                "flags": ",".join(flags) or "ok",
            }
        )

    CACHE.write_text(json.dumps(cache, indent=0))
    header = [
        "instance_id",
        "ss_type",
        "gene",
        "uniprot",
        "uniprot_name",
        "pfam",
        "uniprot_len",
        "expected_aa",
        "uniprot_organism",
        "row_organism",
        "flags",
    ]
    out.sort(key=lambda x: (x["flags"] == "ok", x["instance_id"]))
    write_tsv(OUT, header, out)

    from collections import Counter

    flagged = [x for x in out if x["flags"] != "ok"]
    print(f"domain_audit.tsv: {len(out)} rows with a UniProt accession ({len(no_acc)} rows have NO accession)")
    cats = Counter(f.split("(")[0].split(":")[0] for x in flagged for f in x["flags"].split(","))
    print(f"  flag categories: {dict(cats)}")
    if no_acc:
        print(f"\n  NO UNIPROT ACCESSION ({len(no_acc)}) -- never identity-checked:")
        for iid, ss, gene, loc in no_acc:
            print(f"    {iid:9} {ss:6} {gene:16} locus={loc}")
    if flagged:
        print(f"\n  FLAGGED ({len(flagged)}):")
        for x in flagged:
            print(f"    {x['instance_id']:9} {x['gene']:14} [{x['flags']}]")
            print(f"        name={x['uniprot_name']!r} pfam={x['pfam']!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
