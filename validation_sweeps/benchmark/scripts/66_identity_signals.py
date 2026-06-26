#!/usr/bin/env python3
"""Phase A second pass: deterministic identity signals for every gold-list row.

These are the mechanical checks that should NOT be left to an LLM -- they are computed here and fed into
the agent sweep batches so agents adjudicate the flags instead of re-deriving (and possibly hallucinating)
the numbers. Signals per row:

  len_ratio   - UniProt sequence length / expected aa from the gold CDS coordinates ((stop-start+1)/3 - 1).
                ~1.0 = the accession matches the gene at those coords; <<1 = the accession is a FRAGMENT or a
                different (shorter) gene; >>1 = the coords are a fragment of a longer protein, or wrong gene.
                This is the check that caught T4SS_05 (251 bp vs 684 aa) and T5SS_03 (496 aa Fragment vs
                1815 aa) by hand; run over all rows it flags every such mismatch.
  reviewed    - is the cited accession a Swiss-Prot (Reviewed) entry? TrEMBL rows where a Reviewed entry
                for the same gene+organism exists are candidates to upgrade (higher-confidence identity).
  fragment    - UniProt flags the entry itself as a Fragment.
  dup         - the accession or the effector_locus appears on more than one gold row.
  no_acc      - row has no UniProt accession; len_ratio can't be computed from UniProt (needs CDS translation).

Network: one batched call to the UniProt /accessions endpoint, cached to _identity_cache.json so re-runs are
offline and reproducible. Delete the cache to refetch.

Inputs : data/phase2/verification_phase_a/gold_list_final.tsv
Outputs: data/phase2/verification_phase_a/identity_signals.tsv
Run    : .venv/bin/python scripts/66_identity_signals.py
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
OUT = VD / "identity_signals.tsv"
CACHE = VD / "_identity_cache.json"
API = "https://rest.uniprot.org/uniprotkb/accessions"
FIELDS = "accession,reviewed,length,protein_name,gene_names,organism_name,ft_non_ter"

FRAG_LO, FRAG_HI = 0.85, 1.15  # len_ratio outside this band is flagged


def fetch(accs: list[str]) -> dict[str, dict]:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    out: dict[str, dict] = {}
    for i in range(0, len(accs), 80):  # endpoint caps the accession list length
        chunk = accs[i : i + 80]
        url = f"{API}?{urllib.parse.urlencode({'accessions': ','.join(chunk), 'fields': FIELDS, 'format': 'json'})}"
        with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 (trusted host)
            for e in json.load(r).get("results", []):
                acc = e["primaryAccession"]
                out[acc] = {
                    "reviewed": e.get("entryType", ""),
                    "length": e.get("sequence", {}).get("length"),
                    "fragment": any(f.get("type") == "Non-terminal residue" for f in e.get("features", [])),
                    "protein": (e.get("proteinDescription", {}).get("recommendedName", {}) or {})
                    .get("fullName", {})
                    .get("value", ""),
                }
    CACHE.write_text(json.dumps(out, indent=0))
    return out


def main() -> int:
    rows = read_tsv(GOLD)
    accs = sorted({r["uniprot"].strip() for r in rows if r["uniprot"].strip() not in ("", "-")})
    meta = fetch(accs)

    acc_count = Counter(r["uniprot"].strip() for r in rows if r["uniprot"].strip() not in ("", "-"))
    loc_count = Counter(r["effector_locus"].strip() for r in rows if r["effector_locus"].strip())

    out, flagged = [], 0
    for r in rows:
        acc = r["uniprot"].strip()
        try:
            span = int(r["stop"]) - int(r["start"]) + 1
            exp_aa = round(span / 3) - 1
        except (ValueError, KeyError):
            span, exp_aa = None, None
        m = meta.get(acc, {})
        up_len = m.get("length")
        ratio = round(up_len / exp_aa, 3) if (up_len and exp_aa and exp_aa > 0) else ""
        flags = []
        if acc in ("", "-"):
            flags.append("no_acc")
        else:
            if isinstance(ratio, float) and not (FRAG_LO <= ratio <= FRAG_HI):
                flags.append("len_mismatch")
            if m.get("fragment"):
                flags.append("uniprot_fragment")
            if m.get("reviewed") != "UniProtKB reviewed (Swiss-Prot)":
                flags.append("unreviewed")
            if acc_count[acc] > 1:
                flags.append("dup_acc")
        if r["effector_locus"].strip() and loc_count[r["effector_locus"].strip()] > 1:
            flags.append("dup_locus")
        if any(f != "unreviewed" for f in flags):  # unreviewed alone is informational, not a "flag"
            flagged += 1
        out.append(
            {
                "instance_id": r["instance_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "uniprot": acc,
                "reviewed": "yes"
                if m.get("reviewed") == "UniProtKB reviewed (Swiss-Prot)"
                else ("no" if acc not in ("", "-") else ""),
                "uniprot_len": up_len or "",
                "expected_aa": exp_aa if exp_aa is not None else "",
                "len_ratio": ratio,
                "uniprot_protein": m.get("protein", ""),
                "flags": ";".join(flags),
            }
        )

    out.sort(key=lambda x: (x["flags"] == "", x["ss_type"], x["instance_id"]))  # flagged rows first
    write_tsv(OUT, list(out[0].keys()), out)

    fc = Counter(f for r in out for f in r["flags"].split(";") if f)
    print(f"identity_signals.tsv: {len(out)} rows, {flagged} with a non-trivial flag")
    print(f"  flag counts: {dict(fc.most_common())}")
    print("  rows needing a look (any flag except bare 'unreviewed'):")
    for r in out:
        fl = [f for f in r["flags"].split(";") if f]
        if fl and fl != ["unreviewed"]:
            print(
                f"    {r['instance_id']:9} {r['gene']:18} {r['uniprot']:12} "
                f"len {r['uniprot_len']}/{r['expected_aa']} ratio={r['len_ratio']}  -> {';'.join(fl)}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
