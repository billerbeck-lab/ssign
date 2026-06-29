#!/usr/bin/env python3
"""Cross-field consistency audit of the HIGH-CHURN gold rows (>=2 distinct correction types).

Rows touched by several independent correction passes (identity swap + citation fix + re-anchor + found
recompute) are where a *compounding* error hides: one pass updates a field, a dependent field is left
stale, and no single sweep re-checks the combination. This audit re-derives the machine-checkable
relationships on the CURRENT sheet and flags any that no longer hold:

  LEN: does the assigned UniProt's length match the protein the coordinates encode?
       expected_aa = (stop - start + 1)/3 - 1 (drop the stop codon); len_ratio = uniprot_len / expected_aa.
       A re-anchor that moved coords onto the wrong gene, or a uniprot swap to a different-size protein,
       shows up as len_ratio far from 1.0. Tolerance 0.15 is generous on purpose (signal-peptide cleavage,
       passenger processing, mature-vs-precursor isoforms all shift length a little); flagged rows are for
       eyeballing, not auto-correction.
  OVL: do the coords land on a real gene in the tier-2 rerun (n_overlap > 0)? A re-anchor onto an
       intergenic coordinate would overlap nothing. (NleD's genuine Bakta-miss is the one known n_overlap=0.)
  FOUND: is found_by_ssign internally consistent with the overlap recompute (sanity, already enforced by 65).

UniProt lengths are pulled live from the REST API (.fasta), cached to _uniprot_len_cache.json so reruns are
offline. Read-only: prints a table + flag summary, writes highchurn_consistency.tsv. No sheet mutation.

Run: .venv/bin/python scripts/72_highchurn_consistency.py
"""

from __future__ import annotations

import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402
from rerun_coords import RerunIndex  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
OUT = VD / "highchurn_consistency.tsv"
LENCACHE = VD / "_uniprot_len_cache.json"
LEN_TOL = 0.15


def correction_types(corr: str) -> set[str]:
    """Distinct correction *kinds* recorded in a row's correction string (mirrors scripts/71 reporting)."""
    cl, p = corr.lower(), set()
    if any(k in cl for k in ("reanchor", "old locus", "wrong_genome", "swap_up", "->")):
        p.add("identity/anchor")
    if "swap_ref" in cl or "fix_ref" in cl:
        p.add("ref")
    if "citation" in cl or "quote" in cl:
        p.add("quote")
    if "found_recompute" in cl:
        p.add("found")
    if "note" in cl:
        p.add("note")
    return p


def uniprot_len(acc: str, cache: dict) -> int | None:
    acc = acc.strip()
    if not acc:
        return None
    if acc in cache:
        return cache[acc]
    n = None
    try:
        with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=20) as fh:
            seq = "".join(ln.strip() for ln in fh.read().decode().splitlines() if not ln.startswith(">"))
        n = len(seq) or None
    except Exception as e:  # noqa: BLE001 -- network/404 -> don't crash the audit
        print(f"  WARN uniprot fetch failed {acc}: {e}", file=sys.stderr)
    if n is not None:  # cache hits only; a transient miss must re-fetch next run, never freeze a row as
        cache[acc] = n  # length-unknown (that would silently disable its LEN check)
    return n


def main() -> int:
    gold = read_tsv(GOLD)
    cache = json.loads(LENCACHE.read_text()) if LENCACHE.exists() else {}
    ridx = RerunIndex()

    hc = [r for r in gold if len(correction_types(r["correction"])) >= 2]
    out = []
    for r in sorted(hc, key=lambda r: (-len(correction_types(r["correction"])), r["instance_id"])):
        cts = correction_types(r["correction"])
        start, stop = int(r["start"]), int(r["stop"])
        expected_aa = (stop - start + 1) // 3 - 1
        ulen = uniprot_len(r["uniprot"], cache)
        ratio = round(ulen / expected_aa, 3) if (ulen and expected_aa) else None
        eo = ridx.emitted_overlap(r["contig"].strip(), start, stop, strand=r.get("strand"))
        n_ov = eo["n_overlap"] if eo else 0
        flags = []
        if ratio is not None and abs(ratio - 1) > LEN_TOL:
            flags.append("LEN")
        if eo is None:
            flags.append("UNRECONCILED")
        elif n_ov == 0:
            flags.append("NO_OVERLAP")
        if eo and r.get("found_by_ssign", "") != eo["found"]:
            flags.append("FOUND_STALE")
        out.append(
            {
                "instance_id": r["instance_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "uniprot": r["uniprot"],
                "effector_locus": r["effector_locus"],
                "n_types": len(cts),
                "types": ",".join(sorted(cts)),
                "uniprot_len": ulen if ulen is not None else "",
                "expected_aa": expected_aa,
                "len_ratio": ratio if ratio is not None else "",
                "n_overlap": n_ov,
                "overlap_bp": eo["best_overlap_bp"] if eo else "",
                "found_by_ssign": r.get("found_by_ssign", ""),
                "flags": ",".join(flags) or "ok",
            }
        )

    LENCACHE.write_text(json.dumps(cache, indent=0))
    header = [
        "instance_id",
        "ss_type",
        "gene",
        "uniprot",
        "effector_locus",
        "n_types",
        "types",
        "uniprot_len",
        "expected_aa",
        "len_ratio",
        "n_overlap",
        "overlap_bp",
        "found_by_ssign",
        "flags",
    ]
    out.sort(key=lambda x: (x["flags"] == "ok", -x["n_types"], x["instance_id"]))
    write_tsv(OUT, header, out)

    flagged = [x for x in out if x["flags"] != "ok"]
    print(f"highchurn_consistency.tsv: {len(out)} rows (>=2 correction types)")
    print(f"  flag tally: {dict(Counter(x['flags'] for x in out))}")
    if flagged:
        print(f"\n  FLAGGED ({len(flagged)}) -- eyeball these:")
        for x in flagged:
            print(
                f"    {x['instance_id']:9} {x['ss_type']:6} {x['gene']:14} {x['types']:30} "
                f"len {x['uniprot_len']}/{x['expected_aa']} (ratio {x['len_ratio']}) "
                f"n_ov={x['n_overlap']} -> {x['flags']}"
            )
    else:
        print("  all high-churn rows internally consistent (len concordant, coords on a gene, found fresh)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
