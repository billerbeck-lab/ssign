#!/usr/bin/env python3
"""Phase A second-pass review: deterministic locus-IDENTITY audit of every gold row.

Sweep 1 found 7/23 T3SS rows whose effector_locus pointed to the WRONG gene (a transposase / hydrolase /
regulator), not the named effector. The len_ratio precompute (scripts/66) only catches that when the wrong
gene's length differs from the cited accession; a wrong gene that coincidentally length-matches would slip
through. This audit closes that gap mechanically: it asks, for every row that has a UniProt accession,

    does the accession actually sit at the gold effector_locus?

by mapping BOTH the accession's UniProt ordered-locus-name (OLN) and the gold effector_locus through the
gene-order index (which carries old<->new RS locus aliases) and checking they resolve to the SAME gene.

Verdict per row:
  match               - the accession's OLN resolves to the same index locus as the gold effector_locus. CLEAR.
  MISMATCH            - the OLN resolves to a DIFFERENT index locus than the gold locus -> mis-anchored, like
                        the sweep-1 cases. Needs reanchor/drop.
  oln_not_in_index    - OLN present but not in this genome's index (accession likely from another assembly).
  gold_locus_unindexed- the gold effector_locus itself is absent from the index (can't anchor it).
  no_oln              - UniProt has no ordered-locus-name; falls back to the gene-name+length already checked.
  no_acc              - blank accession (none_exists rows); nothing to map.

The match rows are cleared deterministically. Everything else is the residual to resolve (manually or via a
follow-up agent sweep), exactly as sweep 1 handled its flagged subset.

Network: one batched UniProt /accessions call, cached to _locus_audit_cache.json.
Inputs : data/phase2/verification_phase_a/gold_list_final.tsv ; data/phase1/gene_order_index.tsv
Outputs: data/phase2/verification_phase_a/locus_identity_audit.tsv
Run    : .venv/bin/python scripts/70_locus_identity_audit.py
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import locus_index, norm_locus, read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
IDX = BENCH / "data" / "phase1" / "gene_order_index.tsv"
OUT = VD / "locus_identity_audit.tsv"
CACHE = VD / "_locus_audit_cache.json"
API = "https://rest.uniprot.org/uniprotkb/accessions"
FIELDS = "accession,gene_oln,gene_primary"


def fetch(accs: list[str]) -> dict[str, dict]:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    out: dict[str, dict] = {}
    for i in range(0, len(accs), 80):
        url = f"{API}?{urllib.parse.urlencode({'accessions': ','.join(accs[i : i + 80]), 'fields': FIELDS, 'format': 'json'})}"
        with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310 (trusted host)
            for e in json.load(r).get("results", []):
                genes = e.get("genes", [])
                oln = [s["value"] for g in genes for s in g.get("orderedLocusNames", [])]
                primary = next((g["geneName"]["value"] for g in genes if g.get("geneName")), "")
                out[e["primaryAccession"]] = {"oln": oln, "primary": primary}
    CACHE.write_text(json.dumps(out, indent=0))
    return out


def main() -> int:
    rows = read_tsv(GOLD)
    idx = locus_index(IDX)
    accs = sorted({r["uniprot"].strip() for r in rows if r["uniprot"].strip() not in ("", "-")})
    meta = fetch(accs)

    out = []
    for r in rows:
        acc = r["uniprot"].strip()
        gold_loc = r["effector_locus"].strip()
        gold_entry = idx.get(norm_locus(gold_loc))
        m = meta.get(acc, {})
        oln = m.get("oln", [])
        oln_entries = {norm_locus(o): idx.get(norm_locus(o)) for o in oln}

        if acc in ("", "-"):
            verdict, detail = "no_acc", ""
        elif gold_entry is None:
            verdict, detail = "gold_locus_unindexed", f"gold {gold_loc} not in index"
        elif not oln:
            verdict, detail = "no_oln", f"UniProt {acc} has no ordered-locus-name"
        elif any(e == gold_entry for e in oln_entries.values()):
            verdict, detail = "match", f"{'/'.join(oln)} -> {gold_entry[2]}"
        elif any(e is not None for e in oln_entries.values()):
            hit = next((e for e in oln_entries.values() if e is not None))
            verdict, detail = (
                "MISMATCH",
                f"OLN {'/'.join(oln)} -> {hit[2]} (ord {hit[1]}), NOT gold {gold_loc} (ord {gold_entry[1]})",
            )
        else:
            verdict, detail = "oln_not_in_index", f"OLN {'/'.join(oln)} absent from {r['genome']} index"

        # A `match` already proves same genome: the index key is the gold genome's record_acc, so the OLN
        # resolving to the gold locus's index entry means same molecule. (UniProt's own RefSeq nucleotide
        # xref is NOT used as a cross-check -- it points at whatever assembly UniProt mirrors, e.g. NC_002696
        # for the same genome the gold list cites as INSDC AE005673, which would only add false mismatches.)
        out.append(
            {
                "instance_id": r["instance_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "uniprot": acc,
                "gold_locus": gold_loc,
                "uniprot_oln": "/".join(oln),
                "uniprot_gene": m.get("primary", ""),
                "verdict": verdict,
                "detail": detail,
            }
        )

    order = {
        "MISMATCH": 0,
        "oln_not_in_index": 1,
        "gold_locus_unindexed": 2,
        "no_oln": 3,
        "no_acc": 4,
        "match": 5,
    }
    out.sort(key=lambda x: (order.get(x["verdict"], 9), x["instance_id"]))
    write_tsv(OUT, list(out[0].keys()), out)

    tally = Counter(r["verdict"] for r in out)
    print(f"locus_identity_audit.tsv: {len(out)} rows -> {dict(tally)}")
    for r in out:
        if r["verdict"] not in ("match", "no_acc"):
            print(f"  {r['verdict']:20} {r['instance_id']:9} {r['gene']:16} {r['uniprot']:11} {r['detail']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
