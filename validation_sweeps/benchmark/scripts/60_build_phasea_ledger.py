#!/usr/bin/env python3
"""Phase A (benchmark-final-validation 1.1): build the full re-verification ledger.

Dumps every effector in positives_all.tsv (the citation-verified answer key, 337 rows) into one
ledger TSV for a fresh, complete blind re-verification. This supersedes the earlier tasks-74-84 audit,
which covered only 245 rows (the answer key has since grown, mostly T3SS). Each row carries the exact
identity (so 1.6 can locate the source row) plus everything a verifier needs: gene, uniprot, locus,
organism, family, the primary-reference DOI, and the verbatim citation_quote.

Stable row ids PA001.. are assigned in (ss_type, gene, uniprot) order so batches are reproducible.

Output: data/phase2/verification_phase_a/phasea_ledger.tsv
Run   : .venv/bin/python scripts/60_build_phasea_ledger.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
POS = BENCH / "data" / "dataset" / "positives_all.tsv"
OUTDIR = BENCH / "data" / "phase2" / "verification_phase_a"
OUT = OUTDIR / "phasea_ledger.tsv"

SS_ORDER = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS"]
# identity columns copied verbatim (the apply step in 1.6 locates the source row by these) + the
# context a verifier reads. citation_quote is the verbatim experimental sentence already in the key.
FIELDS = [
    "row_id",
    "ss_type",
    "subtype",
    "gene",
    "uniprot",
    "locus_tag",
    "organism",
    "refseq_genome",
    "family",
    "primary_ref",
    "citation_quote",
    "verification_status",
    "audit_tier",
]


def main() -> int:
    rows = read_tsv(POS)
    ss_rank = {s: i for i, s in enumerate(SS_ORDER)}
    rows.sort(key=lambda r: (ss_rank.get(r.get("ss_type", ""), 99), r.get("gene", "").lower(), r.get("uniprot", "")))
    out = []
    for i, r in enumerate(rows, 1):
        rec = {f: r.get(f, "") for f in FIELDS}
        rec["row_id"] = f"PA{i:03d}"
        rec["citation_quote"] = r.get("citation_quote") or r.get("quote", "")  # quote is the older col name
        out.append(rec)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT, FIELDS, out)
    from collections import Counter

    by = Counter(r["ss_type"] for r in out)
    print(f"wrote {OUT.relative_to(BENCH)}: {len(out)} effectors")
    print("  by ss_type: " + ", ".join(f"{s} {by[s]}" for s in SS_ORDER))
    print(f"  missing citation_quote: {sum(1 for r in out if not r['citation_quote'].strip())}")
    print(f"  missing real uniprot:   {sum(1 for r in out if not r['uniprot'].strip() or r['uniprot'] == '-')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
