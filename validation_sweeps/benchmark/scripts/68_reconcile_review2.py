#!/usr/bin/env python3
"""Phase A second-pass review: reconcile the per-sweep agent verdicts into one table per sweep.

Generic over the two sweeps (scripts/67 layout). For each gold row it tallies the K independent agents'
verdicts, picks the plurality, flags anything that is not a unanimous "clean" verdict for human
adjudication, and carries the modal proposed correction (accession for the identity sweep; ref/quote for
the citation sweep) with its support count. Coverage + row_id drift are reported, not assumed.

The output table is the input to the human adjudication step: I read it, decide each non-clean row, and
encode the outcome in scripts/65 DISPOSITIONS (raw gold list is never touched directly).

Usage : .venv/bin/python scripts/68_reconcile_review2.py sweep1_identity
        .venv/bin/python scripts/68_reconcile_review2.py sweep2_citation
Inputs : data/phase2/verification_phase_a/gold_review2/<sweep>/batches/*.tsv
         data/phase2/verification_phase_a/gold_review2/<sweep>/verdicts/agent{k}_<batch>.tsv
Outputs: data/phase2/verification_phase_a/gold_review2/<sweep>/reconciliation.tsv
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import modal, read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
GR2 = BENCH / "data" / "phase2" / "verification_phase_a" / "gold_review2"
N_AGENTS = 3  # independent agents per batch

# Per-sweep config: which column holds the verdict, which verdict values need no change ("clean"),
# the canonical key field carried from the batch, and the proposed-correction columns to take a mode of.
SWEEPS = {
    "sweep1_identity": {
        "verdict_col": "identity_verdict",
        "clean": {"ok", "none_exists"},
        "key_col": "uniprot_in",
        "proposed": {"uniprot": "proposed_uniprot"},
    },
    "sweep2_citation": {
        "verdict_col": "citation_verdict",
        "clean": {"ok"},
        "key_col": "uniprot",
        "proposed": {"primary_ref": "proposed_primary_ref", "quote": "proposed_quote"},
    },
}


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in SWEEPS:
        print(f"usage: {sys.argv[0]} {{{'|'.join(SWEEPS)}}}", file=sys.stderr)
        return 2
    sweep = sys.argv[1]
    cfg = SWEEPS[sweep]
    vcol, clean, key_col, proposed = cfg["verdict_col"], cfg["clean"], cfg["key_col"], cfg["proposed"]
    sdir = GR2 / sweep
    batch_files = sorted((sdir / "batches").glob("*.tsv"))
    if not batch_files:
        print(f"no batch files in {sdir / 'batches'}", file=sys.stderr)
        return 1

    row_meta: dict[str, dict] = {}
    batch_rows: dict[str, list[dict]] = {}
    for bf in batch_files:
        rows = read_tsv(bf)
        batch_rows[bf.stem] = rows
        for r in rows:
            row_meta[r["row_id"]] = r

    verdicts: dict[str, list[dict]] = defaultdict(list)
    coverage, drift = [], []
    for batch, rows in batch_rows.items():
        expected = {r["row_id"] for r in rows}
        contributors = set()  # agents that supplied >=1 expected row (an empty/all-drift file does not count)
        for k in range(1, N_AGENTS + 1):
            vf = sdir / "verdicts" / f"agent{k}_{batch}.tsv"
            try:
                vrows = read_tsv(vf)
            except FileNotFoundError:
                continue
            got = {r["row_id"] for r in vrows}
            if got != expected:
                drift.append(f"agent{k}_{batch}: missing={sorted(expected - got)} extra={sorted(got - expected)}")
            for r in vrows:
                if r["row_id"] in expected:
                    r["_agent_k"] = k
                    verdicts[r["row_id"]].append(r)
                    contributors.add(k)
        present = len(contributors)
        coverage.append(f"{batch}: {present}/{N_AGENTS}" + ("" if present == N_AGENTS else "  <-- INCOMPLETE"))

    fields = (
        [
            "row_id",
            "ss_type",
            "gene",
            "key_id",
            "n_agents",
            "verdict_counts",
            "final_verdict",
            "unanimous",
            "adjudicate",
        ]
        + [f"proposed_{k}" for k in proposed]
        + [f"{k}_support" for k in proposed]
        + ["agent_notes"]
    )
    out = []
    for row_id, meta in row_meta.items():
        vs = verdicts.get(row_id, [])
        n = len(vs)
        counts = Counter(v for x in vs if (v := (x.get(vcol) or "").strip()))
        # deterministic winner: highest count, ties broken by lowest verdict string (never by agent order)
        top = max(counts.values(), default=0)
        winners = sorted(v for v, c in counts.items() if c == top)
        final = winners[0] if winners else "NO_VERDICT"
        if len(winners) > 1:
            final = f"TIE:{'/'.join(winners)}"  # no majority -> name the split, force adjudication
        unanimous = len(counts) == 1
        clean_unanimous = unanimous and final in clean and n == N_AGENTS
        rec = {
            "row_id": row_id,
            "ss_type": meta.get("ss_type", ""),
            "gene": meta.get("gene", ""),
            "key_id": meta.get(key_col, ""),
            "n_agents": n,
            "verdict_counts": ";".join(f"{v}:{c}" for v, c in counts.most_common()),
            "final_verdict": final,
            "unanimous": "yes" if unanimous else "no",
            "adjudicate": "" if clean_unanimous else "yes",
        }
        actors = [x for x in vs if (x.get(vcol) or "").strip() not in clean]
        for k, col in proposed.items():
            val, sup = modal([x.get(col, "") for x in actors])
            rec[f"proposed_{k}"] = val
            rec[f"{k}_support"] = f"{sup}/{len(actors)}" if val else ""
        rec["agent_notes"] = " || ".join(
            f"a{x['_agent_k']}:{(x.get(vcol) or '?')} {(x.get('notes') or '')[:140]}" for x in actors
        )
        out.append(rec)

    out.sort(key=lambda r: (r["adjudicate"] == "", r["ss_type"], r["row_id"]))  # adjudicate rows first
    write_tsv(sdir / "reconciliation.tsv", fields, out)

    tally = Counter(r["final_verdict"] for r in out)
    adj = [r["row_id"] for r in out if r["adjudicate"] == "yes"]
    print(f"[{sweep}] coverage:")
    for line in coverage:
        print("  " + line)
    if drift:
        print("ROW_ID DRIFT (re-run these agents):")
        for x in drift:
            print("  " + x)
    print(f"\nreconciliation.tsv: {len(out)} rows -> verdict tally {dict(tally)}")
    print(f"  needs adjudication (not unanimous-clean): {len(adj)}")
    for r in out:
        if r["adjudicate"] == "yes":
            extra = " ".join(f"{k}={r[f'proposed_{k}']}({r[f'{k}_support']})" for k in proposed if r[f"proposed_{k}"])
            print(f"    {r['row_id']:10} {r['gene']:18} [{r['verdict_counts']}] {extra}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
