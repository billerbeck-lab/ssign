#!/usr/bin/env python3
"""Phase A task 1.5: reconcile the 4 blind-agent verdicts per gold-list row into one final verdict.

Each gold-list instance was audited by 4 independent agents (gold_review/verdicts/agent{1..4}_<batch>.tsv).
This merges them: per row_id it tallies the four verdicts, derives a final verdict + an agreement label,
carries the agents' proposed field corrections (modal new_uniprot/new_locus/new_primary_ref among the
agents that voted `fix`), and marks rows that need a human to adjudicate.

Final-verdict rule (n = agents that voted on the row):
  - any `drop`            -> `drop` if unanimous, else ADJUDICATE  (a drop is serious; never auto-resolve a split)
  - else `unclear` plurality -> ADJUDICATE
  - else any `fix`        -> `fix` (effector stays, a field is corrected); agreement unanimous/partial
  - else all `confirmed`  -> `confirmed`

Coverage is reported, not assumed: a batch with <4 verdict files, or an agent whose row_id set drifts from
its batch input, is flagged so it can be re-run before the corrections table (1.6) is built.

Inputs : data/phase2/verification_phase_a/gold_review/batches/<batch>.tsv         (canonical row metadata)
         data/phase2/verification_phase_a/gold_review/verdicts/agent{k}_<batch>.tsv
Outputs: data/phase2/verification_phase_a/gold_review/reconciliation.tsv          (one row per gold-list row)
Run    : .venv/bin/python scripts/64_reconcile_gold_review.py
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
GR = BENCH / "data" / "phase2" / "verification_phase_a" / "gold_review"
BATCHES = GR / "batches"
VERDICTS = GR / "verdicts"
OUT = GR / "reconciliation.tsv"

VERDICT_VALUES = ("confirmed", "fix", "drop", "unclear")
FIELDS = [
    "row_id",
    "ss_type",
    "gene",
    "uniprot",
    "n_agents",
    "n_confirmed",
    "n_fix",
    "n_drop",
    "n_unclear",
    "final_verdict",
    "agreement",
    "adjudicate",
    "defects",
    "proposed_new_uniprot",
    "proposed_new_locus_tag",
    "proposed_new_primary_ref",
    "fix_support",
    "agent_notes",
]


def _modal(values: list[str]) -> tuple[str, int]:
    """Most common non-empty value and how many agents proposed it (blank, 0 if none)."""
    vals = [v.strip() for v in values if (v or "").strip()]
    if not vals:
        return "", 0
    val, n = Counter(vals).most_common(1)[0]
    return val, n


def _final(c: int, f: int, d: int, u: int, n: int) -> tuple[str, str]:
    if d:
        return ("drop", "unanimous") if d == n else ("ADJUDICATE", f"split-drop({d}/{n})")
    if u and u >= max(c, f):
        return ("ADJUDICATE", f"unclear({u}/{n})")
    if f:
        return ("fix", "unanimous" if f == n else f"partial({f}/{n})")
    return ("confirmed", "unanimous" if c == n else f"partial({c}/{n})")


def main() -> int:
    batch_files = sorted(BATCHES.glob("*.tsv"))
    if not batch_files:
        print(f"no batch files in {BATCHES}", file=sys.stderr)
        return 1

    # canonical row metadata + expected row_id set per batch
    batch_rows: dict[str, list[dict]] = {}
    row_meta: dict[str, dict] = {}
    for bf in batch_files:
        rows = read_tsv(bf)
        batch_rows[bf.stem] = rows
        for r in rows:
            row_meta[r["row_id"]] = r

    # collect verdicts per row_id; track coverage + row_id drift
    verdicts: dict[str, list[dict]] = defaultdict(list)
    coverage: list[str] = []
    drift: list[str] = []
    for batch, rows in batch_rows.items():
        expected = {r["row_id"] for r in rows}
        present = 0
        for k in (1, 2, 3, 4):
            vf = VERDICTS / f"agent{k}_{batch}.tsv"
            if not vf.exists():
                continue
            present += 1
            vrows = read_tsv(vf)
            got = {r["row_id"] for r in vrows}
            if got != expected:
                drift.append(f"agent{k}_{batch}: missing={sorted(expected - got)} extra={sorted(got - expected)}")
            for r in vrows:
                if r["row_id"] in expected:
                    verdicts[r["row_id"]].append(r)
        coverage.append(f"{batch}: {present}/4 agent files" + ("" if present == 4 else "  <-- INCOMPLETE"))

    out = []
    for row_id, meta in row_meta.items():
        vs = verdicts.get(row_id, [])
        n = len(vs)
        counts = {v: sum(1 for x in vs if (x.get("verdict") or "").strip() == v) for v in VERDICT_VALUES}
        c, f, d, u = (counts[v] for v in VERDICT_VALUES)
        final, agreement = _final(c, f, d, u, n) if n else ("NO_VERDICT", "missing")

        fixers = [x for x in vs if (x.get("verdict") or "").strip() in ("fix", "drop")]
        nu, nu_s = _modal([x.get("new_uniprot", "") for x in fixers])
        nl, nl_s = _modal([x.get("new_locus_tag", "") for x in fixers])
        nr, nr_s = _modal([x.get("new_primary_ref", "") for x in fixers])
        defects = ";".join(sorted({(x.get("defect") or "").strip() for x in vs if (x.get("defect") or "").strip()}))
        fix_support = ",".join(
            s
            for s in (
                f"uniprot {nu_s}/{len(fixers)}" if nu else "",
                f"locus {nl_s}/{len(fixers)}" if nl else "",
                f"ref {nr_s}/{len(fixers)}" if nr else "",
            )
            if s
        )
        # one compact line per non-confirmed agent so a human can adjudicate from this file alone
        notes = " || ".join(
            f"a{i + 1}:{(x.get('verdict') or '?')}={(x.get('defect') or '').strip()} {(x.get('notes') or '')[:120]}"
            for i, x in enumerate(vs)
            if (x.get("verdict") or "").strip() != "confirmed"
        )

        out.append(
            {
                "row_id": row_id,
                "ss_type": meta.get("ss_type", ""),
                "gene": meta.get("gene", ""),
                "uniprot": meta.get("uniprot", ""),
                "n_agents": n,
                "n_confirmed": c,
                "n_fix": f,
                "n_drop": d,
                "n_unclear": u,
                "final_verdict": final,
                "agreement": agreement,
                "adjudicate": "yes" if (final == "ADJUDICATE" or n < 4) else "",
                "defects": defects,
                "proposed_new_uniprot": nu,
                "proposed_new_locus_tag": nl,
                "proposed_new_primary_ref": nr,
                "fix_support": fix_support,
                "agent_notes": notes,
            }
        )

    out.sort(key=lambda r: (r["ss_type"], r["row_id"]))
    write_tsv(OUT, FIELDS, out)

    tally = Counter(r["final_verdict"] for r in out)
    adj = [r["row_id"] for r in out if r["adjudicate"] == "yes"]
    print("coverage:")
    for line in coverage:
        print("  " + line)
    if drift:
        print("ROW_ID DRIFT (re-run these agents):")
        for x in drift:
            print("  " + x)
    print(f"\nreconciliation.tsv: {len(out)} rows -> {dict(tally)}")
    print(f"  needs adjudication (split/incomplete): {len(adj)} -> {adj}")
    fixes = [r for r in out if r["final_verdict"] == "fix"]
    print(f"  fix (field correction, effector kept): {len(fixes)} -> {[r['row_id'] for r in fixes]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
