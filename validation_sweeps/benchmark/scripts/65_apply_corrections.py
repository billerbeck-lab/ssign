#!/usr/bin/env python3
"""Phase A task 1.6: apply the reconciled blind-review corrections to the gold list.

The dispositions below are the human-adjudicated outcome of the 4-agent blind review (reconciliation.tsv).
Each affected instance gets exactly one disposition:

  relabel    - organism free-text was wrong but genome+locus+coords are already correct -> swap the label only
  swap_up    - cited UniProt accession was deleted/wrong but the locus is the right gene -> swap accession only
  swap_ref   - primary_ref pointed at the wrong paper; the verbatim quote traces to another DOI -> swap ref
  note       - single-agent concern, kept on the 3/4 majority; record the alternative for the curator
  hold       - the answer key is MIS-ANCHORED (locus/genome/coords point at the wrong gene or replicon);
               can't be cosmetically patched, needs coordinate re-derivation or a drop decision. Kept in the
               file but verification_status=hold_reanchor so it is EXCLUDED from scoring until resolved.
  drop       - removed from the gold list (user-confirmed: weak evidence, or deleted accession w/ no successor)

Raw gold_list.tsv is never modified; the corrected list is written to gold_list_final.tsv, and the full
audit trail to gold_review/corrections.tsv. Re-running is idempotent. Edit the DISPOSITIONS table (not the
output) to change a call, then re-run.

Inputs : data/phase2/verification_phase_a/gold_list.tsv
Outputs: data/phase2/verification_phase_a/gold_list_final.tsv
         data/phase2/verification_phase_a/gold_review/corrections.tsv
Run    : .venv/bin/python scripts/65_apply_corrections.py
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VDIR = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VDIR / "gold_list.tsv"
FINAL = VDIR / "gold_list_final.tsv"
CORR = VDIR / "gold_review" / "corrections.tsv"

EPEC = "Escherichia coli O127:H6 (strain E2348/69 / EPEC)"
CROD = "Citrobacter rodentium (strain ICC168)"

# instance_id -> (action, {field: new_value}, agreement, basis)
DISPOSITIONS: dict[str, tuple] = {
    # --- relabel: organism free-text wrong, genome+locus+coords already correct ---
    "T1SS_06": (
        "relabel",
        {"organism": "Photorhabdus laumondii subsp. laumondii TT01"},
        "4/4",
        "genome BX470251 + locus plu0655 are Photorhabdus TT01; organism said Dickeya chrysanthemi",
    ),
    "T1SS_07": (
        "relabel",
        {"organism": "Pseudomonas entomophila (strain L48)"},
        "4/4",
        "genome CT573326 + locus PSEEN1550 are P. entomophila L48; organism said P. aeruginosa PAO1",
    ),
    "T3SS_25": (
        "relabel",
        {"organism": EPEC},
        "4/4",
        "genome NC_011601.1 + locus E2348C_ are EPEC; organism said O157:H7 Sakai",
    ),
    "T3SS_26": (
        "relabel",
        {"organism": EPEC},
        "4/4",
        "genome NC_011601.1 + locus E2348C_ + uniprot B7UI22 are EPEC; organism said O157:H7 Sakai",
    ),
    "T3SS_27": (
        "relabel",
        {"organism": EPEC},
        "4/4",
        "genome NC_011601.1 + locus E2348C_ are EPEC; organism said O157:H7 Sakai",
    ),
    "T3SS_28": (
        "relabel",
        {"organism": CROD},
        "4/4",
        "genome NC_013716.1 + locus ROD_ are C. rodentium ICC168; organism said O157:H7 Sakai",
    ),
    "T3SS_30": (
        "relabel",
        {"organism": CROD},
        "4/4",
        "genome NC_013716.1 + locus ROD_ are C. rodentium ICC168; organism said O157:H7 Sakai",
    ),
    "T3SS_32": (
        "relabel",
        {"organism": CROD},
        "4/4",
        "genome NC_013716.1 + locus ROD_ are C. rodentium ICC168; organism said O157:H7 Sakai",
    ),
    # --- swap_up: deleted/wrong accession, locus is the right gene, coords correct ---
    "T3SS_17": (
        "swap_up",
        {"uniprot": "Q87P32"},
        "4/4",
        "A0AAW8Q5J6 deleted; Q87P32 = VopS, VP1686 exact strain match",
    ),
    "T3SS_21": (
        "swap_up",
        {"uniprot": "Q63K35"},
        "4/4",
        "A0A8A4DTS9 deleted; Q63K35 = BipC/BPSS1531 (chr2 NC_006351.1, matches)",
    ),
    # --- swap_ref: cited paper wrong; verbatim quote traces to another DOI ---
    "T6SS_03": (
        "swap_ref",
        {"primary_ref": "10.1128/mBio.00075-15"},
        "3/4",
        "DOI 10.1101/868539 is the Aux3 MGE preprint; ledger quote is verbatim from Altindis 2015 mBio (PMID 25759499)",
    ),
    # --- note: single-agent concern, kept on 3/4 majority; alt recorded for curator ---
    "T5SS_03": (
        "note",
        {},
        "1/4 fix",
        "kept Q51163 (3/4 confirm, N. meningitidis IgA protease); agent1 alt = Q9K0B4/NMB0700 (clean MC58 entry)",
    ),
    "T4SS_05": (
        "note",
        {},
        "1/4 fix",
        "kept P08061 (3/4 confirm VirE3 C58); agent4 flags P08061 is an 83aa fragment, full-length = A0A2Z2PK47",
    ),
    "T5SS_13": (
        "note",
        {},
        "1/4 fix",
        "kept ref (3/4 confirm it documents LspA1 supernatant secretion); agent1 alt = 10.1128/IAI.72.4.1874-1884.2004",
    ),
    "T5SS_14": (
        "note",
        {},
        "1/4 fix",
        "kept ref (3/4 confirm it covers LspA2 secretion); agent1 alt = 10.1128/IAI.72.4.1874-1884.2004",
    ),
    # --- hold: answer key mis-anchored (wrong gene/genome/replicon); excluded from scoring until re-derived ---
    "T1SS_04": (
        "hold",
        {},
        "4/4 wrong_genome",
        "locus BP0760 + uniprot P0DKX7 are B. pertussis; staged genome BX470248 is B. bronchiseptica RB50",
    ),
    "T2SS_05": (
        "hold",
        {},
        "4/4 wrong_uniprot",
        "coords/locus LPG_RS04595 are lpg0926/ravI; plaA is lpg2837/Q5ZRP3 elsewhere on NC_002942.5",
    ),
    "T2SS_08": (
        "hold",
        {},
        "4/4 wrong_genome",
        "genome NC_007508.1 is Xanthomonas; V. cholerae lipase = P15493/VC_A0221 on NC_002506.1",
    ),
    "T3SS_13": (
        "hold",
        {},
        "4/4 wrong_uniprot",
        "A0ABY6NJB5 deleted; PopC/Q9RBS2 is on the GMI1000 megaplasmid NC_003296.1/RSp0875, not chr locus RS_RS03050",
    ),
    "T3SS_14": (
        "hold",
        {},
        "3/4 wrong_uniprot",
        "A0ABY6NGA6 deleted; only replacement A0ABX7ZTC2 is a different-strain (R. nicotianae) RipJ, not GMI1000",
    ),
    "T4SS_08": (
        "hold",
        {},
        "4/4 wrong_uniprot",
        "Q2YNA0 = BAB1_0782/DUF2069; BtpB = Q2YN91/BAB1_0756, so locus BAB_RS19660 likely the wrong gene too",
    ),
    "T4SS_09": (
        "hold",
        {},
        "1/4 wrong_genome",
        "VceC = BAB1_1058 on chr I NC_007618.1; listed on chr II NC_007624.1 (locus BAB_RS26945)",
    ),
    # --- drop: user-confirmed ---
    "T4SS_03": (
        "drop",
        {},
        "1/4 drop",
        "weak evidence: CBU0937 only listed in a PPI map, no translocation assay; 'CirD' name unsupported",
    ),
    "T2SS_10": ("drop", {}, "4/4 wrong_uniprot", "A1JQR7 deleted, no live successor accession"),
    "T2SS_12": ("drop", {}, "3/3 wrong_uniprot", "A6V5Q1 deleted, no live successor accession"),
    "T6SS_12": ("drop", {}, "4/4 wrong_uniprot", "A0A0H3B0Q3 deleted, no live successor accession"),
}

STATUS = {
    "relabel": "corrected",
    "swap_up": "corrected",
    "swap_ref": "corrected",
    "note": "confirmed_note",
    "hold": "hold_reanchor",
}


def main() -> int:
    gold = read_tsv(GOLD)
    header = list(gold[0].keys())
    if "correction" not in header:
        header = header + ["correction"]

    out, corr_rows = [], []
    for r in gold:
        iid = r["instance_id"]
        disp = DISPOSITIONS.get(iid)
        if not disp:
            r["verification_status"] = "confirmed"
            r["correction"] = ""
            out.append(r)
            continue
        action, overrides, agreement, basis = disp
        if action == "drop":
            corr_rows.append(
                {
                    "instance_id": iid,
                    "ss_type": r["ss_type"],
                    "gene": r["gene"],
                    "action": "drop",
                    "field": "",
                    "old_value": "",
                    "new_value": "",
                    "agreement": agreement,
                    "basis": basis,
                }
            )
            continue
        changed = []
        for field, new in overrides.items():
            old = r.get(field, "")
            if old != new:
                corr_rows.append(
                    {
                        "instance_id": iid,
                        "ss_type": r["ss_type"],
                        "gene": r["gene"],
                        "action": action,
                        "field": field,
                        "old_value": old,
                        "new_value": new,
                        "agreement": agreement,
                        "basis": basis,
                    }
                )
                r[field] = new
                changed.append(f"{field}:{old}->{new}")
        if action in ("note", "hold"):
            corr_rows.append(
                {
                    "instance_id": iid,
                    "ss_type": r["ss_type"],
                    "gene": r["gene"],
                    "action": action,
                    "field": "",
                    "old_value": "",
                    "new_value": "",
                    "agreement": agreement,
                    "basis": basis,
                }
            )
        r["verification_status"] = STATUS[action]
        r["correction"] = f"{action}: {basis}" + (f" [{'; '.join(changed)}]" if changed else "")
        out.append(r)

    write_tsv(FINAL, header, out)
    write_tsv(
        CORR,
        ["instance_id", "ss_type", "gene", "action", "field", "old_value", "new_value", "agreement", "basis"],
        corr_rows,
    )

    dropped = [i for i, d in DISPOSITIONS.items() if d[0] == "drop"]
    held = [i for i, d in DISPOSITIONS.items() if d[0] == "hold"]
    by_status = Counter(r["verification_status"] for r in out)
    by_type = Counter(r["ss_type"] for r in out)
    scorable = [r for r in out if r["verification_status"] != "hold_reanchor"]
    print(f"gold_list_final.tsv: {len(out)} rows kept ({len(gold)} - {len(dropped)} dropped)")
    print(f"  dropped: {dropped}")
    print(f"  held (excluded from scoring until re-anchored): {held}")
    print(f"  by verification_status: {dict(by_status)}")
    print(f"  by ss_type (kept): {dict(sorted(by_type.items()))}")
    print(f"  SCORABLE rows (confirmed + corrected + noted, excludes held): {len(scorable)}")
    print(f"corrections.tsv: {len(corr_rows)} audit rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
