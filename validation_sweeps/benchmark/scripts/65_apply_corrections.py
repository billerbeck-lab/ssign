#!/usr/bin/env python3
"""Phase A task 1.6: apply the reconciled blind-review corrections to the gold list.

The dispositions below are the human-adjudicated outcome of the 4-agent blind review (reconciliation.tsv)
plus the tier-B re-anchor pass (2026-06-26). Each affected instance gets exactly one disposition:

  relabel    - organism free-text was wrong but genome+locus+coords are already correct -> swap the label only
  swap_up    - cited UniProt accession was deleted/wrong but the locus is the right gene -> swap accession only
  swap_ref   - primary_ref pointed at the wrong paper; the verbatim quote traces to another DOI -> swap ref
  note       - single-agent concern, kept on the 3/4 majority; record the alternative for the curator
  reanchor   - the answer key pointed at the WRONG gene/replicon. Supply the correct effector_locus (+ the
               verified UniProt/gene); the geometry (contig, coords, strand, nearest machinery + gene-
               distance, reachability, found-by-ssign) is RECOMPUTED from the same gene-order index and
               machinery table scripts 62/63 use -- no hand-typed numbers. found_by_ssign is read from the
               tier-2 rerun by coordinate join (RerunIndex), exactly as the T5SS rows in 63 are graded.
  drop       - removed: weak evidence, a deleted accession with no live successor, or a mis-anchored row
               that cannot be re-anchored coherently (substrate belongs to a different organism than the
               detected machinery instance, or no clean locus for the cited effector exists in the corpus).

Raw gold_list.tsv is never modified; the corrected list is written to gold_list_final.tsv, and the full
audit trail to gold_review/corrections.tsv. Re-running is idempotent. Edit the DISPOSITIONS table (not the
output) to change a call, then re-run.

Inputs : data/phase2/verification_phase_a/gold_list.tsv
         data/phase1/gene_order_index.tsv, data/machinery/machinery_resolved.tsv  (reanchor geometry)
         rerun/ + rerun_fullasm/                                                  (reanchor found-by-ssign)
Outputs: data/phase2/verification_phase_a/gold_list_final.tsv
         data/phase2/verification_phase_a/gold_review/corrections.tsv
Run    : .venv/bin/python scripts/65_apply_corrections.py
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VDIR = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VDIR / "gold_list.tsv"
FINAL = VDIR / "gold_list_final.tsv"
CORR = VDIR / "gold_review" / "corrections.tsv"
IDX = BENCH / "data" / "phase1" / "gene_order_index.tsv"
MR = BENCH / "data" / "machinery" / "machinery_resolved.tsv"
INF = 10**9

EPEC = "Escherichia coli O127:H6 (strain E2348/69 / EPEC)"
CROD = "Citrobacter rodentium (strain ICC168)"


def _norm(lt: str) -> str:
    return (lt or "").replace("_", "").replace(" ", "").upper()


class Reanchor:
    """Recompute a substrate's coordinates + machinery geometry from the gene-order index and the
    machinery table, so a corrected anchor is internally consistent (same gene-distance ssign's window
    uses; no hand-typed numbers)."""

    def __init__(self) -> None:
        self.idx: dict[str, list[tuple]] = defaultdict(list)  # norm(locus) -> [(rec, ordinal, start, end, strand)]
        for r in read_tsv(IDX):
            entry = (
                r["record_acc"],
                int(r["ordinal"]),
                int(r["start"]) if r.get("start") else None,
                int(r["end"]) if r.get("end") else None,
                r.get("strand", ""),
            )
            keys = {_norm(r["locus_tag"])} | {_norm(a) for a in (r.get("aliases", "") or "").split(";") if a.strip()}
            for k in keys:
                self.idx[k].append(entry)
        self.mach: dict[str, list[tuple]] = defaultdict(list)  # instance -> [(locus_tag, machinery_gene)]
        for m in read_tsv(MR):
            lt = (m.get("locus_tag") or "").strip()
            if lt:
                self.mach[m["instance_id"]].append((lt, m.get("gene", "")))
        try:
            from rerun_coords import RerunIndex  # noqa: PLC0415

            self.ridx = RerunIndex()
        except Exception:
            self.ridx = None

    def _coords(self, locus: str) -> tuple:
        for rec, ordn, s, e, strand in self.idx.get(_norm(locus), []):
            if s is not None:
                return rec, ordn, s, e, strand
        return ("", None, None, None, "")

    def geometry(self, iid: str, new_locus: str) -> dict:
        rec, ordn, s, e, strand = self._coords(new_locus)
        if not rec:
            raise ValueError(f"{iid}: re-anchor locus {new_locus} absent from gene-order index")
        best, bgene, bloc, machrecs = INF, "", "", set()
        for mlt, mgene in self.mach.get(iid, []):
            for mrec, mord, *_ in self.idx.get(_norm(mlt), []):
                machrecs.add(mrec)
                if mrec == rec and mord is not None and ordn is not None and abs(mord - ordn) < best:
                    best, bgene, bloc = abs(mord - ordn), mgene, mlt
        j = self.ridx.join(rec, s, e) if (self.ridx and s and e) else None
        return {
            "genome": rec,
            "contig": rec,
            "start": s,
            "stop": e,
            "strand": strand,
            "stage_replicons": ";".join(sorted({rec} | (machrecs - {""}))),
            "nearest_machinery_gene": bgene if best < INF else "",
            "nearest_machinery_locus": bloc if best < INF else "",
            "distance_to_machinery_genes": best if best < INF else "",
            "reachable_within_3": "yes" if best <= 3 else "no",
            "found_by_ssign": "yes" if (j and j["emitted"]) else "no",
        }


# instance_id -> (action, {field_or_identity_overrides}, agreement, basis)
DISPOSITIONS: dict[str, tuple] = {
    # --- relabel: organism free-text wrong, genome+locus+coords already correct ---
    "T1SS_04": (
        "relabel",
        {"organism": "Bordetella pertussis Tohama I"},
        "4/4 (re-resolved)",
        "BX470248.1 ORGANISM line = B. pertussis Tohama I; cya/BP0760/P0DKX7 + CyaB (BP0761, dist 1, found) "
        "all consistent on it. The hold wrongly assumed BX470248 was B. bronchiseptica -- it is a mislabel, "
        "not a mis-anchor.",
    ),
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
    # --- reanchor: answer key pointed at the wrong gene/replicon; supply the verified locus, recompute geometry ---
    "T2SS_05": (
        "reanchor",
        {"effector_locus": "LPG_RS11775", "uniprot": "Q5ZT22", "gene": "plaA"},
        "4/4 wrong_uniprot",
        "old LPG_RS04595/Q5ZX07 is ravI/lpg0926 (a Dot/Icm T4SS effector). plaA = LPG_RS11775/lpg2343/Q5ZT22 "
        "(lysophospholipase A) on the same NC_002942.5, far from the Lsp machinery -> not reachable.",
    ),
    "T4SS_08": (
        "reanchor",
        {"effector_locus": "BAB_RS19540", "uniprot": "Q2YN91", "gene": "BtpB"},
        "4/4 wrong_uniprot",
        "old BAB_RS19660/Q2YNA0 is BAB1_0782/DUF2069. BtpB = BAB1_0756/BAB_RS19540/Q2YN91 (TIR domain) on chr I "
        "NC_007618.1; the VirB machinery is on chr II -> cross-replicon, not reachable.",
    ),
    "T4SS_09": (
        "reanchor",
        {"effector_locus": "BAB_RS20990", "uniprot": "Q2YQ34", "gene": "VceC"},
        "4/4 wrong_genome",
        "old BAB_RS26945 (chr II) was the wrong gene. VceC = BR1038/BAB1_1058/BAB_RS20990/Q2YQ34 on chr I "
        "NC_007618.1; the VirB machinery is on chr II -> cross-replicon, not reachable.",
    ),
    # --- drop ---
    "T2SS_08": (
        "drop",
        {},
        "4/4 wrong_genome",
        "the detected instance is a Xanthomonas euvesicatoria Xps T2SS (NC_007508.1, XpsE/XpsF...); the cited "
        "substrate is a Vibrio cholerae lipase (P15493/VC_A0221). Substrate organism != machinery organism, so "
        "it cannot be re-anchored coherently to this instance.",
    ),
    "T3SS_13": (
        "drop",
        {},
        "4/4 wrong_uniprot",
        "PopC's cited accession A0ABY6NJB5 was deleted from UniProt; dropped per curator (do not substitute "
        "identifiers on dead-accession rows). A live Reviewed replacement exists if ever reinstated: "
        "PopC = RSp0875/RS_RS21320/Q9RBS2 (POPC_RALN1) on the GMI1000 megaplasmid NC_003296.1, adjacent to HrcC.",
    ),
    "T3SS_14": (
        "drop",
        {},
        "3/4 wrong_uniprot",
        "the listed locus RS_RS21300/RSp0871 is hrpD (a hrp-cluster gene), not RipJ -- the row was anchored onto "
        "the machinery (spurious dist 1). The deleted RipJ accession's only successor (A0ABX7ZTC2) is a different "
        "assembly (R. nicotianae, NZ_CP046674.1), not GMI1000; no clean GMI1000 RipJ locus exists in the corpus.",
    ),
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
    "reanchor": "corrected",
    "note": "confirmed_note",
}


def main() -> int:
    gold = read_tsv(GOLD)
    header = list(gold[0].keys())
    if "correction" not in header:
        header = header + ["correction"]
    engine = Reanchor()

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

        if action == "reanchor":
            old = {
                k: r.get(k, "")
                for k in (
                    "effector_locus",
                    "uniprot",
                    "contig",
                    "start",
                    "stop",
                    "reachable_within_3",
                    "found_by_ssign",
                )
            }
            merged = {**overrides, **engine.geometry(iid, overrides["effector_locus"])}
            for field, new in merged.items():
                r[field] = new
            r["verification_status"] = "corrected"
            r["correction"] = (
                f"reanchor: {basis} [{old['effector_locus']}@{old['contig']}:{old['start']} -> "
                f"{r['effector_locus']}@{r['contig']}:{r['start']}; reach {old['reachable_within_3']}->"
                f"{r['reachable_within_3']}, found {old['found_by_ssign']}->{r['found_by_ssign']}]"
            )
            corr_rows.append(
                {
                    "instance_id": iid,
                    "ss_type": r["ss_type"],
                    "gene": r["gene"],
                    "action": "reanchor",
                    "field": "anchor",
                    "old_value": f"{old['effector_locus']}/{old['uniprot']} @ {old['contig']}:{old['start']}-{old['stop']} "
                    f"(reach {old['reachable_within_3']}, found {old['found_by_ssign']})",
                    "new_value": f"{r['effector_locus']}/{r['uniprot']} @ {r['contig']}:{r['start']}-{r['stop']} "
                    f"(reach {r['reachable_within_3']}, found {r['found_by_ssign']})",
                    "agreement": agreement,
                    "basis": basis,
                }
            )
            out.append(r)
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
        if action == "note":
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
    reanchored = [i for i, d in DISPOSITIONS.items() if d[0] == "reanchor"]
    by_status = Counter(r["verification_status"] for r in out)
    by_type = Counter(r["ss_type"] for r in out)
    print(f"gold_list_final.tsv: {len(out)} rows kept ({len(gold)} - {len(dropped)} dropped)")
    print(f"  dropped: {dropped}")
    print(f"  reanchored (geometry recomputed): {reanchored}")
    print(f"  by verification_status: {dict(by_status)}")
    print(f"  by ss_type (kept): {dict(sorted(by_type.items()))}")
    print(f"  SCORABLE rows (no holds remain): {len(out)}")
    print(f"corrections.tsv: {len(corr_rows)} audit rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
