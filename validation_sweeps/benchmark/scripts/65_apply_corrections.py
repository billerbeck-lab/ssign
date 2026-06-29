#!/usr/bin/env python3
"""Phase A task 1.6: apply the reconciled blind-review corrections to the gold list.

The dispositions below are the human-adjudicated outcome of the 4-agent blind review (reconciliation.tsv)
plus the tier-B re-anchor pass (2026-06-26). Each affected instance gets exactly one disposition:

  relabel    - organism free-text was wrong but genome+locus+coords are already correct -> swap the label only
  swap_up    - cited UniProt accession was deleted/wrong but the locus is the right gene -> swap accession only
  swap_ref   - primary_ref pointed at the wrong paper; the verbatim quote traces to another DOI -> swap ref
  fix_quote  - ref is correct but citation_quote was garbled/non-verbatim/described function not secretion ->
               swap in a verbatim secretion sentence from the same paper (sweep-2 citation pass)
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
# scripts/73 output: recomputed geometry for the RTX-T1SS/T5SS/(reanchored T6SS) rows machinery_resolved.tsv
# can't resolve. Absent -> the geometry-override pass is a no-op. Run order converges as 65 -> 73 -> 65: this
# script first reanchors T6SS_17/18 onto their new coords, scripts/73 then reads those coords to recompute
# their geometry, and a second 65 run applies it. All three steps are idempotent.
GEOM_RC = VDIR / "t1ss_t5ss_geometry_recompute.tsv"
INF = 10**9

EPEC = "Escherichia coli O127:H6 (strain E2348/69 / EPEC)"
CROD = "Citrobacter rodentium (strain ICC168)"
# Verbatim from Ward et al. 2004 Infect Immun abstract (PMID 15039306) -- documents LspA1/LspA2 supernatant secretion.
LSPA_QUOTE = (
    "The LspA1 and LspA2 proteins of Haemophilus ducreyi 35000 are two very large macromolecules that can be "
    "detected in concentrated culture supernatant fluid."
)
# Verbatim from Bowen et al. 2003 Microbiology (PMID 12777498) -- shared by T1SS_06 (P. laumondii TT01) and
# T1SS_R12 (P. luminescens W14), the same conserved Photorhabdus prtA RTX metalloprotease.
PRTA_RTX_QUOTE = (
    "One of the secreted proteins is PrtA, which is shown here to be a repeats-in-toxin (RTX) alkaline zinc "
    "metalloprotease."
)


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
        {
            "organism": "Photorhabdus laumondii subsp. laumondii TT01",
            "primary_ref": "10.1099/mic.0.26171-0",
            "citation_quote": PRTA_RTX_QUOTE,
        },
        "4/4 + citation 3/3",
        "genome BX470251 + locus plu0655 are Photorhabdus TT01 (organism said Dickeya chrysanthemi). Citation sweep: "
        "old ref Akatsuka 1997 is the SERRATIA PrtA heterologous-secretion paper (wrong-organism quote); swapped to "
        "Bowen 2003 Microbiology (PMID 12777498), the Photorhabdus prtA secretion paper.",
    ),
    "T1SS_07": (
        "relabel",
        {"organism": "Pseudomonas entomophila (strain L48)"},
        "4/4",
        "genome CT573326 + locus PSEEN1550 are P. entomophila L48; organism said P. aeruginosa PAO1",
    ),
    "T3SS_25": (
        "relabel",
        {"organism": EPEC, "uniprot": "B7UNX6"},
        "4/4 + identity 3/3",
        "genome NC_011601.1 + locus E2348C_ are EPEC; organism said O157:H7 Sakai. Identity sweep: "
        "blank accession -> B7UNX6 (NleD homolog, E2348/69, 232aa = the CDS span).",
    ),
    "T3SS_26": (
        "relabel",
        {"organism": EPEC},
        "4/4",
        "genome NC_011601.1 + locus E2348C_ + uniprot B7UI22 are EPEC; organism said O157:H7 Sakai",
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
    # --- second-pass upgrades (2026-06-26): on direct review the 4 "note" rows each resolved to a concrete fix ---
    "T5SS_03": (
        "swap_up",
        {
            "uniprot": "Q9K0B4",
            "citation_quote": "The precursor contains three functional domains, the amino-terminal leader which is assumed to initiate the inner membrane transport of the precursor, the protease, and a carboxyl-terminal 'helper' domain apparently required for extracellular secretion (excretion).",
        },
        "agent3 + re-review + citation 2/3",
        "Q51163 is a 496aa Fragment (strain unspecified); Q9K0B4 = full-length 1815aa iga/NMB0700, strain MC58, "
        "RefSeq WP_002225480.1 on NC_003112.2 -- matches the gold genome + the 5447bp span. primary_ref Pohlner "
        "1987 kept (foundational autotransporter-secretion paper); citation sweep swapped to its verbatim "
        "extracellular-secretion sentence.",
    ),
    "T4SS_05": (
        "reanchor",
        {"effector_locus": "ATU_RS23905", "uniprot": "A0A2Z2PK47", "gene": "VirE3"},
        "agent4 + re-review",
        "old ATU_RS23890 = Atu6188/P08061, an 83aa 'VirE 9 kDa' ORF mislabelled virE3 (not the effector). The "
        "VirB/D4-translocated VirE3 effector = Atu6191/ATU_RS23905/A0A2Z2PK47 (684aa) on pTiC58 NC_003065.3; "
        "~6 genes from VirD4 -> not reachable.",
    ),
    "T5SS_13": (
        "swap_ref",
        {"primary_ref": "10.1128/IAI.72.4.1874-1884.2004", "citation_quote": LSPA_QUOTE},
        "agent2 + re-review",
        "kept ref documented antiphagocytic function, not secretion. Ward 2004 IAI (PMID 15039306) directly "
        "documents LspB-mediated TPS secretion of LspA1/LspA2 to the supernatant. Q7VLE8 (35000HP lspA1, 4152aa = "
        "the span) unchanged.",
    ),
    "T5SS_14": (
        "swap_ref",
        {"uniprot": "Q9ZHL0", "primary_ref": "10.1128/IAI.72.4.1874-1884.2004", "citation_quote": LSPA_QUOTE},
        "agent4 + re-review",
        "uniprot was blank; full-length LspA2 = Q9ZHL0 (4919aa = the 14760bp span). Same ref swap to the Ward 2004 "
        "LspB-secretion paper.",
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
        {
            "effector_locus": "BAB_RS20990",
            "uniprot": "Q2YQ34",
            "gene": "VceC",
            "citation_quote": "This analysis resulted in the identification of the proteins encoded by BAB1_1652 (VceA) and BR1038/BAB1_1058 (VceC) as novel protein substrates of the Brucella T4SS.",
        },
        "4/4 wrong_genome + citation 2/2",
        "old BAB_RS26945 (chr II) was the wrong gene. VceC = BR1038/BAB1_1058/BAB_RS20990/Q2YQ34 on chr I "
        "NC_007618.1; the VirB machinery is on chr II -> cross-replicon, not reachable. Citation sweep also fixed the "
        "quote (old had dropped 'Brucella' -> 'theT4SS'; de Jong 2008 Mol Microbiol).",
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
    # ===== second-pass IDENTITY sweep (2026-06-28, gold_review2/sweep1_identity, 3 agents/batch) =====
    # accession ADD: the gold locus IS the named gene (UniProt gene/protein name confirms it); agents found
    # the genome-matched accession whose length equals the CDS span. Deterministically re-validated by scripts/69.
    "T3SS_01": (
        "swap_up",
        {"uniprot": "O84091"},
        "identity 3/3",
        "CopN = CT_089/O84091 (gene lcrE, 421aa = span), Chlamydia trachomatis D.",
    ),
    "T3SS_03": (
        "swap_up",
        {"uniprot": "Q7DB68"},
        "identity 3/3",
        "EspZ = ECs_4571/Q7DB68 (99aa = span), E. coli O157:H7 Sakai.",
    ),
    "T3SS_11": (
        "swap_up",
        {"uniprot": "Q8ZQ60"},
        "identity 3/3",
        "PipA = STM1087/Q8ZQ60 (SPI5 protein, 226aa ~ span), S. Typhimurium LT2.",
    ),
    "T4SS_02": (
        "swap_up",
        {
            "uniprot": "Q5ZYD5",
            "citation_quote": "Legionella pneumophila utilizes the Dot/Icm type IVB secretion system to deliver hundreds of effector proteins inside eukaryotic cells to ensure intracellular replication.",
        },
        "identity 3/3 + citation 2/3",
        "Ceg14 = LPG_RS02190/Q5ZYD5 (Dot/Icm substrate, 666aa = span), L. pneumophila Philadelphia 1. Citation sweep: "
        "old quote was an effector-function list; swapped to the verbatim Dot/Icm-delivery sentence (Patel 2024 Mol Syst Biol).",
    ),
    "T5SS_06": (
        "swap_up",
        {"uniprot": "Q2SV12"},
        "identity 3/3",
        "bcpA = BTH_RS26080/Q2SV12 (filamentous haemagglutinin, 3147aa = span), B. thailandensis E264.",
    ),
    "T5SS_08": (
        "swap_up",
        {"uniprot": "Q9HVG6"},
        "identity 3/3",
        "cdrA = PA4625/Q9HVG6 (TPS FhaB/CdiA family, 2154aa = span), P. aeruginosa PAO1.",
    ),
    "T5SS_10": (
        "swap_up",
        {
            "uniprot": "Q48031",
            "citation_quote": "HMW1B is located in the outer membrane and serves to translocate HMW1 across the outer membrane.",
        },
        "identity 3/3 + citation 2/3",
        "hmw1A = R2846_RS03620/Q48031 (HMW1 adhesin, 1536aa = span), H. influenzae R2846. Citation sweep swapped to the "
        "verbatim HMW1-translocation sentence (St Geme & Grass 1998 Mol Microbiol).",
    ),
    "T6SS_03": (
        "swap_up",
        {
            "uniprot": "Q9KMN9",
            "primary_ref": "10.1128/mBio.00075-15",
            "citation_quote": "We show that the product of the gene VCA0285 is likely a new peptidoglycan hydrolase that is secreted by T6SS and that its cognate immunity protein is encoded by the gene that is immediately downstream (VCA0286).",
        },
        "identity 3/3 + citation 3/3",
        "TseH = VC_RS14690/Q9KMN9 (T6SS effector TseH, 223aa = span), V. cholerae N16961. Citation sweep: cited "
        "10.1101/868539 is the Aux3 MGE preprint; the quote traces to Altindis 2015 mBio (10.1128/mBio.00075-15). This "
        "ref swap had been LOST when the sweep-1 swap_up entry shadowed the first-pass swap_ref (dup key) -- now consolidated.",
    ),
    "T6SS_08": (
        "swap_up",
        {
            "uniprot": "Q7CUP8",
            "citation_quote": "Although Tde1 and Tae were readily secreted when bacteria were grown on acidic AB-MES agar plate",
        },
        "identity 3/3 + citation (tie -> documented)",
        "Tae = ATU_RS20360/Q7CUP8 (T6SS amidase effector, 166aa = span), A. tumefaciens C58. Citation sweep: old quote "
        "was a Tae-Tai classification statement; swapped to the verbatim secretion observation from the same ref "
        "(Ma 2014 Cell Host Microbe).",
    ),
    "T6SS_09": (
        "swap_up",
        {"uniprot": "Q2T422", "organism": "Burkholderia thailandensis E264"},
        "identity 3/3",
        "TseM = BTH_RS09640/Q2T422 (141aa = span); organism field had carried the contig id NC_007650.",
    ),
    "T6SS_10": (
        "swap_up",
        {"uniprot": "Q2SV36", "organism": "Burkholderia thailandensis E264"},
        "identity 3/3",
        "Tle1 = BTH_RS32260/Q2SV36 (T6SS phospholipase Tle1, 800aa = span); organism field had the contig id NC_007651.",
    ),
    "T6SS_17": (
        "swap_up",
        {"uniprot": "A0A0F6AX88"},
        "identity 3/3",
        "Tae4 = STM14_RS02020/A0A0F6AX88 (262aa = span), S. Typhimurium 14028s.",
    ),
    "T6SS_18": (
        "swap_up",
        {"uniprot": "D3GUV9"},
        "identity 3/3",
        "Tle1 = EC042_RS24190/D3GUV9 (T6SS Tle1, 576aa = span), E. coli 042.",
    ),
    # organism-only fix: no genome-matched UniProt entry exists for this strain protein (accession kept blank),
    # but the organism field had carried the contig id instead of the species name.
    "T6SS_13": (
        "relabel",
        {"organism": "Yersinia pseudotuberculosis YPIII"},
        "identity 3/3 none_exists",
        "YezP: no genome-matched UniProt accession exists (kept blank); organism field had the contig id NC_010465.1.",
    ),
    # MIS-ANCHORED (coord/gene mismatch): the gold effector_locus encodes a DIFFERENT gene (transposase /
    # hydrolase / regulator); reanchor to the true effector locus (length confirmed against the index + UniProt).
    "T3SS_16": (
        "reanchor",
        {"effector_locus": "PSPTO_RS07380", "uniprot": "G3XDB3", "gene": "HrpK1"},
        "identity 3/3 (coord/gene mismatch)",
        "old locus PSPTO_RS07250 is hrpS (a T3SS transcriptional regulator), NOT HrpK1. True HrpK1 = "
        "PSPTO_1405/PSPTO_RS07380/G3XDB3 (780aa) on NC_004578.1.",
    ),
    "T3SS_27": (
        "reanchor",
        {"effector_locus": "E2348C_RS07680", "uniprot": "B7UR63", "gene": "NleF", "organism": EPEC},
        "identity 3/3 (coord/gene mismatch)",
        "old locus E2348C_RS09145 is ydgJ (a predicted oxidoreductase), NOT NleF. True NleF = "
        "E2348C_1445/E2348C_RS07680/B7UR63 (189aa) on NC_011601.1; organism also said O157:H7 Sakai.",
    ),
    "T3SS_28": (
        "reanchor",
        {"effector_locus": "ROD_RS14675", "uniprot": "D2TKE3", "gene": "EspA", "organism": CROD},
        "identity 3/3 (coord/gene mismatch)",
        "old locus ROD_RS14860 is an HTH-domain protein, NOT EspA. True EspA = ROD_29761/ROD_RS14675/"
        "D2TKE3 (192aa) on NC_013716.1; organism also said O157:H7 Sakai.",
    ),
    "T3SS_29": (
        "reanchor",
        {"effector_locus": "ROD_RS01735", "uniprot": "D2TJZ4", "gene": "NleA_EspI", "organism": CROD},
        "identity 3/3 (coord/gene mismatch)",
        "old locus ROD_RS06085 is an IS110 transposase, NOT NleA. True NleA/EspI = ROD_03401/ROD_RS01735/"
        "D2TJZ4 (430aa) on NC_013716.1; replaces the length-mismatched Q6T6T6.",
    ),
    "T3SS_30": (
        "reanchor",
        {"effector_locus": "ROD_RS02910", "uniprot": "D2TML3", "gene": "NleD", "organism": CROD},
        "identity 3/3 (coord/gene mismatch)",
        "old locus ROD_RS03135 is an alpha/beta hydrolase, NOT NleD. True NleD = ROD_05791/ROD_RS02910/"
        "D2TML3 (235aa) on NC_013716.1; organism also said O157:H7 Sakai.",
    ),
    "T3SS_31": (
        "reanchor",
        {"effector_locus": "ROD_RS05475", "uniprot": "D2TT37", "gene": "NleB", "organism": CROD},
        "identity 3/3 (coord/gene mismatch)",
        "old locus ROD_RS14635 is an IS3 transposase (spuriously dist 5 from EscF, found=yes). True NleB1 = "
        "ROD_10831/ROD_RS05475/D2TT37 (329aa) on NC_013716.1.",
    ),
    "T3SS_32": (
        "reanchor",
        {"effector_locus": "ROD_RS24095", "uniprot": "D2TRY1", "gene": "EspJ", "organism": CROD},
        "identity 3/3 (coord/gene mismatch)",
        "old locus ROD_RS12795 is a lytic transglycosylase, NOT EspJ. True EspJ = ROD_48901/ROD_RS24095/"
        "D2TRY1 (217aa) on NC_013716.1; organism also said O157:H7 Sakai.",
    ),
    # ===== second-pass CITATION sweep (2026-06-28, gold_review2/sweep2_citation, 3 agents/batch) =====
    # quote-only: ref is correct but citation_quote was garbled / non-verbatim / described function not secretion.
    # Agents supplied a verbatim secretion sentence from the SAME ref (modal support noted; >=2/3 unless flagged).
    "T1SS_R12": (
        "fix_quote",
        {"citation_quote": PRTA_RTX_QUOTE},
        "citation 3/3",
        "Bowen 2003 Microbiology (PMID 12777498), Photorhabdus prtA1; old quote described the prt locus structure, not secretion.",
    ),
    "T3SS_20": (
        "fix_quote",
        {
            "citation_quote": "Shigella deliver a subset of effectors into the host cell via the type III secretion system, that stimulate host cell signal pathways to modulate the actin dynamics required for invasion of epithelial cells."
        },
        "citation 2/3",
        "Yoshida 2002 EMBO J (PMID 12065406); old quote was a VirA effector-function sentence, not the T3SS-delivery statement.",
    ),
    "T4SS_07": (
        "fix_quote",
        {
            "citation_quote": "These Bartonella-translocated effector proteins (Beps) A-G are encoded together with the VirB system and the T4S coupling protein VirD4 on a Bartonella-specific pathogenicity island."
        },
        "citation 3/3",
        "Schulein 2005 PNAS (PMID 15642951); old quote carried a '[Bartonella]' placeholder -> restored verbatim.",
    ),
    "T4SS_10": (
        "fix_quote",
        {
            "citation_quote": "A hypothetical protein was identified as a putative effector, hereby named Anaplasma translocated substrate 1 (Ats-1)."
        },
        "citation 2/3",
        "Niu 2010 PLoS Pathog (PMID 20174550); old quote was garbled ('namedranslocatedubstrate 1') -> restored verbatim.",
    ),
    "T4SS_11": (
        "fix_quote",
        {
            "citation_quote": "ECH0825 was translocated from the bacterium into the host-cell cytoplasm and localized to mitochondria."
        },
        "citation 3/3",
        "Liu 2012 Cell Microbiol (PMID 22348527), Etf-1=ECH0825; citation_quote was EMPTY -> filled with the verbatim translocation sentence.",
    ),
    "T5SS_04": (
        "fix_quote",
        {
            "citation_quote": "We have identified and characterized a secreted protein, designated Pic, which is encoded on the chromosomes of enteroaggregative Escherichia coli (EAEC) 042 and Shigella flexneri 2457T."
        },
        "citation 3/3",
        "Henderson 1999 Infect Immun (PMC96930); old quote had dropped 'Escherichia coli'/'Shigella flexneri' -> restored verbatim.",
    ),
    "T5SS_05": (
        "fix_quote",
        {
            "citation_quote": "Sat, the secreted autotransporter toxin of uropathogenic Escherichia coli, is a vacuolating cytotoxin for bladder and kidney epithelial cells."
        },
        "citation 2/3",
        "Guyer 2002 Infect Immun (PMC128167); old quote had the typo 'uropathogenicCFT073' (missing 'E. coli') -> verbatim title sentence.",
    ),
    "T5SS_09": (
        "fix_quote",
        {
            "citation_quote": "Although partly surface associated, it is also very efficiently secreted into the extracellular milieu. Its secretion depends on the outer membrane accessory protein FhaC."
        },
        "citation 2/3",
        "Jacob-Dubuisson 1996 Mol Microbiol; old quote was the FhaB preproprotein/adhesin-family sentence, not the secretion statement.",
    ),
    "T5SS_11": (
        "fix_quote",
        {"citation_quote": "HpmB is necessary for the extracellular secretion and activation of HpmA"},
        "citation (corrects cross-wired quote)",
        "Uphoff & Welch 1990 J Bacteriol (PMC208585); old quote referenced FhaC (a Bordetella TpsB) on this Proteus HpmA "
        "row -> swapped to the verbatim HpmB-secretes-HpmA TPS sentence.",
    ),
    # ref swaps: cited paper did not document secretion (or wrong organism). Swapped to a verified, organism-matched
    # secretion paper + a verbatim sentence pulled from its abstract (Europe PMC, confirmed this session).
    "T3SS_08": (
        "swap_ref",
        {
            "primary_ref": "10.1371/journal.ppat.1000067",
            "citation_quote": "This is due to the activity of a Type III secreted effector protein, designated YopJ in Y. pseudotuberculosis and Y. pestis, and YopP in the closely related Y. enterocolitica.",
        },
        "citation fix_ref:2/3",
        "old Orth 2000 Science documents YopJ enzymatic mechanism + homology, no secretion assay. Brodsky 2008 PLoS "
        "Pathog (DOI 10.1371/journal.ppat.1000067) directly addresses YopJ Type-III secretion; Y. pestis matches.",
    ),
    "T1SS_R06": (
        "swap_ref",
        {
            "primary_ref": "10.1002/j.1460-2075.1991.tb04890.x",
            "citation_quote": "The alternative secretion pathway which exports hemolysin across both Escherichia coli membranes into the surrounding medium is directed by an uncleaved C-terminal targeting signal and the membrane translocator proteins HlyD and HlyB.",
        },
        "citation fix_ref:3/3",
        "old Felmlee 1985 is the hemolysin DNA-sequence/transcription paper (quote was a molecular-mass list). Koronakis "
        "1991 EMBO J directly documents HlyB/HlyD-dependent HlyA secretion across both E. coli membranes.",
    ),
    "T1SS_R11": (
        "swap_ref",
        {
            "primary_ref": "10.1007/BF00279652",
            "citation_quote": "Erwinia chrysanthemi, a phytopathogenic enterobacterium, secretes three proteases (PrtA, PrtB and PrtC) into the extracellular medium.",
        },
        "citation fix_ref:1/3 (verified)",
        "old Akatsuka 1997 documents SERRATIA PrtA secreted via the E. chrysanthemi exporter (wrong-organism quote for "
        "this Erwinia/Dickeya PrtA row). Ghigo 1992 MGG is the organism-matched E. chrysanthemi PrtA secretion paper; "
        "quote verified verbatim via Europe PMC.",
    ),
    # weak-but-valid: protein is a bona fide substrate; the cited quote describes structure/function, not a direct
    # secretion assay, and no better organism-matched secretion sentence was found. Kept -> confirmed_note.
    "T1SS_R15": (
        "note",
        {},
        "citation weak:3/3 (kept)",
        "RtxA is a bona fide RTX/T1SS toxin (Kingella kingae); Kehl-Fie & St Geme 2007 establishes the rtx locus + "
        "cytotoxicity. Quote describes the mutagenesis screen, not a supernatant assay -- kept as a valid T1SS substrate.",
    ),
    "T5SS_02": (
        "note",
        {},
        "citation fix_quote:3/3 (current quote kept)",
        "Ag43 (flu) is a canonical T5aSS autotransporter; the current quote ('surface-displayed autotransporter protein "
        "that mediates...') is a valid autotransport statement, and the agents' proposed replacement was a weaker "
        "structure/homology sentence -- so the current quote is kept.",
    ),
    "T5SS_16": (
        "note",
        {},
        "citation weak:3/3 (kept)",
        "NadA is a canonical trimeric autotransporter adhesin (T5cSS); its passenger is translocated across the OM by its "
        "own C-terminal barrel. Comanducci 2002 quote describes oligomeric structure -- kept as a valid T5cSS substrate.",
    ),
    "T5SS_17": (
        "note",
        {},
        "citation weak:3/3 (kept)",
        "YadA is the prototypical trimeric autotransporter adhesin (T5cSS); Hoiczyk 2000 quote describes the head/stalk/"
        "anchor architecture -- kept as a valid T5cSS substrate.",
    ),
    # ===== third-pass HIGH-CHURN re-review (2026-06-29, blind 3-agent + my deterministic Pfam check) =====
    # Two rows were anchored onto the WRONG gene by an earlier swap_up. Confirmed by Pfam domain (authoritative,
    # not the agent vote): the cited accession's Pfam is a structural/regulatory domain, not the effector's.
    # Re-anchored to the real effector-immunity gene (located in the rerun + RefSeq product + UniProt) so geometry
    # and found_by_ssign are recomputed against the correct CDS.
    "T6SS_17": (
        "reanchor",
        {"effector_locus": "STM14_RS01975", "uniprot": "A0A0F6AX79", "gene": "Tae4_Stm"},
        "id re-review 1/3 + Pfam (authoritative)",
        "old STM14_RS02020/A0A0F6AX88 (262aa) carries Pfam PF09867 TagF_N -- a T6SS regulator, NOT the Tae4 amidase "
        "effector. Real Tae4 = STM14_RS01975/STM14_0325/A0A0F6AX79 (161aa), RefSeq 'T6SS amidase effector Tae4', "
        "beside its Tai4 immunity gene (STM14_RS01980) in SPI-6, NC_016856.1:318368-318853.",
    ),
    "T6SS_18": (
        "reanchor",
        {"effector_locus": "EC042_RS24220", "uniprot": "B7LG63", "gene": "Tle1_Sci1"},
        "id re-review 3/3 + Pfam (authoritative)",
        "old EC042_RS24190/D3GUV9 (576aa) carries Pfam PF00691 OmpA -- the TagL PG-binding accessory protein (PDB "
        "5M38/7BBA), NOT the Tle1 phospholipase. Real Tle1 = EC042_RS24220/EC042_4534/B7LG63 (560aa, Pfam PF09994 "
        "Tle1-like_cat), RefSeq 'T6SS effector phospholipase Tle1-EAEC', beside its two Tli1 immunity genes in sci-1, "
        "NC_017626.1:4865605-4867287.",
    ),
}

STATUS = {
    "relabel": "corrected",
    "swap_up": "corrected",
    "swap_ref": "corrected",
    "fix_quote": "corrected",
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

    # found_by_ssign is a MEASURED quantity, not a curated label: recompute it for every kept row from
    # its current coordinates via the canonical any-overlap rule (RerunIndex.emitted_overlap), so the
    # column can't go stale (e.g. the 3 RTX toxins the rerun_fullasm re-run only later emitted) or stay
    # over-credited (a T6SS effector whose ORF Bakta called but ssign never emitted). The contig is
    # molecule-reconciled (RefSeq<->INSDC) inside emitted_overlap; None = un-reconcilable -> keep stored.
    n_found_fixed = 0
    for r in out:
        if engine.ridx is None:
            break
        try:
            eo = engine.ridx.emitted_overlap(
                r["contig"].strip(), int(r["start"]), int(r["stop"]), strand=r.get("strand")
            )
        except (ValueError, KeyError):
            eo = None
        if eo is None or eo["found"] == r.get("found_by_ssign", ""):
            continue
        corr_rows.append(
            {
                "instance_id": r["instance_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "action": "found_recompute",
                "field": "found_by_ssign",
                "old_value": r.get("found_by_ssign", ""),
                "new_value": eo["found"],
                "agreement": "any-overlap (scripts/71)",
                "basis": f"{eo['reason']}: n_overlap={eo['n_overlap']} bp={eo['best_overlap_bp']} emit={eo['emitted_locus'] or '-'}",
            }
        )
        r["found_by_ssign"] = eo["found"]
        r["correction"] = (
            f"{r['correction']}; found_recompute->{eo['found']}"
            if r.get("correction")
            else f"found_recompute->{eo['found']}"
        )
        n_found_fixed += 1

    # geometry override for the RTX-T1SS + T5SS rows (scripts/73): machinery_resolved.tsv has no instance for
    # these types, so the Reanchor engine can't set their nearest/distance/reachable. scripts/73 derives them
    # from the secretion-system components TXSScan actually DETECTED in the rerun (gene-order distance to the
    # nearest T1SS component / the T5bSS TpsB translocator; self-secreting autotransporters get distance 0,
    # reachable "self"). This corrects the prior T5bSS mislabel ("(self-secreting)" -> real TpsB distance, so
    # fhaB/lspA1 flip reachable yes->no). It does NOT touch found_by_ssign (recall is unchanged).
    GEOM_FIELDS = {
        "nearest_machinery_gene": "rc_nearest",
        "nearest_machinery_locus": "rc_locus",
        "distance_to_machinery_genes": "rc_dist",
        "reachable_within_3": "rc_reach",
    }
    n_geom_fixed = 0
    geom_rc = {x["instance_id"]: x for x in read_tsv(GEOM_RC)} if GEOM_RC.exists() else {}
    APPLICABLE = ("recomputed", "self_secreting", "self_no_model")  # allow-list: every other status
    for r in out:  # (contig_absent / effector_not_in_order /
        g = geom_rc.get(r["instance_id"])  # no_*_component_on_contig) blanks the recompute,
        if not g or g["status"] not in APPLICABLE:  # so keep the stored value instead of overwriting it
            continue
        changed = [(f, r.get(f, ""), g[src]) for f, src in GEOM_FIELDS.items() if r.get(f, "") != g[src]]
        if not changed:
            continue
        for f, old, new in changed:
            r[f] = new
            corr_rows.append(
                {
                    "instance_id": r["instance_id"],
                    "ss_type": r["ss_type"],
                    "gene": r["gene"],
                    "action": "geometry_recompute",
                    "field": f,
                    "old_value": old,
                    "new_value": new,
                    "agreement": "rerun-detected components (scripts/73)",
                    "basis": g["status"],
                }
            )
        tag = f"geometry_recompute({g['status']})"
        r["correction"] = f"{r['correction']}; {tag}" if r.get("correction") else tag
        n_geom_fixed += 1

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
    print(f"  found_by_ssign recomputed (any-overlap): {n_found_fixed} value(s) corrected vs raw")
    print(f"  geometry recomputed (RTX-T1SS/T5SS, scripts/73): {n_geom_fixed} row(s) corrected")
    print(f"  SCORABLE rows (no holds remain): {len(out)}")
    print(f"corrections.tsv: {len(corr_rows)} audit rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
