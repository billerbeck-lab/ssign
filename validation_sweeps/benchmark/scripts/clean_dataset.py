#!/usr/bin/env python3
"""Shared filter so the benchmark figures reflect the citation-verified answer key + the T1SS fix.

After the full-table citation audit, positives_all.tsv holds the 337 citation-verified effectors, but
the runtime per-effector tables (actual_per_effector, ceiling_per_effector) still carry every effector
from the original (largely fabricated) answer key. This helper keeps only the verified effectors.

It also (optionally) flips the RTX T1SS toxins whose detection failed solely because they were staged
on a plasmid / partial WGS contig that lacked the chromosomal TolC (hlyA, apxIA, ltxA, lktA). The flag
is now SUPERSEDED by the CX3 full-assembly rerun (2026-06-25): re-running hlyA/apxIA/ltxA on their full
assemblies emits all three as secreted (DLP+DSE, DLP 0.87-0.98, nearby T1SS), so 4.4 writes the real
emission straight into actual_per_effector and the flag is left OFF. lktA P55117 (NZ_SMAM01000003.1)
was NOT in the fullasm rerun, but the other two lktA instances (P0C082, Q9ETX2) already emit from
complete genomes, so every RTX toxin type is shown to emit on a complete assembly. The flag stays as a
documented fallback only. Flipped rows are marked recovered_by_staging_fix=yes for transparency.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv  # noqa: E402

DATASET = Path(__file__).resolve().parents[1] / "data" / "dataset"
GOLD_BUILD = Path(__file__).resolve().parents[1] / "data" / "gold_build"
# effectors quarantined by the answer-key audit: round 1 (citation/dup/misclass) + round 2 (re-audit
# strict same-species-evidence bar). Both figure scripts (52, 54) must filter on the SAME set or the
# plotted recall and recall_figure_proteins.tsv silently disagree -- so it lives here, once.
AUDIT_REMOVED = [
    GOLD_BUILD / "effector_gold_set.removed_audit.tsv",
    GOLD_BUILD / "effector_gold_set.removed_reaudit2.tsv",
]
# hlyA, apxIA, ltxA, lktA: reachable@3 T1SS, full-assembly T1SS detection confirmed (script 50).
# DISABLED 2026-06-24 per Teo (recall reports STRICT emission, no draft-fragmentation rescue), and now
# RETIRED 2026-06-25: the CX3 fullasm rerun emits hlyA/apxIA/ltxA for real, so 4.4 writes their true
# emission into actual_per_effector instead. Kept only as a documented fallback (set =True to restore).
T1SS_STAGING_FIX = {"P08715", "P55128", "P16462", "P55117"}
APPLY_T1SS_STAGING_FIX = False


def _real(u: str) -> str:
    u = (u or "").strip()
    return u if u and u != "-" else ""


def _keep_keys():
    pos = read_tsv(DATASET / "positives_all.tsv")
    uni = {_real(r.get("uniprot")) for r in pos if _real(r.get("uniprot"))}
    loc = {(r.get("locus_tag") or "").strip() for r in pos if (r.get("locus_tag") or "").strip()}
    gss = {((r.get("gene") or "").strip().lower(), (r.get("ss_type") or "").strip()) for r in pos}
    return uni, loc, gss


def _verified(row: dict, keys) -> bool:
    uni, loc, gss = keys
    return (
        _real(row.get("uniprot")) in uni
        or (row.get("effector_locus") or "").strip() in loc
        or ((row.get("gene") or "").strip().lower(), (row.get("ss_type") or "").strip()) in gss
    )


def load_clean_actual(path) -> list[dict]:
    """actual_per_effector restricted to verified effectors, the 4 confirmed T1SS flipped to found."""
    keys = _keep_keys()
    rows = [r for r in read_tsv(path) if _verified(r, keys)]
    for r in rows:
        if APPLY_T1SS_STAGING_FIX and _real(r.get("uniprot")) in T1SS_STAGING_FIX:
            r["ssign_call"] = "emitted_secreted"
            r["recovered_by_staging_fix"] = "yes"
    return rows


def clean_ceiling(path) -> list[dict]:
    """ceiling_per_effector restricted to verified effectors."""
    keys = _keep_keys()
    return [r for r in read_tsv(path) if _verified(r, keys)]


def dropped_id() -> set:
    """(gene, uniprot) + (gene, locus_tag) pairs quarantined by the answer-key audit (both rounds).

    Keyed by genome identity, NOT (ss_type, gene): a gene can have several rows differing only by
    organism (cya in B. bronchiseptica vs B. pertussis; lktA across 4 hosts) where only some are
    dropped. Membership checks assume these pairs are disjoint from kept rows (true for the current
    gold set); a future audit dropping one host of a gene whose sibling shares an accession would
    need a finer key.
    """
    out = set()
    for f in AUDIT_REMOVED:
        if f.exists():
            for r in read_tsv(f):
                if r["uniprot"] and r["uniprot"] != "-":
                    out.add((r["gene"], r["uniprot"]))
                if r["locus_tag"]:
                    out.add((r["gene"], r["locus_tag"]))
    return out
