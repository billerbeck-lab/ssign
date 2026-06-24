#!/usr/bin/env python3
"""Side-by-side sheet: known experimental annotation vs ssign's predicted annotation.

Scope is recall-gated. ssign only runs the annotation tools (InterProScan, EggNOG,
pLM-BLAST/ECOD) on proteins it EMITS as secreted, so the only corpus effectors we
can grade are the ones ssign both found and emitted: the 51 `emitted_secreted` rows
in the Phase-2 panel. For each we line up:

  known  = what the experimental literature says it is (family + gene + a verbatim
           citation quote, from the curated corpus positives_all.tsv)
  ssign  = each annotation tool's free-text call (interpro / eggnog / pLM-BLAST-ECOD
           / Pfam), read from that genome's results_raw row at the effector's locus.

No keyword/LLM scoring: the output is a TSV for manual eyeballing. We also print a
yield summary (how often each tool produced ANY annotation).

    .venv/bin/python annotation_accuracy_sheet.py

Genome resolution: the panel carries ssign_locus but not which fleet genome it came
from, so we map locus_tag -> genome by scanning all 67 results_raw. PAO1 appears
twice in the fleet (AE004091 = INSDC, NC_002516.2 = RefSeq, same genome), so a
handful of PA loci resolve to both; we break the tie toward the corpus's stated
refseq_genome and fall back deterministically (the annotations are identical anyway).
"""

from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import clean_dataset  # noqa: E402
from bench_index import accession_base, load_tsv  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.abspath(os.path.join(HERE, ".."))
PANEL = os.path.join(BENCH, "data", "phase2", "actual_per_effector.panel_genbank_t3ss.tsv")
POS = os.path.join(BENCH, "data", "dataset", "positives_all.tsv")
FLEET = "/tmp/ssign_fleet_67"
OUT = os.path.join(HERE, "annotation_accuracy_sheet.tsv")
UNIPROT_CACHE = os.path.join(HERE, "annotation_uniprot_cache.json")

# annotation free-text columns in results_raw, one per tool
ANN_COLS = {
    "ssign_interpro": "interpro_descriptions",
    "ssign_eggnog": "eggnog_description",
    "ssign_plmblast_ecod": "ecod_top1_description",
    "ssign_pfam": "pfam_ids",
}
OUT_FIELDS = [
    "uniprot",
    "gene",
    "organism",
    "ss_type",
    # known/ground-truth annotation: UniProt is the richest (covers ~all entries);
    # corpus family/quote are sparse fallbacks (family only ~11/51).
    "known_uniprot_name",
    "known_uniprot_families",
    "known_uniprot_function",
    "known_family",
    "known_quote",
    *ANN_COLS.keys(),
    "genome",
    "ssign_locus",
]


def fetch_uniprot_functions(accessions: set[str]) -> dict:
    """uniprot accession -> {name, families, function}, cached on disk.

    Queries the UniProt REST accessions endpoint once for any not already cached, so
    reruns are offline. Accessions UniProt doesn't return (obsolete/demerged) are
    cached as blanks so they aren't re-queried."""
    import json
    import urllib.parse
    import urllib.request

    cache = json.load(open(UNIPROT_CACHE)) if os.path.exists(UNIPROT_CACHE) else {}
    missing = sorted(a for a in accessions if a and a != "-" and a not in cache)
    for i in range(0, len(missing), 100):
        chunk = missing[i : i + 100]
        params = urllib.parse.urlencode(
            {
                "accessions": ",".join(chunk),
                "fields": "accession,protein_name,protein_families,cc_function",
                "format": "tsv",
            }
        )
        url = f"https://rest.uniprot.org/uniprotkb/accessions?{params}"
        with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310 (fixed UniProt host)
            lines = resp.read().decode().splitlines()
        for r in csv.DictReader(lines, delimiter="\t"):
            func = (r.get("Function [CC]", "") or "").replace("FUNCTION: ", "").strip()
            cache[r.get("Entry", "")] = {
                "name": r.get("Protein names", ""),
                "families": r.get("Protein families", ""),
                "function": func[:300] + ("..." if len(func) > 300 else ""),
            }
        for a in chunk:
            cache.setdefault(a, {"name": "", "families": "", "function": ""})
    json.dump(cache, open(UNIPROT_CACHE, "w"), indent=0)
    return cache


def build_locus_index(loci: set[str]):
    """locus_tag -> {genome: results_raw_row} for every fleet genome carrying it."""
    idx: dict[str, dict[str, dict]] = {}
    for d in sorted(os.listdir(FLEET)):
        raw = os.path.join(FLEET, d, "results", f"{d}_results_raw.csv")
        if not os.path.exists(raw):
            continue
        with open(raw) as fh:
            for row in csv.DictReader(fh):
                lt = row.get("locus_tag", "")
                if lt in loci:
                    idx.setdefault(lt, {})[d] = row
    return idx


def pick_genome(candidates: dict[str, dict], refseq_genome: str) -> str:
    """Choose one genome among the candidates carrying this locus.

    Prefer an exact refseq_genome match, then an accession-base match (folds
    version/prefix drift), then the first sorted name (deterministic; the PAO1
    duplicate has identical annotations either way)."""
    names = sorted(candidates)
    if refseq_genome:
        if refseq_genome in candidates:
            return refseq_genome
        base = accession_base(refseq_genome)
        for n in names:
            if accession_base(n) == base:
                return n
    return names[0]


def match_corpus(e: dict, pos: list[dict]) -> dict:
    """Find this panel effector's corpus row (its known experimental annotation).

    The panel is built on GenBank assemblies, so its locus tags live in a
    different namespace than the corpus's RefSeq placement — we can't join on
    locus. Match on UniProt accession when the effector has one (unique), else
    on gene name. (The old code keyed a dict by uniprot, which collapsed all
    117 accession-less corpus rows onto a single entry and broadcast one wrong
    annotation across every UniProt-less effector.) Disambiguate multi-gene
    hits by ss_type, then by genome.
    """
    u = e.get("uniprot", "").strip()
    cands: list[dict] = []
    if u and u != "-":
        cands = [r for r in pos if r.get("uniprot", "").strip() == u]
    if not cands:
        g = e.get("gene", "").strip().lower()
        if g:
            cands = [r for r in pos if r.get("gene", "").strip().lower() == g]
    if not cands:
        return {}
    same_ss = [r for r in cands if r.get("ss_type", "") == e.get("ss_type", "")]
    pool = same_ss or cands
    gb = accession_base(e.get("unit_id", ""))
    same_genome = [r for r in pool if gb and accession_base(r.get("refseq_genome", "")) == gb]
    return (same_genome or pool)[0]


# SecReT6 external-DB entries that leaked into the emitted set and are NOT backed by the audited
# corpus. Two are gene-name collisions with the corpus's T3SS Bop effectors (a different
# B. thailandensis T6SS-locus protein matched the corpus BopA/BopE by gene name); two have no
# primary reference (Tle3 paralogs). Graded ground truth must come from the audited corpus, so
# these are dropped. (EFF00142, a duplicate of TseA, is removed separately by dedup_by_locus.)
EXCLUDE_EXTERNAL = {
    ("BopA", "BTH_RS04545"),
    ("BopE", "BTH_RS04535"),
    ("EFF00136", "BTH_RS00455"),
    ("EFF00150", "BTH_RS28610"),
    # Broken ground truth: Q3BPB1 is a deleted UniProt entry and the corpus 'Serralysin' family
    # conflicts with the unanimous S8-subtilisin tool calls, so it cannot be graded.
    ("prtA", "IS_RS19530"),
}


def dedup_by_locus(rows: list[dict]) -> list[dict]:
    """One emitted protein is graded once. Two corpus/external entries can point at the same
    ssign_locus (e.g. SecReT6 'EFF00142' and corpus 'TseA_T6SS1' are both BTH_RS25940); keep the
    one with a UniProt accession / non-EFF gene name so the named ground truth wins."""

    def rank(r):
        has_uni = bool((r.get("uniprot") or "").strip() and r["uniprot"] != "-")
        return (has_uni, not (r.get("gene") or "").startswith("EFF"))

    best, order = {}, []
    for r in rows:
        k = (r.get("ssign_locus") or r.get("effector_locus") or "").strip() or f"_{len(order)}"
        if k not in best:
            best[k] = r
            order.append(k)
        elif rank(r) > rank(best[k]):
            best[k] = r
    return [best[k] for k in order]


def main():
    panel = dedup_by_locus([r for r in load_tsv(PANEL) if r["ssign_call"] == "emitted_secreted"])
    panel = [r for r in panel if (r.get("gene", ""), r.get("ssign_locus", "")) not in EXCLUDE_EXTERNAL]
    # Drop effectors quarantined by the answer-key audit (e.g. redundant strain-duplicates cya Q57506,
    # lktA Q9ETX2) so the annotation grade reflects the same verified set as the recall figures.
    dropped = clean_dataset.dropped_id()
    panel = [
        r
        for r in panel
        if (r.get("gene", ""), r.get("uniprot", "")) not in dropped
        and (r.get("gene", ""), r.get("effector_locus", "")) not in dropped
    ]
    pos = load_tsv(POS)
    loci = {r["ssign_locus"] for r in panel if r["ssign_locus"]}
    idx = build_locus_index(loci)
    uniprot = fetch_uniprot_functions({e["uniprot"] for e in panel})

    out_rows = []
    for e in panel:
        u = e["uniprot"]
        p = match_corpus(e, pos)
        up = uniprot.get(u, {})
        locus = e["ssign_locus"]
        cand = idx.get(locus, {})
        genome = pick_genome(cand, p.get("refseq_genome", "")) if cand else ""
        raw = cand.get(genome, {})

        row = {
            "uniprot": u,
            "gene": e.get("gene", "") or p.get("gene", ""),
            "organism": p.get("organism", ""),
            "ss_type": e.get("ss_type", ""),
            "known_uniprot_name": up.get("name", ""),
            "known_uniprot_families": up.get("families", ""),
            "known_uniprot_function": up.get("function", ""),
            "known_family": p.get("family", ""),
            # citation_quote is the verbatim experimental sentence; fall back to `quote`
            "known_quote": p.get("citation_quote") or p.get("quote", ""),
            "genome": genome,
            "ssign_locus": locus,
        }
        for out_col, raw_col in ANN_COLS.items():
            row[out_col] = (raw.get(raw_col, "") or "").strip()
        out_rows.append(row)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    # --- yield summary ---
    n = len(out_rows)
    print(f"{n} emitted_secreted effectors -> {OUT}\n")
    print("ground-truth coverage:")
    print(f"  known_uniprot_name present:     {sum(1 for r in out_rows if r['known_uniprot_name'])}/{n}")
    print(f"  known_uniprot_function present: {sum(1 for r in out_rows if r['known_uniprot_function'])}/{n}")
    print(f"  known_family (corpus) present:  {sum(1 for r in out_rows if r['known_family'])}/{n}")
    print(f"  known_quote (corpus) present:   {sum(1 for r in out_rows if r['known_quote'])}/{n}")
    print("\nssign annotation yield (produced ANY call):")
    for c in ANN_COLS:
        got = sum(1 for r in out_rows if r[c])
        print(f"  {c:22s}: {got}/{n}  ({got / n:.0%})")
    any_ann = sum(1 for r in out_rows if any(r[c] for c in ANN_COLS))
    print(f"\n  at least one tool annotated: {any_ann}/{n}  ({any_ann / n:.0%})")
    none_ann = [r for r in out_rows if not any(r[c] for c in ANN_COLS)]
    if none_ann:
        print(f"  NO annotation from any tool ({len(none_ann)}): " + ", ".join(r["uniprot"] for r in none_ann))


if __name__ == "__main__":
    main()
