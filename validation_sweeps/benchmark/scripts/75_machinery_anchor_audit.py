#!/usr/bin/env python3
"""Re-check every gold row for anchoring to MACHINERY instead of the effector (Teo, 2026-06-30).

Four rows so far were found anchored to a secretion-system component, not the substrate: T6SS_17 (TagF
regulator), T6SS_18 (TagL/OmpA accessory), T3SS_23 (YscL apparatus), T6SS_21 (TagJ accessory). Those were
caught one at a time. This sweeps all 88 rows at once by asking, per row, what the protein AT THE ROW'S
LOCUS is actually ANNOTATED as -- using two independent sources:
  - the RefSeq gene name + aliases from data/phase1/gene_order_index.tsv (effector_locus -> gene), and
  - the Bakta product of the rerun ORF that overlaps the row's coordinates.
A row is flagged MACHINERY if either annotation matches an apparatus/accessory/structural/regulatory class
(the things a SECRETED SUBSTRATE is never annotated as). High-precision term list -- effector-y enzyme names
(toxin, amidase, phospholipase, nuclease, ADP-ribosyl...) are NOT machinery and are left alone; VgrG/Hcp are
NOT flagged because evolved VgrG/Hcp can be genuine effectors (eyeball those by hand if they appear).

Also flags PLACEHOLDER_CLUSTER: >=2 rows on one genome whose locus ordinals are tightly sequential (the
pattern that exposed the source-corpus Ysa block YE3556-3564), which suggests made-up consecutive loci.

Read-only; prints flagged rows + writes machinery_anchor_audit.tsv.
Run: .venv/bin/python scripts/75_machinery_anchor_audit.py
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402
from rerun_coords import RerunIndex, _raw_path, _unit_dirs  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
VD = BENCH / "data" / "phase2" / "verification_phase_a"
GOLD = VD / "gold_list_final.tsv"
IDX = BENCH / "data" / "phase1" / "gene_order_index.tsv"
OUT = VD / "machinery_anchor_audit.tsv"

# Annotation substrings that mark a NON-substrate (machinery / structural / accessory / regulatory). A secreted
# effector is never annotated with these. Kept high-precision; matched case-insensitively as substrings.
MACHINERY = (
    "apparatus",
    "accessory protein",
    "-associated protein",
    "associated protein",
    "baseplate",
    "membrane subunit",
    "tube protein",
    "spike protein",
    "contractile sheath",
    "sheath protein",
    "pilotin",
    "sorting platform",
    "secretion system protein",
    "secretion protein",
    "core component",
    "structural",
    "inner membrane",
    "translocase",
    "usher",
    "outer membrane protein assembly",
    "bama",
    "omp85",
    "ompa",
    "hrpe",
    "yscl",
    "yscn",
    "sctn",
    "sctl",
    "sctj",
    "sctv",
    "sctc",
    "sctu",
    "sctr",
    "sctt",
    "sccs",
    "tssa",
    "tssb",
    "tssc",
    "tsse",
    "tssf",
    "tssg",
    "tssj",
    "tssk",
    "tssl",
    "tssm",
    "clpv",
    "tagf",
    "tagj",
    "tagl",
    "tagh",
    "icmf",
    "dotg",
    "doth",
    "dotu",
    "virb",
    "vird4",
    "trbi",
    "pilus",
    "pilin",
    "fimbri",
    "flagell",
    "tad ",
    "response regulator",
    "sensor kinase",
    "sensor histidine",
    "two-component",
    "transcriptional regulator",
    "anti-sigma",
    "chaperone",
    "atp synthase",
    "atpase",
    "lytic transglycosylase",
)
# never flag these even if a machinery term substring-matches: evolved spike/tube effectors + the autotransporter
# barrel that is part of the substrate itself.
EFFECTOR_SAFE = ("vgrg", "hcp", "autotransporter", "rhs", "wapa")


def machinery_hits(text: str) -> list[str]:
    t = text.lower()
    if any(s in t for s in EFFECTOR_SAFE):
        return []
    return [m for m in MACHINERY if m in t]


def main() -> int:
    gold = read_tsv(GOLD)
    # RefSeq gene/aliases per locus_tag
    idx_gene: dict[str, str] = {}
    for r in read_tsv(IDX):
        idx_gene[r["locus_tag"]] = " ".join([r.get("gene", ""), r.get("aliases", "")])
    ridx = RerunIndex()
    unit_path = {u.name: u for u in _unit_dirs()}  # unit name -> Path (handles rerun/ and rerun_fullasm/)
    raw_cache: dict[str, dict] = {}

    def rerun_product(contig: str, start: int, stop: int, strand: str) -> tuple[str, str]:
        eo = ridx.emitted_overlap(contig, start, stop, strand=strand)
        if not eo or not eo["best_locus"]:
            return "", ""
        unit = eo["unit"]
        if unit not in raw_cache:
            up = unit_path.get(unit)
            raw_cache[unit] = {r["locus_tag"]: r for r in csv.DictReader(open(_raw_path(up)))} if up else {}
        row = raw_cache[unit].get(eo["best_locus"], {})
        return row.get("product", ""), eo["best_locus"]

    out = []
    by_genome: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for r in gold:
        loc = r["effector_locus"].strip()
        refseq = idx_gene.get(loc, "")
        product, rerun_locus = rerun_product(r["contig"].strip(), int(r["start"]), int(r["stop"]), r.get("strand", ""))
        hits = machinery_hits(refseq) + machinery_hits(product)
        out.append(
            {
                "instance_id": r["instance_id"],
                "ss_type": r["ss_type"],
                "gene": r["gene"],
                "effector_locus": loc,
                "refseq_gene": refseq.strip(),
                "rerun_product": product,
                "flags": "MACHINERY:" + "/".join(sorted(set(hits))) if hits else "ok",
            }
        )
        m = re.search(r"(\d+)", loc)
        if m:
            by_genome[r["contig"].strip()].append((int(m.group(1)), r["instance_id"]))

    # placeholder cluster: same genome, >=2 rows with locus numbers within 6 of each other
    cluster = []
    for contig, items in by_genome.items():
        items.sort()
        for i in range(len(items) - 1):
            if 0 < items[i + 1][0] - items[i][0] <= 6:
                cluster.append((contig, items[i][1], items[i + 1][1], items[i][0], items[i + 1][0]))

    write_tsv(
        OUT,
        ["instance_id", "ss_type", "gene", "effector_locus", "refseq_gene", "rerun_product", "flags"],
        sorted(out, key=lambda x: (x["flags"] == "ok", x["instance_id"])),
    )
    flagged = [x for x in out if x["flags"] != "ok"]
    print(f"machinery_anchor_audit.tsv: {len(out)} rows, {len(flagged)} flagged MACHINERY")
    for x in flagged:
        print(f"  {x['instance_id']:9} {x['gene']:14} {x['flags']}")
        print(f"        refseq={x['refseq_gene']!r}  rerun_product={x['rerun_product']!r}")
    print(f"\nplaceholder clusters (same genome, near-sequential loci): {len(cluster)}")
    for c in cluster:
        print(f"  {c[1]} & {c[2]} on {c[0]} (loci ...{c[3]} & ...{c[4]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
