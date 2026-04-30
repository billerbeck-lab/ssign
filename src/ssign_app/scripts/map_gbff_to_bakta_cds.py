#!/usr/bin/env python3
"""Map original GenBank CDS annotations onto Bakta's fresh CDS predictions.

Phase 3.3 re-annotates every input with Bakta by default, including
GenBank files. That gives uniform Bakta-produced annotations across a
cohort, but the user's original GenBank annotations (from whatever
pipeline produced their file) would otherwise be lost — Bakta assigns
its own locus tags and may shift CDS boundaries a few bp in either
direction, so a naive locus-tag join loses the mapping.

This script bridges the two: for each Bakta CDS, find the GenBank CDS
with the largest reciprocal coordinate overlap (above a threshold),
and carry the GenBank's `product` across as `gbff_annotation`. Bakta
CDS with no sufficiently-overlapping GenBank CDS get an empty
`gbff_annotation`. The result feeds annotation-consensus voting as a
third opinion alongside Bakta's own product and EggNOG's description.

Input TSVs expected:
    --bakta-gene-info       Bakta's gene_info.tsv (from run_bakta.py)
    --genbank-gene-info     Original GenBank's gene_info.tsv (from
                            extract_proteins.py's --out-gene-info)

Both have columns locus_tag, contig, start (0-based), end, strand, product.

Output: the Bakta gene_info.tsv with an added `gbff_annotation` column.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Default minimum reciprocal overlap. For two CDS (a, b) with lengths
# La, Lb and overlap O, reciprocal overlap is min(O/La, O/Lb).
# 0.8 = "at least 80% of each CDS falls inside the other" — tight enough
# to reject chance alignments, loose enough to absorb the few-bp shifts
# common across gene callers.
_DEFAULT_MIN_OVERLAP = 0.8


def _read_gene_info(path: str):
    """Read a gene_info.tsv into a dict[contig] -> list[full row dicts].

    Preserves every column in the input TSV so Bakta's rich annotation
    fields (ec_numbers, cog_ids, go_terms, kegg_ko, refseq_ids, pfam_ids,
    gene, protein_id, ...) survive the mapping step. start/end are
    coerced to int for the overlap arithmetic.
    """
    by_contig = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for raw in reader:
            row = dict(raw)
            contig = (row.get("contig") or "").strip()
            if not contig:
                continue
            try:
                row["start"] = int(row.get("start", 0))
                row["end"] = int(row.get("end", 0))
            except (ValueError, TypeError):
                continue
            row["contig"] = contig
            row["strand"] = (row.get("strand") or "+").strip()
            by_contig[contig].append(row)
    # Sort each contig's CDS by start so overlap search can short-circuit
    for contig in by_contig:
        by_contig[contig].sort(key=lambda r: r["start"])
    return by_contig


def _reciprocal_overlap(a_start, a_end, b_start, b_end) -> float:
    """Min of (overlap / len(a), overlap / len(b)). 0.0 if no overlap."""
    o_start = max(a_start, b_start)
    o_end = min(a_end, b_end)
    overlap = o_end - o_start
    if overlap <= 0:
        return 0.0
    la = a_end - a_start
    lb = b_end - b_start
    if la <= 0 or lb <= 0:
        return 0.0
    return min(overlap / la, overlap / lb)


def best_gbff_match(
    bakta_cds: dict,
    candidates: list,
    min_overlap: float = _DEFAULT_MIN_OVERLAP,
    strand_must_match: bool = True,
) -> dict | None:
    """Return the best-matching GenBank CDS on the same contig, or None.

    "Best" = highest reciprocal overlap above `min_overlap`. Strand
    mismatches are rejected by default (same bp range but antisense is
    a different gene).
    """
    best = None
    best_score = 0.0
    for gbff in candidates:
        if strand_must_match and gbff["strand"] != bakta_cds["strand"]:
            continue
        # Candidates are sorted by start; once start exceeds bakta's end
        # there can be no more overlap
        if gbff["start"] >= bakta_cds["end"]:
            break
        score = _reciprocal_overlap(
            bakta_cds["start"], bakta_cds["end"], gbff["start"], gbff["end"]
        )
        if score >= min_overlap and score > best_score:
            best = gbff
            best_score = score
    return best


def map_gbff_to_bakta(
    bakta_by_contig: dict,
    gbff_by_contig: dict,
    min_overlap: float = _DEFAULT_MIN_OVERLAP,
):
    """Yield Bakta rows augmented with `gbff_annotation`.

    Iterates Bakta's CDS; for each, looks up the best GenBank match on
    the same contig by reciprocal overlap. Emits rows in Bakta's order
    so the output preserves the locus-tag assignment downstream code
    keys on. All Bakta columns (ec_numbers, cog_ids, etc.) are
    preserved unchanged — only `gbff_annotation` is added.
    """
    for contig, bakta_rows in bakta_by_contig.items():
        gbff_candidates = gbff_by_contig.get(contig, [])
        for bakta in bakta_rows:
            match = best_gbff_match(bakta, gbff_candidates, min_overlap=min_overlap)
            row = dict(bakta)
            row["contig"] = contig
            row["gbff_annotation"] = (match or {}).get("product", "")
            yield row


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Map original GenBank CDS annotations onto Bakta's CDS via "
            "reciprocal coordinate overlap."
        )
    )
    parser.add_argument(
        "--bakta-gene-info",
        required=True,
        help="Bakta's gene_info.tsv (canonical CDS set)",
    )
    parser.add_argument(
        "--genbank-gene-info",
        required=True,
        help="Original GenBank-derived gene_info.tsv (source of gbff_annotation)",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output gene_info.tsv — Bakta's rows + gbff_annotation column",
    )
    parser.add_argument(
        "--min-overlap",
        type=float,
        default=_DEFAULT_MIN_OVERLAP,
        help=f"Minimum reciprocal overlap to accept a match (default: {_DEFAULT_MIN_OVERLAP})",
    )
    args = parser.parse_args()

    bakta = _read_gene_info(args.bakta_gene_info)
    gbff = _read_gene_info(args.genbank_gene_info)

    # Preserve every column the Bakta TSV had; tack on `gbff_annotation`
    # at the end. We re-read Bakta's header directly so any future Bakta
    # columns flow through without code changes here.
    with open(args.bakta_gene_info) as f:
        bakta_header = next(csv.reader(f, delimiter="\t"))
    fieldnames = list(bakta_header)
    if "gbff_annotation" not in fieldnames:
        fieldnames.append("gbff_annotation")

    n_mapped = 0
    n_total = 0
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in map_gbff_to_bakta(bakta, gbff, min_overlap=args.min_overlap):
            writer.writerow({k: row.get(k, "") for k in fieldnames})
            n_total += 1
            if row["gbff_annotation"]:
                n_mapped += 1

    pct = (n_mapped / n_total * 100.0) if n_total else 0.0
    logger.info(
        f"Mapped {n_mapped}/{n_total} Bakta CDS to GenBank annotations "
        f"({pct:.1f}%) at min_overlap={args.min_overlap}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
