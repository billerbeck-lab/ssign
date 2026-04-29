#!/usr/bin/env python3
"""Integrate all annotation sources into master CSV.

Left-join merges on locus_tag, adds GBFF annotation from gene_info,
protein sequences, and computes annotation tool hit counts.

Annotation tools counted: BLASTp, HHpred (Pfam), HHpred (PDB),
InterProScan, ProtParam, GBFF (original genome annotation).

SignalP is a secretion prediction tool (not an annotation tool) and
is already included in the substrates table via cross_validate.
"""

import argparse
import logging
import os

import pandas as pd
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Annotation tool columns — used for n_tools_with_hits counting.
# SignalP is NOT here (it's a prediction tool, not annotation).
# Bakta and EggNOG are new in Phase 3.2.c; their columns are populated
# by runner.py once 3.2.d wires those steps into the master-CSV merge.
# Missing columns here are skipped at consensus time (see _compute_consensus).
TOOL_HIT_COLUMNS = {
    "BLASTp": "blastp_hit_description",
    "HHpred_Pfam": "pfam_top1_description",
    "HHpred_PDB": "pdb_top1_description",
    "InterProScan": "interpro_descriptions",
    "pLM-BLAST": "ecod70_top1_description",
    "Bakta": "bakta_product",
    "EggNOG": "eggnog_description",
    "GBFF": "gbff_annotation",
}


def _compute_consensus(df):
    """Compute annotation consensus and tool counts for each protein."""
    # Relative import (script directory is on sys.path when invoked as
    # a subprocess by runner.py, AND when invoked directly via
    # `python integrate_annotations.py`). The previous absolute path
    # `ssign_app.scripts.annotation_consensus` only worked if the
    # package was pip-installed.
    import os as _os
    import sys as _sys

    _here = _os.path.dirname(_os.path.abspath(__file__))
    if _here not in _sys.path:
        _sys.path.insert(0, _here)
    from annotation_consensus import compute_consensus

    consensus_rows = []
    for _, row in df.iterrows():
        # Collect descriptions from each tool
        tool_descs = {}
        for tool, col in TOOL_HIT_COLUMNS.items():
            if col not in df.columns:
                continue
            val = row.get(col, "")
            if pd.isna(val) or not str(val).strip():
                continue
            # Skip generic GBFF annotations
            if tool == "GBFF" and str(val).strip().lower() in (
                "hypothetical protein",
                "uncharacterized protein",
                "",
            ):
                continue
            tool_descs[tool] = str(val).strip()

        consensus = compute_consensus(tool_descs)
        consensus_rows.append(consensus)

    # Add consensus columns to dataframe
    consensus_df = pd.DataFrame(consensus_rows)
    for col in consensus_df.columns:
        df[col] = consensus_df[col].values

    # Also add the simple annotation_tools list
    df["annotation_tools"] = [
        ",".join(sorted(td.keys())) if td else ""
        for td in [
            {
                t: row.get(c, "")
                for t, c in TOOL_HIT_COLUMNS.items()
                if c in df.columns
                and pd.notna(row.get(c, ""))
                and str(row.get(c, "")).strip()
                and not (
                    t == "GBFF"
                    and str(row.get(c, "")).strip().lower()
                    in ("hypothetical protein", "uncharacterized protein", "")
                )
            }
            for _, row in df.iterrows()
        ]
    ]

    return df


def main():
    parser = argparse.ArgumentParser(description="Integrate annotations")
    parser.add_argument("--substrates-filtered", required=True)
    parser.add_argument("--substrates-all", required=True)
    parser.add_argument("--annotations", nargs="*", default=[])
    parser.add_argument(
        "--gene-info",
        default="",
        help="Gene info TSV from extract_proteins (for GBFF annotation)",
    )
    parser.add_argument(
        "--proteins", default="", help="Proteins FASTA (to include sequences in output)"
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load base substrate table
    df = pd.read_csv(args.substrates_filtered, sep="\t")
    logger.info(f"Base: {len(df)} filtered substrates for {args.sample}")

    # Add GBFF annotation from gene_info (original genome annotations).
    # If Phase 3.3.c re-annotation ran, gene_info already has a
    # `gbff_annotation` column carrying the original GenBank product —
    # use it directly. Otherwise (use_input_annotations=True, FASTA
    # contigs, etc.) fall back to renaming `product`.
    if args.gene_info and os.path.exists(args.gene_info):
        try:
            gi = pd.read_csv(args.gene_info, sep="\t")
            if "locus_tag" in gi.columns:
                if "gbff_annotation" in gi.columns:
                    gi_ann = gi[["locus_tag", "gbff_annotation"]].copy()
                elif "product" in gi.columns:
                    gi_ann = gi[["locus_tag", "product"]].copy()
                    gi_ann = gi_ann.rename(columns={"product": "gbff_annotation"})
                else:
                    gi_ann = None
                if gi_ann is not None:
                    gi_ann = gi_ann.drop_duplicates(subset="locus_tag")
                    before = len(df)
                    df = df.merge(gi_ann, on="locus_tag", how="left")
                    assert len(df) == before
                    logger.info("Added GBFF annotations from gene_info")
        except Exception as e:
            logger.warning(f"Failed to add GBFF annotations: {e}")

    # Add protein sequences
    if args.proteins and os.path.exists(args.proteins):
        try:
            seqs = {}
            for rec in SeqIO.parse(args.proteins, "fasta"):
                seqs[rec.id] = str(rec.seq)
            if seqs and "locus_tag" in df.columns:
                df["sequence"] = df["locus_tag"].map(seqs)
                df["aa_length"] = df["sequence"].apply(
                    lambda s: len(s) if pd.notna(s) else 0
                )
                logger.info(
                    f"Added sequences for {df['sequence'].notna().sum()} proteins"
                )
        except Exception as e:
            logger.warning(f"Failed to add sequences: {e}")

    # Left-join each annotation file
    for ann_file in args.annotations:
        if not os.path.exists(ann_file):
            continue

        try:
            ann_df = pd.read_csv(ann_file, sep=None, engine="python")

            # Determine join column
            join_col = "locus_tag" if "locus_tag" in ann_df.columns else None
            if not join_col:
                for col in ann_df.columns:
                    if "protein" in col.lower() or "id" == col.lower():
                        join_col = col
                        break

            if join_col and join_col in df.columns:
                # Avoid duplicate columns
                overlap = set(ann_df.columns) & set(df.columns) - {join_col}
                if overlap:
                    ann_df = ann_df.drop(columns=list(overlap))

                before = len(df)
                df = df.merge(ann_df, on=join_col, how="left")
                assert len(df) == before, (
                    f"Row count changed after merging {ann_file}: {before} -> {len(df)}"
                )
                logger.info(f"Merged {ann_file}: +{len(ann_df.columns) - 1} columns")
            else:
                logger.warning(f"No join column found in {ann_file}, skipping")

        except Exception as e:
            logger.warning(f"Failed to merge {ann_file}: {e}")

    # Compute annotation consensus across tools
    df = _compute_consensus(df)

    # Write output
    df.to_csv(args.output, index=False)
    logger.info(f"Wrote {len(df)} rows x {len(df.columns)} columns to {args.output}")


if __name__ == "__main__":
    main()
