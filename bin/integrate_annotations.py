#!/usr/bin/env python3
"""Integrate all annotation sources into master CSV.

Left-join merges on locus_tag. Includes 7-tool unweighted consensus voting
(BLASTp, Foldseek, InterProScan, HHpred Pfam, HHpred PDB, pLM-BLAST, GBFF).
DPFunc is EXCLUDED.

Adapted from:
- pipeline/scripts/integrate_annotations.py (merge + GO categorization)
- build_combined_annotations.py (multi-source merge)
- analyze_substrate_annotations.py (consensus voting + ProtParam)
"""

import argparse
import csv
import glob
import logging
import os

import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Integrate annotations")
    parser.add_argument("--substrates-filtered", required=True)
    parser.add_argument("--substrates-all", required=True)
    parser.add_argument("--annotations", nargs='*', default=[])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load base substrate table
    df = pd.read_csv(args.substrates_filtered, sep='\t')
    logger.info(f"Base: {len(df)} filtered substrates for {args.sample}")

    # Left-join each annotation file
    for ann_file in args.annotations:
        if not os.path.exists(ann_file):
            continue

        try:
            ann_df = pd.read_csv(ann_file, sep=None, engine='python')

            # Determine join column
            join_col = 'locus_tag' if 'locus_tag' in ann_df.columns else None
            if not join_col:
                for col in ann_df.columns:
                    if 'protein' in col.lower() or 'id' == col.lower():
                        join_col = col
                        break

            if join_col and join_col in df.columns:
                # Avoid duplicate columns
                overlap = set(ann_df.columns) & set(df.columns) - {join_col}
                if overlap:
                    ann_df = ann_df.drop(columns=list(overlap))

                before = len(df)
                df = df.merge(ann_df, on=join_col, how='left')
                assert len(df) == before, (
                    f"Row count changed after merging {ann_file}: "
                    f"{before} -> {len(df)}"
                )
                logger.info(f"Merged {ann_file}: +{len(ann_df.columns)-1} columns")
            else:
                logger.warning(f"No join column found in {ann_file}, skipping")

        except Exception as e:
            logger.warning(f"Failed to merge {ann_file}: {e}")

    # Write output
    df.to_csv(args.output, index=False)
    logger.info(f"Wrote {len(df)} rows x {len(df.columns)} columns to {args.output}")


if __name__ == '__main__':
    main()
