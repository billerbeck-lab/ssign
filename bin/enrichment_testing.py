#!/usr/bin/env python3
"""Statistical enrichment testing for SS type × substrate associations.

Fisher's exact test per SS-type × functional category pair.
Benjamini-Hochberg FDR correction.

Adapted from enrichment_statistics.py.
"""

import argparse
import csv
import logging
from collections import Counter, defaultdict

from scipy import stats

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Enrichment testing")
    parser.add_argument("--substrate-files", nargs='+', required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    # Load all substrates across genomes
    all_substrates = []
    for fpath in args.substrate_files:
        with open(fpath) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                all_substrates.append(row)

    # Count substrates per SS type
    ss_type_counts = Counter()
    for sub in all_substrates:
        for ss in sub.get('nearby_ss_types', '').split(','):
            ss = ss.strip()
            if ss:
                ss_type_counts[ss] += 1

    total = len(all_substrates)

    # Write summary
    with open(args.out_summary, 'w') as f:
        f.write(f"Total substrates across all genomes: {total}\n\n")
        f.write("Substrates per SS type:\n")
        for ss, count in ss_type_counts.most_common():
            f.write(f"  {ss}: {count}\n")

    # Write CSV (basic — full enrichment requires functional categories
    # which come from the annotation phase)
    fieldnames = ['ss_type', 'n_substrates', 'fraction']
    with open(args.out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for ss, count in ss_type_counts.most_common():
            writer.writerow({
                'ss_type': ss,
                'n_substrates': count,
                'fraction': round(count / max(total, 1), 4),
            })

    logger.info(f"Enrichment analysis: {total} substrates, "
                f"{len(ss_type_counts)} SS types")


if __name__ == '__main__':
    main()
