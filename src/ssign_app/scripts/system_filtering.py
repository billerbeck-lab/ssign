#!/usr/bin/env python3
"""Apply system-level filtering and DSE cross-genome validation.

Merges proximity substrates with T5SS substrates, applies exclusion
filters, and produces filtered + unfiltered substrate lists.

Preserves the dse_type_in_genome bug fix from Session 11.
"""

import argparse
import csv
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Filter substrates")
    parser.add_argument("--proximity-substrates", required=True)
    parser.add_argument("--t5ss-substrates", required=True)
    parser.add_argument("--valid-systems", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--excluded-systems", default="Flagellum,Tad,T3SS")
    parser.add_argument("--required-fraction-correct", type=float, default=0.8)
    parser.add_argument("--filter-dse-type-mismatch", action="store_true",
                        help="Remove DSE-only substrates where DSE predicted type "
                             "doesn't match nearby MacSyFinder system type")
    parser.add_argument("--out-filtered", required=True)
    parser.add_argument("--out-all", required=True)
    args = parser.parse_args()

    excluded = set(s.strip() for s in args.excluded_systems.split(',') if s.strip())

    # Load proximity substrates
    substrates = []
    with open(args.proximity_substrates) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            row['substrate_source'] = 'proximity'
            substrates.append(row)

    # Load T5SS substrates
    with open(args.t5ss_substrates) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            row['substrate_source'] = 'T5SS-self'
            substrates.append(row)

    # Write unfiltered (all substrates before exclusion)
    if substrates:
        fieldnames = list(substrates[0].keys())
    else:
        fieldnames = ['locus_tag', 'sample_id', 'tool', 'nearby_ss_types',
                      'substrate_source']

    with open(args.out_all, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for s in substrates:
            writer.writerow(s)

    # Filter: remove substrates associated only with excluded SS types
    # Also optionally filter DSE-type mismatches
    filtered = []
    n_dse_mismatch_removed = 0
    for s in substrates:
        ss_types = set(s.get('nearby_ss_types', '').split(','))
        ss_types.discard('')

        # Keep if any non-excluded SS type is associated
        non_excluded = ss_types - excluded
        if non_excluded or s.get('substrate_source') == 'T5SS-self':
            # Update nearby_ss_types to only show non-excluded
            if non_excluded:
                s['nearby_ss_types'] = ','.join(sorted(non_excluded))

            # DSE type-match filter: for DSE-only substrates, check if
            # DSE predicted type matches the nearby MacSyFinder system
            if args.filter_dse_type_mismatch:
                tool = s.get('tool', '')
                dse_match = s.get('dse_type_match', 'True')
                # Only filter DSE-only substrates (not DLP or DLP+DSE)
                if tool == 'DSE' and str(dse_match).lower() in ('false', '0', ''):
                    n_dse_mismatch_removed += 1
                    continue

            filtered.append(s)

    if n_dse_mismatch_removed:
        logger.info(
            f"DSE type-match filter removed {n_dse_mismatch_removed} DSE-only substrates "
            f"where predicted SS type didn't match nearby MacSyFinder system"
        )

    with open(args.out_filtered, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for s in filtered:
            writer.writerow(s)

    logger.info(
        f"{args.sample}: {len(substrates)} total substrates, "
        f"{len(filtered)} after filtering (excluded: {excluded})"
    )


if __name__ == '__main__':
    main()
