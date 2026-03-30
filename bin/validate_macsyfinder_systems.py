#!/usr/bin/env python3
"""Validate MacSyFinder v2 results and extract SS component mappings.

Parses MacSyFinder output, applies wholeness threshold, maps components
to locus_tags, and optionally excludes specified system types.

Adapted from extract_substrate_sequences.py system validation logic.
"""

import argparse
import csv
import logging
import os
import re
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def parse_sys_id(sys_id: str) -> str:
    """Extract SS type from MacSyFinder sys_id string.

    Examples:
        'TXSS__T2SS' -> 'T2SS'
        'TXSS__T5aSS' -> 'T5aSS'
        'TXSS__Flagellum' -> 'Flagellum'
    """
    parts = sys_id.split('__')
    return parts[-1] if parts else sys_id


def parse_macsyfinder_results(msf_dir: str):
    """Parse MacSyFinder all_systems.tsv output.

    Returns list of system dicts with: sys_id, ss_type, wholeness,
    components (list of {gene, gene_name, ...}).
    """
    systems_file = os.path.join(msf_dir, "all_systems.tsv")

    if not os.path.exists(systems_file):
        # Try alternative location
        for fname in os.listdir(msf_dir):
            if fname.endswith("_all_systems.tsv") or fname == "all_systems.tsv":
                systems_file = os.path.join(msf_dir, fname)
                break

    if not os.path.exists(systems_file):
        logger.warning(f"No all_systems.tsv found in {msf_dir}")
        return []

    systems = {}
    with open(systems_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sys_id = row.get('sys_id', '')
            if not sys_id:
                continue

            if sys_id not in systems:
                # Parse wholeness from the system-level fields
                try:
                    wholeness = float(row.get('sys_wholeness', 0))
                except (ValueError, TypeError):
                    wholeness = 0.0

                systems[sys_id] = {
                    'sys_id': sys_id,
                    'ss_type': parse_sys_id(sys_id),
                    'wholeness': wholeness,
                    'components': [],
                }

            # Add this component
            systems[sys_id]['components'].append({
                'gene': row.get('hit_id', ''),
                'gene_name': row.get('gene_name', ''),
                'sys_status': row.get('sys_status', ''),
                'gene_status': row.get('gene_status', ''),
            })

    return list(systems.values())


def main():
    parser = argparse.ArgumentParser(description="Validate MacSyFinder systems")
    parser.add_argument("--msf-dir", required=True, help="MacSyFinder output directory")
    parser.add_argument("--gene-info", required=True, help="Gene info TSV")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--wholeness-threshold", type=float, default=0.8)
    parser.add_argument("--excluded-systems", default="Flagellum,Tad,T3SS",
                        help="Comma-separated system types to exclude")
    parser.add_argument("--out-components", required=True)
    parser.add_argument("--out-systems", required=True)
    args = parser.parse_args()

    excluded = set(s.strip() for s in args.excluded_systems.split(',') if s.strip())

    # Parse MacSyFinder results
    all_systems = parse_macsyfinder_results(args.msf_dir)
    logger.info(f"Found {len(all_systems)} systems in MacSyFinder output")

    # Filter by wholeness
    valid_systems = [s for s in all_systems if s['wholeness'] >= args.wholeness_threshold]
    logger.info(f"{len(valid_systems)} systems pass wholeness >= {args.wholeness_threshold}")

    # Separate excluded from included
    included_systems = [s for s in valid_systems if s['ss_type'] not in excluded]
    excluded_count = len(valid_systems) - len(included_systems)
    if excluded_count:
        logger.info(f"Excluded {excluded_count} systems of types: {excluded}")

    # Write valid systems summary
    with open(args.out_systems, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample_id', 'sys_id', 'ss_type', 'wholeness', 'n_components', 'excluded'
        ], delimiter='\t')
        writer.writeheader()
        for s in valid_systems:
            writer.writerow({
                'sample_id': args.sample,
                'sys_id': s['sys_id'],
                'ss_type': s['ss_type'],
                'wholeness': s['wholeness'],
                'n_components': len(s['components']),
                'excluded': s['ss_type'] in excluded,
            })

    # Write component-level detail (all valid systems, including excluded for reference)
    with open(args.out_components, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample_id', 'sys_id', 'ss_type', 'locus_tag', 'gene_name',
            'gene_status', 'wholeness', 'excluded'
        ], delimiter='\t')
        writer.writeheader()
        for s in valid_systems:
            for comp in s['components']:
                writer.writerow({
                    'sample_id': args.sample,
                    'sys_id': s['sys_id'],
                    'ss_type': s['ss_type'],
                    'locus_tag': comp['gene'],
                    'gene_name': comp['gene_name'],
                    'gene_status': comp['gene_status'],
                    'wholeness': s['wholeness'],
                    'excluded': s['ss_type'] in excluded,
                })


if __name__ == '__main__':
    main()
