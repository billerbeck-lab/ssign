#!/usr/bin/env python3
"""Validate MacSyFinder v2 results and extract SS component mappings.

Parses MacSyFinder output, applies wholeness threshold, maps components
to locus_tags, and optionally excludes specified system types.

Adapted from extract_substrate_sequences.py system validation logic.
"""

import argparse
import csv
import io
import logging
import os
import re

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def parse_sys_id(sys_id: str, model_fqn: str = "") -> str:
    """Extract SS type from MacSyFinder sys_id or model_fqn.

    MacSyFinder v2 sys_id format: 'replicon_SStype_N'
    model_fqn is more reliable: 'TXSScan/bacteria/diderm/T5aSS' -> 'T5aSS'

    Falls back to regex extraction from sys_id if model_fqn not available.
    """
    # Prefer model_fqn (most reliable)
    if model_fqn:
        return model_fqn.rstrip('/').split('/')[-1]

    # Fallback: try double-underscore format (older MacSyFinder)
    if '__' in sys_id:
        return sys_id.split('__')[-1]

    # MacSyFinder v2 format: repliconname_SStype_N
    # Match known SS type patterns
    m = re.search(r'(T[1-9][a-z]*SS[a-z]*|Flagellum|Tad|Com|MSH|T4P|pT4SS[it])',
                  sys_id, re.IGNORECASE)
    if m:
        return m.group(1)

    return sys_id


def parse_macsyfinder_results(msf_dir: str):
    """Parse MacSyFinder best_solution.tsv output.

    `best_solution.tsv` is the canonical authoritative output: highest-
    scoring non-overlapping combination of systems. `all_systems.tsv`
    contains overlapping candidates (the same components can appear in
    multiple competing system calls, e.g. T1SS_1 and T1SS_2 both built
    from the same hit_ids), which downstream proximity analysis would
    then double-count. See MacSyFinder docs (Néron et al. 2023).

    Returns list of system dicts with: sys_id, ss_type, wholeness,
    components (list of {gene, gene_name, ...}).

    Handles MacSyFinder v2 TSV format with # comment headers and blank lines.
    """
    systems_file = os.path.join(msf_dir, "best_solution.tsv")

    if not os.path.exists(systems_file):
        # Some MacSyFinder versions / sample-prefixed runs use {sample}_best_solution.tsv
        for fname in os.listdir(msf_dir):
            if fname.endswith("_best_solution.tsv") or fname == "best_solution.tsv":
                systems_file = os.path.join(msf_dir, fname)
                break

    if not os.path.exists(systems_file):
        logger.warning(f"No best_solution.tsv found in {msf_dir}")
        return []

    # Read file, skip comment lines (starting with #) and blank lines
    data_lines = []
    with open(systems_file) as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith('#'):
                data_lines.append(stripped)

    if not data_lines:
        logger.warning("best_solution.tsv is empty (no data lines)")
        return []

    # First non-comment line is the header
    reader = csv.DictReader(io.StringIO('\n'.join(data_lines)), delimiter='\t')

    systems = {}
    for row in reader:
        sys_id = row.get('sys_id', '')
        if not sys_id:
            continue

        model_fqn = row.get('model_fqn', '')

        if sys_id not in systems:
            try:
                wholeness = float(row.get('sys_wholeness', 0))
            except (ValueError, TypeError):
                wholeness = 0.0

            systems[sys_id] = {
                'sys_id': sys_id,
                'ss_type': parse_sys_id(sys_id, model_fqn),
                'wholeness': wholeness,
                'components': [],
            }

        # Add this component
        systems[sys_id]['components'].append({
            'gene': row.get('hit_id', ''),
            'gene_name': row.get('gene_name', ''),
            'sys_status': row.get('hit_status', ''),
            'gene_status': row.get('hit_status', ''),
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
