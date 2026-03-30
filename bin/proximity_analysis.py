#!/usr/bin/env python3
"""Per-component proximity analysis for substrate identification.

Finds extracellular proteins within +/- N genes of each secretion system
component on the same contig. NEVER spans contigs.

CRITICAL: Uses per-COMPONENT proximity, NOT system-boundary proximity.
This was a bug fix — system-boundary proximity (finding proteins within
the full span of a secretion system) is less biologically correct.

Adapted from extract_substrate_sequences.py::get_near_ss_positions()
"""

import argparse
import csv
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

from ssign_lib.constants import DSE_TO_MACSYFINDER


def dse_type_in_genome(dse_ss_type: str, genome_ss_types: set) -> bool:
    """Check if the DSE-predicted SS type has a validated system in this genome.

    This prevents cross-genome leakage where a protein gets classified as a
    substrate because DeepSecE predicts a SS type that exists in other genomes
    but not in this protein's own genome.

    Bug fix from Session 11 — removed 26 false positives.
    """
    dse_str = str(dse_ss_type)
    macsyfinder_names = DSE_TO_MACSYFINDER.get(dse_str, [dse_str])
    return any(
        mf_name in gt
        for mf_name in macsyfinder_names
        for gt in genome_ss_types
    )


def main():
    parser = argparse.ArgumentParser(description="Per-component proximity analysis")
    parser.add_argument("--gene-order", required=True)
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--window", type=int, default=3)
    parser.add_argument("--conf-threshold", type=float, default=0.8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load gene order: {contig: [{gene_index, locus_tag, ...}]}
    genes_by_contig = {}
    locus_to_info = {}
    with open(args.gene_order) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            contig = row['contig']
            idx = int(row['gene_index'])
            if contig not in genes_by_contig:
                genes_by_contig[contig] = []
            genes_by_contig[contig].append(row)
            locus_to_info[row['locus_tag']] = {
                'contig': contig,
                'gene_index': idx,
            }

    # Load SS components (non-excluded only)
    ss_component_loci = set()
    genome_ss_types = set()
    component_ss_types = {}  # locus_tag -> ss_type
    with open(args.ss_components) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('excluded', 'False').lower() == 'true':
                continue
            locus = row['locus_tag']
            ss_component_loci.add(locus)
            genome_ss_types.add(row['ss_type'])
            component_ss_types[locus] = row['ss_type']

    # Load predictions
    predictions = {}
    with open(args.predictions) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            predictions[row['locus_tag']] = row

    # Build index: for each contig, map gene_index -> locus_tag
    contig_index = {}
    for contig, genes in genes_by_contig.items():
        idx_map = {}
        for g in genes:
            idx_map[int(g['gene_index'])] = g['locus_tag']
        contig_index[contig] = idx_map

    # Per-component proximity: for each SS component, find proteins within +/- window
    substrates = {}  # locus_tag -> {metadata}
    for comp_locus in ss_component_loci:
        info = locus_to_info.get(comp_locus)
        if not info:
            continue

        contig = info['contig']
        comp_idx = info['gene_index']
        idx_map = contig_index.get(contig, {})
        max_idx = max(idx_map.keys()) if idx_map else 0

        # Check proteins within window on same contig
        for offset in range(-args.window, args.window + 1):
            neighbor_idx = comp_idx + offset
            if neighbor_idx < 0 or neighbor_idx > max_idx:
                continue

            neighbor_locus = idx_map.get(neighbor_idx)
            if not neighbor_locus:
                continue

            # Skip SS components themselves (they're not substrates, except T5SS)
            if neighbor_locus in ss_component_loci:
                continue

            pred = predictions.get(neighbor_locus, {})

            # Check DeepLocPro extracellular
            try:
                dlp_ext_prob = float(pred.get('extracellular_prob', 0))
            except (ValueError, TypeError):
                dlp_ext_prob = 0.0

            is_dlp = dlp_ext_prob >= args.conf_threshold

            # Check DeepSecE
            dse_type = pred.get('dse_ss_type', 'Non-secreted')
            try:
                dse_max_prob = float(pred.get('dse_max_prob', 0))
            except (ValueError, TypeError):
                dse_max_prob = 0.0

            is_dse = (
                dse_type not in ('Non-secreted', 'T3SS', '')
                and dse_max_prob >= args.conf_threshold
            )

            # DSE cross-genome leakage fix: verify SS type exists in this genome
            if is_dse and not is_dlp:
                if not dse_type_in_genome(dse_type, genome_ss_types):
                    is_dse = False

            if is_dlp or is_dse:
                tool = []
                if is_dlp:
                    tool.append('DLP')
                if is_dse:
                    tool.append('DSE')

                ss_type = component_ss_types.get(comp_locus, '')

                if neighbor_locus not in substrates:
                    substrates[neighbor_locus] = {
                        'locus_tag': neighbor_locus,
                        'sample_id': args.sample,
                        'tool': '+'.join(tool),
                        'nearby_ss_types': {ss_type},
                        'dlp_extracellular_prob': dlp_ext_prob,
                        'dse_ss_type': dse_type,
                        'dse_max_prob': dse_max_prob,
                        'product': pred.get('product', ''),
                    }
                else:
                    substrates[neighbor_locus]['nearby_ss_types'].add(ss_type)
                    # Upgrade tool if both detected
                    existing_tools = set(substrates[neighbor_locus]['tool'].split('+'))
                    existing_tools.update(tool)
                    substrates[neighbor_locus]['tool'] = '+'.join(sorted(existing_tools))

    # Write output
    fieldnames = ['locus_tag', 'sample_id', 'tool', 'nearby_ss_types',
                  'dlp_extracellular_prob', 'dse_ss_type', 'dse_max_prob', 'product']
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for sub in substrates.values():
            sub['nearby_ss_types'] = ','.join(sorted(sub['nearby_ss_types']))
            writer.writerow(sub)

    logger.info(f"Found {len(substrates)} putative substrates in {args.sample}")


if __name__ == '__main__':
    main()
