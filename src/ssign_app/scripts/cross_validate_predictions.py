#!/usr/bin/env python3
"""Cross-validate DeepLocPro and DeepSecE predictions.

Merges prediction outputs and flags unreliable predictions:
- DeepSecE T3SS predictions without MacSyFinder T3SS support
  (MacSyFinder found 0 T3SS across all 74 genomes; DeepSecE predicts
  1808 T3SS — mostly hypothetical proteins and flagellar misclassification)
"""

import argparse
import csv
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Cross-validate predictions")
    parser.add_argument("--deeplocpro", required=True)
    parser.add_argument("--deepsece", default="")
    parser.add_argument("--valid-systems", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--conf-threshold", type=float, default=0.8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Check which SS types exist in this genome
    genome_ss_types = set()
    with open(args.valid_systems) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('excluded', 'False').lower() != 'true':
                genome_ss_types.add(row['ss_type'])

    has_t3ss = any('T3SS' in st for st in genome_ss_types)

    # Load DeepLocPro predictions
    dlp_data = {}
    with open(args.deeplocpro) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            dlp_data[row['locus_tag']] = row

    # Load DeepSecE predictions (optional)
    dse_data = {}
    if args.deepsece and os.path.exists(args.deepsece):
        with open(args.deepsece) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                dse_data[row['locus_tag']] = row
    else:
        logger.info("DeepSecE not available — running with DeepLocPro only")

    # Merge and cross-validate
    all_loci = sorted(set(dlp_data.keys()) | set(dse_data.keys()))

    fieldnames = [
        'locus_tag', 'sample_id',
        # DeepLocPro
        'predicted_localization', 'extracellular_prob',
        'periplasmic_prob', 'outer_membrane_prob', 'cytoplasmic_prob',
        # DeepSecE
        'dse_ss_type', 'dse_max_prob', 'dse_T3SS_flagged',
        # Agreement
        'is_extracellular', 'product',
    ]

    n_flagged_t3ss = 0
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for locus in all_loci:
            dlp = dlp_data.get(locus, {})
            dse = dse_data.get(locus, {})

            try:
                ext_prob = float(dlp.get('extracellular_prob', 0))
            except (ValueError, TypeError):
                ext_prob = 0.0

            dse_type = dse.get('dse_ss_type', 'Non-secreted')
            try:
                dse_max = float(dse.get('dse_max_prob', 0))
            except (ValueError, TypeError):
                dse_max = 0.0

            # Flag DeepSecE T3SS predictions without MacSyFinder support
            t3ss_flagged = (dse_type == 'T3SS' and not has_t3ss)
            if t3ss_flagged:
                n_flagged_t3ss += 1

            writer.writerow({
                'locus_tag': locus,
                'sample_id': args.sample,
                'predicted_localization': dlp.get('predicted_localization', ''),
                'extracellular_prob': ext_prob,
                'periplasmic_prob': dlp.get('periplasmic_prob', ''),
                'outer_membrane_prob': dlp.get('outer_membrane_prob', ''),
                'cytoplasmic_prob': dlp.get('cytoplasmic_prob', ''),
                'dse_ss_type': dse_type,
                'dse_max_prob': dse_max,
                'dse_T3SS_flagged': t3ss_flagged,
                'is_extracellular': ext_prob >= args.conf_threshold,
                'product': dlp.get('product', dse.get('product', '')),
            })

    if n_flagged_t3ss:
        logger.warning(
            f"Flagged {n_flagged_t3ss} DeepSecE T3SS predictions "
            f"(no MacSyFinder T3SS found in {args.sample})"
        )

    logger.info(f"Cross-validated {len(all_loci)} proteins for {args.sample}")


if __name__ == '__main__':
    main()
