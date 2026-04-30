#!/usr/bin/env python3
"""Handle T5SS (Type V Secretion System) self-secreting autotransporters.

T5aSS proteins ARE their own substrates (self-secreting). PF03797 is the
definitive marker — DeepLocPro thresholds are NOT used for T5aSS (median
outer membrane probability is only 0.47, making DLP unreliable for these).

Domain classification (adapted from split_t5a_domains.py):
- Classical AT (has passenger): PF03797 barrel + passenger >= 100aa
- Minimal passenger: PF03797 barrel + passenger 1-100aa
- Barrel-only: PF03797 barrel, no passenger
- OMP/Porin: PF13505 (not PF03797), 231-264aa
"""

import argparse
import csv
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

import os as _os
import sys as _sys

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)


def main():
    parser = argparse.ArgumentParser(description="Handle T5SS autotransporters")
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-substrates", required=True)
    parser.add_argument("--out-domains", required=True)
    args = parser.parse_args()

    # Identify T5SS components
    t5ss_components = []
    with open(args.ss_components) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            ss_type = row.get('ss_type', '')
            if ss_type.startswith('T5'):
                t5ss_components.append(row)

    # Load predictions for product info
    predictions = {}
    with open(args.predictions) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            predictions[row['locus_tag']] = row

    # T5aSS components are self-substrates
    substrates = []
    for comp in t5ss_components:
        locus = comp['locus_tag']
        pred = predictions.get(locus, {})

        try:
            dlp_prob = float(pred.get('dlp_extracellular_prob',
                                            pred.get('extracellular_prob', 0)))
        except (ValueError, TypeError):
            dlp_prob = 0.0

        substrates.append({
            'locus_tag': locus,
            'sample_id': args.sample,
            'tool': 'T5SS-self',
            'nearby_ss_types': comp.get('ss_type', 'T5aSS'),
            'dlp_extracellular_prob': dlp_prob,
            'predicted_localization': pred.get('predicted_localization', ''),
            'dlp_max_localization': pred.get('dlp_max_localization', ''),
            'dlp_max_probability': pred.get('dlp_max_probability', ''),
            'dse_ss_type': pred.get('dse_ss_type', ''),
            'dse_max_prob': pred.get('dse_max_prob', ''),
            'signalp_prediction': pred.get('signalp_prediction', ''),
            'signalp_probability': pred.get('signalp_probability', ''),
            'signalp_cs_position': pred.get('signalp_cs_position', ''),
            'product': pred.get('product', ''),
        })

    # Write substrates
    sub_fields = ['locus_tag', 'sample_id', 'tool', 'nearby_ss_types',
                  'dlp_extracellular_prob', 'predicted_localization',
                  'dlp_max_localization', 'dlp_max_probability',
                  'dse_ss_type', 'dse_max_prob',
                  'signalp_prediction', 'signalp_probability', 'signalp_cs_position',
                  'product']
    with open(args.out_substrates, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sub_fields, delimiter='\t')
        writer.writeheader()
        for s in substrates:
            writer.writerow(s)

    # Write domain classification (basic — full InterProScan-based classification
    # requires InterProScan results which may not be available yet)
    domain_fields = ['locus_tag', 'sample_id', 'ss_type', 'domain_group']
    with open(args.out_domains, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=domain_fields, delimiter='\t')
        writer.writeheader()
        for comp in t5ss_components:
            writer.writerow({
                'locus_tag': comp['locus_tag'],
                'sample_id': args.sample,
                'ss_type': comp.get('ss_type', ''),
                'domain_group': 'T5SS-component',
            })

    logger.info(f"Found {len(substrates)} T5SS self-substrates in {args.sample}")


if __name__ == '__main__':
    main()
