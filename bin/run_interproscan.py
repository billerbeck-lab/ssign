#!/usr/bin/env python3
"""Run InterProScan locally.

Preserves the TSV column-index parsing + GO term extraction from the
original pipeline/scripts/parse_interproscan.py.
"""

import argparse
import csv
import logging
import os
import re
import subprocess
import sys
import tempfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.fasta_io import read_fasta
from ssign_lib.substrates import load_substrate_ids
from dedup_sequences import deduplicate_dict, expand_results_dict

# InterProScan TSV column indices (0-based, no header)
_COL_PROTEIN_ID = 0
_COL_SIG_ACC = 4
_COL_SIG_DESC = 5
_COL_IPR_ACC = 11
_COL_IPR_DESC = 12
_COL_GO_TERMS = 13
_MISSING = "-"
_GO_ID_RE = re.compile(r"(GO:\d+)")

# Bacteria-relevant member DBs. PANTHER (eukaryote-leaning, slowest IPS
# member) is excluded by default. Pass --applications "" to run all DBs.
DEFAULT_IPS_APPLICATIONS = (
    "Pfam",
    "TIGRFAM",
    "HAMAP",
    "SMART",
    "PIRSF",
    "SUPERFAMILY",
    "Gene3D",
    "ProSiteProfiles",
    "ProSitePatterns",
    "CDD",
)


def run_local_interproscan(
    query_fasta,
    db_path,
    output_dir,
    applications="",
    offline=False,
):
    """Run InterProScan locally and return path to the TSV output.

    By default, IPS queries the EBI precalculated-match lookup service
    (5-10× speedup on previously-seen sequences; falls back to local
    HMMs after a ~30 s connection timeout if unreachable). Pass
    `offline=True` to add `-dp` and skip the lookup outright — useful
    on isolated networks where the timeout would dominate runtime.

    `applications` restricts the member-DB scan; pass "" to use IPS
    defaults.
    """
    output_file = os.path.join(output_dir, "results.tsv")
    cmd = [
        "interproscan.sh",
        "-i",
        query_fasta,
        "-o",
        output_file,
        "-f",
        "tsv",
        "-goterms",
        "-pathways",
    ]
    if offline:
        cmd.append("-dp")
    if applications:
        cmd.extend(["-appl", applications])
    if db_path:
        cmd.extend(["-d", db_path])

    logger.info("Running local InterProScan...")
    # FRAGILE: subprocess call requires interproscan.sh on PATH
    # If this breaks: install InterProScan from https://www.ebi.ac.uk/interpro/download/
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"InterProScan (interproscan.sh) not found: {e}\n"
            f"  Common causes:\n"
            f"    - InterProScan is not installed or not on PATH\n"
            f"  How to fix:\n"
            f"    - Download: https://www.ebi.ac.uk/interpro/download/\n"
            f"    - See docs/how-to/install.md for full setup steps"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"InterProScan timed out after 4 hours: {e}\n  How to fix:\n    - Reduce the number of input sequences"
        ) from e

    if result.returncode != 0:
        logger.error(f"InterProScan failed: {result.stderr[:500]}")
        raise RuntimeError(f"InterProScan exit code {result.returncode}")

    return output_file


def parse_interproscan_tsv(tsv_path, target_ids=None):
    """Parse InterProScan TSV output and aggregate per protein."""
    per_protein = {}

    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue

            protein_id = parts[_COL_PROTEIN_ID]
            if target_ids and protein_id not in target_ids:
                continue

            if protein_id not in per_protein:
                per_protein[protein_id] = {
                    "domains": set(),
                    "go_terms": set(),
                    "pfam_ids": set(),
                    "descriptions": set(),
                }

            sig_acc = parts[_COL_SIG_ACC] if len(parts) > _COL_SIG_ACC else _MISSING
            ipr_acc = parts[_COL_IPR_ACC] if len(parts) > _COL_IPR_ACC else _MISSING
            ipr_desc = parts[_COL_IPR_DESC] if len(parts) > _COL_IPR_DESC else _MISSING
            go_raw = parts[_COL_GO_TERMS] if len(parts) > _COL_GO_TERMS else _MISSING

            if ipr_acc != _MISSING:
                per_protein[protein_id]["domains"].add(ipr_acc)
            if ipr_desc != _MISSING:
                per_protein[protein_id]["descriptions"].add(ipr_desc)
            if sig_acc != _MISSING and sig_acc.startswith("PF"):
                per_protein[protein_id]["pfam_ids"].add(sig_acc)

            if go_raw != _MISSING:
                for match in _GO_ID_RE.finditer(go_raw):
                    per_protein[protein_id]["go_terms"].add(match.group(1))

    results = {}
    for pid, data in per_protein.items():
        results[pid] = {
            "locus_tag": pid,
            "interpro_domains": ";".join(sorted(data["domains"])),
            "interpro_go_terms": ";".join(sorted(data["go_terms"])),
            "interpro_pfam_ids": ";".join(sorted(data["pfam_ids"])),
            "interpro_descriptions": ";".join(sorted(data["descriptions"])),
        }

    return results


def main():
    parser = argparse.ArgumentParser(description="Run InterProScan locally")
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--db",
        default="",
        help="Optional -d path to custom InterProScan data directory",
    )
    parser.add_argument(
        "--applications",
        default=",".join(DEFAULT_IPS_APPLICATIONS),
        help=(
            "Comma-separated InterProScan member DBs to run (-appl). "
            "Default skips PANTHER (eukaryote-leaning, slowest IPS member). "
            "Pass empty string to run all DBs."
        ),
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help=(
            "Add -dp to skip the EBI precalc-match lookup service. "
            "Use on isolated networks where the lookup-timeout (~30 s/run) "
            "would dominate runtime."
        ),
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    substrate_ids = load_substrate_ids(args.substrates)
    all_seqs = read_fasta(args.proteins)
    sub_seqs = {k: v for k, v in all_seqs.items() if k in substrate_ids}

    # Deduplicate to avoid redundant InterProScan work on identical substrates
    unique_seqs, seq_groups = deduplicate_dict(sub_seqs)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = os.path.join(tmpdir, "substrates.fasta")
        with open(tmp_fasta, "w") as f:
            for pid, seq in unique_seqs.items():
                f.write(f">{pid}\n{seq}\n")

        tsv_path = run_local_interproscan(
            tmp_fasta,
            args.db,
            tmpdir,
            applications=args.applications,
            offline=args.offline,
        )
        results_unique = parse_interproscan_tsv(tsv_path, set(unique_seqs.keys()))

    results = expand_results_dict(results_unique, seq_groups)

    logger.info(f"Annotated {len(results)}/{len(substrate_ids)} substrates for {args.sample}")

    fieldnames = [
        "locus_tag",
        "interpro_domains",
        "interpro_go_terms",
        "interpro_pfam_ids",
        "interpro_descriptions",
    ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results.values():
            writer.writerow(r)


if __name__ == "__main__":
    main()
