#!/usr/bin/env python3
"""Merge PLM-Effector's per-type TSV outputs into one combined TSV.

PLM-Effector predicts bacterial secreted effectors per secretion-system
type (T1SE / T2SE / T3SE / T4SE / T6SE). Each run writes one TSV. The
ssign pipeline only cares whether a protein is secreted by **any** SS
type — so this script ORs `passes_threshold` across the per-type
outputs and records which types flagged the protein.

Usage:
    merge_plm_effector_outputs.py \\
        --inputs t1se.tsv t2se.tsv t3se.tsv t4se.tsv t6se.tsv \\
        --out plm_effector_merged.tsv

The merged TSV has one row per unique `seq_id` across all inputs. Each
row carries the max `stacking` probability across types that flagged
the protein, the comma-separated list of flagging types, and a single
`passes_threshold` that is 1 if any input flagged it.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def _coerce_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _passes(value) -> bool:
    """Interpret `passes_threshold` column (per-type TSV) as a bool.

    Per-type TSVs use "1"/"0" (or "True"/"False") per run_plm_effector's
    output format. Anything else is treated as 0.
    """
    return str(value).strip() in ("1", "True", "true")


def merge_per_type_outputs(per_type_paths):
    """Read per-type TSVs and yield merged rows, one per unique seq_id.

    Args:
        per_type_paths: iterable of paths to `run_plm_effector.py` outputs.
            Missing paths are silently skipped — the caller is expected to
            have filtered them already if desired.

    Yields:
        dicts with keys: seq_id, passes_threshold (int 0/1),
        flagging_types (comma-separated), max_stacking (float),
        effector_type (always "merged" to distinguish from per-type files).
    """
    # seq_id -> dict(passes_by_type: {type: bool}, max_stacking: float)
    merged: dict = {}

    for path in per_type_paths:
        if not path or not os.path.exists(path):
            logger.info(f"PLM-Effector per-type file missing, skipping: {path}")
            continue
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                seq_id = (row.get("seq_id") or "").strip()
                if not seq_id:
                    continue
                eff_type = (row.get("effector_type") or "").strip() or "?"
                passes = _passes(row.get("passes_threshold", "0"))
                stacking = _coerce_float(row.get("stacking", 0))

                entry = merged.setdefault(
                    seq_id,
                    {"passes_by_type": {}, "max_stacking": 0.0},
                )
                entry["passes_by_type"][eff_type] = passes
                if passes and stacking > entry["max_stacking"]:
                    entry["max_stacking"] = stacking

    for seq_id, entry in merged.items():
        flagging = sorted(t for t, p in entry["passes_by_type"].items() if p)
        yield {
            "seq_id": seq_id,
            "locus_tag": seq_id,
            "passes_threshold": "1" if flagging else "0",
            "flagging_types": ",".join(flagging),
            "max_stacking": f"{entry['max_stacking']:.4f}",
            "effector_type": "merged",
        }


def write_merged_tsv(rows, out_path: str) -> int:
    """Write merged rows to a TSV. Returns the count written."""
    fieldnames = [
        "seq_id",
        "locus_tag",
        "passes_threshold",
        "flagging_types",
        "max_stacking",
        "effector_type",
    ]
    n = 0
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            n += 1
    return n


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge PLM-Effector per-type TSVs into one combined file"
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Per-type PLM-Effector TSVs (one per SS type)",
    )
    parser.add_argument("--out", required=True, help="Merged output TSV")
    args = parser.parse_args()

    n = write_merged_tsv(merge_per_type_outputs(args.inputs), args.out)
    logger.info(f"Merged {len(args.inputs)} per-type files into {n} rows at {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
