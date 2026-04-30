#!/usr/bin/env python3
"""Run HH-suite locally: hhblits for MSA generation, then hhsearch against
Pfam-A and PDB70."""

import argparse
import csv
import logging
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import (
    HHBLITS_ITERATIONS,
    HHBLITS_TIMEOUT_S,
    HHSEARCH_TIMEOUT_S,
    HHSUITE_MIN_PROB,
)
from ssign_lib.fasta_io import read_fasta
from ssign_lib.substrates import load_substrate_ids

# Regex for HHR hit table lines — works for both PDB (desc ends with ;) and Pfam (multiple ;)
# Fields (from right): (TLen) TRange QRange Cols SS Score Pval Eval Prob
_HHR_HIT_RE = re.compile(
    r"^\s*\d+\s+(\S+)\s+(.*?)\s+"
    r"(\d[\d.E+-]+)\s+"  # Prob
    r"([\d.E+-]+)\s+"  # E-value
    r"[\d.E+-]+\s+"  # P-value
    r"([\d.E+-]+)\s+"  # Score
    r"[\d.E+-]+\s+"  # SS
    r"\d+\s+"  # Cols
    r"\d+-\d+\s+"  # Query HMM range
    r"\d+-\d+\s+"  # Template HMM range
    r"\(\d+\)\s*$"  # Template length
)


def _run_one(
    pid,
    seq,
    pfam_db,
    pdb70_db,
    uniclust_db,
    output_dir,
    cpu_per_job=2,
    min_prob=HHSUITE_MIN_PROB,
):
    """Run hhblits + hhsearch for a single protein. Returns {} on per-protein failure.

    Raises RuntimeError only for binary-not-found (whole-pipeline fatal).
    """
    query_file = os.path.join(output_dir, f"{pid}.fasta")
    with open(query_file, "w") as f:
        f.write(f">{pid}\n{seq}\n")

    a3m_file = os.path.join(output_dir, f"{pid}.a3m")
    cmd_hhblits = [
        "hhblits",
        "-i",
        query_file,
        "-d",
        uniclust_db,
        "-n",
        str(HHBLITS_ITERATIONS),
        "-cpu",
        str(cpu_per_job),
        # -o /dev/null: discard hhblits' text report; we only need the .a3m
        "-o",
        "/dev/null",
        "-oa3m",
        a3m_file,
    ]

    # FRAGILE: subprocess requires hhblits on PATH.
    try:
        subprocess.run(
            cmd_hhblits,
            capture_output=True,
            text=True,
            timeout=HHBLITS_TIMEOUT_S,
            check=True,
        )
    except FileNotFoundError as e:
        raise RuntimeError(
            f"hhblits binary not found: {e}\n"
            f"  How to fix:\n"
            f"    - Debian/Ubuntu: sudo apt install hhsuite\n"
            f"    - Conda:         conda install -c bioconda hhsuite"
        ) from e
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        # type name disambiguates "DB too slow" (TimeoutExpired) from
        # "DB corrupt / missing" (CalledProcessError) when debugging.
        logger.warning(f"hhblits {type(e).__name__} for {pid}: {e}")
        return {}

    entry = {}
    for db_name, db_path in [("pfam", pfam_db), ("pdb", pdb70_db)]:
        if not db_path:
            continue
        hhr_file = os.path.join(output_dir, f"{pid}_{db_name}.hhr")
        cmd_hhsearch = [
            "hhsearch",
            "-i",
            a3m_file,
            "-d",
            db_path,
            "-cpu",
            str(cpu_per_job),
            "-o",
            hhr_file,
        ]
        try:
            subprocess.run(
                cmd_hhsearch,
                capture_output=True,
                text=True,
                timeout=HHSEARCH_TIMEOUT_S,
                check=True,
            )
            entry.update(parse_hhr(hhr_file, db_name, min_prob=min_prob))
        except FileNotFoundError as e:
            raise RuntimeError(
                f"hhsearch binary not found: {e}\n"
                f"  How to fix:\n"
                f"    - Debian/Ubuntu: sudo apt install hhsuite\n"
                f"    - Conda:         conda install -c bioconda hhsuite"
            ) from e
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.warning(f"hhsearch {db_name} {type(e).__name__} for {pid}: {e}")

    if entry:
        entry["locus_tag"] = pid
    return entry


def run_hhsuite_parallel(
    sequences,
    pfam_db,
    pdb70_db,
    uniclust_db,
    output_dir,
    max_workers=4,
    cpu_per_job=2,
    min_prob=HHSUITE_MIN_PROB,
):
    """Run hhblits + hhsearch across all proteins in a thread pool.

    ThreadPool (not Process): hhblits/hhsearch are subprocess calls, so the GIL
    isn't contended — threads give the same parallelism as processes with less
    overhead. max_workers × cpu_per_job should not exceed available cores.

    `min_prob` is the HHR-Prob cutoff for keeping the top-1 hit per database.
    """
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(
                _run_one,
                pid,
                seq,
                pfam_db,
                pdb70_db,
                uniclust_db,
                output_dir,
                cpu_per_job,
                min_prob,
            ): pid
            for pid, seq in sequences.items()
        }
        for fut in as_completed(futures):
            pid = futures[fut]
            entry = fut.result()  # RuntimeError (binary missing) propagates
            if entry:
                results[pid] = entry
    return results


def parse_hhr(hhr_path, db_prefix, min_prob=HHSUITE_MIN_PROB):
    """Parse HHR result file and extract top hit if Prob >= min_prob.

    Returns {} when the file is missing, the table is empty, or the top
    hit's Prob falls below the cutoff. Filtering on Prob (not E-value)
    matches Söding-lab guidance for remote homology.
    """
    if not os.path.exists(hhr_path):
        return {}

    with open(hhr_path) as f:
        in_hits = False
        for line in f:
            if line.startswith(" No Hit"):
                in_hits = True
                continue
            if in_hits and line.strip():
                m = _HHR_HIT_RE.match(line)
                if m:
                    hit_id, desc, prob, evalue, score = m.group(1, 2, 3, 4, 5)
                    prob_f = float(prob)
                    if prob_f < min_prob:
                        return {}
                    return {
                        f"{db_prefix}_top1_id": hit_id,
                        f"{db_prefix}_top1_description": desc.strip("; ")[:200],
                        f"{db_prefix}_top1_probability": prob_f,
                        f"{db_prefix}_top1_evalue": float(evalue),
                        f"{db_prefix}_top1_score": float(score),
                    }
                break
    return {}


def main():
    parser = argparse.ArgumentParser(description="Run HH-suite locally")
    parser.add_argument("--substrates", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pfam-db", default="")
    parser.add_argument("--pdb70-db", default="")
    parser.add_argument(
        "--uniclust-db",
        required=True,
        help="Path to UniRef30/UniClust30 HH-suite database (required for hhblits MSA)",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--min-prob",
        type=float,
        default=HHSUITE_MIN_PROB,
        help=(
            f"Minimum HHR Prob (0-100) to keep the top-1 hit per DB. "
            f"Defaults to {HHSUITE_MIN_PROB} (Söding-lab guidance: ≥95 "
            f"near-certain, ≥50 worth considering)."
        ),
    )
    args = parser.parse_args()

    if not args.pfam_db and not args.pdb70_db:
        parser.error("At least one of --pfam-db or --pdb70-db must be provided")

    substrate_ids = load_substrate_ids(args.substrates)
    all_seqs = read_fasta(args.proteins)
    sub_seqs = {k: v for k, v in all_seqs.items() if k in substrate_ids}

    logger.info(f"Processing {len(sub_seqs)} substrate proteins for {args.sample}")

    with tempfile.TemporaryDirectory() as tmpdir:
        results = run_hhsuite_parallel(
            sub_seqs,
            args.pfam_db,
            args.pdb70_db,
            args.uniclust_db,
            tmpdir,
            min_prob=args.min_prob,
        )

    # Determine output columns from results
    all_cols = set()
    for r in results.values():
        all_cols.update(r.keys())
    fieldnames = ["locus_tag"] + sorted(c for c in all_cols if c != "locus_tag")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results.values():
            writer.writerow(r)

    logger.info(f"Wrote {len(results)} HH-suite annotations for {args.sample}")


if __name__ == "__main__":
    main()
