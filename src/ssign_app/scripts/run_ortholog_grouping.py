#!/usr/bin/env python3
"""Assign ortholog groups to substrate proteins using all-vs-all BLASTp.

Algorithm:
1. Collect all filtered substrate sequences across genomes
2. Run all-vs-all BLASTp (local makeblastdb + blastp)
3. Filter hits: >= min_pident AND >= min_qcov
4. Single-linkage clustering via Union-Find
5. Output: CSV mapping locus_tag → ortholog_group + group stats

Requires NCBI BLAST+ (makeblastdb, blastp) on PATH. If not installed the
step is skipped gracefully (every substrate becomes its own singleton
group). All other failure modes (DB build, search, parse) raise.
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import TOOL_TIMEOUT_S
from ssign_lib.fasta_io import read_fasta

# BLAST outfmt 6 column indices, mirroring the order in BLAST_OUTFMT below.
BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
_COL_QSEQID = 0
_COL_SSEQID = 1
_COL_PIDENT = 2
_COL_ALN_LEN = 3
_COL_QLEN = 12
_BLAST_MIN_FIELDS = 14


# Soft-skip sentinel — distinguishes "BLAST+ not installed" (caller writes
# singleton-only output) from "BLAST+ ran but failed" (caller raises).
class BlastpUnavailableError(RuntimeError):
    """Raised when blastp/makeblastdb cannot be located on PATH."""


def _find_blast_binary(name: str) -> str:
    """Locate a BLAST+ binary on PATH.

    Raises BlastpUnavailableError if missing — caller decides whether to
    soft-skip (no installation present) or fail.
    """
    # FRAGILE: requires NCBI BLAST+ on PATH (blastp + makeblastdb).
    # If this breaks: install BLAST+:
    #   - sudo apt install ncbi-blast+   (Debian/Ubuntu)
    #   - conda install -c bioconda blast (conda)
    #   - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    try:
        subprocess.run([name, "-version"], capture_output=True, timeout=10, check=True)
    except FileNotFoundError as e:
        raise BlastpUnavailableError(f"{name} not found on PATH") from e
    except subprocess.CalledProcessError as e:
        # Binary exists but `-version` failed — corrupted install or
        # incompatible version. Fail loudly rather than masking as "not
        # installed."
        raise RuntimeError(
            f"{name} is on PATH but `{name} -version` exited with "
            f"code {e.returncode}. Likely a corrupted or incompatible "
            f"BLAST+ install."
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(f"{name} did not respond to `-version` within 10s — likely hung. Check the install.") from e
    return name


# Cap on BLAST+ -max_target_seqs. High enough that all-vs-all hit recovery
# is complete for any realistic substrate cohort (BLAST clamps to actual
# hits available, so this is not a per-query cost). See Shah 2019,
# Bioinformatics 35:1613 for why max_target_seqs is a footgun.
_MAX_TARGET_SEQS = 10000


def run_local_blast(
    fasta_path: str,
    min_pident: float,
    min_qcov: float,
    evalue: float = 1e-5,
    num_threads: int = 4,
):
    """Run all-vs-all BLASTp locally.

    Returns list of (query, subject, pident, qcov) tuples.

    Raises:
        BlastpUnavailableError: BLAST+ not installed (caller soft-skips).
        RuntimeError: BLAST+ ran but failed (caller fails the step).
    """
    blastp = _find_blast_binary("blastp")
    makeblastdb = _find_blast_binary("makeblastdb")

    with tempfile.TemporaryDirectory(prefix="ortho_blast_") as tmpdir:
        db_path = os.path.join(tmpdir, "substrates_db")
        out_path = os.path.join(tmpdir, "blast_results.txt")

        # FRAGILE: subprocess call requires BLAST+ (makeblastdb) on PATH.
        # If this breaks: see _find_blast_binary docstring.
        cmd_db = [makeblastdb, "-in", fasta_path, "-dbtype", "prot", "-out", db_path]
        result = subprocess.run(cmd_db, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(f"makeblastdb failed (rc={result.returncode}): {result.stderr[:500]}")

        cmd_blast = [
            blastp,
            "-query",
            fasta_path,
            "-db",
            db_path,
            "-out",
            out_path,
            "-outfmt",
            BLAST_OUTFMT,
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(_MAX_TARGET_SEQS),
            "-num_threads",
            str(num_threads),
        ]
        logger.info("Running all-vs-all BLASTp...")
        # FRAGILE: subprocess call requires BLAST+ (blastp) on PATH.
        # If this breaks: see _find_blast_binary docstring.
        result = subprocess.run(cmd_blast, capture_output=True, text=True, timeout=TOOL_TIMEOUT_S)
        if result.returncode != 0:
            raise RuntimeError(f"BLASTp failed (rc={result.returncode}): {result.stderr[:500]}")

        hits = []
        with open(out_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < _BLAST_MIN_FIELDS:
                    continue
                query, subject = parts[_COL_QSEQID], parts[_COL_SSEQID]
                if query == subject:
                    continue
                pident = float(parts[_COL_PIDENT])
                aln_len = int(parts[_COL_ALN_LEN])
                qlen = int(parts[_COL_QLEN])
                qcov = 100.0 * aln_len / max(qlen, 1)

                if pident >= min_pident and qcov >= min_qcov:
                    hits.append((query, subject, pident, qcov))

        logger.info(f"Found {len(hits)} ortholog-quality hits (>={min_pident}% id, >={min_qcov}% qcov)")
        return hits


def cluster_union_find(hits, all_protein_ids):
    """Single-linkage clustering via Union-Find. Returns dict of group_id → set of members."""
    parent = {pid: pid for pid in all_protein_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for query, subject, _, _ in hits:
        if query in parent and subject in parent:
            union(query, subject)

    groups = defaultdict(set)
    for pid in all_protein_ids:
        groups[find(pid)].add(pid)

    return groups


def compute_group_stats(groups, hits, all_protein_ids):
    """Compute per-group statistics. Returns list of dicts."""
    hit_identities = defaultdict(list)
    for query, subject, pident, _ in hits:
        hit_identities[(query, subject)].append(pident)

    stats = []
    for idx, (rep, members) in enumerate(sorted(groups.items(), key=lambda x: -len(x[1])), 1):
        within_ids = []
        member_list = sorted(members)
        for i, m1 in enumerate(member_list):
            for m2 in member_list[i + 1 :]:
                for key in [(m1, m2), (m2, m1)]:
                    if key in hit_identities:
                        within_ids.extend(hit_identities[key])

        mean_id = sum(within_ids) / len(within_ids) if within_ids else 100.0

        stats.append(
            {
                "ortholog_group": f"OG_{idx:03d}",
                "n_members": len(members),
                "members": ";".join(sorted(members)),
                "mean_pident": round(mean_id, 1),
            }
        )

    return stats


_OUTPUT_FIELDS = ["locus_tag", "ortholog_group", "og_n_members", "og_mean_pident"]


def _write_singleton_output(seqs, output_path):
    """Write each protein as its own singleton group (used when BLAST is unavailable)."""
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDS)
        writer.writeheader()
        for idx, pid in enumerate(sorted(seqs.keys()), 1):
            writer.writerow(
                {
                    "locus_tag": pid,
                    "ortholog_group": f"OG_{idx:03d}",
                    "og_n_members": 1,
                    "og_mean_pident": 100.0,
                }
            )


def main():
    parser = argparse.ArgumentParser(description="Assign ortholog groups via all-vs-all BLASTp")
    parser.add_argument(
        "--substrates-fasta", required=True, help="FASTA of all substrate proteins (combined across genomes)"
    )
    parser.add_argument(
        "--min-pident", type=float, default=40.0, help="Minimum %% identity for ortholog assignment (default: 40)"
    )
    parser.add_argument(
        "--min-qcov", type=float, default=70.0, help="Minimum query coverage %% for ortholog assignment (default: 70)"
    )
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for BLASTp (default: 1e-5)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for BLASTp -num_threads (default: 4)")
    parser.add_argument("--output", required=True, help="Output CSV with ortholog group assignments")
    parser.add_argument("--output-groups", default="", help="Optional: output CSV with group-level statistics")
    args = parser.parse_args()

    seqs = read_fasta(args.substrates_fasta)
    if not seqs:
        logger.warning("No sequences in input FASTA — writing empty output")
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDS)
            writer.writeheader()
        return 0

    logger.info(f"Loaded {len(seqs)} substrate sequences")

    if len(seqs) < 2:
        logger.info("Only 1 substrate — assigning singleton group OG_001")
        _write_singleton_output(seqs, args.output)
        return 0

    try:
        hits = run_local_blast(
            args.substrates_fasta,
            args.min_pident,
            args.min_qcov,
            args.evalue,
            args.threads,
        )
    except BlastpUnavailableError as e:
        logger.warning(
            "BLAST+ not installed — skipping ortholog grouping (every substrate becomes its own singleton group).\n"
            "  To enable ortholog grouping: install NCBI BLAST+:\n"
            "    - sudo apt install ncbi-blast+   (Debian/Ubuntu)\n"
            "    - conda install -c bioconda blast (conda)\n"
            f"  Detail: {e}"
        )
        _write_singleton_output(seqs, args.output)
        return 0

    all_ids = set(seqs.keys())
    groups = cluster_union_find(hits, all_ids)
    logger.info(f"Clustered into {len(groups)} ortholog groups")

    group_stats = compute_group_stats(groups, hits, all_ids)

    pid_to_group = {}
    pid_to_stats = {}
    for gs in group_stats:
        for member in gs["members"].split(";"):
            pid_to_group[member] = gs["ortholog_group"]
            pid_to_stats[member] = gs

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDS)
        writer.writeheader()
        for pid in sorted(seqs.keys()):
            og = pid_to_group.get(pid, "")
            stats = pid_to_stats.get(pid, {})
            writer.writerow(
                {
                    "locus_tag": pid,
                    "ortholog_group": og,
                    "og_n_members": stats.get("n_members", 1),
                    "og_mean_pident": stats.get("mean_pident", 100.0),
                }
            )

    logger.info(f"Wrote {len(seqs)} protein assignments to {args.output}")

    if args.output_groups:
        with open(args.output_groups, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["ortholog_group", "n_members", "members", "mean_pident"])
            writer.writeheader()
            for gs in group_stats:
                writer.writerow(gs)
        logger.info(f"Wrote {len(group_stats)} group stats to {args.output_groups}")

    sizes = [gs["n_members"] for gs in group_stats]
    n_singleton = sum(1 for s in sizes if s == 1)
    n_multi = sum(1 for s in sizes if s > 1)
    logger.info(f"Groups: {n_singleton} singletons, {n_multi} multi-member (largest: {max(sizes)} members)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
