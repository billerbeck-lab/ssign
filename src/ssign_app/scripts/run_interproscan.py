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
from dedup_sequences import deduplicate_dict, expand_results_dict
from ssign_lib.constants import TOOL_TIMEOUT_S
from ssign_lib.fasta_io import read_fasta
from ssign_lib.resources import resolve_threads
from ssign_lib.subprocess_diag import dump_failure_log
from ssign_lib.substrates import load_substrate_ids

# InterProScan TSV column indices (0-based, no header)
_COL_PROTEIN_ID = 0
_COL_SIG_ACC = 4
_COL_SIG_DESC = 5
_COL_IPR_ACC = 11
_COL_IPR_DESC = 12
_COL_GO_TERMS = 13
_MISSING = "-"
_GO_ID_RE = re.compile(r"(GO:\d+)")

_OUTPUT_FIELDNAMES = (
    "locus_tag",
    "interpro_domains",
    "interpro_go_terms",
    "interpro_pfam_ids",
    "interpro_descriptions",
)

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


def _resolve_interproscan_binary(install_dir):
    """Return an absolute path to interproscan.sh, or the bare name for PATH lookup.

    Resolution order:
      1. `install_dir` argument (from --interproscan-db CLI flag)
      2. SSIGN_INTERPROSCAN_PATH environment variable
      3. Fall back to "interproscan.sh" (PATH lookup)

    When invoked via the ssign runner, the env-var branch is unreachable —
    PipelineConfig.__post_init__ already copies SSIGN_INTERPROSCAN_PATH into
    `interproscan_db` and the runner passes that as `install_dir`. The branch
    only fires when this wrapper is invoked standalone (e.g. unit tests,
    manual debugging). Keep the precedence in sync with runner.py if either
    side changes.
    """
    for candidate in (install_dir, os.environ.get("SSIGN_INTERPROSCAN_PATH", "")):
        if candidate:
            binary = os.path.join(candidate, "interproscan.sh")
            if os.path.isfile(binary):
                return binary
    return "interproscan.sh"


def run_local_interproscan(
    query_fasta,
    install_dir,
    output_dir,
    applications="",
    offline=False,
    cpu=None,
    timeout=TOOL_TIMEOUT_S,
):
    """Run InterProScan locally and return path to the TSV output.

    `install_dir` is the InterProScan install directory containing
    `interproscan.sh`. Empty string falls back to SSIGN_INTERPROSCAN_PATH
    env var, then to a PATH lookup.

    By default, IPS queries the EBI precalculated-match lookup service
    (5-10× speedup on previously-seen sequences; falls back to local
    HMMs after a ~30 s connection timeout if unreachable). Pass
    `offline=True` to add `-dp` and skip the lookup outright — useful
    on isolated networks where the timeout would dominate runtime.

    `applications` restricts the member-DB scan; pass "" to use IPS
    defaults.

    `cpu` controls IPS's `-cpu N` flag (worker count). `None` (default)
    autodetects from `os.cpu_count()` so HPC nodes use all available
    cores; pass an explicit integer to override. Without this flag IPS
    falls back to its `interproscan.properties` defaults (typically 6
    workers regardless of node CPU count), which silently throttles
    runs on 16+ core nodes.
    """
    cpu = resolve_threads(cpu)
    output_file = os.path.join(output_dir, "results.tsv")
    binary = _resolve_interproscan_binary(install_dir)
    cmd = [
        binary,
        "-i",
        query_fasta,
        "-o",
        output_file,
        "-f",
        "tsv",
        "-goterms",
        "-pathways",
        "-cpu",
        str(cpu),
    ]
    # IPS writes per-analysis scratch (hmmsearch output files, etc.) under its temp
    # dir. Its host-configured temporary.file.directory can be an absolute path
    # (e.g. $HOME/temp) that is NOT mounted inside an --containall container, so
    # hmmsearch fails at "open output file for writing" and aborts (exit 134).
    # Point IPS at a container-writable local temp dir under the real-disk scratch
    # (TMPDIR); harmless natively (mkdtemp falls back to the system temp).
    ips_tempdir = tempfile.mkdtemp(prefix="ips_temp_", dir=(os.environ.get("TMPDIR") or None))
    cmd.extend(["-T", ips_tempdir])
    if offline:
        cmd.append("-dp")
    if applications:
        cmd.extend(["-appl", applications])

    logger.info("Running local InterProScan...")
    # FRAGILE: needs interproscan.sh — either at $install_dir/, via
    # $SSIGN_INTERPROSCAN_PATH, or on PATH.
    # If this breaks: install InterProScan from https://www.ebi.ac.uk/interpro/download/

    # Java heap auto-scale: IPS's interproscan.sh hard-codes -Xmx
    # (currently -Xmx15G) and does NOT read $JAVA_OPTS — confirmed
    # against the upstream launcher (ebi-pf-team/interproscan, the
    # bundled support-mini-x86-32/interproscan.sh). _JAVA_OPTIONS is the
    # one env var the JVM ALWAYS reads regardless of launcher script; it
    # prints a "Picked up _JAVA_OPTIONS" stderr notice and the last -Xmx
    # wins, so our value overrides the launcher's hard-coded one.
    #
    # Formula: (effective_ram_gb - 2) * 0.5. Reserves 2 GB for OS / IPS
    # Perl helpers / file cache, then halves the rest (the JVM also
    # needs metaspace and Code Cache on top of the heap). Clamped to
    # [4, 64] GB so a 256 GB box doesn't ask for 192 GB and a 16 GB
    # laptop doesn't crash because IPS asked for 12.
    from ssign_lib.resources import effective_ram_gb

    env = os.environ.copy()
    heap_gb = max(4, min(64, int((effective_ram_gb() - 2) * 0.5)))
    java_opt = f"-Xmx{heap_gb}g"
    existing = env.get("_JAVA_OPTIONS", "")
    env["_JAVA_OPTIONS"] = f"{existing} {java_opt}".strip() if existing else java_opt
    logger.info("InterProScan _JAVA_OPTIONS=%s", env["_JAVA_OPTIONS"])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"InterProScan (interproscan.sh) not found: {e}\n"
            f"  Common causes:\n"
            f"    - InterProScan is not installed\n"
            f"    - Wrong --interproscan-db path / SSIGN_INTERPROSCAN_PATH\n"
            f"    - interproscan.sh not on PATH\n"
            f"  How to fix:\n"
            f"    - Run: scripts/fetch_databases.sh --tier extended\n"
            f"    - Then pass --interproscan-db <install_dir> or export SSIGN_INTERPROSCAN_PATH=<install_dir>"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"InterProScan timed out after 4 hours: {e}\n  How to fix:\n    - Reduce the number of input sequences"
        ) from e

    if result.returncode != 0:
        raise dump_failure_log("InterProScan", result, cmd, output_dir)

    return output_file


def parse_interproscan_tsv(tsv_path, target_ids=None):
    """Parse InterProScan TSV output and aggregate per protein."""
    per_protein: dict[str, dict[str, set[str]]] = {}

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
        help=(
            "InterProScan install directory (contains interproscan.sh). "
            "Falls back to $SSIGN_INTERPROSCAN_PATH, then to PATH lookup."
        ),
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
    parser.add_argument(
        "--cpu",
        type=int,
        default=None,
        help=(
            "Worker count for InterProScan (-cpu N). When omitted, resolves "
            "to the scheduler-aware effective CPU count, divided by the "
            "parallel-group size if the runner launched this wrapper inside "
            "one. Pass an explicit value to override on shared machines."
        ),
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--local-cache-dir",
        default=None,
        help=(
            "Stage the InterProScan install directory (~24 GB) to this "
            "directory before running. No-op when on local filesystem; "
            "otherwise rsyncs from gpfs/nfs/lustre. Pass PBS/SLURM job-"
            "local TMPDIR. Mixed effect — HMMER members are mostly "
            "sequential, but ProSiteProfiles + CDD do random reads."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=TOOL_TIMEOUT_S,
        help="Subprocess timeout in seconds (the runner passes a size-scaled value; default is the flat floor)",
    )
    args = parser.parse_args()

    # Short-circuit on empty substrate input BEFORE the DB stage. Reduced-
    # genome endosymbionts (e.g. Buchnera) yield zero substrates; staging
    # ~10 GB of IPS data to local SSD only to throw it away is pure waste.
    # Previous bug here: empty FASTA -> IPS produced no results.tsv -> the
    # wrapper crashed with FileNotFoundError on parse.
    substrate_ids = load_substrate_ids(args.substrates)
    if not substrate_ids:
        logger.info(f"No substrates for {args.sample}; writing empty InterProScan output.")
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES)
            writer.writeheader()
        return

    if args.local_cache_dir and args.db:
        from ssign_lib.resources import stage_directory_tree_to_local_ssd_if_remote

        args.db = stage_directory_tree_to_local_ssd_if_remote(args.db, args.local_cache_dir)

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
            cpu=args.cpu,
            timeout=args.timeout,
        )
        results_unique = parse_interproscan_tsv(tsv_path, set(unique_seqs.keys()))

    results = expand_results_dict(results_unique, seq_groups)

    logger.info(f"Annotated {len(results)}/{len(substrate_ids)} substrates for {args.sample}")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES)
        writer.writeheader()
        for r in results.values():
            writer.writerow(r)


if __name__ == "__main__":
    main()
