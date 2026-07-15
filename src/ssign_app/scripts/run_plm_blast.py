#!/usr/bin/env python3
"""Run pLM-BLAST remote-homology search and convert output to ssign format.

pLM-BLAST (Kaminski et al. 2023, https://github.com/labstructbioinf/pLM-BLAST)
detects remote homologues by comparing ProtT5 residue embeddings. We use
it with the precomputed ECOD30 database to annotate secreted proteins with
ECOD structural domains, complementing the sequence-homology searches
BLAST and HH-suite provide.

Default DB is ECOD30 (10 GB, 30%-identity clustered representatives) and
default `--cpc` is 90 — the setting the Kaminski 2023 paper benchmarks at.
ECOD30 still has ≥1 representative per F-group, so annotation labels are
unchanged from ECOD70; it just drops within-family near-duplicates that
the embedding step doesn't need.

Usage:
    run_plm_blast.py --input proteins.faa --ecod-db /path/to/ecod30 \\
        --out plm_blast.tsv [--threads 4] [--cpc 90]

The `--ecod-db` argument points at a pre-built pLM-BLAST database directory
fetched from
    http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs
Output is a TSV mapping each query protein to its top pLM-BLAST hit with
score and alignment coordinates (one row per query, top-1 by score).
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import PLMBLAST_TIMEOUT_S  # noqa: E402
from ssign_lib.resources import resolve_threads  # noqa: E402
from ssign_lib.subprocess_diag import dump_failure_log  # noqa: E402
from ssign_lib.substrates import (  # noqa: E402
    load_substrate_ids,
    write_substrates_only_fasta,
)

# Column schema emitted to ssign. Matches the contract that
# integrate_annotations.TOOL_HIT_COLUMNS["pLM-BLAST"] expects
# (`ecod_top1_description`) and uses `locus_tag` as the join key so the
# merge against the substrates table works. Pre-fix the wrapper emitted
# `protein_id, target_id, score, ...` and produced multiple rows per
# query, which caused integrate_annotations to silently drop the file
# (no `locus_tag` column, and join-column heuristic resolved to
# `protein_id` which the substrate table did not have) — see task #80.
_OUTPUT_FIELDNAMES = [
    "locus_tag",
    "ecod_top1_id",
    "ecod_top1_description",
    "ecod_top1_score",
    "ecod_top1_qstart",
    "ecod_top1_qend",
    "ecod_top1_tstart",
    "ecod_top1_tend",
]


def _write_empty_output(out_path):
    """Write a header-only TSV when there are no substrates to search."""
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()


# pLM-BLAST output CSV columns (v1.x). Exact header names verified against
# the upstream `scripts/plmblast.py` result-writer on first integration
# run — see the TODO below if upstream renames columns. `sdesc` is added
# by upstream `alntools/postprocess/format.py:prepare_output` when the DB
# CSV has a `description` column (true for the public ECOD DB).
# Sibling-file extensions plmblast.py probes for next to the DB directory.
# `.csv/.p/.pkl/.fas/.fasta` are the metadata index (alntools/settings.py:EXTENSIONS,
# `DataObject.from_dir` picks the first match); `.pt` is the optional pooled
# embedding file `DataObject._find_datatype` probes independently. The
# metadata file is required; `.pt` is optional but speeds up loading when
# present. Both can co-exist, so we copy every match we find.
_PLMBLAST_DB_SIDECAR_EXTS = (".csv", ".p", ".pkl", ".fas", ".fasta", ".pt")

_COL_QUERY = "qid"
_COL_TARGET = "sid"
_COL_SDESC = "sdesc"
_COL_SCORE = "score"
_COL_QSTART = "qstart"
_COL_QEND = "qend"
_COL_TSTART = "tstart"
_COL_TEND = "tend"


def find_plmblast_script() -> str | None:
    """Return the resolved path to pLM-BLAST's `plmblast.py`, or None.

    Resolution order (same as `_resolve_plmblast_script` below, but
    non-raising so callers can use it as a pre-flight availability check):
      1. `SSIGN_PLMBLAST_SCRIPT` env var.
      2. `plmblast.py` / `plmblast` on PATH.
    """
    override = os.environ.get("SSIGN_PLMBLAST_SCRIPT")
    if override and os.path.isfile(override):
        return override
    for candidate in ("plmblast.py", "plmblast"):
        found = shutil.which(candidate)
        if found:
            return found
    return None


def _resolve_plmblast_script() -> str:
    """Locate the `plmblast.py` script pLM-BLAST ships.

    Raises RuntimeError with install hints if not found. For a non-raising
    availability check use `find_plmblast_script()`.
    """
    override = os.environ.get("SSIGN_PLMBLAST_SCRIPT")
    if override and not os.path.exists(override):
        raise RuntimeError(f"SSIGN_PLMBLAST_SCRIPT points at a missing file: {override}")
    found = find_plmblast_script()
    if found:
        return found
    return "plmblast.py"


def _use_cuda_for_embedding() -> bool:
    """Whether to pass `--cuda` to pLM-BLAST's `embeddings.py`.

    Auto-detects CUDA via torch, falling back to CPU if torch is missing
    or no GPU is visible. `SSIGN_PLMBLAST_FORCE_CPU=1|true|yes` overrides
    to CPU. Without `--cuda`, ProtT5 embedding runs ~100x slower than on
    GPU — silently — so the flag is worth getting right.
    """
    if os.environ.get("SSIGN_PLMBLAST_FORCE_CPU", "").lower() in ("1", "true", "yes"):
        return False
    try:
        import torch

        return torch.cuda.is_available()
    except ImportError:
        return False


def _resolve_embeddings_script(plmblast_script: str) -> str:
    """Locate pLM-BLAST's `embeddings.py` (one level above `scripts/plmblast.py`)."""
    # `_resolve_plmblast_script` returns the literal "plmblast.py" when
    # neither SSIGN_PLMBLAST_SCRIPT nor a PATH entry resolves. Catch that
    # here so we emit a clear "not installed" error instead of a
    # misleading "embeddings.py not found at <cwd>/.." message.
    if plmblast_script == "plmblast.py" and not os.path.isabs(plmblast_script):
        raise RuntimeError(
            "pLM-BLAST is not installed: plmblast.py is not on PATH and "
            "SSIGN_PLMBLAST_SCRIPT is not set.\n"
            "  Install:\n"
            "    git clone https://github.com/labstructbioinf/pLM-BLAST.git ~/tools/pLM-BLAST\n"
            "    export SSIGN_PLMBLAST_SCRIPT=~/tools/pLM-BLAST/scripts/plmblast.py\n"
            "  Or pass --skip-plmblast to disable this step."
        )
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(plmblast_script)))
    candidate = os.path.join(repo_root, "embeddings.py")
    if os.path.exists(candidate):
        return candidate
    raise RuntimeError(
        f"embeddings.py not found next to plmblast.py at {repo_root}. "
        f"pLM-BLAST install appears incomplete; re-clone from "
        f"https://github.com/labstructbioinf/pLM-BLAST"
    )


def _embed_query_fasta(
    embed_script: str,
    proteins_fasta: str,
    out_base: str,
    embedder: str = "pt",
    failure_log_dir: str = "",
    timeout: int = PLMBLAST_TIMEOUT_S,
) -> None:
    """Run pLM-BLAST's embeddings.py to embed a FASTA, producing a
    single pooled `<base>.pt` file plus a sibling `<base>.csv` index.

    `embedder` defaults to `pt` (ProtT5-XL UniRef50), matching the
    embedder used to build the public ECOD databases. ESM-2 etc. will
    only work against a database embedded with that same model.

    NB: embeddings.py writes the index to `<base>.pt.csv`, but
    pLM-BLAST's search-side `_find_datatype()` looks for `<base>.csv`.
    We rename it after embedding to bridge the gap.
    """
    # embeddings.py only accepts .csv/.p/.pkl/.fas/.fasta extensions —
    # NOT .faa, which is what ssign uses internally and what most
    # protein-FASTA tools emit. Symlink to a .fasta path so the upstream
    # extension check passes; the file content is identical.
    if proteins_fasta.endswith((".fasta", ".fas")):
        embed_input = proteins_fasta
    else:
        embed_input = os.path.join(os.path.dirname(out_base), "_proteins_for_embed.fasta")
        os.symlink(os.path.abspath(proteins_fasta), embed_input)

    pt_path = f"{out_base}.pt"
    # Route through _plm_blast_embed_bootstrap.py so torch.distributed /
    # torch.multiprocessing get stubbed if absent — see that script's
    # docstring for the why. Without the bootstrap, stripped torch wheels
    # (some HPC modules, some CPU-only builds) crash at module load before
    # any embedding work happens.
    bootstrap = os.path.join(_scripts_dir, "_plm_blast_embed_bootstrap.py")
    cmd = [
        sys.executable,
        bootstrap,
        embed_script,
        "start",
        embed_input,
        pt_path,
        "-embedder",
        embedder,
    ]
    use_cuda = _use_cuda_for_embedding()
    if use_cuda:
        cmd.append("--cuda")
        logger.info("pLM-BLAST embedding: CUDA detected, passing --cuda to embeddings.py")
    else:
        logger.info("pLM-BLAST embedding: running on CPU (~100x slower than GPU)")
    logger.info(f"Embedding query: embeddings.py start {proteins_fasta} {pt_path} -embedder {embedder}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"pLM-BLAST embeddings.py not runnable: {e}\n"
            f"  How to fix:\n"
            f"    - pip install git+https://github.com/labstructbioinf/pLM-BLAST.git\n"
            f"    - Or clone the repo and set SSIGN_PLMBLAST_SCRIPT to its scripts/plmblast.py"
        ) from e
    except subprocess.TimeoutExpired as e:
        if use_cuda:
            gpu_hint = (
                "    - Embedding ran on GPU but still timed out — check VRAM is "
                "sufficient for the batch size and reduce query size if needed"
            )
        else:
            gpu_hint = (
                "    - GPU not detected — install CUDA-enabled torch and verify "
                "with `python -c 'import torch; print(torch.cuda.is_available())'`\n"
                "    - Or unset SSIGN_PLMBLAST_FORCE_CPU if it's overriding auto-detect"
            )
        raise RuntimeError(
            f"pLM-BLAST embedding timed out after {timeout}s: {e}\n  How to fix:\n    - Reduce query size\n{gpu_hint}"
        ) from e
    if result.returncode != 0:
        # out_base lives in a TemporaryDirectory the parent will nuke
        # before the user sees the RuntimeError; prefer the persistent
        # dir threaded in from main().
        log_dir = failure_log_dir or os.path.dirname(out_base) or "."
        raise dump_failure_log("pLM-BLAST embedding", result, cmd, log_dir)

    # Rename the index file (see docstring NB).
    src_csv = f"{pt_path}.csv"
    dst_csv = f"{out_base}.csv"
    if not os.path.exists(src_csv):
        raise RuntimeError(
            f"embeddings.py exited 0 but did not produce {src_csv}; pLM-BLAST output layout may have changed upstream"
        )
    os.rename(src_csv, dst_csv)


def run_plmblast(
    proteins_fasta: str,
    ecod_db: str,
    out_csv: str,
    cpc: int = 70,
    threads: int | None = None,
    failure_log_dir: str = "",
    timeout: int = PLMBLAST_TIMEOUT_S,
) -> str:
    """Run pLM-BLAST against the ECOD DB and return the CSV output path.

    Two-step pipeline: ProtT5-embed the query FASTA into a single
    pooled .pt file, then search the embedding DB. Embedding dominates
    wall time on CPU (~5-10 sec per 500-aa protein).
    """
    threads = resolve_threads(threads)
    script = _resolve_plmblast_script()
    embed_script = _resolve_embeddings_script(script)

    with tempfile.TemporaryDirectory() as tmpdir:
        query_emb_base = os.path.join(tmpdir, "query")
        _embed_query_fasta(
            embed_script,
            proteins_fasta,
            query_emb_base,
            failure_log_dir=failure_log_dir,
            timeout=timeout,
        )

        # plmblast.py uses single-dash long flags (`-cpc`, `-workers`)
        # per its argparse setup in alntools/parser.py.
        cmd = [
            sys.executable,
            script,
            ecod_db,
            query_emb_base,
            out_csv,
            "-cpc",
            str(cpc),
            "-workers",
            str(threads),
        ]

        # scripts/plmblast.py does `import alntools` but alntools is a sibling
        # directory at the pLM-BLAST repo root, not pip-installable. Python's
        # subprocess only auto-adds the script's own directory to sys.path, so
        # without this env override the import fails. Setting PYTHONPATH to
        # the repo root lets the subprocess resolve alntools.
        env = os.environ.copy()
        plm_blast_root = os.path.dirname(os.path.dirname(os.path.abspath(script)))
        existing_pp = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = f"{plm_blast_root}{os.pathsep}{existing_pp}" if existing_pp else plm_blast_root

        logger.info(f"Running pLM-BLAST: {' '.join(cmd[:3])} <db> <query_emb> {out_csv} -cpc {cpc}")
        # FRAGILE: subprocess call requires pLM-BLAST's scripts/plmblast.py
        # on PATH, or set SSIGN_PLMBLAST_SCRIPT to its absolute path.
        # If this breaks: pip install git+https://github.com/labstructbioinf/pLM-BLAST.git
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        except FileNotFoundError as e:
            raise RuntimeError(
                f"pLM-BLAST script not found: {e}\n"
                f"  Common causes:\n"
                f"    - pLM-BLAST is not installed\n"
                f"    - plmblast.py is not on PATH\n"
                f"  How to fix:\n"
                f"    - pip install git+https://github.com/labstructbioinf/pLM-BLAST.git\n"
                f"    - Or set SSIGN_PLMBLAST_SCRIPT=/path/to/pLM-BLAST/scripts/plmblast.py"
            ) from e
        except subprocess.TimeoutExpired as e:
            raise RuntimeError(
                f"pLM-BLAST timed out after 4 hours: {e}\n"
                f"  How to fix:\n"
                f"    - Reduce query size\n"
                f"    - Increase --threads if more CPUs are available\n"
                f"    - Use a smaller ECOD subset (e.g. ECOD30 instead of ECOD50/70)"
            ) from e

        if result.returncode != 0:
            # out_csv is inside the tmpdir contextmanager we're still in;
            # the caller's failure_log_dir (args.out's parent) survives.
            log_dir = failure_log_dir or os.path.dirname(out_csv) or "."
            raise dump_failure_log("pLM-BLAST", result, cmd, log_dir)

        if not os.path.exists(out_csv):
            raise FileNotFoundError(f"pLM-BLAST output not found: {out_csv}")

    return out_csv


def parse_plmblast_csv(csv_path: str):
    """Parse pLM-BLAST CSV output into a list of hit dicts.

    Returns a list keyed per-hit (a protein can have multiple hits) with
    fields: protein_id, target_id, description, score, qstart, qend,
    tstart, tend. `description` is sourced from pLM-BLAST's `sdesc`
    column when present (added by upstream when the DB has a description
    column — true for the public ECOD DB); otherwise empty.

    Top-1-per-query reduction happens in `_reduce_to_top1` before the
    ssign-facing output is written. Callers that want the full hit list
    (e.g. the per-tool integration test) use this function directly.
    """
    entries = []
    with open(csv_path) as f:
        # pLM-BLAST emits ``;``-separated rows (verified upstream 2026-06-04
        # against ECOD30; header is ``qid;score;ident;similarity;sid;...``).
        # Reading with the csv default (``,``) silently collapses every row
        # into one field and drops every hit — the failure mode that hit
        # task #80 in spirit, masked there by an earlier staging crash.
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            query_id = row.get(_COL_QUERY, "").strip()
            if not query_id:
                continue
            entries.append(
                {
                    "protein_id": query_id,
                    "target_id": row.get(_COL_TARGET, "").strip(),
                    "description": row.get(_COL_SDESC, "").strip(),
                    "score": row.get(_COL_SCORE, "").strip(),
                    "qstart": row.get(_COL_QSTART, "").strip(),
                    "qend": row.get(_COL_QEND, "").strip(),
                    "tstart": row.get(_COL_TSTART, "").strip(),
                    "tend": row.get(_COL_TEND, "").strip(),
                }
            )
    return entries


def _entry_score(entry: dict) -> float:
    """Parse an entry's score to float; NaN/inf and unparseable values
    become 0.0 so they can't latch as "best" via NaN's weird comparison
    semantics (NaN > anything is False, which would freeze the first
    NaN-scored hit in place even when real-scored hits arrive later)."""

    try:
        score = float(entry.get("score", "") or 0)
    except ValueError:
        return 0.0
    return score if math.isfinite(score) else 0.0


def _reduce_to_top1(entries):
    """Reduce a per-hit entry list to one row per query, picking the
    highest-scoring hit. pLM-BLAST emits ranked-but-unsorted hits and
    can list several entries per query; integrate_annotations.py needs
    one row per locus_tag for the left-join to preserve row count."""
    best: dict[str, dict] = {}
    for e in entries:
        pid = e["protein_id"]
        prev = best.get(pid)
        if prev is None or _entry_score(e) > _entry_score(prev):
            best[pid] = e
    return list(best.values())


def _stage_db_sidecar(src_db: str, staged_db: str) -> None:
    """Copy every pLM-BLAST DB sibling file from src to staged.

    pLM-BLAST stores its DB as a directory (`<DB>/`) plus one or more
    sibling files (`<DB>.csv`, `<DB>.pt`, ...) at the same level. When
    `stage_directory_tree_to_local_ssd_if_remote` copies the directory
    to local SSD, the siblings are left behind on the network mount and
    `DataObject.from_dir(staged_db)` crashes trying to read the metadata.
    """
    src_parent = os.path.dirname(src_db.rstrip("/"))
    staged_parent = os.path.dirname(staged_db.rstrip("/"))
    basename = os.path.basename(src_db.rstrip("/"))
    copied = 0
    for ext in _PLMBLAST_DB_SIDECAR_EXTS:
        src_sidecar = os.path.join(src_parent, basename + ext)
        if os.path.exists(src_sidecar):
            dst_sidecar = os.path.join(staged_parent, basename + ext)
            shutil.copy2(src_sidecar, dst_sidecar)
            logger.info(f"Staged DB sidecar: {basename + ext} -> {staged_parent}")
            copied += 1
    if copied == 0:
        exts = "/".join(_PLMBLAST_DB_SIDECAR_EXTS)
        logger.warning(
            f"No DB sidecar (expected one of {basename}{{{exts}}} in {src_parent}); "
            f"plmblast.py will likely fail loading the DB index"
        )


def _to_output_row(entry: dict) -> dict:
    """Rename a parse_plmblast_csv entry to the ssign output schema."""
    return {
        "locus_tag": entry["protein_id"],
        "ecod_top1_id": entry["target_id"],
        "ecod_top1_description": entry.get("description", ""),
        "ecod_top1_score": entry["score"],
        "ecod_top1_qstart": entry["qstart"],
        "ecod_top1_qend": entry["qend"],
        "ecod_top1_tstart": entry["tstart"],
        "ecod_top1_tend": entry["tend"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run pLM-BLAST against the ECOD database (default: ECOD30) on "
            "the filtered substrates and convert output to ssign format. "
            "pLM-BLAST is an annotation-tier tool — it runs only on the "
            "~50 filtered substrates, not the full ~5000-protein genome."
        )
    )
    parser.add_argument(
        "--substrates",
        required=True,
        help="Filtered substrates TSV (columns include locus_tag)",
    )
    parser.add_argument(
        "--proteins",
        required=True,
        help="Full protein FASTA (pre-called by Bakta / Pyrodigal)",
    )
    parser.add_argument(
        "--ecod-db",
        required=True,
        help=(
            "Path to a pLM-BLAST ECOD database directory (any clustering "
            "level — ECOD30 default, ECOD50/70/90 also available from "
            "ftp.tuebingen.mpg.de)"
        ),
    )
    parser.add_argument("--out", required=True, help="Output hits TSV (ssign format)")
    parser.add_argument(
        "--local-cache-dir",
        default=None,
        help=(
            "Stage the ECOD DB tree (~11 GB for ECOD30, up to ~25 GB for "
            "ECOD90) to this directory before running. No-op when the DB "
            "is on local filesystem; otherwise rsyncs from gpfs/nfs/"
            "lustre. Pass PBS/SLURM job-local TMPDIR. ProtT5 embedding "
            "does heavy random I/O on the embedding shards; warm page "
            "cache cuts wallclock substantially on shared FS."
        ),
    )
    parser.add_argument(
        "--cpc",
        type=int,
        default=90,
        help=(
            "pLM-BLAST -cpc cosine percentile cutoff (default: 90). The "
            "Kaminski 2023 paper benchmarks at 90; the upstream argparse "
            "default of 70 is more conservative but ~3x slower with little "
            "annotation gain for ssign's use case."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "CPU threads (passed to plmblast.py -workers). When omitted, "
            "resolves to the scheduler-aware effective CPU count, divided "
            "by the parallel-group size if the runner launched this "
            "wrapper inside one."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=PLMBLAST_TIMEOUT_S,
        help="Subprocess timeout in seconds (the runner passes a size-scaled "
        "value; default is the pLM-BLAST 24h floor)",
    )
    args = parser.parse_args()

    if not os.path.exists(args.substrates):
        print(f"ERROR: Substrates TSV not found: {args.substrates}", file=sys.stderr)
        return 2
    if not os.path.exists(args.proteins):
        print(f"ERROR: Protein FASTA not found: {args.proteins}", file=sys.stderr)
        return 2
    if not os.path.isdir(args.ecod_db):
        print(
            f"ERROR: ECOD database directory not found: {args.ecod_db}\n"
            f"  Fetch from http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs",
            file=sys.stderr,
        )
        return 2

    # 30 GB free floor covers all four ECOD variants (ECOD30 ~11 GB
    # through ECOD90 ~25 GB) with headroom; the helper logs a warning
    # and skips staging if the cache dir doesn't have it.
    if args.local_cache_dir:
        from ssign_lib.resources import stage_directory_tree_to_local_ssd_if_remote

        src_db = args.ecod_db
        args.ecod_db = stage_directory_tree_to_local_ssd_if_remote(src_db, args.local_cache_dir, min_free_gb=30.0)
        if args.ecod_db != src_db:
            _stage_db_sidecar(src_db, args.ecod_db)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)

    substrate_ids = load_substrate_ids(args.substrates)
    if not substrate_ids:
        logger.info("No substrates to search — writing empty output")
        _write_empty_output(args.out)
        return 0

    with tempfile.TemporaryDirectory() as tmpdir:
        filtered_fasta = os.path.join(tmpdir, "substrates.faa")
        n = write_substrates_only_fasta(args.proteins, substrate_ids, filtered_fasta)
        if n == 0:
            logger.warning(
                f"No substrate proteins found in {args.proteins} (expected {len(substrate_ids)}); writing empty output"
            )
            _write_empty_output(args.out)
            return 0
        logger.info(f"pLM-BLAST will search {n} substrate proteins (of {len(substrate_ids)} listed)")

        raw_csv = os.path.join(tmpdir, "plm_blast_raw.csv")
        try:
            run_plmblast(
                proteins_fasta=filtered_fasta,
                ecod_db=args.ecod_db,
                out_csv=raw_csv,
                cpc=args.cpc,
                threads=args.threads,
                failure_log_dir=os.path.dirname(os.path.abspath(args.out)),
                timeout=args.timeout,
            )
        except RuntimeError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1

        entries = parse_plmblast_csv(raw_csv)

    logger.info(f"Parsed {len(entries)} pLM-BLAST hits")
    top1 = _reduce_to_top1(entries)
    logger.info(f"Reduced to {len(top1)} top-1 hits across queries")

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for entry in top1:
            writer.writerow(_to_output_row(entry))

    logger.info(f"Done: wrote {len(top1)} top-1 hits to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
