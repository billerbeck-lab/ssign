#!/usr/bin/env python3
"""Run pLM-BLAST remote-homology search and convert output to ssign format.

pLM-BLAST (Kaminski et al. 2023, https://github.com/labstructbioinf/pLM-BLAST)
detects remote homologues by comparing ProtT5 residue embeddings. We use
it with the precomputed ECOD70 database to annotate secreted proteins with
ECOD structural domains, complementing the sequence-homology searches
BLAST and HH-suite provide.

Usage:
    run_plm_blast.py --input proteins.faa --ecod-db /path/to/ecod70 \\
        --out plm_blast.tsv [--threads 4] [--cpc 70]

The `--ecod-db` argument points at a pre-built pLM-BLAST database directory
fetched from
    http://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/plmblast_dbs
Output is a TSV mapping each query protein to its top pLM-BLAST hits with
score and alignment coordinates.
"""

from __future__ import annotations

import argparse
import csv
import logging
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
from ssign_lib.constants import TOOL_TIMEOUT_S  # noqa: E402
from ssign_lib.substrates import (  # noqa: E402
    load_substrate_ids,
    write_substrates_only_fasta,
)


_OUTPUT_FIELDNAMES = [
    "protein_id",
    "target_id",
    "score",
    "qstart",
    "qend",
    "tstart",
    "tend",
]


def _write_empty_output(out_path):
    """Write a header-only TSV when there are no substrates to search."""
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()


# pLM-BLAST output CSV columns (v1.x). Exact header names verified against
# the upstream `scripts/plmblast.py` result-writer on first integration
# run — see the TODO below if upstream renames columns.
_COL_QUERY = "qid"
_COL_TARGET = "sid"
_COL_SCORE = "score"
_COL_QSTART = "qstart"
_COL_QEND = "qend"
_COL_TSTART = "tstart"
_COL_TEND = "tend"


def _resolve_plmblast_script() -> str:
    """Locate the `plmblast.py` script pLM-BLAST ships.

    Resolution order:
      1. `SSIGN_PLMBLAST_SCRIPT` env var (explicit override).
      2. `plmblast.py` / `plmblast` on PATH (entry point installed by pip).
      3. Fall back to assuming `plmblast.py` — the subprocess call will
         raise with install hints if it's not found.
    """
    override = os.environ.get("SSIGN_PLMBLAST_SCRIPT")
    if override:
        if not os.path.exists(override):
            raise RuntimeError(f"SSIGN_PLMBLAST_SCRIPT points at a missing file: {override}")
        return override

    for candidate in ("plmblast.py", "plmblast"):
        found = shutil.which(candidate)
        if found:
            return found

    return "plmblast.py"


def _resolve_embeddings_script(plmblast_script: str) -> str:
    """Locate pLM-BLAST's `embeddings.py` (one level above `scripts/plmblast.py`)."""
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
) -> None:
    """Run pLM-BLAST's embeddings.py to embed a FASTA, producing a
    single pooled `<base>.pt` file plus a sibling `<base>.csv` index.

    `embedder` defaults to `pt` (ProtT5-XL UniRef50), matching the
    embedder used to build the public ECOD70 database. ESM-2 etc. will
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
    cmd = [
        "python",
        embed_script,
        "start",
        embed_input,
        pt_path,
        "-embedder",
        embedder,
    ]
    logger.info(f"Embedding query: embeddings.py start {proteins_fasta} {pt_path} -embedder {embedder}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TOOL_TIMEOUT_S)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"pLM-BLAST embeddings.py not runnable: {e}\n"
            f"  How to fix:\n"
            f"    - pip install git+https://github.com/labstructbioinf/pLM-BLAST.git\n"
            f"    - Or clone the repo and set SSIGN_PLMBLAST_SCRIPT to its scripts/plmblast.py"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"pLM-BLAST embedding timed out after 4 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce query size\n"
            f"    - Use a GPU (--cuda) — embedding is ~100x faster"
        ) from e
    if result.returncode != 0:
        logger.error(f"pLM-BLAST embedding failed:\n{result.stderr[:1000]}")
        raise RuntimeError(f"pLM-BLAST embedding exit code {result.returncode}: {result.stderr[:300]}")

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
    threads: int = 4,
) -> str:
    """Run pLM-BLAST against ECOD70 and return the CSV output path.

    Two-step pipeline: ProtT5-embed the query FASTA into a single
    pooled .pt file, then search the embedding DB. Embedding dominates
    wall time on CPU (~5-10 sec per 500-aa protein).
    """
    script = _resolve_plmblast_script()
    embed_script = _resolve_embeddings_script(script)

    with tempfile.TemporaryDirectory() as tmpdir:
        query_emb_base = os.path.join(tmpdir, "query")
        _embed_query_fasta(embed_script, proteins_fasta, query_emb_base)

        # plmblast.py uses single-dash long flags (`-cpc`, `-workers`)
        # per its argparse setup in alntools/parser.py.
        cmd = [
            "python",
            script,
            ecod_db,
            query_emb_base,
            out_csv,
            "-cpc",
            str(cpc),
            "-workers",
            str(threads),
        ]

        logger.info(f"Running pLM-BLAST: {' '.join(cmd[:3])} <db> <query_emb> {out_csv} -cpc {cpc}")
        # FRAGILE: subprocess call requires pLM-BLAST's scripts/plmblast.py
        # on PATH, or set SSIGN_PLMBLAST_SCRIPT to its absolute path.
        # If this breaks: pip install git+https://github.com/labstructbioinf/pLM-BLAST.git
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=TOOL_TIMEOUT_S)
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
                f"    - Use a smaller ECOD subset (e.g. ECOD50 instead of ECOD70)"
            ) from e

        if result.returncode != 0:
            logger.error(f"pLM-BLAST failed:\n{result.stderr[:1000]}")
            raise RuntimeError(f"pLM-BLAST exit code {result.returncode}")

        if not os.path.exists(out_csv):
            raise FileNotFoundError(f"pLM-BLAST output not found: {out_csv}")

    return out_csv


def parse_plmblast_csv(csv_path: str):
    """Parse pLM-BLAST CSV output into a list of hit dicts.

    Returns a list keyed per-hit (a protein can have multiple hits) with
    fields: protein_id, target_id, score, qstart, qend, tstart, tend.

    TODO (Phase 3.2.a): surface the ECOD domain classification (T-group /
    F-group IDs embedded in target IDs) once annotation_consensus.py is
    extended to vote on structural domains. Currently the raw target ID
    passes through as-is.
    """
    entries = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            query_id = row.get(_COL_QUERY, "").strip()
            if not query_id:
                continue
            entries.append(
                {
                    "protein_id": query_id,
                    "target_id": row.get(_COL_TARGET, "").strip(),
                    "score": row.get(_COL_SCORE, "").strip(),
                    "qstart": row.get(_COL_QSTART, "").strip(),
                    "qend": row.get(_COL_QEND, "").strip(),
                    "tstart": row.get(_COL_TSTART, "").strip(),
                    "tend": row.get(_COL_TEND, "").strip(),
                }
            )
    return entries


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run pLM-BLAST against ECOD70 on the filtered substrates and "
            "convert output to ssign format. pLM-BLAST is an annotation-"
            "tier tool — it runs only on the ~50 filtered substrates, not "
            "the full ~5000-protein genome."
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
        help="Path to ECOD70 database directory (from ftp.tuebingen.mpg.de)",
    )
    parser.add_argument("--out", required=True, help="Output hits TSV (ssign format)")
    parser.add_argument(
        "--cpc",
        type=int,
        default=70,
        help="pLM-BLAST -cpc cluster percent cutoff (default: 70)",
    )
    parser.add_argument("--threads", type=int, default=4, help="CPU threads")
    args = parser.parse_args()

    if not os.path.exists(args.substrates):
        print(f"ERROR: Substrates TSV not found: {args.substrates}", file=sys.stderr)
        return 2
    if not os.path.exists(args.proteins):
        print(f"ERROR: Protein FASTA not found: {args.proteins}", file=sys.stderr)
        return 2
    if not os.path.isdir(args.ecod_db):
        print(
            f"ERROR: ECOD70 database directory not found: {args.ecod_db}\n"
            f"  Fetch from http://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/plmblast_dbs",
            file=sys.stderr,
        )
        return 2

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
            )
        except RuntimeError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1

        entries = parse_plmblast_csv(raw_csv)

    logger.info(f"Parsed {len(entries)} pLM-BLAST hits")

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for e in entries:
            writer.writerow(e)

    logger.info(f"Done: wrote {len(entries)} hits to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
