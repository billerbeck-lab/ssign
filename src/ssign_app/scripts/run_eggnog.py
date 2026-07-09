#!/usr/bin/env python3
"""Run EggNOG-mapper functional annotation and convert output to ssign format.

EggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper) provides ortholog-
based functional annotation via DIAMOND searches against the EggNOG database.
Requires a database download (~50 GB for the full set; smaller taxonomy-
restricted subsets are also supported).

Usage:
    emapper.py -i proteins.faa --output-dir outdir -o sample --data_dir /path/to/eggnog_db

Output file consumed:
    {prefix}.emapper.annotations  — TSV with COG / KEGG / EC / GO / PFAM fields

Database download:
    download_eggnog_data.py --data_dir /path/to/eggnog_db
"""

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
from ssign_lib.resources import resolve_threads  # noqa: E402
from ssign_lib.subprocess_diag import dump_failure_log  # noqa: E402
from ssign_lib.substrates import (  # noqa: E402
    load_substrate_ids,
    write_substrates_only_fasta,
)

# emapper .annotations column names (v2.1+). The header line starts with
# '#query'; any line starting with '##' is a comment and must be skipped.
_COL_QUERY = "#query"
_COL_SEED_ORTHOLOG = "seed_ortholog"
_COL_EVALUE = "evalue"
_COL_DESCRIPTION = "Description"
_COL_PREFERRED_NAME = "Preferred_name"

# Rich-annotation columns surfaced for annotation-consensus voting
# (Phase 3.2.c). Fields come through as comma- or slash-separated lists
# in the raw TSV; we re-split and re-join on semicolon for our output
# TSV so downstream code sees a consistent separator across tools.
_COL_COG_CATEGORY = "COG_category"
_COL_EC = "EC"
_COL_KEGG_KO = "KEGG_ko"
_COL_GOS = "GOs"
_COL_PFAMS = "PFAMs"

# emapper uses "-" to mean "no annotation" in any rich field.
_EMAPPER_MISSING = "-"


# Column names follow the ssign convention used by every other annotation
# wrapper (run_blastp, run_interproscan, run_hhsuite, ...):
#   - join key is named `locus_tag` so integrate_annotations.py merges it
#     against the substrate master table without falling through to its
#     "no join column" branch and silently dropping the file.
#   - description-style fields are tool-namespaced (`eggnog_description`,
#     same shape as `bakta_product`, `blastp_hit_description`,
#     `interpro_descriptions`) so the consensus voting in
#     integrate_annotations.py `TOOL_HIT_COLUMNS` resolves them.
_OUTPUT_FIELDNAMES = [
    "locus_tag",
    "seed_ortholog",
    "evalue",
    "eggnog_description",
    "preferred_name",
    "cog_category",
    "ec_numbers",
    "kegg_ko",
    "go_terms",
    "pfam_ids",
]
_LIST_FIELDS = {"ec_numbers", "kegg_ko", "go_terms", "pfam_ids"}


def _write_empty_output(out_path):
    """Write a header-only TSV when there are no substrates to annotate."""
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()


# --dbmem loads ~44 GB resident; gate it on total RAM. Picked 50 GB so a
# 64 GB CPU node (Imperial CX3 high-mem) keeps the speedup but a 32 GB GPU
# node falls back to on-disk SQLite. <50 substrates pays only seconds extra
# on-disk vs in-RAM (the lookups are point queries, not table scans).
_EGGNOG_DBMEM_MIN_GB = 50

# --block_size and --index_chunks are emapper.py passthroughs to DIAMOND.
# Both are speed-vs-RAM tradeoffs: bigger block_size processes more queries
# in parallel, smaller index_chunks holds the DB index in fewer (bigger)
# passes. DIAMOND defaults of 2.0 / 4 are tuned for ~16 GB machines; on
# 64-128 GB HPC nodes we can run much faster.
_EGGNOG_BLOCKSIZE_DEFAULT = 2
_EGGNOG_BLOCKSIZE_CAP = 8  # DIAMOND wiki: gains plateau past ~8
# DIAMOND --block_size N ≈ 6×N GB resident; --dbmem reserves ~44 GB on top.
_EGGNOG_GB_PER_BLOCKSIZE_UNIT = 6
_EGGNOG_DBMEM_RESERVED_GB = 44
# Headroom we leave for OS / other processes when --dbmem is off.
_EGGNOG_NONDB_RESERVED_GB = 16
# At ≥ 60 GB the search index fits in 1 chunk (DIAMOND -c1) instead of 4.
_EGGNOG_INDEX_CHUNKS_MIN_RAM_GB = 60


def _autodetect_dbmem() -> bool:
    """True only when the job has headroom for emapper's in-RAM SQLite copy.

    Uses `ssign_lib.resources.parallel_share_ram_gb` — this subprocess's RAM
    *share*, not the whole node. It respects PBS `Resource_List.mem`, SLURM
    `SLURM_MEM_PER_NODE`, and cgroup limits, then divides by the parallel-group
    size when EggNOG runs inside the annotation group. Without the share divisor,
    --dbmem + a node-sized --block_size would overcommit and OOM-kill DIAMOND
    when the other annotation tools run concurrently (seen live at full tier).

    Inside the full-tier group the share is usually below the threshold, so
    --dbmem stays off. The shared-filesystem mmap stall that --dbmem also
    guarded against is then covered by the runner staging eggnog.db to local
    $TMPDIR (--local-cache-dir); a bare shared-FS run with $TMPDIR unset is the
    only gap (env misconfig, and run_emapper warns when staging is skipped).
    """
    from ssign_lib.resources import parallel_share_ram_gb

    return parallel_share_ram_gb() >= _EGGNOG_DBMEM_MIN_GB


def _autodetect_block_size(dbmem: bool) -> int:
    """Pick DIAMOND --block_size to use the RAM the job actually has.

    Reserves either ~44 GB (when --dbmem is loading the SQLite into RAM)
    or ~16 GB (system headroom otherwise), then sizes the DIAMOND working
    set against the remainder. Caps at 8 (DIAMOND gains plateau there).
    Sizes against this subprocess's RAM *share* (parallel-group aware), so a
    concurrent annotation group doesn't collectively overcommit the node.
    """
    from ssign_lib.resources import parallel_share_ram_gb

    reserved = _EGGNOG_DBMEM_RESERVED_GB if dbmem else _EGGNOG_NONDB_RESERVED_GB
    free = parallel_share_ram_gb() - reserved
    if free <= _EGGNOG_BLOCKSIZE_DEFAULT * _EGGNOG_GB_PER_BLOCKSIZE_UNIT:
        return _EGGNOG_BLOCKSIZE_DEFAULT
    sized = int(free // _EGGNOG_GB_PER_BLOCKSIZE_UNIT)
    return min(sized, _EGGNOG_BLOCKSIZE_CAP)


def _autodetect_index_chunks() -> int:
    """DIAMOND --index_chunks. Default 4 → 4 sequential DB-index passes.
    On a 60+ GB job the index fits in 1 chunk (~25% faster overall). Gauged on
    this subprocess's RAM share (parallel-group aware), matching --block_size."""
    from ssign_lib.resources import parallel_share_ram_gb

    return 1 if parallel_share_ram_gb() >= _EGGNOG_INDEX_CHUNKS_MIN_RAM_GB else 4


# Eggnog DB files emapper actually reads at run time. Listed exhaustively
# so we copy only what's needed (eggnog/ also holds optional HMM data and
# per-taxon downloads we don't touch). Filenames stable across v2.1.x.
# Files emapper requires at run time (raises if missing) vs nice-to-have.
_EGGNOG_REQUIRED_FILES = ("eggnog.db", "eggnog_proteins.dmnd")
_EGGNOG_OPTIONAL_FILES = ("eggnog.taxa.db", "eggnog.taxa.db.traverse.pkl")

# Headroom above the on-disk DB footprint: ~50 GB of files + temp space.
_EGGNOG_LOCAL_CACHE_MIN_FREE_GB = 60


def _stage_eggnog_db_to_local(src_dir: str, cache_dir: str) -> str:
    """Copy emapper runtime files from src_dir to cache_dir/eggnog/ and return the new path.

    Network filesystems (gpfs, nfs, lustre) make random-access mmap on the
    41 GB eggnog.db pathological — emapper hangs at near-zero CPU for hours
    instead of running. Copying to node-local SSD (~2 min on CX3 gpfs read)
    converts every subsequent lookup to a local page fault. Re-runs in the
    same job reuse the existing copy. Writes are atomic (tmp + os.replace)
    so a concurrent staging from another job can't corrupt the cache.
    """
    local_eggnog_dir = os.path.join(cache_dir, "eggnog")
    os.makedirs(local_eggnog_dir, exist_ok=True)
    for name in _EGGNOG_REQUIRED_FILES:
        if not os.path.exists(os.path.join(src_dir, name)):
            raise FileNotFoundError(
                f"Required eggnog DB file not found in {src_dir}: {name}. "
                f"Run download_eggnog_data.py --data_dir {src_dir} to populate it."
            )
    for name in _EGGNOG_REQUIRED_FILES + _EGGNOG_OPTIONAL_FILES:
        src = os.path.join(src_dir, name)
        if not os.path.exists(src):
            continue
        dst = os.path.join(local_eggnog_dir, name)
        if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src):
            logger.info(f"EggNOG DB already cached: {dst}")
            continue
        logger.info(f"Caching {name} ({os.path.getsize(src) / 2**30:.1f} GB) -> {local_eggnog_dir}")
        tmp = f"{dst}.tmp.{os.getpid()}"
        shutil.copy2(src, tmp)
        os.replace(tmp, dst)  # atomic; safe against concurrent stagers
    return local_eggnog_dir


def _build_emapper_cmd(
    proteins_fasta,
    db_path,
    sample_id,
    output_dir,
    threads=4,
    tax_scope="2",
    sensmode="sensitive",
    dbmem=None,
    block_size=None,
    index_chunks=None,
):
    """Build the emapper.py argv list. Exposed for unit testing.

    `dbmem=None` auto-detects via host RAM (`_autodetect_dbmem`); pass
    True/False to force.

    `block_size=None` and `index_chunks=None` auto-detect from RAM via
    `_autodetect_block_size` / `_autodetect_index_chunks`. Pass an int
    to force a specific value.
    """
    if dbmem is None:
        dbmem = _autodetect_dbmem()
    if block_size is None:
        block_size = _autodetect_block_size(dbmem)
    if index_chunks is None:
        index_chunks = _autodetect_index_chunks()
    cmd = [
        "emapper.py",
        "-i",
        proteins_fasta,
        # Explicit --itype proteins: EggNOG-mapper can also do ab-initio gene
        # prediction (wrapping Prodigal) when given DNA contigs. ssign gives
        # it Bakta-predicted proteins, so we force `proteins` to guarantee no
        # double CDS-calling even if EggNOG's default ever changes.
        "--itype",
        "proteins",
        # NB: emapper uses --output_dir / --data_dir (underscore, not dash).
        # Confirmed against eggnog-mapper 2.1.13.
        "--output_dir",
        output_dir,
        "-o",
        sample_id,
        "--data_dir",
        db_path,
        "--cpu",
        str(threads),
        "--tax_scope",
        tax_scope,
        "--sensmode",
        sensmode,
        "--override",
    ]
    if dbmem:
        # `--dbmem` loads the ~39 GB eggnog.db SQLite into RAM (~44 GB
        # resident). Without it, emapper mmaps the file — pathological on
        # NFS-backed shared scratch (Imperial CX3 RDS, similar clusters),
        # where startup hangs silently for hours instead of running.
        # eggnog-mapper issue #267 / #277.
        cmd.append("--dbmem")
    # Always emit so the runner log records the values actually used,
    # not the DIAMOND defaults that may or may not match autodetect.
    cmd.extend(["--block_size", str(block_size)])
    cmd.extend(["--index_chunks", str(index_chunks)])
    return cmd


def run_emapper(
    proteins_fasta,
    db_path,
    sample_id,
    output_dir,
    threads=None,
    tax_scope="2",
    sensmode="sensitive",
    dbmem=None,
    local_cache_dir=None,
    timeout=TOOL_TIMEOUT_S,
):
    """Run emapper.py on a protein FASTA file.

    `threads=None` (default) auto-detects from `resolve_threads()`.
    Argparse default already used this; the Python-API signature pinned
    4 and starved direct callers on multi-core machines.

    `tax_scope` restricts which orthologous groups are reported (post-
    filter on DIAMOND seed-ortholog hits) — ssign is gram-negative-
    bacteria-only, so NCBI taxid "2" (Bacteria) keeps eukaryotic OGs
    out of the output.

    `sensmode` is the DIAMOND sensitivity preset. "sensitive" is ~10× the
    DIAMOND default; "more-sensitive" is another ~2× slower and rarely
    rescues additional substrate-set orthologs.

    `dbmem=None` (default) auto-decides via `_autodetect_dbmem()` — on a
    32 GB host the load OOM-kills emapper mid-startup with no clean
    error, so the default fall back to on-disk SQLite there.

    `local_cache_dir`, when set and large enough, copies the eggnog runtime
    files to it before invoking emapper and points `--data_dir` there.
    Required on shared filesystems (gpfs/nfs/lustre): random-access mmap on
    the 41 GB SQLite hangs at near-zero CPU for hours otherwise. Pass the
    PBS/SLURM job-local TMPDIR.

    Returns:
        str: path to the `.emapper.annotations` file written by emapper.
    """
    threads = resolve_threads(threads)
    effective_db_path = db_path
    if local_cache_dir:
        free_gb = shutil.disk_usage(local_cache_dir).free / 2**30
        if free_gb >= _EGGNOG_LOCAL_CACHE_MIN_FREE_GB:
            effective_db_path = _stage_eggnog_db_to_local(db_path, local_cache_dir)
        else:
            logger.warning(
                f"Skipping local DB cache: {local_cache_dir} has only {free_gb:.1f} GB free "
                f"(need {_EGGNOG_LOCAL_CACHE_MIN_FREE_GB}). Emapper will mmap from {db_path} — "
                f"likely to hang on a shared filesystem."
            )

    cmd = _build_emapper_cmd(
        proteins_fasta,
        effective_db_path,
        sample_id,
        output_dir,
        threads=threads,
        tax_scope=tax_scope,
        sensmode=sensmode,
        dbmem=dbmem,
    )

    logger.info(f"Running EggNOG-mapper: emapper.py -i {proteins_fasta} -o {sample_id} --cpu {threads}")
    # FRAGILE: subprocess call requires `emapper.py` on PATH.
    # Eggnog-mapper is not in the ssign[extended] extras: its hard pin
    # biopython==1.76 conflicts with bakta>=1.78. Users install it
    # separately. See docs/how-to/install.md § EggNOG-mapper.
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"emapper.py binary not found: {e}\n"
            f"  How to fix:\n"
            f"    - conda install -c bioconda eggnog-mapper   (recommended)\n"
            f"    - Or:  pip install --no-deps eggnog-mapper\n"
            f"  See docs/how-to/install.md § EggNOG-mapper for details."
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"EggNOG-mapper timed out after 4 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input sequences\n"
            f"    - Use a taxonomy-restricted subset of the EggNOG DB"
        ) from e

    if result.returncode != 0:
        raise dump_failure_log("EggNOG-mapper", result, cmd, output_dir)

    annotations_path = os.path.join(output_dir, f"{sample_id}.emapper.annotations")
    if not os.path.exists(annotations_path):
        raise FileNotFoundError(f"EggNOG annotations not found: {annotations_path}")

    return annotations_path


def _split_rich_field(value):
    """Split an emapper multi-value cell into a clean list.

    emapper uses comma separators for most multi-value columns (GOs,
    EC, PFAMs, KEGG_ko), but uses "-" as a sentinel for "no annotation"
    and emits empty strings on unannotated rows. This helper normalises
    all three edge cases.
    """
    if value is None:
        return []
    stripped = value.strip()
    if not stripped or stripped == _EMAPPER_MISSING:
        return []
    return [part.strip() for part in stripped.split(",") if part.strip()]


def parse_eggnog_annotations(annotations_path):
    """Parse an EggNOG-mapper `.emapper.annotations` TSV.

    Real emapper output has '##' comment lines before the '#query' header
    and more '##' lines after the data. We strip them, then use DictReader
    on the remaining tab-separated rows.

    Returns list of dicts with columns:
        locus_tag, seed_ortholog, evalue, eggnog_description, preferred_name,
        cog_category, ec_numbers, kegg_ko, go_terms, pfam_ids

    Multi-value fields (ec_numbers, kegg_ko, go_terms, pfam_ids) are
    returned as lists; callers joining them for TSV output should
    semicolon-separate to stay consistent with Bakta's gene_info output.
    cog_category is a short single-letter-or-string code and stays as a
    string.
    """
    entries = []
    with open(annotations_path) as f:
        lines = [line for line in f if not line.startswith("##")]

    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        locus_tag = row.get(_COL_QUERY, "").strip()
        if not locus_tag or locus_tag == "-":
            continue

        cog_raw = row.get(_COL_COG_CATEGORY, "").strip()
        cog_category = "" if cog_raw in ("", _EMAPPER_MISSING) else cog_raw

        entries.append(
            {
                "locus_tag": locus_tag,
                "seed_ortholog": row.get(_COL_SEED_ORTHOLOG, "").strip(),
                "evalue": row.get(_COL_EVALUE, "").strip(),
                "eggnog_description": row.get(_COL_DESCRIPTION, "").strip() or "-",
                "preferred_name": row.get(_COL_PREFERRED_NAME, "").strip(),
                "cog_category": cog_category,
                "ec_numbers": _split_rich_field(row.get(_COL_EC)),
                "kegg_ko": _split_rich_field(row.get(_COL_KEGG_KO)),
                "go_terms": _split_rich_field(row.get(_COL_GOS)),
                "pfam_ids": _split_rich_field(row.get(_COL_PFAMS)),
            }
        )

    return entries


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run EggNOG-mapper on the filtered secreted-protein substrates "
            "and convert output to ssign format. Filters the full protein "
            "FASTA down to just the substrates before calling emapper — "
            "EggNOG is an annotation-tier tool, not a whole-genome step."
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
    parser.add_argument("--db", required=True, help="Path to EggNOG database directory")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "CPU threads. When omitted, resolves to the scheduler-aware "
            "effective CPU count, divided by the parallel-group size if "
            "the runner launched this wrapper inside one."
        ),
    )
    parser.add_argument(
        "--tax-scope",
        default="2",
        help=(
            "EggNOG --tax_scope: NCBI taxid restricting which orthologous "
            "groups are reported. Defaults to '2' (Bacteria) since ssign "
            "targets gram-negative bacteria. Pass 'auto' for unrestricted, "
            "or another taxid (e.g. '1224' for Proteobacteria)."
        ),
    )
    parser.add_argument(
        "--sensmode",
        default="sensitive",
        choices=[
            "default",
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ],
        help="DIAMOND sensitivity preset (default: sensitive).",
    )
    parser.add_argument(
        "--dbmem",
        dest="dbmem",
        action="store_const",
        const=True,
        default=None,
        help="Force --dbmem on (~44 GB RAM). Default auto-detects via host RAM.",
    )
    parser.add_argument(
        "--no-dbmem",
        dest="dbmem",
        action="store_const",
        const=False,
        help="Force --dbmem off (use on-disk SQLite).",
    )
    parser.add_argument(
        "--local-cache-dir",
        default=None,
        help=(
            "Copy the eggnog runtime DB (~50 GB) to this directory before "
            "calling emapper. Required on shared filesystems (gpfs/nfs): "
            "random-access mmap on the 41 GB SQLite stalls for tens of "
            "minutes otherwise. Typical: PBS/SLURM job-local TMPDIR. "
            "Skipped if directory has <60 GB free."
        ),
    )
    parser.add_argument("--out", required=True, help="Output annotations TSV (ssign format)")
    parser.add_argument(
        "--timeout",
        type=int,
        default=TOOL_TIMEOUT_S,
        help="Subprocess timeout in seconds (the runner passes a size-scaled value; default is the flat floor)",
    )
    args = parser.parse_args()

    substrate_ids = load_substrate_ids(args.substrates)
    if not substrate_ids:
        logger.info("No substrates to annotate — writing empty output")
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
        logger.info(f"EggNOG will annotate {n} substrate proteins (of {len(substrate_ids)} listed)")

        annotations_path = run_emapper(
            filtered_fasta,
            args.db,
            args.sample,
            tmpdir,
            threads=args.threads,
            tax_scope=args.tax_scope,
            sensmode=args.sensmode,
            dbmem=args.dbmem,
            local_cache_dir=args.local_cache_dir,
            timeout=args.timeout,
        )
        entries = parse_eggnog_annotations(annotations_path)

    logger.info(f"Parsed {len(entries)} annotations from EggNOG-mapper")

    # Multi-value fields are lists in `entries`; join with semicolons for
    # a single TSV cell (same convention as run_bakta.py's gene_info).
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for e in entries:
            row = {k: (";".join(e[k]) if k in _LIST_FIELDS else e[k]) for k in _OUTPUT_FIELDNAMES}
            writer.writerow(row)

    logger.info(f"Done: wrote {len(entries)} annotations to {args.out}")
    return 0


if __name__ == "__main__":
    main()
