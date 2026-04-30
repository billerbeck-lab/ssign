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
import subprocess
import sys
import tempfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
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


_OUTPUT_FIELDNAMES = [
    "protein_id",
    "seed_ortholog",
    "evalue",
    "description",
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


def run_emapper(
    proteins_fasta,
    db_path,
    sample_id,
    output_dir,
    threads=4,
    tax_scope="2",
    sensmode="more-sensitive",
):
    """Run emapper.py on a protein FASTA file.

    `tax_scope` restricts which orthologous groups are reported (post-
    filter on DIAMOND seed-ortholog hits) — ssign is gram-negative-
    bacteria-only, so NCBI taxid "2" (Bacteria) keeps eukaryotic OGs
    out of the output. Numeric taxids are version-stable; named scopes
    like "bacteria" only resolve if they match an eggnog tax-tree node
    label and have shifted across emapper releases.

    `sensmode` selects the DIAMOND sensitivity preset; "more-sensitive"
    is ~3× slower than DIAMOND default but recovers more remote
    orthologs (acceptable for the ~30 substrate proteins ssign typically
    annotates).

    Returns:
        str: path to the `.emapper.annotations` file written by emapper.
    """
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

    logger.info(f"Running EggNOG-mapper: emapper.py -i {proteins_fasta} -o {sample_id} --cpu {threads}")
    # FRAGILE: subprocess call requires `emapper.py` on PATH
    # If this breaks: pip install ssign[extended] (or pip install eggnog-mapper)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"emapper.py binary not found: {e}\n"
            f"  Common causes:\n"
            f"    - EggNOG-mapper is not installed or not on PATH\n"
            f"  How to fix:\n"
            f"    - pip install ssign[extended]   # installs eggnog-mapper\n"
            f"    - Or:  pip install eggnog-mapper\n"
            f"    - Or:  conda install -c bioconda eggnog-mapper"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"EggNOG-mapper timed out after 4 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce the number of input sequences\n"
            f"    - Use a taxonomy-restricted subset of the EggNOG DB"
        ) from e

    if result.returncode != 0:
        logger.error(f"EggNOG-mapper failed:\n{result.stderr[:1000]}")
        raise RuntimeError(f"EggNOG-mapper exit code {result.returncode}")

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
        protein_id, seed_ortholog, evalue, description, preferred_name,
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
        protein_id = row.get(_COL_QUERY, "").strip()
        if not protein_id or protein_id == "-":
            continue

        cog_raw = row.get(_COL_COG_CATEGORY, "").strip()
        cog_category = "" if cog_raw in ("", _EMAPPER_MISSING) else cog_raw

        entries.append(
            {
                "protein_id": protein_id,
                "seed_ortholog": row.get(_COL_SEED_ORTHOLOG, "").strip(),
                "evalue": row.get(_COL_EVALUE, "").strip(),
                "description": row.get(_COL_DESCRIPTION, "").strip() or "-",
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
    parser.add_argument("--threads", type=int, default=4, help="CPU threads")
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
        default="more-sensitive",
        choices=[
            "default",
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ],
        help="DIAMOND sensitivity preset (default: more-sensitive).",
    )
    parser.add_argument("--out", required=True, help="Output annotations TSV (ssign format)")
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
