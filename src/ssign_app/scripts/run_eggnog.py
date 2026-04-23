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
import tempfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# emapper .annotations column names (v2.1+). The header line starts with
# '#query'; any line starting with '##' is a comment and must be skipped.
_COL_QUERY = "#query"
_COL_SEED_ORTHOLOG = "seed_ortholog"
_COL_EVALUE = "evalue"
_COL_DESCRIPTION = "Description"
_COL_PREFERRED_NAME = "Preferred_name"


def run_emapper(proteins_fasta, db_path, sample_id, output_dir, threads=4):
    """Run emapper.py on a protein FASTA file.

    Returns:
        str: path to the `.emapper.annotations` file written by emapper.
    """
    cmd = [
        "emapper.py",
        "-i",
        proteins_fasta,
        "--output-dir",
        output_dir,
        "-o",
        sample_id,
        "--data_dir",
        db_path,
        "--cpu",
        str(threads),
        "--override",
    ]

    logger.info(
        f"Running EggNOG-mapper: emapper.py -i {proteins_fasta} "
        f"-o {sample_id} --cpu {threads}"
    )
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


def parse_eggnog_annotations(annotations_path):
    """Parse an EggNOG-mapper `.emapper.annotations` TSV.

    Real emapper output has '##' comment lines before the '#query' header
    and more '##' lines after the data. We strip them, then use DictReader
    on the remaining tab-separated rows.

    Returns list of dicts with columns:
        protein_id, seed_ortholog, evalue, description, preferred_name

    TODO (Phase 3.2.a): surface COG_category / EC / KEGG_ko / GOs / PFAMs
    columns once annotation_consensus.py is extended to vote with them.
    Currently these richer fields are discarded.
    """
    entries = []
    with open(annotations_path) as f:
        lines = [line for line in f if not line.startswith("##")]

    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        protein_id = row.get(_COL_QUERY, "").strip()
        if not protein_id or protein_id == "-":
            continue

        entries.append(
            {
                "protein_id": protein_id,
                "seed_ortholog": row.get(_COL_SEED_ORTHOLOG, "").strip(),
                "evalue": row.get(_COL_EVALUE, "").strip(),
                "description": row.get(_COL_DESCRIPTION, "").strip() or "-",
                "preferred_name": row.get(_COL_PREFERRED_NAME, "").strip(),
            }
        )

    return entries


def main():
    parser = argparse.ArgumentParser(
        description="Run EggNOG-mapper and convert output to ssign format"
    )
    parser.add_argument("--input", required=True, help="Input protein FASTA")
    parser.add_argument("--db", required=True, help="Path to EggNOG database directory")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--threads", type=int, default=4, help="CPU threads")
    parser.add_argument(
        "--out", required=True, help="Output annotations TSV (ssign format)"
    )
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        annotations_path = run_emapper(
            args.input, args.db, args.sample, tmpdir, args.threads
        )
        entries = parse_eggnog_annotations(annotations_path)

    logger.info(f"Parsed {len(entries)} annotations from EggNOG-mapper")

    fieldnames = [
        "protein_id",
        "seed_ortholog",
        "evalue",
        "description",
        "preferred_name",
    ]
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for e in entries:
            writer.writerow(e)

    logger.info(f"Done: wrote {len(entries)} annotations to {args.out}")


if __name__ == "__main__":
    main()
