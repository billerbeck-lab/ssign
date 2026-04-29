#!/usr/bin/env python3
"""Run Bakta genome annotation and convert output to ssign standard format.

Bakta (https://github.com/oschwengers/bakta) provides rich annotation
compared to Prodigal — includes functional annotations, gene names, and
database cross-references. Requires a database download (~30GB full,
~2GB light).

Usage:
    bakta --db /path/to/db --output outdir --prefix sample contigs.fasta

Output files used:
    {prefix}.faa  — protein sequences
    {prefix}.tsv  — annotation table with locus tags, coordinates, products

Database download:
    bakta_db download --output /path/to/db --type light   # ~2GB
    bakta_db download --output /path/to/db --type full    # ~30GB
"""

import argparse
import csv
import logging
import os
import subprocess
import tempfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Bakta TSV column names (tab-separated, with header)
# Sequence Id | Type | Start | Stop | Strand | Locus Tag | Gene | Product | DbXrefs
_COL_SEQ_ID = "Sequence Id"
_COL_TYPE = "Type"
_COL_START = "Start"
_COL_STOP = "Stop"
_COL_STRAND = "Strand"
_COL_LOCUS_TAG = "Locus Tag"
_COL_GENE = "Gene"
_COL_PRODUCT = "Product"
_COL_DBXREFS = "DbXrefs"

# Only these feature types contain proteins. Bakta 1.12+ writes them
# lowercase (`cds`, `sorf`); compare case-insensitively to stay robust
# across Bakta versions.
_PROTEIN_TYPES = {"cds", "sorf"}

# Prefixes we surface from the DbXrefs column. EC/COG/GO/KEGG/RefSeq/Pfam
# feed annotation_consensus voting (Phase 3.2.c); UniRef/UniParc are
# carried for provenance (most-populated fields in light-DB output).
# Other Bakta prefixes (e.g. SO sequence-ontology terms) are dropped.
#
# Field names align with run_eggnog.py so annotation_consensus can treat
# Bakta and EggNOG rows uniformly.
_DBXREF_PREFIX_TO_OUTPUT_KEY = {
    "EC": "ec_numbers",
    "COG": "cog_ids",
    "GO": "go_terms",
    "KEGG": "kegg_ko",
    "RefSeq": "refseq_ids",
    "Pfam": "pfam_ids",
    "UniRef": "uniref_ids",
    "UniParc": "uniparc_ids",
}
_DBXREF_OUTPUT_KEYS = tuple(_DBXREF_PREFIX_TO_OUTPUT_KEY.values())


def run_bakta(contigs_fasta, db_path, sample_id, output_dir, threads=4):
    """Run Bakta on a genome FASTA file.

    Returns:
        tuple: (proteins_faa_path, tsv_path)
    """
    cmd = [
        "bakta",
        "--db",
        db_path,
        "--output",
        output_dir,
        "--prefix",
        sample_id,
        "--threads",
        str(threads),
        "--force",  # Overwrite if exists
        "--skip-plot",  # Skip plot (not needed in pipeline)
        contigs_fasta,
    ]

    logger.info(f"Running Bakta: bakta --db {db_path} --prefix {sample_id} ...")
    # FRAGILE: subprocess call requires the `bakta` binary on PATH
    # If this breaks: pip install ssign[bakta] (or pip install bakta>=1.5)
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=14400,
        )
    except FileNotFoundError as e:
        raise RuntimeError(
            f"Bakta binary not found: {e}\n"
            f"  Common causes:\n"
            f"    - Bakta is not installed or not on PATH\n"
            f"  How to fix:\n"
            f"    - pip install ssign[bakta]      # installs bakta>=1.5\n"
            f"    - Or:  pip install bakta\n"
            f"    - Or:  conda install -c bioconda bakta"
        ) from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"Bakta timed out after 4 hours: {e}\n"
            f"  How to fix:\n"
            f"    - Reduce genome size or split into contigs\n"
            f"    - Increase --threads if more CPUs are available"
        ) from e

    if result.returncode != 0:
        logger.error(f"Bakta failed:\n{result.stderr[:1000]}")
        raise RuntimeError(f"Bakta exit code {result.returncode}")

    proteins_faa = os.path.join(output_dir, f"{sample_id}.faa")
    tsv_path = os.path.join(output_dir, f"{sample_id}.tsv")

    if not os.path.exists(proteins_faa):
        raise FileNotFoundError(f"Bakta proteins not found: {proteins_faa}")
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Bakta TSV not found: {tsv_path}")

    return proteins_faa, tsv_path


def parse_dbxrefs(dbxrefs_field):
    """Split Bakta's DbXrefs column into a per-prefix dict of ID lists.

    Bakta emits entries as a comma-separated string like:
        "EC:3.6.5.n1, COG:COG0481, GO:GO:0001234, RefSeq:WP_123.1"

    GO IDs are themselves prefixed (`GO:GO:...`), so we split on the
    first colon only. Entries whose prefix isn't in
    `_DBXREF_PREFIX_TO_OUTPUT_KEY` are discarded.

    Returns a dict like:
        {"ec_numbers": ["3.6.5.n1"], "cog_ids": ["COG0481"], "go_terms": [...], ...}
    Missing prefixes map to empty lists.
    """
    result = {key: [] for key in _DBXREF_OUTPUT_KEYS}
    if not dbxrefs_field or not dbxrefs_field.strip():
        return result

    for raw_entry in dbxrefs_field.split(","):
        entry = raw_entry.strip()
        if not entry or ":" not in entry:
            continue
        prefix, _, value = entry.partition(":")
        prefix = prefix.strip()
        value = value.strip()
        if not value or prefix not in _DBXREF_PREFIX_TO_OUTPUT_KEY:
            continue
        result[_DBXREF_PREFIX_TO_OUTPUT_KEY[prefix]].append(value)

    return result


def parse_bakta_tsv(tsv_path):
    """Parse Bakta TSV annotation table.

    Returns list of dicts in gene_info.tsv shape:
        locus_tag, protein_id, gene, product, contig, start, end, strand,
        ec_numbers, cog_ids, go_terms, kegg_ko, refseq_ids, pfam_ids

    The six cross-reference fields are lists; they're joined with
    semicolons when we write to TSV (an empty list becomes the empty
    string). Annotation-consensus voting (Phase 3.2.c) uses these to
    map each protein to broad functional categories.
    """
    # Bakta TSV starts with "# " comment lines (software version, DB
    # version, DOI, URL) followed by a header row that *also* begins with
    # "#" but is the column header itself (e.g. "#Sequence Id\tType\t...").
    # Drop the comment lines and strip the leading "#" off the header so
    # DictReader sees the column names.
    with open(tsv_path) as f:
        lines = []
        for line in f:
            if line.startswith("# "):
                continue
            if line.startswith("#"):
                line = line[1:]
            lines.append(line)

    entries = []
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        feat_type = row.get(_COL_TYPE, "").strip().lower()
        if feat_type not in _PROTEIN_TYPES:
            continue

        locus_tag = row.get(_COL_LOCUS_TAG, "").strip()
        if not locus_tag:
            continue

        try:
            start = int(row.get(_COL_START, 0)) - 1  # Convert to 0-based
            end = int(row.get(_COL_STOP, 0))
        except (ValueError, TypeError):
            start, end = 0, 0

        xrefs = parse_dbxrefs(row.get(_COL_DBXREFS, ""))

        entries.append(
            {
                "locus_tag": locus_tag,
                "protein_id": "",
                "gene": row.get(_COL_GENE, "").strip(),
                "product": row.get(_COL_PRODUCT, "hypothetical protein").strip()
                or "hypothetical protein",
                "contig": row.get(_COL_SEQ_ID, "").strip(),
                "start": start,
                "end": end,
                "strand": row.get(_COL_STRAND, "+").strip(),
                **xrefs,
            }
        )

    return entries


def write_proteins_fasta(bakta_faa_path, entries, out_fasta):
    """Write cleaned proteins FASTA using Bakta locus tags.

    Bakta's .faa headers use locus tags directly:
    >{locus_tag} {product} [{contig}]
    """
    # Build set of locus tags we want
    wanted = {e["locus_tag"] for e in entries}

    # Parse Bakta FASTA, keeping only CDS/sORF proteins
    seqs = {}
    current_tag = None
    current_seq = []

    with open(bakta_faa_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_tag and current_tag in wanted:
                    seqs[current_tag] = "".join(current_seq).rstrip("*")
                # Bakta header: >{locus_tag} {description}
                current_tag = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_tag and current_tag in wanted:
            seqs[current_tag] = "".join(current_seq).rstrip("*")

    n_written = 0
    with open(out_fasta, "w") as f:
        for e in entries:
            tag = e["locus_tag"]
            seq = seqs.get(tag, "")
            if not seq:
                continue
            f.write(f">{tag}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
            n_written += 1

    logger.info(f"Wrote {n_written} protein sequences")
    return n_written


def main():
    parser = argparse.ArgumentParser(
        description="Run Bakta and convert output to ssign format"
    )
    parser.add_argument("--input", required=True, help="Input genome FASTA (contigs)")
    parser.add_argument("--db", required=True, help="Path to Bakta database")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--threads", type=int, default=4, help="CPU threads")
    parser.add_argument("--out-proteins", required=True, help="Output proteins FASTA")
    parser.add_argument("--out-gene-info", required=True, help="Output gene_info.tsv")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        # Run Bakta
        bakta_faa, bakta_tsv = run_bakta(
            args.input, args.db, args.sample, tmpdir, args.threads
        )

        # Parse annotation table
        entries = parse_bakta_tsv(bakta_tsv)
        logger.info(f"Parsed {len(entries)} CDS/sORF features from Bakta TSV")

        # Write cleaned proteins FASTA
        write_proteins_fasta(bakta_faa, entries, args.out_proteins)

    # Write gene_info TSV. List-valued cross-reference fields are joined
    # with semicolons for a single-line TSV cell — annotation_consensus.py
    # (Phase 3.2.c) splits them again on the consuming side.
    # List fields and their column ordering are derived from
    # _DBXREF_OUTPUT_KEYS so adding a new prefix above can't desync the writer.
    fieldnames = [
        "locus_tag",
        "protein_id",
        "gene",
        "product",
        "contig",
        "start",
        "end",
        "strand",
        *_DBXREF_OUTPUT_KEYS,
    ]
    _LIST_FIELDS = set(_DBXREF_OUTPUT_KEYS)
    with open(args.out_gene_info, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for e in entries:
            row = {
                k: (";".join(e[k]) if k in _LIST_FIELDS else e[k]) for k in fieldnames
            }
            writer.writerow(row)

    logger.info(f"Done: {len(entries)} proteins annotated for {args.sample}")


if __name__ == "__main__":
    main()
