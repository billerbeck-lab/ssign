"""Shared helpers for ssign unit tests — bare functions, no fixtures.

Importable from any test file in tests/unit/. conftest.py also imports from
here so the auto-discovered fixtures stay in sync with what tests build by
hand. The leading underscore signals "internal to the test suite" so this
file isn't picked up as a test module by pytest collection.
"""

import csv
import sys

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def write_tsv(path, fieldnames, rows, delimiter="\t"):
    """Write rows as a delimited table with the given fieldnames; returns the path.

    Default delimiter is tab (the ssign convention). Pass delimiter="," when a
    test needs to mimic a tool whose output is comma-separated (e.g. DeepSecE's
    raw CSV).
    """
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def read_tsv_rows(path, delimiter="\t"):
    """Read a TSV/CSV into a list of dict rows."""
    with open(path) as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def run_script_main(monkeypatch, main_fn, argv):
    """Invoke a script's main() in-process with argv monkey-patched.

    Standardises the pattern used to drive each pipeline script's CLI
    entry point from a test: assemble argv, swap sys.argv, call main(),
    rely on monkeypatch to restore on teardown.
    """
    monkeypatch.setattr(sys, "argv", argv)
    main_fn()


# NOTE: FASTA writing lives in ssign_app.scripts.ssign_lib.fasta_io.write_fasta;
# tests should import that production helper directly rather than re-rolling.


# ---------------------------------------------------------------------------
# Field constants — match the canonical schemas the production scripts emit
# ---------------------------------------------------------------------------
# These mirror the inline `csv.DictWriter(fieldnames=...)` literals in
# extract_gene_order.py, prodigal_to_gene_info.py, validate_macsyfinder_systems.py,
# and cross_validate_predictions.py. Centralising the canonical column lists
# in a "ssign_lib/io_schemas.py" module is a future refactor.

GENE_INFO_FIELDS = [
    "locus_tag",
    "protein_id",
    "gene",
    "product",
    "contig",
    "start",
    "end",
    "strand",
]

GENE_ORDER_FIELDS = [
    "contig",
    "gene_index",
    "locus_tag",
    "start",
    "end",
    "strand",
]

# validate_macsyfinder_systems.py:168-171
SS_COMPONENT_FIELDS = [
    "sample_id",
    "sys_id",
    "ss_type",
    "locus_tag",
    "gene_name",
    "gene_status",
    "wholeness",
    "excluded",
]

PREDICTIONS_FIELDS = [
    "locus_tag",
    "dlp_extracellular_prob",
    "predicted_localization",
    "dlp_max_localization",
    "dlp_max_probability",
    "dse_ss_type",
    "dse_max_prob",
    "signalp_prediction",
    "signalp_probability",
    "signalp_cs_position",
    "product",
]


# ---------------------------------------------------------------------------
# Row builders
# ---------------------------------------------------------------------------


def two_contig_genes():
    """15-gene 2-contig layout: contig_A genes 0-9, contig_B genes 0-4."""
    rows = []
    for i in range(10):
        rows.append(
            {
                "locus_tag": f"GENE_{i:04d}",
                "protein_id": f"PROT_A_{i}",
                "gene": "",
                "product": "hypothetical protein",
                "contig": "contig_A",
                "start": str(i * 1000),
                "end": str(i * 1000 + 999),
                "strand": "+",
            }
        )
    for i in range(5):
        rows.append(
            {
                "locus_tag": f"GENEB_{i:04d}",
                "protein_id": f"PROT_B_{i}",
                "gene": "",
                "product": "hypothetical protein",
                "contig": "contig_B",
                "start": str(i * 1000),
                "end": str(i * 1000 + 999),
                "strand": "+",
            }
        )
    return rows


def gene_info_to_gene_order(gene_info_rows):
    """Project gene_info rows onto the gene_order schema (per-contig 0-based index).

    Assumes locus_tag ends in a numeric suffix that doubles as the contig-local
    gene_index. The conftest fixture genes follow this convention; tests that
    deviate must build gene_order rows themselves.
    """
    rows = []
    for r in gene_info_rows:
        gene_index = int(r["locus_tag"].split("_")[-1])
        rows.append(
            {
                "contig": r["contig"],
                "gene_index": str(gene_index),
                "locus_tag": r["locus_tag"],
                "start": r["start"],
                "end": r["end"],
                "strand": r["strand"],
            }
        )
    return rows


def make_ss_component_row(
    locus_tag,
    ss_type,
    gene_name="generic",
    *,
    excluded="False",
    sys_id=None,
):
    """Build a single row matching SS_COMPONENT_FIELDS with sensible defaults."""
    return {
        "sample_id": "test_sample",
        "sys_id": sys_id or f"contig_A_{ss_type}_1",
        "ss_type": ss_type,
        "locus_tag": locus_tag,
        "gene_name": gene_name,
        "gene_status": "mandatory",
        "wholeness": "1.0",
        "excluded": excluded,
    }


def make_prediction_row(
    locus_tag,
    dlp_ext=0.0,
    dse_type="Non-secreted",
    dse_prob=0.0,
    signalp_pred="OTHER",
    signalp_prob=0.0,
    product="hypothetical protein",
):
    """Builder for a single predictions row — defaults are non-secreted."""
    return {
        "locus_tag": locus_tag,
        "dlp_extracellular_prob": str(dlp_ext),
        "predicted_localization": "Extracellular" if dlp_ext >= 0.8 else "Cytoplasmic",
        "dlp_max_localization": "Extracellular" if dlp_ext >= 0.8 else "Cytoplasmic",
        "dlp_max_probability": str(dlp_ext),
        "dse_ss_type": dse_type,
        "dse_max_prob": str(dse_prob),
        "signalp_prediction": signalp_pred,
        "signalp_probability": str(signalp_prob),
        "signalp_cs_position": "",
        "product": product,
    }


# ---------------------------------------------------------------------------
# Tool-output row builders
# ---------------------------------------------------------------------------
# These mirror the upstream tool output formats so wrapper-parser tests can
# build realistic input rows without re-rolling the field layout per file.

# BLAST outfmt 6 column order — must stay in lockstep with run_blastp.BLAST_OUTFMT.
BLAST_OUTFMT_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
    "qlen",
    "slen",
]


def make_blast_outfmt_row(
    qseqid="GENE_0001",
    sseqid="WP_000000001.1",
    pident=95.0,
    aln_len=200,
    mismatch=5,
    gapopen=0,
    qstart=1,
    qend=200,
    sstart=1,
    send=200,
    evalue=1e-50,
    bitscore=400.0,
    stitle="hemolysin",
    qlen=200,
    slen=200,
):
    """Single tab-separated BLAST outfmt-6 row matching BLAST_OUTFMT_COLS."""
    return "\t".join(
        str(x)
        for x in [
            qseqid,
            sseqid,
            pident,
            aln_len,
            mismatch,
            gapopen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
            stitle,
            qlen,
            slen,
        ]
    )


# DeepSecE raw output columns (pre-rename) — must stay in lockstep with the
# keys of run_deepsece._COLUMN_MAP.
DSE_RAW_FIELDS = [
    "protein_id",
    "deepsece_prediction",
    "deepsece_ss_type",
    "max_prob",
    "nonsec_prob",
    "T1_prob",
    "T2_prob",
    "T3_prob",
    "T4_prob",
    "T6_prob",
]


def make_dse_raw_row(
    protein_id,
    ss_type="Non-secreted",
    max_prob=0.5,
    nonsec_prob=0.5,
    T1_prob=0.05,
    T2_prob=0.05,
    T3_prob=0.05,
    T4_prob=0.05,
    T6_prob=0.05,
):
    """Single DeepSecE raw output row (pre-rename), keyed by DSE_RAW_FIELDS."""
    return {
        "protein_id": protein_id,
        "deepsece_prediction": ss_type,
        "deepsece_ss_type": ss_type,
        "max_prob": str(max_prob),
        "nonsec_prob": str(nonsec_prob),
        "T1_prob": str(T1_prob),
        "T2_prob": str(T2_prob),
        "T3_prob": str(T3_prob),
        "T4_prob": str(T4_prob),
        "T6_prob": str(T6_prob),
    }
