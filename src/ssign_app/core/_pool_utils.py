"""Pool and split helpers for multi-genome batched runs.

When ssign processes N>1 genomes in one invocation, prediction and
annotation tools see a single pooled FASTA / substrates TSV instead
of N per-genome files. Each protein's id is rewritten as
``<sample_id>__<original_id>`` at pool time and stripped back at
split time, so wrappers never need to know about genome provenance.

The double-underscore separator is FASTA-safe and distinct from
anything Bakta or Prodigal emit. Always go through ``make_prefixed_id``
/ ``split_prefixed_id`` rather than hardcoding ``"__"``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import IO

from ssign_app.scripts.ssign_lib.fasta_io import read_fasta_records

logger = logging.getLogger(__name__)

SEPARATOR = "__"

# Map file extension → csv delimiter. ssign wrappers emit CSV
# (interproscan, blastp, protparam — csv.DictWriter default ',') or TSV
# (eggnog, plm_blast, signalp, etc. — explicit '\t'); pool_tsvs and
# split_tsv_by_source preserve whatever the wrapper chose so downstream
# integrate_annotations sees the per-genome file in the same format the
# single-genome path produces.
_DELIMITER_BY_EXT = {".csv": ",", ".tsv": "\t"}


def _delimiter_for(path: Path | str) -> str:
    """Return the delimiter implied by a path's extension.

    Defaults to tab when the extension isn't ``.csv`` or ``.tsv`` (ssign's
    convention is TSV for ad-hoc intermediate files). Callers should pass
    paths whose suffix matches the wrapper's emitted format.
    """
    return _DELIMITER_BY_EXT.get(Path(path).suffix.lower(), "\t")


def make_prefixed_id(sample_id: str, locus_tag: str) -> str:
    return f"{sample_id}{SEPARATOR}{locus_tag}"


def split_prefixed_id(prefixed: str) -> tuple[str, str]:
    """Return ``(sample_id, locus_tag)``; splits on the FIRST ``SEPARATOR``.

    Locus tags themselves are allowed to contain ``SEPARATOR`` — only
    the leading ``<sample_id>__`` is stripped. sample_id is guaranteed
    not to contain ``SEPARATOR`` by ``validate_sample_id``.
    """
    sample_id, sep, locus_tag = prefixed.partition(SEPARATOR)
    if not sep:
        raise ValueError(f"Not a prefixed id (missing {SEPARATOR!r}): {prefixed!r}")
    return sample_id, locus_tag


def validate_sample_id(sample_id: str) -> None:
    if not sample_id:
        raise ValueError("sample_id must be non-empty")
    if SEPARATOR in sample_id:
        raise ValueError(
            f"sample_id {sample_id!r} contains the pool separator {SEPARATOR!r}; "
            f"choose a sample_id without {SEPARATOR!r}"
        )


def pool_fastas(sources: list[tuple[str, Path]], dest: Path) -> int:
    """Concatenate per-genome FASTAs into one pooled file with prefixed ids.

    Each input record's first-token id is rewritten to
    ``<sample_id>__<original_id>``; any header metadata after the first
    whitespace is preserved.

    Returns the total number of sequences written.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    n_total = 0
    with open(dest, "w") as out:
        for sample_id, src in sources:
            validate_sample_id(sample_id)
            n = 0
            for header, seq in read_fasta_records(src):
                parts = header.split(None, 1)
                original_id = parts[0]
                rest = f" {parts[1]}" if len(parts) > 1 else ""
                out.write(f">{make_prefixed_id(sample_id, original_id)}{rest}\n{seq}\n")
                n += 1
            n_total += n
            logger.debug("pool_fastas: %s contributed %d records from %s", sample_id, n, src)
    return n_total


def split_fasta_by_source(pooled_fasta: Path, out_dir: Path) -> dict[str, Path]:
    """Partition a pooled FASTA back into per-genome ``<sample_id>.faa`` files.

    Records whose first-token id is missing the ``<sample_id>__`` prefix
    are skipped with a warning.

    Returns ``{sample_id: output_path}`` for every genome present in the input.
    """
    pooled_fasta = Path(pooled_fasta)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    handles: dict[str, IO[str]] = {}
    paths: dict[str, Path] = {}

    try:
        for header, seq in read_fasta_records(pooled_fasta):
            parts = header.split(None, 1)
            first_token = parts[0]
            rest = f" {parts[1]}" if len(parts) > 1 else ""
            if SEPARATOR not in first_token:
                logger.warning("split_fasta_by_source: header missing prefix: %r", first_token)
                continue
            sample_id, original_id = split_prefixed_id(first_token)
            if sample_id not in handles:
                p = out_dir / f"{sample_id}.faa"
                handles[sample_id] = open(p, "w")
                paths[sample_id] = p
            handles[sample_id].write(f">{original_id}{rest}\n{seq}\n")
    finally:
        for h in handles.values():
            h.close()

    return paths


def pool_tsvs(
    sources: list[tuple[str, Path]],
    dest: Path,
    id_column: str = "locus_tag",
) -> int:
    """Concatenate per-genome TSVs, prefixing each row's id column.

    Output columns are the union of source columns in first-seen order;
    cells missing from a source are written empty. Rows with a blank
    ``id_column`` value are skipped (matches ``load_substrate_ids``
    tolerance for blank-tag header artefacts).

    Returns the total number of rows written. If no source has any rows
    and no source ever provided a header, writes an empty file and
    returns 0.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    all_fields: list[str] = []
    seen: set[str] = set()
    pooled_rows: list[dict[str, str]] = []

    for sample_id, src in sources:
        validate_sample_id(sample_id)
        with open(src, newline="") as f:
            reader = csv.DictReader(f, delimiter=_delimiter_for(src))
            source_fields = list(reader.fieldnames or [])
            if id_column not in source_fields:
                raise ValueError(f"pool_tsvs: id_column {id_column!r} not in {src} columns {source_fields}")
            for col in source_fields:
                if col not in seen:
                    all_fields.append(col)
                    seen.add(col)
            for row in reader:
                tag = (row.get(id_column) or "").strip()
                if not tag:
                    continue
                row[id_column] = make_prefixed_id(sample_id, tag)
                pooled_rows.append(row)

    if not all_fields:
        dest.write_text("")
        return 0

    with open(dest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_fields, delimiter=_delimiter_for(dest))
        writer.writeheader()
        writer.writerows(pooled_rows)
    return len(pooled_rows)


def split_tsv_by_source(
    pooled_tsv: Path,
    out_dir: Path,
    id_column: str = "locus_tag",
) -> dict[str, Path]:
    """Partition a pooled TSV back into per-genome ``<sample_id>.tsv`` files.

    Strips the ``<sample_id>__`` prefix from the id column on the way out.
    Rows whose id-column value is unprefixed are skipped with a warning.

    Returns ``{sample_id: output_path}`` for every genome present in the input.
    """
    pooled_tsv = Path(pooled_tsv)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    delimiter = _delimiter_for(pooled_tsv)
    # Mirror the input extension on output so downstream readers (e.g.
    # integrate_annotations) consume per-genome split files exactly as
    # they would the single-genome wrapper output of the same tool.
    out_ext = pooled_tsv.suffix or ".tsv"

    buffers: dict[str, list[dict[str, str]]] = {}
    fieldnames: list[str] = []

    with open(pooled_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            return {}
        if id_column not in fieldnames:
            raise ValueError(f"split_tsv_by_source: id_column {id_column!r} not in {pooled_tsv} columns {fieldnames}")
        for row in reader:
            raw = (row.get(id_column) or "").strip()
            if not raw or SEPARATOR not in raw:
                logger.warning("split_tsv_by_source: row missing %r prefix: %r", SEPARATOR, raw)
                continue
            sample_id, locus_tag = split_prefixed_id(raw)
            row[id_column] = locus_tag
            buffers.setdefault(sample_id, []).append(row)

    paths: dict[str, Path] = {}
    for sample_id, rows in buffers.items():
        out_path = out_dir / f"{sample_id}{out_ext}"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
        paths[sample_id] = out_path

    return paths
