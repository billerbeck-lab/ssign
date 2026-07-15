"""Build a passenger-substituted FASTA for T5aSS annotation routing.

The annotation tools that default to passenger-only for T5aSS substrates
(EggNOG / BLASTp / pLM-BLAST / HHsearch / ProtParam) read their input
sequences from a single FASTA keyed by locus_tag. To route the passenger
domain through those tools without changing each wrapper's contract, we
build one FASTA where T5aSS substrates carry their passenger subsequence
under the same locus_tag, and every other entry carries the full protein
unchanged. The wrappers see a normal locus_tag-keyed FASTA.

Routing rules (per T5aSS substrate):
- t5_quality_flag is non-empty (barrel_only, no_signalp, no_sec_signal,
  omp_porin_no_at, unclassified) → fall back to full protein.
- passenger_length < MIN_PASSENGER_FOR_ANNOTATION → fall back to full.
- Otherwise → write the passenger subsequence
  seq[sp_end : barrel_start - LINKER_LENGTH]  (0-indexed slice).

Non-T5aSS substrates and proteins absent from the classifications TSV
are written with their full sequence.

Side-car TSV records, per locus_tag, whether the FASTA entry is the
passenger or the full protein; integrate_annotations.py uses this to
populate the t5_annotation_source column in the master CSV.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Literal

from .constants import LINKER_LENGTH, MIN_PASSENGER_FOR_ANNOTATION
from .fasta_io import read_fasta, write_fasta
from .parsing import parse_int_or_none

logger = logging.getLogger(__name__)

AnnotationSource = Literal["passenger", "full"]
PASSENGER: AnnotationSource = "passenger"
FULL: AnnotationSource = "full"


def load_t5_classifications(classifications_tsv: str | Path) -> dict[str, dict]:
    """Read t5ss_handler's per-substrate domain classifications TSV.

    Returns a dict keyed by locus_tag with the geometry fields needed to
    extract the passenger subsequence:
      ss_type, t5_quality_flag, passenger_length, sp_end, barrel_start.
    """
    classifications_tsv = Path(classifications_tsv)
    if not classifications_tsv.exists():
        return {}

    out: dict[str, dict] = {}
    with open(classifications_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            locus = row.get("locus_tag")
            if not locus:
                continue
            out[locus] = {
                "ss_type": row.get("ss_type", ""),
                "t5_quality_flag": row.get("t5_quality_flag", ""),
                "passenger_length": parse_int_or_none(row.get("passenger_length", "")) or 0,
                "sp_end": parse_int_or_none(row.get("sp_end", "")),
                "barrel_start": parse_int_or_none(row.get("barrel_start", "")),
            }
    return out


def _passenger_slice(seq: str, sp_end: int, barrel_start: int) -> str:
    """Extract the passenger subsequence from a full protein.

    sp_end and barrel_start are the 1-indexed positions emitted by
    t5ss_handler. Passenger spans 1-indexed [sp_end+1, barrel_start-LINKER_LENGTH]
    inclusive, which is seq[sp_end : barrel_start - LINKER_LENGTH] in
    Python's 0-indexed half-open form.
    """
    start = sp_end
    end = barrel_start - LINKER_LENGTH
    if start < 0 or end <= start or end > len(seq):
        return ""
    return seq[start:end]


def _routing_decision(
    locus: str,
    classification: dict | None,
    min_passenger_length: int,
) -> AnnotationSource:
    """Decide passenger vs full for one locus_tag. Pure function, easy to test."""
    if classification is None or classification.get("ss_type") != "T5aSS":
        return FULL
    if classification.get("t5_quality_flag"):
        return FULL
    sp_end = classification.get("sp_end")
    barrel_start = classification.get("barrel_start")
    if sp_end is None or barrel_start is None:
        return FULL
    passenger_length = classification.get("passenger_length", 0)
    if passenger_length < min_passenger_length:
        return FULL
    return PASSENGER


def build_passenger_substituted_fasta(
    proteins_fasta: str | Path,
    classifications_tsv: str | Path,
    out_fasta: str | Path,
    out_source_tsv: str | Path,
    min_passenger_length: int = MIN_PASSENGER_FOR_ANNOTATION,
) -> dict[str, AnnotationSource]:
    """Write a FASTA where T5aSS substrates carry their passenger subsequence.

    Every protein in ``proteins_fasta`` is written to ``out_fasta`` with the
    same locus_tag. T5aSS substrates with a clean t5_quality_flag and
    passenger_length >= ``min_passenger_length`` carry their passenger
    sequence; all other entries (including T5bSS/T5cSS, every non-T5SS
    substrate, and T5aSS substrates that failed the quality gate) carry the
    full protein sequence.

    Writes a side-car TSV (``out_source_tsv``) with columns
    ``locus_tag``, ``t5_annotation_source`` recording the choice for each
    T5aSS substrate. Non-T5aSS entries are omitted from the side-car —
    downstream readers treat absent rows as full-protein.

    Returns the {locus_tag: source} mapping (T5aSS rows only).
    """
    sequences = read_fasta(proteins_fasta)
    classifications = load_t5_classifications(classifications_tsv)

    routing: dict[str, AnnotationSource] = {}
    out_sequences: dict[str, str] = {}
    fallback_quality = 0
    fallback_short = 0
    passenger_count = 0

    for locus, seq in sequences.items():
        classification = classifications.get(locus)
        decision = _routing_decision(locus, classification, min_passenger_length)

        if decision == PASSENGER:
            # _routing_decision only returns PASSENGER when classification is present.
            assert classification is not None
            passenger_seq = _passenger_slice(
                seq,
                classification["sp_end"],
                classification["barrel_start"],
            )
            if not passenger_seq:
                logger.warning(
                    "Passenger slice empty for %s (sp_end=%s, barrel_start=%s,"
                    " seq_len=%d); falling back to full protein",
                    locus,
                    classification["sp_end"],
                    classification["barrel_start"],
                    len(seq),
                )
                out_sequences[locus] = seq
                routing[locus] = FULL
            else:
                out_sequences[locus] = passenger_seq
                routing[locus] = PASSENGER
                passenger_count += 1
        else:
            out_sequences[locus] = seq
            if classification is not None and classification.get("ss_type") == "T5aSS":
                routing[locus] = FULL
                if classification.get("t5_quality_flag"):
                    fallback_quality += 1
                else:
                    fallback_short += 1

    write_fasta(out_sequences, out_fasta)

    out_source_tsv = Path(out_source_tsv)
    out_source_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_source_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["locus_tag", "t5_annotation_source"])
        for locus in sorted(routing):
            writer.writerow([locus, routing[locus]])

    logger.info(
        "T5aSS annotation routing: %d passenger, %d full (quality-flag), %d full (passenger<%daa)",
        passenger_count,
        fallback_quality,
        fallback_short,
        min_passenger_length,
    )
    return routing
