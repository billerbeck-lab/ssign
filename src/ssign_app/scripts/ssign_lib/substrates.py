"""Read locus_tags from a substrates TSV; write a substrates-only FASTA.

Used by all downstream-annotation wrappers (run_blastp, run_eggnog,
run_hhsuite, run_interproscan, run_plm_blast, compute_protparam).
The defensive ``(row.get("locus_tag") or "").strip()`` read tolerates
blank-tag header artefacts that have appeared in upstream
cross_validate output — older inlined copies of this function used
``row["locus_tag"]`` directly and would KeyError on the same input.
"""

import csv

from .fasta_io import read_fasta


def load_substrate_ids(substrates_path: str) -> set[str]:
    """Return the set of locus_tags listed in a filtered substrates TSV.

    Tolerates rows with a missing or whitespace-only ``locus_tag`` column
    by silently skipping them — defensive because upstream cross_validate
    output has, on rare occasions, included blank-tag header artefacts.
    """
    ids: set[str] = set()
    with open(substrates_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tag = (row.get("locus_tag") or "").strip()
            if tag:
                ids.add(tag)
    return ids


def write_substrates_only_fasta(
    proteins_path: str,
    substrate_ids: set[str],
    out_path: str,
) -> int:
    """Write a FASTA containing only proteins whose ID is in ``substrate_ids``.

    Returns the number of sequences written. If zero, the caller should
    short-circuit — the downstream tool has nothing to annotate.
    """
    all_seqs = read_fasta(proteins_path)
    n = 0
    with open(out_path, "w") as f:
        for locus_tag in substrate_ids:
            seq = all_seqs.get(locus_tag)
            if not seq:
                continue
            f.write(f">{locus_tag}\n{seq}\n")
            n += 1
    return n
