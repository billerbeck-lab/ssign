"""Single FASTA read/write implementation for the entire pipeline.

This is the ONLY FASTA I/O module -- all pipeline scripts import from here.
No duplicate ``read_fasta()`` or ``write_fasta()`` functions elsewhere.
"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def read_fasta(fasta_path: str | Path) -> dict[str, str]:
    """Read a FASTA file into an ``{id: sequence}`` dictionary.

    Header parsing takes the first whitespace-delimited token after ``>``.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dictionary mapping sequence IDs to their sequences.

    Raises:
        FileNotFoundError: If *fasta_path* does not exist.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif current_id is not None:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    logger.debug("Read %d sequences from %s", len(sequences), fasta_path)
    return sequences


def read_fasta_records(fasta_path: str | Path):
    """Stream ``(full_header, sequence)`` tuples from a FASTA file.

    Use when the caller needs the raw header line (e.g. tool output where
    metadata is encoded after the first whitespace token). Plain
    ``read_fasta`` keeps only the first token; this preserves everything
    after the leading ``>``.

    Lazy: yields one record at a time so memory usage is bounded by the
    longest single sequence.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    current_header: str | None = None
    current_seq: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            elif current_header is not None:
                current_seq.append(line.strip())

    if current_header is not None:
        yield current_header, "".join(current_seq)


def count_sequences(source: str | Path | bytes) -> int:
    """Count records in a FASTA without loading sequences into memory.

    Accepts either a path (read line-by-line) or already-loaded bytes
    (split-and-count). String paths are treated as paths; pass bytes
    explicitly when the content is in memory and reading from disk
    would duplicate it.
    """
    if isinstance(source, bytes):
        # bytes.count is C-implemented and fast; handle the file-starts-
        # with-> case explicitly so the count matches the line-by-line path.
        n = source.count(b"\n>")
        if source.startswith(b">"):
            n += 1
        return n

    with open(source) as f:
        return sum(1 for line in f if line.startswith(">"))


def write_fasta(sequences: dict[str, str], output_path: str | Path, line_width: int = 80) -> int:
    """Write a ``{id: sequence}`` dictionary to a FASTA file.

    Creates parent directories if they do not exist.  Skips entries with
    empty sequences (logs a warning for each).

    Args:
        sequences: Mapping of sequence IDs to amino-acid / nucleotide strings.
        output_path: Destination file path.
        line_width: Maximum characters per sequence line (default 80).

    Returns:
        Number of sequences actually written (excluding skipped empties).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with open(output_path, "w") as f:
        for seq_id, sequence in sequences.items():
            if not sequence:
                logger.warning("Skipping %s: empty sequence", seq_id)
                continue
            f.write(f">{seq_id}\n")
            for i in range(0, len(sequence), line_width):
                f.write(sequence[i : i + line_width] + "\n")
            count += 1

    logger.info("Wrote %d sequences to %s", count, output_path)
    return count
