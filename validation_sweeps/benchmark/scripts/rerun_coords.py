#!/usr/bin/env python3
"""Coordinate-join corpus effectors to the Bakta-reannotated tier-2 rerun.

Bakta renames every locus_tag on re-annotation, so a corpus effector's RefSeq locus
cannot be matched to the rerun by name. We bridge by genome COORDINATES instead: the
rerun preserves the input contig accessions (Bakta keeps the sequence IDs from the
input GenBank), so we join an effector's (contig, start, stop) to the max-overlap
rerun protein on the same contig. Emission is then read from the per-genome
`<unit>_results.csv` (the filtered, emitted-secreted set).

Used by 4.2 (T5SS annotation grading) and 4.4 (recall reconciliation). The 4 RTX/T1SS
genomes that the main rerun ran on fragment inputs were re-run into `rerun_fullasm/`;
this index prefers `rerun_fullasm/<unit>` over `rerun/<unit>` automatically once those
results land, so callers fold them in without code changes.
"""

from __future__ import annotations

import csv
import glob
import io
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
RERUN = BENCH / "rerun"
RERUN_FULLASM = BENCH / "rerun_fullasm"

# Corpus contig -> rerun internal contig, for genomes whose rerun unit uses a different
# replicon accession for the SAME assembly/coordinate system. MC58 is staged under its
# INSDC accession (AE002098.2) but the corpus places its effectors on the RefSeq
# accession (NC_003112.2); the two are the identical sequence, so coords join 1:1.
CONTIG_ALIAS = {"NC_003112.2": "AE002098.2"}


def _unit_dirs() -> list[Path]:
    """Every rerun unit dir, preferring rerun_fullasm/<unit> when both exist."""
    dirs: dict[str, Path] = {}
    for root in (RERUN, RERUN_FULLASM):  # fullasm second so it overrides
        if not root.is_dir():
            continue
        for d in root.iterdir():
            if d.is_dir() and glob.glob(str(d / "results" / "*_results_raw.csv")):
                dirs[d.name] = d
    return list(dirs.values())


def _raw_path(unit: Path) -> str:
    return glob.glob(str(unit / "results" / "*_results_raw.csv"))[0]


def _emitted_path(unit: Path) -> str | None:
    hits = [p for p in glob.glob(str(unit / "results" / "*_results.csv")) if not p.endswith("_raw.csv")]
    return hits[0] if hits else None


def _read_emitted_loci(path: str) -> set[str]:
    """Locus tags in the leading `# Secreted Proteins` section of a results.csv.

    The file holds further `#`-prefixed sections (e.g. `# Secretion Systems`) after a blank
    line, so we read only up to that first blank line; a naive `#`-line filter would map those
    later rows onto this section's header (cf. bench_runout._read_secreted_chunk)."""
    with open(path, newline="") as fh:
        lines = fh.read().split("\n")
    start = next((i for i, ln in enumerate(lines) if ln.strip().lower().startswith("# secreted")), None)
    if start is None:
        return set()
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == ""), len(lines))
    block = lines[start + 1 : end]
    return {r["locus_tag"] for r in csv.DictReader(io.StringIO("\n".join(block)))} if block else set()


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


class RerunIndex:
    """Lazily-loaded contig->unit map + per-unit row/emission cache for coordinate joins."""

    def __init__(self) -> None:
        self._contig2unit: dict[str, Path] = {}
        self._rows: dict[str, list[dict]] = {}
        self._emitted: dict[str, set[str]] = {}
        # Cheap pre-scan: read only the contig column (csv.reader, no per-row dict) to map every
        # contig to its unit. Full rows are read once, on demand, in _load() for joined units only.
        for unit in _unit_dirs():
            with open(_raw_path(unit), newline="") as fh:
                rd = csv.reader(fh)
                header = next(rd, None)
                if not header or "contig" not in header:
                    continue
                ci = header.index("contig")
                for row in rd:
                    if len(row) > ci:
                        self._contig2unit.setdefault(row[ci], unit)

    def _load(self, unit: Path) -> tuple[list[dict], set[str]]:
        key = str(unit)
        if key not in self._rows:
            with open(_raw_path(unit), newline="") as fh:
                self._rows[key] = list(csv.DictReader(fh))
            ep = _emitted_path(unit)
            self._emitted[key] = _read_emitted_loci(ep) if ep else set()
        return self._rows[key], self._emitted[key]

    def join(self, contig: str, start: int, stop: int) -> dict | None:
        """Best-overlap rerun protein for a corpus locus at (contig, start, stop).

        Returns {unit, contig, locus_tag, overlap_frac, emitted, raw} or None if the
        contig is absent from the rerun. overlap_frac is the intersection over the
        longer of the two spans (near 1.0 for a clean hit)."""
        ic = contig if contig in self._contig2unit else CONTIG_ALIAS.get(contig)
        if ic not in self._contig2unit:
            return None
        unit = self._contig2unit[ic]
        rows, emitted = self._load(unit)
        best, best_ovl = None, -1
        for row in rows:
            if row["contig"] != ic or not row["start"] or not row["end"]:
                continue
            o = _overlap(start, stop, int(row["start"]), int(row["end"]))
            if o > best_ovl:
                best_ovl, best = o, row
        if best is None:
            return None
        span = max(stop - start, int(best["end"]) - int(best["start"])) or 1
        return {
            "unit": unit.name,
            "contig": ic,
            "locus_tag": best["locus_tag"],
            "overlap_frac": best_ovl / span,
            "emitted": best["locus_tag"] in emitted,
            "raw": best,
        }
