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

from bench_runout import _strand_norm, _three_prime  # exact scripts/24 coord semantics (strand, 3'-stop)

BENCH = Path(__file__).resolve().parents[1]
RERUN = BENCH / "rerun"
RERUN_FULLASM = BENCH / "rerun_fullasm"

# Corpus contig -> rerun internal contig, for genomes whose rerun unit uses a different
# replicon accession for the SAME assembly/coordinate system. RefSeq (NC_/NZ_) and INSDC
# (AE/BX/CP...) are two NCBI IDs for the identical DNA at identical coordinates, so coords
# join 1:1. Each pair below was verified by coordinate+length concordance (scripts/71): the
# gold CDS span overlaps a rerun protein of the same length at the same position (lenratio
# 1.000). gold contig (absent from rerun) -> the rerun-staged accession for the same molecule.
CONTIG_ALIAS = {
    "NC_003112.2": "AE002098.2",  # N. meningitidis MC58 chromosome
    "AE004091.2": "NC_002516.2",  # P. aeruginosa PAO1 chromosome
    "NC_002929.2": "BX470248.1",  # B. pertussis Tohama I chromosome
}


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

    def emitted_overlap(self, contig: str, start: int, stop: int, strand=None) -> dict | None:
        """Canonical found-by-ssign rule (Teo, 2026-06-29): does ssign EMIT a secreted protein that
        overlaps the gold span by >=1 bp on the reconciled molecule?

          found=yes : at least one emitted protein overlaps  -> ssign recovered the effector.
          found=no, n_overlap>0, not emitted : Bakta called the ORF but ssign didn't emit it.
          found=no, n_overlap==0             : Bakta missed the gene entirely (still a genuine ssign fail).

        Unlike `join` (max-overlap, which can report a 0-overlap protein's emission since best_ovl
        starts at -1), this requires real overlap AND tests every overlapper for emission, so a small
        emitted overlapper isn't masked by a larger non-emitted one. Returns None if the contig is
        absent from the rerun (molecule un-reconcilable -> data fix, not a fail).

        When `strand` is given, also reports the scripts/24 bridge verdict (3'-stop within 3 bp, same
        strand) so the audit can show where that stricter rule diverges from any-overlap."""
        ic = contig if contig in self._contig2unit else CONTIG_ALIAS.get(contig)
        if ic not in self._contig2unit:
            return None
        unit = self._contig2unit[ic]
        rows, emitted = self._load(unit)
        st = _strand_norm(strand) if strand is not None else 0
        gtp = _three_prime(start, stop, st) if st else None
        hits, tp_hit = [], None
        for row in rows:
            if row["contig"] != ic or not row["start"] or not row["end"]:
                continue
            rs, re_ = int(row["start"]), int(row["end"])
            # _overlap is half-open arithmetic on 1-based inclusive coords, so it under-counts the shared
            # length by 1 bp and misses a hypothetical single-shared-base boundary case. Inert here: real
            # effector<->ORF overlaps span the whole gene body (min 293 bp across the 90 rows), nowhere near
            # the boundary, and Bakta only shifts the start by ~1 bp (stop identical). `> 0` is correct as used.
            ov = _overlap(start, stop, rs, re_)
            if ov > 0:
                hits.append((ov, row["locus_tag"], row["locus_tag"] in emitted))
            if (
                gtp is not None
                and _strand_norm(row.get("strand", "")) == st
                and abs(_three_prime(rs, re_, st) - gtp) <= 3
            ):
                tp_hit = row["locus_tag"]
        hits.sort(reverse=True)
        emitted_hits = [h for h in hits if h[2]]
        return {
            "contig": ic,
            "unit": unit.name,
            "n_overlap": len(hits),
            "best_locus": hits[0][1] if hits else "",
            "best_overlap_bp": hits[0][0] if hits else 0,
            "emitted_locus": emitted_hits[0][1] if emitted_hits else "",
            "found": "yes" if emitted_hits else "no",
            "reason": "emitted_overlap"
            if emitted_hits
            else ("overlap_not_emitted" if hits else "no_overlap_bakta_miss"),
            "three_prime_locus": tp_hit or "",
            "three_prime_found": "yes" if (tp_hit and tp_hit in emitted) else "no",
        }
