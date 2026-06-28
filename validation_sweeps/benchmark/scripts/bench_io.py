#!/usr/bin/env python3
"""Shared TSV I/O + small reporting helpers for the dataset-build scripts (30-35).

Single source for the read/write/normalize boilerplate that was otherwise copied per
script. write_tsv matches 01_build_gold_set.write_tsv (path, header, rows; blanks fill
missing keys) so the two are interchangeable.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in header})


def norm_doi(d: str) -> str:
    """Strip a leading 'doi:' and surrounding whitespace; matches doi_resolves' cache key."""
    return (d or "").strip().removeprefix("doi:").strip()


def by_type(rows: list[dict], field: str = "ss_type") -> dict:
    """Sorted per-value counts for a column, for run-summary logging."""
    return dict(sorted(Counter((r.get(field) or "?").strip() for r in rows).items()))


def modal(values: list[str]) -> tuple[str, int]:
    """Most common non-empty value + how many proposed it. Ties broken deterministically
    (lowest value wins) so the same inputs always yield the same pick across runs."""
    vals = [v.strip() for v in values if (v or "").strip()]
    if not vals:
        return "", 0
    counts = Counter(vals)
    top = max(counts.values())
    return min(v for v, c in counts.items() if c == top), top
