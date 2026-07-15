#!/usr/bin/env python3
"""Audit on-disk size of each ssign database, grouped by install tier.

Walks the entries in ``ssign_lib.dependency_manifest`` against a live
install root, sums actual bytes per database, and produces a Markdown
table plus per-tier rollups. Use to refresh the size claims in the
README and ``docs/how-to/install.md`` after a database upgrade.

Pure-Python (no ``du`` dependency) so it works on any OS the user
might want to run it on.

Usage::

    scripts/audit_database_sizes.py --db-root ~/.ssign/databases

Optional flags::

    --weights-root PATH    where DeepSecE checkpoint sits (default ~/.ssign)
    --json                 emit machine-readable JSON instead of Markdown
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Iterable

# Make the manifest importable when running this script directly without an
# editable install. ``src/ssign_app/scripts/ssign_lib`` is the canonical
# location; tests stub this path via PYTHONPATH.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from ssign_app.scripts.ssign_lib.dependency_manifest import (  # noqa: E402
    _TIER_ORDER,
    DATABASE_PATHS,
    MODEL_WEIGHTS,
)

_GB = 1024**3


def path_size_bytes(path: str) -> int:
    """Total on-disk bytes for ``path`` (file or directory tree).

    Uses ``lstat`` so a single symlink to a 50 GB file outside the tree
    counts only as the link entry itself, not 50 GB. Hardlinked files
    are deduped by ``(st_dev, st_ino)`` because tar extractions of large
    bundles (HH-suite UniRef, BFD) occasionally include hardlinked
    duplicates and counting both would overstate on-disk usage.
    """
    if not os.path.exists(path):
        return 0
    if os.path.isfile(path):
        return os.lstat(path).st_size

    total = 0
    seen: set[tuple[int, int]] = set()
    stack = [path]
    while stack:
        cur = stack.pop()
        try:
            it = os.scandir(cur)
        except OSError:
            continue
        with it:
            for entry in it:
                try:
                    st = entry.stat(follow_symlinks=False)
                except OSError:
                    continue
                if entry.is_dir(follow_symlinks=False):
                    stack.append(entry.path)
                    continue
                key = (st.st_dev, st.st_ino)
                if key in seen:
                    continue
                seen.add(key)
                total += st.st_size
    return total


# Old name kept so tests/external callers don't have to change.
dir_size_bytes = path_size_bytes


def _format_gb(n_bytes: int | None) -> str:
    if n_bytes is None:
        return "—"
    return f"{n_bytes / _GB:.1f} GB"


def audit(db_root: str, weights_root: str) -> list[dict]:
    """Return one row per database/weights entry from the manifest.

    NOTE: ``DatabasePath.resolve_path`` consults the entry's env var first.
    If ``$SSIGN_BAKTA_DB`` etc. are set in the shell when this script runs,
    sizes will be reported against the env-pointed dir, not ``--db-root``.
    Unset those before auditing if you want strict per-tree numbers.
    """
    rows: list[dict] = []
    for db in DATABASE_PATHS:
        path = db.resolve_path(db_root)
        rows.append(
            {
                "kind": "database",
                "name": db.name,
                "tier": db.tier,
                "path": path,
                "bytes": path_size_bytes(path) if path else None,
            }
        )
    for w in MODEL_WEIGHTS:
        base = db_root if w.under_db_root else weights_root
        candidate = os.path.join(base, w.default_subpath)
        # MODEL_WEIGHTS entries can be either a directory (e.g. a weights
        # tree) or a single file (DeepSecE checkpoint .pt). Accept both via
        # os.path.exists; path_size_bytes handles each.
        path = candidate if os.path.exists(candidate) else None
        rows.append(
            {
                "kind": "weights",
                "name": w.name,
                "tier": w.tier,
                "path": path,
                "bytes": path_size_bytes(path) if path else None,
            }
        )
    return rows


def cumulative_tier_totals(rows: Iterable[dict]) -> dict[str, int]:
    """Sum bytes for each tier *cumulatively* (extended ⊇ base, full ⊇ extended).

    Returns {"base": ..., "extended": ..., "full": ...}. Missing entries
    contribute 0. The expected user-facing claim is "to install at tier
    X you need this many GB" — that's a cumulative sum.
    """
    totals = {t: 0 for t in _TIER_ORDER}
    for row in rows:
        size = row["bytes"] or 0
        idx = _TIER_ORDER.index(row["tier"])
        for higher_tier in _TIER_ORDER[idx:]:
            totals[higher_tier] += size
    return totals


def render_markdown(rows: list[dict], totals: dict[str, int]) -> str:
    lines = ["# ssign database size audit", ""]
    lines.append("Generated by `scripts/audit_database_sizes.py`. Re-run after a database upgrade.")
    lines.append("")
    lines.append("## Per-database")
    lines.append("")
    lines.append("| Name | Tier | Path | Size |")
    lines.append("| --- | --- | --- | ---: |")
    for row in rows:
        path_str = row["path"] or "—"
        lines.append(f"| {row['name']} | {row['tier']} | `{path_str}` | {_format_gb(row['bytes'])} |")
    lines.append("")
    lines.append("## Cumulative per tier")
    lines.append("")
    lines.append("(each row sums everything at that tier or below)")
    lines.append("")
    lines.append("| Tier | Total |")
    lines.append("| --- | ---: |")
    for tier in _TIER_ORDER:
        lines.append(f"| {tier} | {_format_gb(totals[tier])} |")
    lines.append("")
    return "\n".join(lines)


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--db-root",
        default=os.path.expanduser("~/.ssign/databases"),
        help="Root passed to fetch_databases.sh (default: ~/.ssign/databases)",
    )
    parser.add_argument(
        "--weights-root",
        default=os.path.expanduser("~/.ssign"),
        help="Root where auto-downloaded weights land (default: ~/.ssign)",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of Markdown")
    args = parser.parse_args(argv)

    rows = audit(args.db_root, args.weights_root)
    totals = cumulative_tier_totals(rows)

    if args.json:
        out = {"rows": rows, "cumulative_bytes": totals}
        print(json.dumps(out, indent=2))
    else:
        print(render_markdown(rows, totals))
    return 0


if __name__ == "__main__":
    sys.exit(main())
