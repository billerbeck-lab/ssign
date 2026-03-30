"""TSV-based per-protein status tracking.

The :class:`Manifest` class maintains an in-memory dictionary of protein
processing records that can be saved to / loaded from a tab-separated file.
It supports incremental updates, resume-on-restart, and summary reporting.
"""

import csv
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class Manifest:
    """Track per-protein processing status as a TSV file.

    Supports incremental updates, resume-on-restart, and summary reporting.

    Args:
        path: File path for the TSV manifest.
        columns: List of column names (``protein_id`` is always first and
            is added automatically if not included).
    """

    def __init__(self, path: str, columns: list[str]) -> None:
        self.path = Path(path)
        self.columns = ["protein_id"] + [c for c in columns if c != "protein_id"]
        self._entries: dict[str, dict] = {}
        if self.path.exists():
            self._load()

    def _load(self) -> None:
        """Load existing manifest from disk."""
        with open(self.path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                pid = row.get("protein_id", "")
                if pid:
                    self._entries[pid] = dict(row)
        logger.debug("Loaded manifest with %d entries from %s", len(self._entries), self.path)

    def set(self, protein_id: str, **kwargs) -> None:
        """Set or update a protein's record.

        Args:
            protein_id: Unique protein identifier.
            **kwargs: Column values to store (e.g. ``status='success'``).
        """
        self._entries[protein_id] = {"protein_id": protein_id, **kwargs}

    def get(self, protein_id: str) -> Optional[dict]:
        """Return the record for *protein_id*, or ``None`` if absent."""
        return self._entries.get(protein_id)

    def get_successful(self) -> list[str]:
        """Return protein IDs with ``status == 'success'``."""
        return [
            pid
            for pid, entry in self._entries.items()
            if entry.get("status") == "success"
        ]

    def get_pending(self, all_ids: list[str]) -> list[str]:
        """Return IDs from *all_ids* that are not yet successfully processed.

        A protein is considered pending if it has no entry or its status is
        neither ``'success'`` nor ``'skipped'``.
        """
        return [
            pid
            for pid in all_ids
            if pid not in self._entries
            or self._entries[pid].get("status") not in ("success", "skipped")
        ]

    def save(self) -> None:
        """Write the manifest to disk as a TSV file.

        Creates parent directories if they do not exist.
        """
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=self.columns,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            for entry in self._entries.values():
                writer.writerow(entry)
        logger.debug("Saved manifest with %d entries to %s", len(self._entries), self.path)

    def summary(self) -> dict[str, int]:
        """Return a dictionary of ``{status: count}``."""
        statuses: dict[str, int] = {}
        for entry in self._entries.values():
            s = entry.get("status", "unknown")
            statuses[s] = statuses.get(s, 0) + 1
        return statuses

    def __len__(self) -> int:
        return len(self._entries)
