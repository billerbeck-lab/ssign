"""Tests for run_ortholog_grouping.py.

Pure-Python surfaces here:

1. `cluster_union_find` — single-linkage clustering via path-compressed
   union-find. Hits between known IDs unify; transitivity must hold;
   IDs not mentioned in any hit must end up as their own singleton.
2. `compute_group_stats` — per-group `OG_NNN` IDs by descending size,
   mean within-group %identity (bidirectional hit lookup), singleton
   default identity = 100.
3. `_find_blast_binary` — must distinguish "not installed"
   (BlastpUnavailableError) from "broken install" (RuntimeError) so
   main() can soft-skip vs hard-fail correctly.
4. `main()` end-to-end via subprocess stubs — empty / singleton /
   blast-unavailable / blast-failure paths.

FASTA reading is provided by ssign_lib.fasta_io.read_fasta — covered in
its own test module; not re-tested here.
"""

import csv
import os
import subprocess
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import run_ortholog_grouping  # noqa: E402
from run_ortholog_grouping import (  # noqa: E402
    BlastpUnavailableError,
    _find_blast_binary,
    cluster_union_find,
    compute_group_stats,
)

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta  # noqa: E402

# ---------------------------------------------------------------------------
# cluster_union_find
# ---------------------------------------------------------------------------


class TestClusterUnionFind:
    def test_no_hits_yields_all_singletons(self):
        groups = cluster_union_find(hits=[], all_protein_ids={"A", "B", "C"})
        # Every protein its own group
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 1, 1]

    def test_pairwise_hit_unifies_two(self):
        groups = cluster_union_find(
            hits=[("A", "B", 80.0, 80.0)],
            all_protein_ids={"A", "B", "C"},
        )
        # 2 groups: {A,B} and {C}
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 2]

    def test_transitive_closure(self):
        # A-B and B-C → all three in one group
        groups = cluster_union_find(
            hits=[("A", "B", 80.0, 80.0), ("B", "C", 80.0, 80.0)],
            all_protein_ids={"A", "B", "C"},
        )
        assert len(groups) == 1
        assert next(iter(groups.values())) == {"A", "B", "C"}

    def test_two_disconnected_components(self):
        groups = cluster_union_find(
            hits=[
                ("A", "B", 80.0, 80.0),
                ("C", "D", 80.0, 80.0),
            ],
            all_protein_ids={"A", "B", "C", "D"},
        )
        member_sets = sorted([frozenset(s) for s in groups.values()], key=lambda s: sorted(s)[0])
        assert member_sets == [frozenset({"A", "B"}), frozenset({"C", "D"})]

    def test_hits_referencing_unknown_ids_skipped(self):
        # "X" is not in all_protein_ids → the hit should silently no-op
        groups = cluster_union_find(
            hits=[("A", "X", 80.0, 80.0)],
            all_protein_ids={"A", "B"},
        )
        # Both A and B remain as singletons
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 1]

    def test_self_hit_idempotent(self):
        # A union(A, A) shouldn't change cluster structure
        groups = cluster_union_find(
            hits=[("A", "A", 100.0, 100.0)],
            all_protein_ids={"A", "B"},
        )
        sizes = sorted(len(s) for s in groups.values())
        assert sizes == [1, 1]

    def test_chain_of_five(self):
        # A-B-C-D-E chain → one group of 5
        hits = [
            ("A", "B", 80.0, 80.0),
            ("B", "C", 80.0, 80.0),
            ("C", "D", 80.0, 80.0),
            ("D", "E", 80.0, 80.0),
        ]
        groups = cluster_union_find(hits=hits, all_protein_ids={"A", "B", "C", "D", "E"})
        assert len(groups) == 1
        assert next(iter(groups.values())) == {"A", "B", "C", "D", "E"}


# ---------------------------------------------------------------------------
# compute_group_stats
# ---------------------------------------------------------------------------


class TestComputeGroupStats:
    def test_groups_ordered_by_size_desc(self):
        # 3-member group first (OG_001), then two singletons
        groups = {
            "rep1": {"A", "B", "C"},
            "rep2": {"D"},
            "rep3": {"E"},
        }
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 90.0, 90.0), ("B", "C", 95.0, 95.0)],
            all_protein_ids={"A", "B", "C", "D", "E"},
        )
        # First entry is the 3-member group
        assert stats[0]["ortholog_group"] == "OG_001"
        assert stats[0]["n_members"] == 3
        # Subsequent entries are singletons (in some order, both n_members=1)
        assert all(s["n_members"] == 1 for s in stats[1:])

    def test_singleton_default_identity_100(self):
        groups = {"rep1": {"A"}}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={"A"},
        )
        assert stats[0]["mean_pident"] == 100.0
        assert stats[0]["n_members"] == 1

    def test_mean_identity_within_group(self):
        # Group {A, B, C} with within-group hits at 80 and 90% → mean 85.0
        groups = {"rep1": {"A", "B", "C"}}
        stats = compute_group_stats(
            groups,
            hits=[
                ("A", "B", 80.0, 80.0),
                ("B", "C", 90.0, 90.0),
                # Cross-group hit (irrelevant — A and B both in same group)
            ],
            all_protein_ids={"A", "B", "C"},
        )
        assert stats[0]["mean_pident"] == 85.0

    def test_bidirectional_hit_lookup(self):
        # Hit recorded only as (A, B); compute_group_stats must also find (B, A)
        groups = {"rep1": {"A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 75.0, 80.0)],
            all_protein_ids={"A", "B"},
        )
        assert stats[0]["mean_pident"] == 75.0

    def test_members_field_sorted_and_semicolon_joined(self):
        groups = {"rep1": {"C", "A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={"A", "B", "C"},
        )
        assert stats[0]["members"] == "A;B;C"

    def test_og_id_zero_padded(self):
        # 3-digit zero-padding contract pinned: OG_001, OG_002, ...
        groups = {f"rep{i}": {f"P{i}"} for i in range(1, 6)}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={f"P{i}" for i in range(1, 6)},
        )
        og_ids = [s["ortholog_group"] for s in stats]
        assert og_ids == ["OG_001", "OG_002", "OG_003", "OG_004", "OG_005"]

    def test_no_hits_in_multi_member_group_falls_back_to_100(self):
        # group {A, B} but no hits between them — defensively defaults to 100
        # rather than crashing on division-by-zero
        groups = {"rep1": {"A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[],
            all_protein_ids={"A", "B"},
        )
        assert stats[0]["mean_pident"] == 100.0


# ---------------------------------------------------------------------------
# _find_blast_binary
# ---------------------------------------------------------------------------


class TestFindBlastBinary:
    """Probe must distinguish 'not installed' (BlastpUnavailableError) from
    'broken install' (RuntimeError) so main() can soft-skip vs hard-fail."""

    def test_not_on_path_raises_unavailable(self, monkeypatch):
        def _raise_fnf(*a, **k):
            raise FileNotFoundError("blastp")

        monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", _raise_fnf)
        with pytest.raises(BlastpUnavailableError):
            _find_blast_binary("blastp")

    def test_corrupt_install_raises_runtime_error(self, monkeypatch):
        def _raise_called_process(*a, **k):
            raise subprocess.CalledProcessError(127, ["blastp", "-version"])

        monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", _raise_called_process)
        with pytest.raises(RuntimeError, match="corrupted or incompatible"):
            _find_blast_binary("blastp")

    def test_hung_install_raises_runtime_error(self, monkeypatch):
        def _raise_timeout(*a, **k):
            raise subprocess.TimeoutExpired(["blastp", "-version"], 10)

        monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", _raise_timeout)
        with pytest.raises(RuntimeError, match="hung"):
            _find_blast_binary("blastp")

    def test_returns_name_when_present(self, monkeypatch):
        monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", lambda *a, **k: None)
        assert _find_blast_binary("blastp") == "blastp"


# ---------------------------------------------------------------------------
# main() end-to-end via subprocess stub
# ---------------------------------------------------------------------------


def _read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


class TestMainEndToEnd:
    """Black-box main() coverage: empty input, singleton, BLAST unavailable,
    BLAST hard-fail. The happy path with hits is covered by integration."""

    def test_empty_fasta_writes_header_only(self, monkeypatch, tmp_dir):
        fasta = os.path.join(tmp_dir, "empty.fasta")
        open(fasta, "w").close()
        out = os.path.join(tmp_dir, "out.csv")

        monkeypatch.setattr(sys, "argv", ["x", "--substrates-fasta", fasta, "--output", out])
        assert run_ortholog_grouping.main() == 0

        rows = _read_csv(out)
        assert rows == []

    def test_single_substrate_writes_OG_001(self, monkeypatch, tmp_dir):
        fasta = os.path.join(tmp_dir, "one.fasta")
        write_fasta({"P1": "MKT"}, fasta)
        out = os.path.join(tmp_dir, "out.csv")

        monkeypatch.setattr(sys, "argv", ["x", "--substrates-fasta", fasta, "--output", out])
        assert run_ortholog_grouping.main() == 0

        rows = _read_csv(out)
        assert len(rows) == 1
        assert rows[0]["locus_tag"] == "P1"
        assert rows[0]["ortholog_group"] == "OG_001"

    def test_blast_unavailable_writes_singleton_csv_and_exits_0(self, monkeypatch, tmp_dir):
        # ≥2 substrates so we reach the BLAST branch, then the binary probe
        # raises BlastpUnavailableError → main() must soft-skip with rc=0.
        fasta = os.path.join(tmp_dir, "two.fasta")
        write_fasta({"P1": "MKT", "P2": "GGG"}, fasta)
        out = os.path.join(tmp_dir, "out.csv")

        def _missing(*a, **k):
            raise FileNotFoundError("blastp")

        monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", _missing)
        monkeypatch.setattr(sys, "argv", ["x", "--substrates-fasta", fasta, "--output", out])
        assert run_ortholog_grouping.main() == 0

        rows = _read_csv(out)
        assert {r["locus_tag"] for r in rows} == {"P1", "P2"}
        # Singleton output: each protein in its own group.
        assert {r["ortholog_group"] for r in rows} == {"OG_001", "OG_002"}
        assert all(r["og_n_members"] == "1" for r in rows)

    def test_blast_hard_failure_propagates(self, monkeypatch, tmp_dir):
        # BLAST is "installed" (probe passes) but makeblastdb returncode != 0.
        # Must raise RuntimeError, NOT silently write empty output.
        fasta = os.path.join(tmp_dir, "two.fasta")
        write_fasta({"P1": "MKT", "P2": "GGG"}, fasta)
        out = os.path.join(tmp_dir, "out.csv")

        class _Result:
            def __init__(self, rc, stderr=""):
                self.returncode = rc
                self.stderr = stderr
                self.stdout = ""

        call_count = {"n": 0}

        def _stub_run(cmd, *a, **k):
            call_count["n"] += 1
            # First two calls are -version probes for blastp + makeblastdb
            if "-version" in cmd:
                return _Result(0)
            # Next call is makeblastdb → fail
            return _Result(1, stderr="bad input fasta")

        monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", _stub_run)
        monkeypatch.setattr(sys, "argv", ["x", "--substrates-fasta", fasta, "--output", out])
        with pytest.raises(RuntimeError, match="makeblastdb failed"):
            run_ortholog_grouping.main()
