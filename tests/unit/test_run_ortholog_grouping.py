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
import run_ortholog_grouping
from run_ortholog_grouping import (
    BlastpUnavailableError,
    _find_blast_binary,
    cluster_union_find,
    compute_group_stats,
)

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta

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

    def test_equal_sized_groups_get_stable_ids(self):
        # Three 2-member groups all tie on size. `groups` is keyed by a Union-Find
        # representative and built by iterating a SET of ids, so its order follows string
        # hashing and differs between interpreter runs; sorting on size alone let the
        # OG_NNN labels shuffle. Ties break on the smallest member, so the assignment is
        # a function of the data alone.
        groups = {"r1": {"p004", "p005"}, "r2": {"p000", "p001"}, "r3": {"p002", "p003"}}
        hits = [("p000", "p001", 90.0, 90.0), ("p002", "p003", 90.0, 90.0), ("p004", "p005", 90.0, 90.0)]
        stats = compute_group_stats(groups, hits, {f"p{i:03d}" for i in range(6)})
        assert [(s["ortholog_group"], s["members"]) for s in stats] == [
            ("OG_001", "p000;p001"),
            ("OG_002", "p002;p003"),
            ("OG_003", "p004;p005"),
        ]

    def test_size_still_outranks_the_tiebreak(self):
        # The tiebreak must only apply within a size, never reorder across sizes: the
        # 2-member group wins OG_001 even though its members sort after the singleton's.
        groups = {"r1": {"z1"}, "r2": {"z8", "z9"}}
        stats = compute_group_stats(groups, hits=[("z8", "z9", 90.0, 90.0)], all_protein_ids={"z1", "z8", "z9"})
        assert [(s["ortholog_group"], s["n_members"]) for s in stats] == [("OG_001", 2), ("OG_002", 1)]

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


class TestEdgeDensity:
    """Single linkage needs only a spanning tree, so a group can be a chain rather
    than a family. edge_density is the only column that tells them apart."""

    def test_every_pair_hit_is_100(self):
        groups = {"rep1": {"A", "B", "C"}}
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 90.0, 90.0), ("B", "C", 90.0, 90.0), ("A", "C", 90.0, 90.0)],
            all_protein_ids={"A", "B", "C"},
        )
        assert stats[0]["edge_density"] == 100.0

    def test_chain_is_below_100(self):
        # A-B-C with no A-C hit: 2 of 3 possible pairs. Same members and same
        # mean_pident as the clique above, so only edge_density separates them.
        groups = {"rep1": {"A", "B", "C"}}
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 90.0, 90.0), ("B", "C", 90.0, 90.0)],
            all_protein_ids={"A", "B", "C"},
        )
        assert stats[0]["edge_density"] == 66.7
        assert stats[0]["mean_pident"] == 90.0

    def test_reciprocal_hits_counted_once(self):
        # BLAST reports each pair in both directions; that is one edge, not two,
        # so a fully-linked pair must not exceed 100.
        groups = {"rep1": {"A", "B"}}
        stats = compute_group_stats(
            groups,
            hits=[("A", "B", 90.0, 90.0), ("B", "A", 88.0, 90.0)],
            all_protein_ids={"A", "B"},
        )
        assert stats[0]["edge_density"] == 100.0

    def test_singleton_is_100(self):
        # No pairs at all. Reported as 100 to match mean_pident's convention for
        # the same case, rather than 0, which would read as "maximally chained".
        stats = compute_group_stats({"rep1": {"A"}}, hits=[], all_protein_ids={"A"})
        assert stats[0]["edge_density"] == 100.0


# ---------------------------------------------------------------------------
# run_local_blast — hit filtering
# ---------------------------------------------------------------------------


def _blast_row(qseqid, sseqid, pident, aln_len, qlen, slen, q_aligned=None, s_aligned=None):
    """One BLAST outfmt-6 line in the column order of BLAST_OUTFMT.

    `q_aligned`/`s_aligned` are how many RESIDUES of each sequence the alignment spans,
    which is what coverage is measured on. They default to `aln_len`, i.e. a gapless
    alignment where columns and residues coincide. Pass them explicitly to model a gapped
    alignment, where `length` counts gap columns that neither sequence contributes a
    residue to.
    """
    q_aligned = aln_len if q_aligned is None else q_aligned
    s_aligned = aln_len if s_aligned is None else s_aligned
    return "\t".join(
        str(v)
        for v in (
            qseqid,
            sseqid,
            pident,
            aln_len,
            0,
            0,
            1,
            q_aligned,
            1,
            s_aligned,
            "1e-40",
            100.0,
            qlen,
            slen,
        )
    )


def _stub_blast(monkeypatch, rows):
    """Make run_local_blast see `rows` as the BLASTp output, without BLAST+."""

    class _Result:
        returncode = 0
        stderr = ""
        stdout = ""

    def _run(cmd, *a, **k):
        if "-out" in cmd and "-query" in cmd:  # the blastp call
            with open(cmd[cmd.index("-out") + 1], "w") as fh:
                fh.write("\n".join(rows) + ("\n" if rows else ""))
        return _Result()

    monkeypatch.setattr(run_ortholog_grouping.subprocess, "run", _run)


class TestRunLocalBlastCoverage:
    """Coverage must be required of BOTH sequences.

    BLAST reports every pair in both directions and a hit passing either way is
    kept, so checking the query alone collapses to "covers 70% of the SHORTER
    sequence" and leaves the longer one unconstrained. The short-vs-long pair here
    is the real one that produced a 118-member chained group on a 74-genome
    Xanthobacter run: 177 aa aligning over its whole length to 8.6% of 2069 aa.
    """

    def test_short_fully_covered_by_long_fragment_is_rejected(self, monkeypatch, tmp_dir):
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"short": "MKT", "long": "GGG"}, fasta)
        _stub_blast(monkeypatch, [_blast_row("short", "long", 45.0, 177, 177, 2069)])

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert hits == []

    def test_same_pair_reported_the_other_way_is_also_rejected(self, monkeypatch, tmp_dir):
        # Query coverage is now the failing side rather than the passing one. The
        # verdict must not depend on which way round BLAST happened to report it.
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"short": "MKT", "long": "GGG"}, fasta)
        _stub_blast(monkeypatch, [_blast_row("long", "short", 45.0, 177, 2069, 177)])

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert hits == []

    def test_genuine_ortholog_still_passes(self, monkeypatch, tmp_dir):
        # Similar lengths, aligned nearly end to end: 96% of both.
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(monkeypatch, [_blast_row("A", "B", 80.0, 480, 500, 500)])

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert len(hits) == 1
        query, subject, pident, coverage = hits[0]
        assert (query, subject, pident) == ("A", "B", 80.0)
        assert coverage == pytest.approx(96.0)

    def test_coverage_reported_is_the_lower_of_the_two(self, monkeypatch, tmp_dir):
        # 400/500 = 80% of the query, 400/450 = 88.9% of the subject. Both clear
        # 70, and the tuple must carry the binding one so a caller cannot read the
        # flattering side.
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(monkeypatch, [_blast_row("A", "B", 75.0, 400, 500, 450)])

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert hits[0][3] == pytest.approx(80.0)

    def test_identity_threshold_still_applies(self, monkeypatch, tmp_dir):
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(monkeypatch, [_blast_row("A", "B", 30.0, 480, 500, 500)])

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert hits == []


class TestCoverageCountsResiduesNotColumns:
    """Coverage is residues aligned, not outfmt-6 `length`.

    `length` counts alignment COLUMNS including gap columns, which no sequence
    contributes a residue to. Dividing it by a sequence length therefore credits gaps as
    coverage: the gappier the alignment the more it over-credits, which is worst for
    exactly the weak pairs min_qcov exists to exclude.
    """

    def test_gapped_alignment_no_longer_over_credits(self, monkeypatch, tmp_dir):
        # 400 aligned columns, but 100 of them are gaps in the query, so the query
        # contributes only 300 of its 500 residues.
        #   by columns:  400/500 = 80%  -> passes 70 and should not
        #   by residues: 300/500 = 60%  -> correctly rejected
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(
            monkeypatch,
            [_blast_row("A", "B", 45.0, 400, 500, 400, q_aligned=300, s_aligned=400)],
        )

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert hits == []

    def test_coverage_cannot_exceed_100_percent(self, monkeypatch, tmp_dir):
        # 600 aligned columns spanning a 500 aa query and a 550 aa subject -- the extra
        # columns are gaps. By columns BOTH sides read over 100% (120.0 and 109.1), so
        # even taking the lower of the two leaves an impossible coverage; by residues
        # each is exactly fully covered.
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(
            monkeypatch,
            [_blast_row("A", "B", 85.0, 600, 500, 550, q_aligned=500, s_aligned=550)],
        )

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert len(hits) == 1
        assert hits[0][3] == pytest.approx(100.0)

    def test_gapless_alignment_is_unaffected(self, monkeypatch, tmp_dir):
        # Where there are no gaps, columns and residues coincide and the answer is the
        # one it always was. This is the case every other test in this file models.
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(monkeypatch, [_blast_row("A", "B", 80.0, 480, 500, 500)])

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert hits[0][3] == pytest.approx(96.0)

    def test_self_hits_dropped(self, monkeypatch, tmp_dir):
        fasta = os.path.join(tmp_dir, "in.faa")
        write_fasta({"A": "MKT", "B": "GGG"}, fasta)
        _stub_blast(
            monkeypatch,
            [_blast_row("A", "A", 100.0, 500, 500, 500), _blast_row("A", "B", 85.0, 490, 500, 500)],
        )

        hits = run_ortholog_grouping.run_local_blast(fasta, min_pident=40.0, min_qcov=70.0)
        assert [(h[0], h[1]) for h in hits] == [("A", "B")]


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
