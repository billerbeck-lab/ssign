"""Tests for dedup_sequences.py.

Five public surfaces:

- `deduplicate_dict` — in-memory dedup (used by run_blastp / run_eggnog
  before tool calls).
- `deduplicate_fasta` — file-on-disk variant.
- `expand_results_dict` / `expand_results_tsv` / `expand_results_csv` —
  broadcast tool output from representatives to all duplicate members.

The load-bearing invariant for the rest of the pipeline: representative
selection is first-seen (order-preserving), and every duplicate member
ends up with the same annotation as its representative.
"""

import csv
import os
import sys

from hypothesis import given
from hypothesis import strategies as st

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from _helpers import write_tsv  # noqa: E402
from dedup_sequences import (  # noqa: E402
    deduplicate_dict,
    deduplicate_fasta,
    expand_results_csv,
    expand_results_dict,
    expand_results_tsv,
)

# ---------------------------------------------------------------------------
# deduplicate_dict
# ---------------------------------------------------------------------------


class TestDeduplicateDict:
    def test_no_duplicates_is_identity(self):
        seqs = {"A": "MKT", "B": "MFV", "C": "MGA"}
        unique, groups = deduplicate_dict(seqs)
        assert unique == seqs
        assert all(len(group) == 1 for group in groups.values())

    def test_duplicate_sequences_collapse_to_first_seen(self):
        seqs = {"A": "MKT", "B": "MKT", "C": "MFV"}
        unique, groups = deduplicate_dict(seqs)
        # Representative is whichever ID was inserted first for that seq
        assert "A" in unique
        assert "B" not in unique
        assert "C" in unique
        # Group for the MKT-rep contains both A and B (in insertion order)
        assert groups["A"] == ["A", "B"]
        assert groups["C"] == ["C"]

    def test_three_way_duplicate(self):
        seqs = {"FIRST": "MKT", "SECOND": "MKT", "THIRD": "MKT"}
        unique, groups = deduplicate_dict(seqs)
        assert len(unique) == 1
        assert "FIRST" in unique
        assert sorted(groups["FIRST"]) == ["FIRST", "SECOND", "THIRD"]

    def test_empty_input(self):
        assert deduplicate_dict({}) == ({}, {})

    def test_groups_partition_input(self):
        # Every input ID lands in exactly one group; the union of groups
        # equals the input ID set.
        seqs = {"A": "MKT", "B": "MKT", "C": "MFV", "D": "MGA", "E": "MFV"}
        _, groups = deduplicate_dict(seqs)
        all_members = sorted(m for group in groups.values() for m in group)
        assert all_members == sorted(seqs.keys())


# ---------------------------------------------------------------------------
# expand_results_dict — broadcast back to all duplicates
# ---------------------------------------------------------------------------


class TestExpandResultsDict:
    def test_broadcasts_annotation_to_every_member(self):
        seq_groups = {"REP": ["REP", "DUP1", "DUP2"]}
        results = {"REP": {"locus_tag": "REP", "score": 99.7}}
        expanded = expand_results_dict(results, seq_groups, id_key="locus_tag")
        assert set(expanded.keys()) == {"REP", "DUP1", "DUP2"}
        assert all(e["score"] == 99.7 for e in expanded.values())

    def test_id_key_updated_per_member(self):
        seq_groups = {"REP": ["REP", "DUP1"]}
        results = {"REP": {"locus_tag": "REP", "score": 1.0}}
        expanded = expand_results_dict(results, seq_groups, id_key="locus_tag")
        assert expanded["DUP1"]["locus_tag"] == "DUP1"
        assert expanded["REP"]["locus_tag"] == "REP"

    def test_id_not_in_seq_groups_passes_through(self):
        # Singleton (never had duplicates) → not in seq_groups → kept as-is
        seq_groups = {"REP": ["REP"]}
        results = {
            "REP": {"locus_tag": "REP", "score": 1.0},
            "ORPHAN": {"locus_tag": "ORPHAN", "score": 2.0},
        }
        expanded = expand_results_dict(results, seq_groups, id_key="locus_tag")
        assert expanded["ORPHAN"]["score"] == 2.0


# ---------------------------------------------------------------------------
# Property: dedup → annotate-rep-only → expand → all duplicates carry it
# ---------------------------------------------------------------------------


@given(
    n_unique=st.integers(min_value=1, max_value=8),
    duplications=st.lists(st.integers(min_value=1, max_value=4), min_size=1, max_size=8),
)
def test_dedup_expand_round_trip_preserves_annotations(n_unique, duplications):
    """For any input where N unique sequences have varying duplication
    counts, the dedup → expand round-trip must end with every input ID
    carrying its sequence's representative annotation."""
    # Construct: n_unique distinct sequences, each duplicated some number of times
    seqs = {}
    expected_groups = {}  # seq_index → list of input IDs
    next_id = 0
    for seq_idx in range(n_unique):
        # synthesise a sequence unique to seq_idx; AAs are arbitrary
        seq = "M" + "K" * (seq_idx + 1)
        n_copies = duplications[seq_idx % len(duplications)]
        for copy_idx in range(n_copies):
            pid = f"P{next_id:04d}"
            seqs[pid] = seq
            expected_groups.setdefault(seq_idx, []).append(pid)
            next_id += 1

    unique, groups = deduplicate_dict(seqs)
    # Annotate representatives only (mimic running a tool on unique seqs)
    annotated = {rep: {"locus_tag": rep, "score": float(rep[1:])} for rep in unique}
    expanded = expand_results_dict(annotated, groups, id_key="locus_tag")

    # Every input ID has an entry in expanded, and its score matches its rep's score
    assert set(expanded.keys()) == set(seqs.keys())
    for seq_idx, member_ids in expected_groups.items():
        # All members in this group share the same score (the rep's)
        scores = {expanded[mid]["score"] for mid in member_ids}
        assert len(scores) == 1


# ---------------------------------------------------------------------------
# deduplicate_fasta — file round-trip
# ---------------------------------------------------------------------------


class TestDeduplicateFasta:
    def test_writes_only_unique_sequences(self, tmp_dir):
        input_fasta = os.path.join(tmp_dir, "in.fasta")
        with open(input_fasta, "w") as f:
            f.write(">A\nMKT\n>B\nMKT\n>C\nMFV\n")
        output_fasta = os.path.join(tmp_dir, "out.fasta")
        groups = deduplicate_fasta(input_fasta, output_fasta)
        # Two unique sequences; A is the rep for MKT, C is the rep for MFV
        with open(output_fasta) as f:
            written = f.read()
        assert ">A\n" in written
        assert ">B\n" not in written
        assert ">C\n" in written
        assert sorted(groups["A"]) == ["A", "B"]
        assert groups["C"] == ["C"]

    def test_no_duplicates_writes_all(self, tmp_dir):
        input_fasta = os.path.join(tmp_dir, "in.fasta")
        with open(input_fasta, "w") as f:
            f.write(">A\nMKT\n>B\nMFV\n")
        output_fasta = os.path.join(tmp_dir, "out.fasta")
        groups = deduplicate_fasta(input_fasta, output_fasta)
        with open(output_fasta) as f:
            written = f.read()
        assert ">A\n" in written
        assert ">B\n" in written
        assert all(len(g) == 1 for g in groups.values())


# ---------------------------------------------------------------------------
# expand_results_tsv / _csv — file expansion
# ---------------------------------------------------------------------------


_RESULTS_FIELDS = ["locus_tag", "score"]


class TestExpandResultsFile:
    def test_tsv_expansion_broadcasts_rows(self, tmp_dir):
        seq_groups = {"REP": ["REP", "DUP1", "DUP2"]}
        in_path = write_tsv(
            os.path.join(tmp_dir, "in.tsv"),
            _RESULTS_FIELDS,
            [{"locus_tag": "REP", "score": "99.7"}],
        )
        out_path = os.path.join(tmp_dir, "out.tsv")
        expand_results_tsv(in_path, out_path, seq_groups)
        with open(out_path) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        assert {r["locus_tag"] for r in rows} == {"REP", "DUP1", "DUP2"}
        assert all(r["score"] == "99.7" for r in rows)

    def test_tsv_expansion_handles_missing_input_silently(self, tmp_dir):
        # Pre-condition: many wrappers skip annotation when the substrate set
        # is empty, leaving expand_results_tsv pointed at a missing file.
        # Must be a silent no-op rather than a crash.
        out_path = os.path.join(tmp_dir, "out.tsv")
        expand_results_tsv(
            os.path.join(tmp_dir, "does_not_exist.tsv"),
            out_path,
            {"REP": ["REP"]},
        )
        assert not os.path.exists(out_path)

    def test_csv_expansion_broadcasts_rows(self, tmp_dir):
        seq_groups = {"REP": ["REP", "DUP1"]}
        in_path = os.path.join(tmp_dir, "in.csv")
        out_path = os.path.join(tmp_dir, "out.csv")
        with open(in_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["locus_tag", "score"])
            writer.writeheader()
            writer.writerow({"locus_tag": "REP", "score": "1.0"})
        expand_results_csv(in_path, out_path, seq_groups)
        with open(out_path) as f:
            rows = list(csv.DictReader(f))
        assert {r["locus_tag"] for r in rows} == {"REP", "DUP1"}

    def test_orphan_id_passes_through_unchanged(self, tmp_dir):
        seq_groups = {"REP": ["REP", "DUP1"]}
        in_path = write_tsv(
            os.path.join(tmp_dir, "in.tsv"),
            _RESULTS_FIELDS,
            [
                {"locus_tag": "REP", "score": "1.0"},
                {"locus_tag": "ORPHAN", "score": "2.0"},
            ],
        )
        out_path = os.path.join(tmp_dir, "out.tsv")
        expand_results_tsv(in_path, out_path, seq_groups)
        with open(out_path) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        by_locus = {r["locus_tag"]: r for r in rows}
        assert by_locus["ORPHAN"]["score"] == "2.0"
        assert by_locus["DUP1"]["score"] == "1.0"


# ---------------------------------------------------------------------------
# Hash determinism — same sequence always picks same representative across runs
# ---------------------------------------------------------------------------


@given(
    seqs=st.lists(
        st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=3, max_size=20),
        min_size=2,
        max_size=10,
    )
)
def test_dedup_is_order_preserving_first_seen(seqs):
    """The representative for each sequence-group is always the first-seen
    ID. Re-ordering the input dict reorders groups but never changes the
    representative-selection rule."""
    d1 = {f"P{i}": s for i, s in enumerate(seqs)}
    _, groups1 = deduplicate_dict(d1)
    d2 = {f"P{i}": s for i, s in reversed(list(enumerate(seqs)))}
    _, groups2 = deduplicate_dict(d2)

    for d, groups in [(d1, groups1), (d2, groups2)]:
        for rep, _members in groups.items():
            first_seen = next(pid for pid in d if d[pid] == d[rep])
            assert rep == first_seen
