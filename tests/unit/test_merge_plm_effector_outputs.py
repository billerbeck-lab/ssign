"""Unit tests for merge_plm_effector_outputs.py.

Exercises the pure-Python `merge_per_type_outputs()` generator with
synthetic per-type TSVs. Verifies the OR across secretion-system types,
the `flagging_types` list, the `max_stacking` pick, and the skip-missing
behaviour.
"""

import os
import sys

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from merge_plm_effector_outputs import (  # noqa: E402, F401
    merge_per_type_outputs,
    write_merged_tsv,
)

_HEADER = (
    "seq_id\tmodel1\tmodel2\tmodel3\tmodel4\tmodel5\tmodel6\t"
    "stacking\tpasses_threshold\teffector_type\n"
)


def _write_per_type(path, effector_type, rows):
    """Write a per-type TSV (shape matches run_plm_effector.py output).

    `rows` is a list of (seq_id, stacking, passes) tuples.
    """
    with open(path, "w") as f:
        f.write(_HEADER)
        for seq_id, stacking, passes in rows:
            f.write(
                f"{seq_id}\t0\t0\t0\t0\t0\t0\t{stacking}\t"
                f"{'1' if passes else '0'}\t{effector_type}\n"
            )


class TestMergePerTypeOutputs:
    def test_protein_flagged_only_by_one_type(self, tmp_dir):
        t1 = os.path.join(tmp_dir, "t1se.tsv")
        t2 = os.path.join(tmp_dir, "t2se.tsv")
        _write_per_type(t1, "T1SE", [("G1", 0.91, True), ("G2", 0.3, False)])
        _write_per_type(t2, "T2SE", [("G1", 0.4, False), ("G2", 0.2, False)])

        rows = list(merge_per_type_outputs([t1, t2]))
        by_id = {r["seq_id"]: r for r in rows}

        assert by_id["G1"]["passes_threshold"] == "1"
        assert by_id["G1"]["flagging_types"] == "T1SE"
        assert by_id["G2"]["passes_threshold"] == "0"
        assert by_id["G2"]["flagging_types"] == ""

    def test_protein_flagged_by_multiple_types(self, tmp_dir):
        t1 = os.path.join(tmp_dir, "t1se.tsv")
        t4 = os.path.join(tmp_dir, "t4se.tsv")
        _write_per_type(t1, "T1SE", [("G1", 0.85, True)])
        _write_per_type(t4, "T4SE", [("G1", 0.92, True)])

        rows = list(merge_per_type_outputs([t1, t4]))
        assert len(rows) == 1
        r = rows[0]
        assert r["passes_threshold"] == "1"
        assert r["flagging_types"] == "T1SE,T4SE"  # sorted
        # max_stacking should pick the higher of the two
        assert float(r["max_stacking"]) == 0.92

    def test_max_stacking_ignores_non_flagging_types(self, tmp_dir):
        """If a type doesn't flag a protein, its stacking score is not eligible
        for max_stacking even if it's numerically higher."""
        t1 = os.path.join(tmp_dir, "t1se.tsv")
        t3 = os.path.join(tmp_dir, "t3se.tsv")
        _write_per_type(t1, "T1SE", [("G1", 0.55, True)])
        _write_per_type(t3, "T3SE", [("G1", 0.99, False)])  # higher score, not flagged

        rows = list(merge_per_type_outputs([t1, t3]))
        assert rows[0]["flagging_types"] == "T1SE"
        assert float(rows[0]["max_stacking"]) == 0.55

    def test_protein_in_some_files_not_others(self, tmp_dir):
        """Proteins appear in per-type files only if PLM-Effector ran on them
        for that type; the merge must take the union of seq_ids."""
        t1 = os.path.join(tmp_dir, "t1se.tsv")
        t6 = os.path.join(tmp_dir, "t6se.tsv")
        _write_per_type(t1, "T1SE", [("G_common", 0.7, True), ("G_t1_only", 0.8, True)])
        _write_per_type(
            t6, "T6SE", [("G_common", 0.1, False), ("G_t6_only", 0.6, True)]
        )

        rows = list(merge_per_type_outputs([t1, t6]))
        ids = {r["seq_id"] for r in rows}
        assert ids == {"G_common", "G_t1_only", "G_t6_only"}

    def test_missing_input_file_silently_skipped(self, tmp_dir):
        existing = os.path.join(tmp_dir, "t1se.tsv")
        missing = os.path.join(tmp_dir, "nonexistent_t2se.tsv")
        _write_per_type(existing, "T1SE", [("G1", 0.8, True)])

        rows = list(merge_per_type_outputs([existing, missing]))
        assert len(rows) == 1
        assert rows[0]["seq_id"] == "G1"

    def test_empty_input_list_yields_nothing(self):
        assert list(merge_per_type_outputs([])) == []

    def test_locus_tag_mirrors_seq_id(self, tmp_dir):
        """Downstream steps key on `locus_tag`; the merge emits both columns
        for convenience so the merged file slots into the master CSV merge."""
        t1 = os.path.join(tmp_dir, "t1se.tsv")
        _write_per_type(t1, "T1SE", [("G1", 0.9, True)])

        rows = list(merge_per_type_outputs([t1]))
        assert rows[0]["locus_tag"] == rows[0]["seq_id"]


class TestWriteMergedTsv:
    def test_writes_header_and_rows(self, tmp_dir):
        t1 = os.path.join(tmp_dir, "t1se.tsv")
        _write_per_type(t1, "T1SE", [("G1", 0.9, True), ("G2", 0.3, False)])

        out = os.path.join(tmp_dir, "merged.tsv")
        n = write_merged_tsv(merge_per_type_outputs([t1]), out)
        assert n == 2
        assert os.path.exists(out)

        with open(out) as f:
            header = f.readline().strip().split("\t")
        assert header == [
            "seq_id",
            "locus_tag",
            "passes_threshold",
            "flagging_types",
            "max_stacking",
            "effector_type",
        ]

    def test_creates_parent_dirs(self, tmp_dir):
        out = os.path.join(tmp_dir, "sub", "dir", "merged.tsv")
        write_merged_tsv([], out)
        assert os.path.exists(out)
