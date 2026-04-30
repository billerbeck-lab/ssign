"""Unit tests for run_plm_blast.py — CSV parser.

Does not exercise the `plmblast.py` subprocess itself (that requires the
pLM-BLAST install plus a ~20 GB ECOD70 database). Tests target the
pure-Python helper that parses pLM-BLAST's CSV output.
"""

import os
import sys

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from run_plm_blast import (  # noqa: E402, F401
    load_substrate_ids,
    parse_plmblast_csv,
    write_substrates_only_fasta,
)

# Representative pLM-BLAST CSV fixture. Real output columns verified on
# first integration run — swap in the actual captured output as a fixture
# once we run pLM-BLAST against the T1SS fixture on CX3.
_PLMBLAST_FIXTURE = (
    "qid,sid,score,qstart,qend,tstart,tend\n"
    "GENE_00001,ecod_1a2bA1,0.912,10,145,5,140\n"
    "GENE_00001,ecod_3c4dB2,0.871,20,150,10,140\n"
    "GENE_00003,ecod_5e6fC1,0.789,1,200,1,199\n"
)


class TestParsePlmblastCsv:
    def test_extracts_all_hits(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        assert len(entries) == 3

    def test_multiple_hits_per_protein_preserved(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        hits_for_gene1 = [e for e in entries if e["protein_id"] == "GENE_00001"]
        assert len(hits_for_gene1) == 2

    def test_captures_all_expected_fields(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        first = next(e for e in entries if e["protein_id"] == "GENE_00001")
        assert first["target_id"] == "ecod_1a2bA1"
        assert first["score"] == "0.912"
        assert first["qstart"] == "10"
        assert first["qend"] == "145"
        assert first["tstart"] == "5"
        assert first["tend"] == "140"

    def test_empty_qid_row_skipped(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(
                "qid,sid,score,qstart,qend,tstart,tend\n"
                ",ecod_empty,0.5,1,10,1,10\n"
                "GENE_X,ecod_real,0.9,1,100,1,100\n"
            )

        entries = parse_plmblast_csv(path)
        assert len(entries) == 1
        assert entries[0]["protein_id"] == "GENE_X"

    def test_empty_csv_returns_empty_list(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write("qid,sid,score,qstart,qend,tstart,tend\n")

        assert parse_plmblast_csv(path) == []

    def test_output_shape_is_stable(self, tmp_dir):
        """Every entry has the six expected fields for downstream consumers."""
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        required = {
            "protein_id",
            "target_id",
            "score",
            "qstart",
            "qend",
            "tstart",
            "tend",
        }
        for e in entries:
            assert required <= set(e.keys())


class TestSubstrateFiltering:
    """Runs before pLM-BLAST itself; verifies we only search substrates, not whole genome."""

    def test_load_substrate_ids_basic(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tsample\n")
            f.write("G1\tsample1\n")
            f.write("G2\tsample1\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_load_substrate_ids_skips_empty_rows(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\n")
            f.write("G1\n")
            f.write("\n")
            f.write("G2\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_write_substrates_only_fasta_filters(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        with open(src, "w") as f:
            f.write(">G1\nMKT\n>G2\nMFV\n>G3\nMQK\n")
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "G3"}, out)
        assert n == 2
        with open(out) as f:
            body = f.read()
        assert ">G1" in body and ">G3" in body
        assert ">G2" not in body

    def test_write_substrates_only_fasta_missing_ids_skipped(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        with open(src, "w") as f:
            f.write(">G1\nMKT\n")
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "NOT_IN_FASTA"}, out)
        assert n == 1
