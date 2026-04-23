"""Unit tests for run_bakta.py — parser and FASTA writer.

Does not exercise the `bakta` subprocess itself (that requires a ~2 GB
database and a working Bakta install). Tests target the pure-Python helpers
that parse Bakta's TSV output and rewrite its FAA into ssign format.
"""

import os
import sys

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from run_bakta import parse_bakta_tsv, write_proteins_fasta  # noqa: E402


# Minimal Bakta TSV fixture. Real Bakta output also contains header lines
# starting with '#' before the column header — we skip them the way
# csv.DictReader would (first non-'#' line is the header).
_BAKTA_TSV_FIXTURE = (
    "Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "contig_1\tCDS\t100\t450\t+\tGENE_00001\tlepA\tGTP-binding protein LepA\tEC:3.6.5.n1, COG:COG0481\n"
    "contig_1\tCDS\t500\t800\t-\tGENE_00002\t\thypothetical protein\t\n"
    "contig_1\ttRNA\t900\t970\t+\tGENE_00003\t\ttRNA-Ala\t\n"
    "contig_2\tsORF\t50\t200\t+\tGENE_00004\t\tsmall hypothetical protein\t\n"
)

_BAKTA_FAA_FIXTURE = (
    ">GENE_00001 GTP-binding protein LepA [contig_1]\n"
    "MKNIRNFSIIAHIDHGKSTLSDRLIQIC\n"
    "GEGDDRLM*\n"
    ">GENE_00002 hypothetical protein [contig_1]\n"
    "MQKALLSAAWLVLLPSTAHA\n"
    ">GENE_00004 small hypothetical protein [contig_2]\n"
    "MKLPLAVLG*\n"
)


class TestParseBaktaTsv:
    def test_extracts_cds_rows(self, tmp_dir):
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(_BAKTA_TSV_FIXTURE)

        entries = parse_bakta_tsv(tsv)
        locus_tags = {e["locus_tag"] for e in entries}
        assert locus_tags == {"GENE_00001", "GENE_00002", "GENE_00004"}

    def test_skips_non_protein_features(self, tmp_dir):
        """tRNAs and other non-CDS/sORF rows must not appear in entries."""
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(_BAKTA_TSV_FIXTURE)

        entries = parse_bakta_tsv(tsv)
        assert all(e["locus_tag"] != "GENE_00003" for e in entries)

    def test_coordinates_converted_to_zero_based_start(self, tmp_dir):
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(_BAKTA_TSV_FIXTURE)

        entries = parse_bakta_tsv(tsv)
        e1 = next(e for e in entries if e["locus_tag"] == "GENE_00001")
        assert e1["start"] == 99  # TSV has 100, parser converts to 0-based
        assert e1["end"] == 450

    def test_empty_product_defaults_to_hypothetical(self, tmp_dir):
        """An empty Product field falls back to 'hypothetical protein'."""
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(
                "Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tCDS\t1\t99\t+\tG1\t\t\t\n"
            )

        entries = parse_bakta_tsv(tsv)
        assert entries[0]["product"] == "hypothetical protein"

    def test_row_without_locus_tag_is_skipped(self, tmp_dir):
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(
                "Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tCDS\t1\t99\t+\t\t\tsomething\t\n"
                "contig_1\tCDS\t200\t299\t+\tG2\t\tsomething else\t\n"
            )

        entries = parse_bakta_tsv(tsv)
        assert len(entries) == 1
        assert entries[0]["locus_tag"] == "G2"

    def test_includes_expected_fields(self, tmp_dir):
        """Output shape stays stable for downstream pipeline consumers."""
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(_BAKTA_TSV_FIXTURE)

        entries = parse_bakta_tsv(tsv)
        required = {
            "locus_tag",
            "protein_id",
            "gene",
            "product",
            "contig",
            "start",
            "end",
            "strand",
        }
        assert required <= set(entries[0].keys())


class TestWriteProteinsFasta:
    def _fixture_files(self, tmp_dir):
        tsv = os.path.join(tmp_dir, "sample.tsv")
        faa = os.path.join(tmp_dir, "sample.faa")
        with open(tsv, "w") as f:
            f.write(_BAKTA_TSV_FIXTURE)
        with open(faa, "w") as f:
            f.write(_BAKTA_FAA_FIXTURE)
        return tsv, faa

    def test_writes_protein_per_cds_entry(self, tmp_dir):
        tsv, faa = self._fixture_files(tmp_dir)
        entries = parse_bakta_tsv(tsv)
        out = os.path.join(tmp_dir, "out.faa")

        n = write_proteins_fasta(faa, entries, out)
        assert n == 3  # GENE_00001, GENE_00002, GENE_00004 (tRNA excluded)

    def test_output_headers_are_locus_tags(self, tmp_dir):
        tsv, faa = self._fixture_files(tmp_dir)
        entries = parse_bakta_tsv(tsv)
        out = os.path.join(tmp_dir, "out.faa")

        write_proteins_fasta(faa, entries, out)

        with open(out) as f:
            content = f.read()
        assert ">GENE_00001\n" in content
        assert ">GENE_00002\n" in content
        assert ">GENE_00004\n" in content
        # No description appended after the locus tag
        assert ">GENE_00001 GTP-binding" not in content

    def test_stop_codon_trimmed(self, tmp_dir):
        tsv, faa = self._fixture_files(tmp_dir)
        entries = parse_bakta_tsv(tsv)
        out = os.path.join(tmp_dir, "out.faa")

        write_proteins_fasta(faa, entries, out)
        with open(out) as f:
            content = f.read()
        # Trailing '*' that Bakta appends must be stripped
        assert "*\n" not in content

    def test_missing_sequence_is_skipped_not_errored(self, tmp_dir):
        """An entry with no sequence in the FAA should be dropped, not crash."""
        tsv, faa = self._fixture_files(tmp_dir)
        entries = parse_bakta_tsv(tsv)
        entries.append(
            {
                "locus_tag": "MISSING",
                "protein_id": "",
                "gene": "",
                "product": "ghost",
                "contig": "contig_9",
                "start": 0,
                "end": 0,
                "strand": "+",
            }
        )
        out = os.path.join(tmp_dir, "out.faa")

        n = write_proteins_fasta(faa, entries, out)
        assert n == 3
        with open(out) as f:
            assert ">MISSING" not in f.read()
