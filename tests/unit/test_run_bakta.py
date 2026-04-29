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

from run_bakta import parse_bakta_tsv, parse_dbxrefs, write_proteins_fasta  # noqa: E402, F401


# Bakta TSV fixture mirroring the real Bakta 1.12 output format:
#   - Comment lines starting with "# " (software / DB / DOI / URL).
#   - A header line that *also* starts with "#" but is the column header
#     itself (e.g. "#Sequence Id\tType\t...").
#   - Type values are lowercase ("cds", "sorf", "trna").
_BAKTA_TSV_FIXTURE = (
    "# Annotated with Bakta\n"
    "# Software: v1.12.0\n"
    "# Database: v6.0, light\n"
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "contig_1\tcds\t100\t450\t+\tGENE_00001\tlepA\tGTP-binding protein LepA\tEC:3.6.5.n1, COG:COG0481\n"
    "contig_1\tcds\t500\t800\t-\tGENE_00002\t\thypothetical protein\t\n"
    "contig_1\ttrna\t900\t970\t+\tGENE_00003\t\ttRNA-Ala\t\n"
    "contig_2\tsorf\t50\t200\t+\tGENE_00004\t\tsmall hypothetical protein\t\n"
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
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tcds\t1\t99\t+\tG1\t\t\t\n"
            )

        entries = parse_bakta_tsv(tsv)
        assert entries[0]["product"] == "hypothetical protein"

    def test_row_without_locus_tag_is_skipped(self, tmp_dir):
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tcds\t1\t99\t+\t\t\tsomething\t\n"
                "contig_1\tcds\t200\t299\t+\tG2\t\tsomething else\t\n"
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
            "ec_numbers",
            "cog_ids",
            "go_terms",
            "kegg_ko",
            "refseq_ids",
            "pfam_ids",
        }
        assert required <= set(entries[0].keys())

    def test_dbxrefs_surfaced_per_entry(self, tmp_dir):
        """The DbXrefs column should populate ec_numbers / cog_ids etc."""
        tsv = os.path.join(tmp_dir, "sample.tsv")
        with open(tsv, "w") as f:
            f.write(_BAKTA_TSV_FIXTURE)

        entries = parse_bakta_tsv(tsv)
        e1 = next(e for e in entries if e["locus_tag"] == "GENE_00001")
        assert e1["ec_numbers"] == ["3.6.5.n1"]
        assert e1["cog_ids"] == ["COG0481"]

        e2 = next(e for e in entries if e["locus_tag"] == "GENE_00002")
        assert e2["ec_numbers"] == []
        assert e2["cog_ids"] == []


class TestParseDbxrefs:
    def test_empty_field_returns_empty_lists(self):
        result = parse_dbxrefs("")
        assert result["ec_numbers"] == []
        assert result["cog_ids"] == []
        assert result["go_terms"] == []

    def test_whitespace_field_returns_empty_lists(self):
        assert parse_dbxrefs("   ")["ec_numbers"] == []

    def test_single_ec_number(self):
        assert parse_dbxrefs("EC:1.2.3.4")["ec_numbers"] == ["1.2.3.4"]

    def test_multiple_prefixes_in_one_field(self):
        result = parse_dbxrefs("EC:1.2.3.4, COG:COG1234, KEGG:K99999")
        assert result["ec_numbers"] == ["1.2.3.4"]
        assert result["cog_ids"] == ["COG1234"]
        assert result["kegg_ko"] == ["K99999"]

    def test_go_id_with_double_prefix_preserved(self):
        """GO IDs embed another 'GO:' prefix; partition on first colon
        should leave the canonical 'GO:0001234' form intact."""
        result = parse_dbxrefs("GO:GO:0001234")
        assert result["go_terms"] == ["GO:0001234"]

    def test_unknown_prefix_ignored(self):
        """SO (sequence-ontology) is a real Bakta prefix that ssign
        doesn't surface — confirm it doesn't bleed into other fields."""
        result = parse_dbxrefs("SO:0001217")
        assert result["ec_numbers"] == []
        assert result["cog_ids"] == []
        assert result["refseq_ids"] == []
        assert result["uniref_ids"] == []
        assert result["uniparc_ids"] == []

    def test_uniref_and_uniparc_surfaced(self):
        """UniRef + UniParc are surfaced for provenance (most-populated
        cross-references in Bakta light-DB output)."""
        result = parse_dbxrefs(
            "UniRef:UniRef50_X1, UniParc:UPI000123"
        )
        assert result["uniref_ids"] == ["UniRef50_X1"]
        assert result["uniparc_ids"] == ["UPI000123"]

    def test_refseq_captured(self):
        assert parse_dbxrefs("RefSeq:WP_123.1")["refseq_ids"] == ["WP_123.1"]

    def test_pfam_captured(self):
        assert parse_dbxrefs("Pfam:PF12345")["pfam_ids"] == ["PF12345"]

    def test_multiple_entries_same_prefix(self):
        result = parse_dbxrefs("EC:1.1.1.1, EC:2.2.2.2")
        assert result["ec_numbers"] == ["1.1.1.1", "2.2.2.2"]

    def test_entry_without_colon_skipped(self):
        """Malformed entries missing the 'prefix:value' colon should not crash."""
        result = parse_dbxrefs("EC:1.2.3.4, malformed")
        assert result["ec_numbers"] == ["1.2.3.4"]


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
