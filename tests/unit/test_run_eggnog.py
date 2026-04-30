"""Unit tests for run_eggnog.py — parser for emapper .annotations output.

Does not exercise the `emapper.py` subprocess itself (that requires a
~50 GB EggNOG database and a working eggnog-mapper install). Tests target
the pure-Python parser that reads the tool's TSV output.
"""

import os
import sys

SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts")
)
sys.path.insert(0, SCRIPTS_DIR)

from run_eggnog import (  # noqa: E402, F401
    _split_rich_field,
    load_substrate_ids,
    parse_eggnog_annotations,
    write_substrates_only_fasta,
)

# Representative emapper 2.1.x .annotations fixture. Lines starting with
# '##' are comments emitted by emapper before and after the data block;
# the parser must skip them. The column header line begins with '#query'.
_EGGNOG_FIXTURE = (
    "## emapper-2.1.12\n"
    "## /path/to/emapper.py -i proteins.faa -o sample --data_dir /db\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
    "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
    "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
    "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
    "GENE_00001\t83333.eco0001\t1e-180\t500.0\tCOG0012@1|root\t2|Bacteria\t"
    "M\tDNA-directed RNA polymerase subunit beta\trpoB\tGO:0003677\t"
    "2.7.7.6\tko:K03040\t-\t-\t-\t-\t-\t-\t-\t-\tRNA_pol_Rpb1\n"
    "GENE_00003\t83333.eco0003\t2e-50\t200.0\tCOG0084@1|root\t2|Bacteria\t"
    "L\tUncharacterized conserved protein\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"
    "-\t-\t-\t-\n"
    "## 2 queries annotated\n"
)


class TestParseEggnogAnnotations:
    def test_skips_comment_lines(self, tmp_dir):
        """Lines starting with '##' before and after the data block are skipped."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        assert len(entries) == 2  # only the two data rows, no ## comments

    def test_extracts_expected_fields(self, tmp_dir):
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e1 = next(e for e in entries if e["protein_id"] == "GENE_00001")
        assert e1["seed_ortholog"] == "83333.eco0001"
        assert e1["evalue"] == "1e-180"
        assert e1["description"] == "DNA-directed RNA polymerase subunit beta"
        assert e1["preferred_name"] == "rpoB"

    def test_dash_description_preserved_as_dash(self, tmp_dir):
        """A literal '-' in the Description column is a real emapper value
        (no description found) and should pass through unchanged."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e3 = next(e for e in entries if e["protein_id"] == "GENE_00003")
        # GENE_00003 has a real description in the fixture
        assert e3["description"] == "Uncharacterized conserved protein"
        # But its Preferred_name is "-" — should be kept as empty string or "-"
        assert e3["preferred_name"] == "-"

    def test_empty_description_defaults_to_dash(self, tmp_dir):
        """An empty Description field falls back to '-' for consistency."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
                "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                "GENE_X\t83333.x\t1e-10\t50\tCOG@1\t2|B\tS\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            )

        entries = parse_eggnog_annotations(path)
        assert entries[0]["description"] == "-"

    def test_empty_annotations_returns_empty_list(self, tmp_dir):
        """A valid emapper output with zero annotated queries returns []."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(
                "## emapper-2.1.12\n"
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
                "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                "## 0 queries annotated\n"
            )

        entries = parse_eggnog_annotations(path)
        assert entries == []

    def test_output_shape_is_stable(self, tmp_dir):
        """Every entry must contain the expected fields for downstream consumers."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        required = {
            "protein_id",
            "seed_ortholog",
            "evalue",
            "description",
            "preferred_name",
            "cog_category",
            "ec_numbers",
            "kegg_ko",
            "go_terms",
            "pfam_ids",
        }
        for e in entries:
            assert required <= set(e.keys())


class TestRichFieldsSurfaced:
    def test_annotated_row_populates_rich_fields(self, tmp_dir):
        """GENE_00001 in the fixture has EC, KEGG_ko, GOs, PFAMs, COG_category."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e1 = next(e for e in entries if e["protein_id"] == "GENE_00001")
        assert e1["cog_category"] == "M"
        assert e1["ec_numbers"] == ["2.7.7.6"]
        assert e1["kegg_ko"] == ["ko:K03040"]
        assert e1["go_terms"] == ["GO:0003677"]
        assert e1["pfam_ids"] == ["RNA_pol_Rpb1"]

    def test_dash_sentinels_become_empty(self, tmp_dir):
        """GENE_00003 has COG_category but '-' for every multi-value field."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e3 = next(e for e in entries if e["protein_id"] == "GENE_00003")
        assert e3["cog_category"] == "L"
        assert e3["ec_numbers"] == []
        assert e3["kegg_ko"] == []
        assert e3["go_terms"] == []
        assert e3["pfam_ids"] == []

    def test_multi_value_row_splits_on_comma(self, tmp_dir):
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
                "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                "GENE_X\t83333.x\t1e-10\t50\tCOG@1\t2|B\tM\tmulti\tname\t"
                "GO:0001,GO:0002\t1.1.1.1,2.2.2.2\tko:K1,ko:K2\t-\t-\t-\t"
                "-\t-\t-\t-\t-\tPF1234,PF5678\n"
            )

        entries = parse_eggnog_annotations(path)
        e = entries[0]
        assert e["go_terms"] == ["GO:0001", "GO:0002"]
        assert e["ec_numbers"] == ["1.1.1.1", "2.2.2.2"]
        assert e["kegg_ko"] == ["ko:K1", "ko:K2"]
        assert e["pfam_ids"] == ["PF1234", "PF5678"]


class TestSplitRichField:
    def test_none_returns_empty_list(self):
        assert _split_rich_field(None) == []

    def test_empty_string_returns_empty_list(self):
        assert _split_rich_field("") == []
        assert _split_rich_field("   ") == []

    def test_dash_sentinel_returns_empty_list(self):
        assert _split_rich_field("-") == []

    def test_single_value_returned_as_one_element_list(self):
        assert _split_rich_field("1.1.1.1") == ["1.1.1.1"]

    def test_comma_separated_values_split(self):
        assert _split_rich_field("a,b,c") == ["a", "b", "c"]

    def test_whitespace_around_entries_stripped(self):
        assert _split_rich_field(" a , b , c ") == ["a", "b", "c"]

    def test_empty_fragments_dropped(self):
        assert _split_rich_field("a,,b") == ["a", "b"]


class TestLoadSubstrateIds:
    """Substrate-TSV locus_tag extraction for filtering the protein FASTA."""

    def test_reads_locus_tags_from_tsv(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tother_col\n")
            f.write("GENE_00001\tx\n")
            f.write("GENE_00002\ty\n")
            f.write("GENE_00003\tz\n")
        assert load_substrate_ids(path) == {"GENE_00001", "GENE_00002", "GENE_00003"}

    def test_empty_locus_tag_rows_skipped(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tproduct\n")
            f.write("G1\thit\n")
            f.write("\tghost\n")  # blank locus_tag — skip
            f.write("G2\thit\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_empty_file_returns_empty_set(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\n")  # header only
        assert load_substrate_ids(path) == set()


class TestWriteSubstratesOnlyFasta:
    def _write_proteins(self, path, seqs):
        with open(path, "w") as f:
            for k, v in seqs.items():
                f.write(f">{k}\n{v}\n")

    def test_writes_only_requested_ids(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        self._write_proteins(
            src,
            {"G1": "MKTLL", "G2": "MFVFL", "G3": "MQKAL", "G4": "MSRLK"},
        )
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G2", "G4"}, out)
        assert n == 2
        with open(out) as f:
            body = f.read()
        assert ">G2" in body
        assert ">G4" in body
        assert ">G1" not in body
        assert ">G3" not in body

    def test_missing_ids_are_silently_dropped(self, tmp_dir):
        """If the substrates TSV references a locus_tag that doesn't exist
        in the protein FASTA (possible if CDS sets drifted), skip rather
        than crash — the output just has fewer rows."""
        src = os.path.join(tmp_dir, "proteins.faa")
        self._write_proteins(src, {"G1": "MKT"})
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "MISSING"}, out)
        assert n == 1

    def test_empty_substrate_set_writes_empty_fasta(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        self._write_proteins(src, {"G1": "MKT"})
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, set(), out)
        assert n == 0
        assert os.path.exists(out)
        assert os.path.getsize(out) == 0
