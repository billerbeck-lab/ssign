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

from run_eggnog import parse_eggnog_annotations  # noqa: E402


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
        """Every entry must contain the five expected fields for downstream consumers."""
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
        }
        for e in entries:
            assert required <= set(e.keys())
