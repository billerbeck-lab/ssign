"""Tests for integrate_annotations.py.

Covers the four merge layers: base substrate load, GBFF annotation
back-fill from gene_info, protein-sequence merge, and the per-file
left-join over arbitrary annotation TSVs. Also covers _compute_consensus —
particularly the "hypothetical protein" filter that prevents generic
GBFF annotations from voting in the consensus.
"""

import os

import pandas as pd
from _helpers import run_script_main, write_tsv
from integrate_annotations import (
    TOOL_HIT_COLUMNS,
    _add_functional_names,
    _compute_consensus,
)
from integrate_annotations import main as integrate_main

SUBSTRATE_FIELDS = [
    "locus_tag",
    "sample_id",
    "tool",
    "nearby_ss_types",
    "product",
]


def _substrate(locus_tag, *, tool="DLP", nearby="T2SS", product=""):
    return {
        "locus_tag": locus_tag,
        "sample_id": "test_sample",
        "tool": tool,
        "nearby_ss_types": nearby,
        "product": product,
    }


def _run_integrate(
    monkeypatch,
    tmp_dir,
    substrates,
    *,
    fields=SUBSTRATE_FIELDS,
    gene_info_rows=None,
    gene_info_fields=None,
    proteins_fasta_content=None,
    annotation_files=(),
):
    """Run integrate_annotations.main(); return the output dataframe."""
    sub_filtered = write_tsv(
        os.path.join(tmp_dir, "substrates_filtered.tsv"),
        fields,
        substrates,
    )
    sub_all = write_tsv(
        os.path.join(tmp_dir, "substrates_all.tsv"),
        fields,
        substrates,
    )
    out = os.path.join(tmp_dir, "master.csv")

    argv = [
        "integrate_annotations",
        "--substrates-filtered",
        sub_filtered,
        "--substrates-all",
        sub_all,
        "--sample",
        "test_sample",
        "--output",
        out,
    ]
    if gene_info_rows is not None:
        gi_path = write_tsv(
            os.path.join(tmp_dir, "gene_info.tsv"),
            gene_info_fields,
            gene_info_rows,
        )
        argv += ["--gene-info", gi_path]
    if proteins_fasta_content is not None:
        prot_path = os.path.join(tmp_dir, "proteins.faa")
        with open(prot_path, "w") as f:
            f.write(proteins_fasta_content)
        argv += ["--proteins", prot_path]
    if annotation_files:
        argv += ["--annotations", *annotation_files]

    run_script_main(monkeypatch, integrate_main, argv)
    return pd.read_csv(out)


class TestBaseLoad:
    def test_substrates_filtered_becomes_master_base(self, monkeypatch, tmp_dir):
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001"), _substrate("GENE_002")],
        )
        assert set(df["locus_tag"]) == {"GENE_001", "GENE_002"}

    def test_dlp_localization_probs_preserved(self, monkeypatch, tmp_dir):
        # The substrate producers carry the DLP localization probs onto each
        # substrate row; integrate's left-joins must not drop them (the
        # autotransporter figure reads dlp_extracellular_prob + outer_membrane_prob).
        sub = _substrate("GENE_001")
        sub.update(dlp_extracellular_prob="0.91", outer_membrane_prob="0.88")
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [sub],
            fields=[*SUBSTRATE_FIELDS, "dlp_extracellular_prob", "outer_membrane_prob"],
        )
        row = df[df["locus_tag"] == "GENE_001"].iloc[0]
        assert float(row["dlp_extracellular_prob"]) == 0.91
        assert float(row["outer_membrane_prob"]) == 0.88


class TestGbffAnnotation:
    def test_gbff_annotation_column_used_directly(self, monkeypatch, tmp_dir):
        # Phase 3.3.c re-annotation flow: gene_info already carries the
        # original GenBank product as `gbff_annotation`. Use as-is.
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001")],
            gene_info_fields=["locus_tag", "gbff_annotation"],
            gene_info_rows=[{"locus_tag": "GENE_001", "gbff_annotation": "real annotation"}],
        )
        assert df.loc[df["locus_tag"] == "GENE_001", "gbff_annotation"].iloc[0] == ("real annotation")

    def test_product_renamed_when_no_gbff_annotation_column(self, monkeypatch, tmp_dir):
        # use_input_annotations=True / FASTA-input flow: only `product` exists.
        # Renamed to gbff_annotation on merge.
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001")],
            gene_info_fields=["locus_tag", "product"],
            gene_info_rows=[{"locus_tag": "GENE_001", "product": "kinase"}],
        )
        assert "gbff_annotation" in df.columns
        assert df.loc[df["locus_tag"] == "GENE_001", "gbff_annotation"].iloc[0] == "kinase"

    def test_no_gene_info_yields_no_gbff_annotation(self, monkeypatch, tmp_dir):
        # Without --gene-info, the merge step is skipped entirely, so the
        # gbff_annotation column should not exist (or be all-NaN) in output.
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001")],
        )
        assert "gbff_annotation" not in df.columns or df["gbff_annotation"].isna().all()


class TestProteinSequences:
    def test_sequence_and_aa_length_added(self, monkeypatch, tmp_dir):
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001")],
            proteins_fasta_content=">GENE_001\nMKTLLLT\n",
        )
        row = df.loc[df["locus_tag"] == "GENE_001"].iloc[0]
        assert row["sequence"] == "MKTLLLT"
        assert row["aa_length"] == 7

    def test_missing_sequence_yields_zero_length(self, monkeypatch, tmp_dir):
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001"), _substrate("GENE_002")],
            proteins_fasta_content=">GENE_001\nMKTLLLT\n",
        )
        row = df.loc[df["locus_tag"] == "GENE_002"].iloc[0]
        assert row["aa_length"] == 0


class TestAnnotationFileMerge:
    def test_left_join_preserves_substrate_row_count(self, monkeypatch, tmp_dir):
        # Annotation file has rows for proteins NOT in substrates → no fan-out
        ann_path = write_tsv(
            os.path.join(tmp_dir, "blastp.tsv"),
            ["locus_tag", "blastp_hit_description"],
            [
                {"locus_tag": "GENE_001", "blastp_hit_description": "hemolysin"},
                {"locus_tag": "GENE_999", "blastp_hit_description": "irrelevant"},
            ],
        )
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001")],
            annotation_files=[ann_path],
        )
        assert len(df) == 1
        assert df["blastp_hit_description"].iloc[0] == "hemolysin"

    def test_missing_protein_yields_nan_in_annotation_columns(self, monkeypatch, tmp_dir):
        ann_path = write_tsv(
            os.path.join(tmp_dir, "blastp.tsv"),
            ["locus_tag", "blastp_hit_description"],
            [{"locus_tag": "GENE_001", "blastp_hit_description": "hemolysin"}],
        )
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001"), _substrate("GENE_002")],
            annotation_files=[ann_path],
        )
        gene2 = df.loc[df["locus_tag"] == "GENE_002"].iloc[0]
        assert pd.isna(gene2["blastp_hit_description"])

    def test_overlapping_columns_dropped_from_annotation(self, monkeypatch, tmp_dir):
        # Annotation file ships its own `product` column that already exists
        # on the substrate base. Base wins; ann's product is dropped.
        ann_path = write_tsv(
            os.path.join(tmp_dir, "blastp.tsv"),
            ["locus_tag", "product", "blastp_hit_description"],
            [
                {
                    "locus_tag": "GENE_001",
                    "product": "ann's product (should be dropped)",
                    "blastp_hit_description": "hemolysin",
                }
            ],
        )
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001", product="base's product")],
            annotation_files=[ann_path],
        )
        assert df["product"].iloc[0] == "base's product"
        assert df["blastp_hit_description"].iloc[0] == "hemolysin"

    def test_missing_annotation_file_skipped(self, monkeypatch, tmp_dir):
        # Non-existent path → silently skip, no crash
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("GENE_001")],
            annotation_files=[os.path.join(tmp_dir, "does_not_exist.tsv")],
        )
        assert len(df) == 1


# ---------------------------------------------------------------------------
# _compute_consensus — direct-call tests (skip the CLI overhead)
# ---------------------------------------------------------------------------


class TestComputeConsensus:
    def test_annotation_tools_lists_each_hit_tool(self):
        df = pd.DataFrame(
            [
                {
                    "locus_tag": "GENE_001",
                    "blastp_hit_description": "hemolysin",
                    "pfam_top1_description": "HSP70",
                    "interpro_descriptions": "",  # empty → not counted
                }
            ]
        )
        result = _compute_consensus(df)
        tools = set(result["annotation_tools"].iloc[0].split(","))
        assert tools == {"BLASTp", "HHpred_Pfam"}

    def test_generic_gbff_excluded_from_voting(self):
        # GBFF "hypothetical protein" must NOT appear in annotation_tools
        df = pd.DataFrame(
            [
                {
                    "locus_tag": "GENE_001",
                    "gbff_annotation": "hypothetical protein",
                    "blastp_hit_description": "hemolysin",
                }
            ]
        )
        result = _compute_consensus(df)
        tools = set(result["annotation_tools"].iloc[0].split(","))
        assert "GBFF" not in tools
        assert "BLASTp" in tools

    def test_specific_gbff_does_vote(self):
        df = pd.DataFrame(
            [
                {
                    "locus_tag": "GENE_001",
                    "gbff_annotation": "type IV secretion system protein VirB4",
                }
            ]
        )
        result = _compute_consensus(df)
        tools = set(result["annotation_tools"].iloc[0].split(","))
        assert "GBFF" in tools

    def test_no_hits_yields_empty_annotation_tools(self):
        df = pd.DataFrame([{"locus_tag": "GENE_001"}])
        result = _compute_consensus(df)
        assert result["annotation_tools"].iloc[0] == ""

    def test_uncharacterized_protein_also_filtered(self):
        df = pd.DataFrame(
            [
                {
                    "locus_tag": "GENE_001",
                    "gbff_annotation": "uncharacterized protein",
                }
            ]
        )
        result = _compute_consensus(df)
        assert "GBFF" not in result["annotation_tools"].iloc[0].split(",")


class TestToolHitColumnsContract:
    """Pin the tool→column mapping. Adding a new annotation tool means
    extending TOOL_HIT_COLUMNS so consensus voting picks it up."""

    def test_signalp_not_in_annotation_tools(self):
        # SignalP is a prediction tool, not an annotation tool — must not
        # appear in TOOL_HIT_COLUMNS.
        assert "SignalP" not in TOOL_HIT_COLUMNS

    def test_canonical_tools_present(self):
        # Pin the full set so adding a new tool requires deliberate update.
        expected = {
            "BLASTp",
            "HHpred_Pfam",
            "HHpred_PDB",
            "InterProScan",
            "pLM-BLAST",
            "Bakta",
            "EggNOG",
            "GBFF",
        }
        assert expected <= set(TOOL_HIT_COLUMNS)


class TestT5QualityFlagSort:
    """Rows with a non-empty t5_quality_flag must sort to the bottom of the
    master CSV. Clean rows come first, then worst-flag-last: 'no_signalp'
    before 'unclassified' before 'barrel_only' before 'omp_porin_no_at'.
    """

    def test_flagged_rows_pushed_to_bottom(self, monkeypatch, tmp_dir):
        substrates_with_flag = [
            {**_substrate("BAD_BARREL"), "t5_quality_flag": "barrel_only"},
            {**_substrate("CLEAN_1"), "t5_quality_flag": ""},
            {**_substrate("BAD_OMP"), "t5_quality_flag": "omp_porin_no_at"},
            {**_substrate("CLEAN_2"), "t5_quality_flag": ""},
            {**_substrate("BAD_NOSIG"), "t5_quality_flag": "no_signalp"},
        ]
        # Need to include the new column in fields to write the TSV.
        sub_filtered = write_tsv(
            os.path.join(tmp_dir, "substrates_filtered.tsv"),
            SUBSTRATE_FIELDS + ["t5_quality_flag"],
            substrates_with_flag,
        )
        sub_all = write_tsv(
            os.path.join(tmp_dir, "substrates_all.tsv"),
            SUBSTRATE_FIELDS + ["t5_quality_flag"],
            substrates_with_flag,
        )
        out = os.path.join(tmp_dir, "master.csv")
        run_script_main(
            monkeypatch,
            integrate_main,
            [
                "integrate_annotations",
                "--substrates-filtered",
                sub_filtered,
                "--substrates-all",
                sub_all,
                "--sample",
                "test_sample",
                "--output",
                out,
            ],
        )
        df = pd.read_csv(out, keep_default_na=False)
        # Clean rows come first in input order
        assert df["locus_tag"].tolist() == [
            "CLEAN_1",
            "CLEAN_2",
            "BAD_NOSIG",
            "BAD_BARREL",
            "BAD_OMP",
        ]

    def test_no_flag_column_leaves_order_untouched(self, monkeypatch, tmp_dir):
        """When no t5_quality_flag column exists (e.g. a proximity-only run),
        the sort is a no-op and input order is preserved."""
        df = _run_integrate(
            monkeypatch,
            tmp_dir,
            [_substrate("B"), _substrate("A"), _substrate("C")],
        )
        assert df["locus_tag"].tolist() == ["B", "A", "C"]


class TestFunctionalNames:
    def test_readable_cog_and_kegg_columns_added(self):
        # Raw COG letters + KEGG KO IDs gain human-readable sibling columns.
        df = pd.DataFrame(
            {
                "locus_tag": ["a", "b", "c"],
                "cog_category": ["U", "MU", ""],
                "kegg_ko": ["ko:K10953", "ko:K15125;ko:K12516", "-"],
            }
        )
        out = _add_functional_names(df)
        assert out.loc[0, "cog_category_name"] == "Intracellular trafficking, secretion, and vesicular transport"
        assert "Cell wall/membrane/envelope biogenesis" in out.loc[1, "cog_category_name"]
        assert out.loc[2, "cog_category_name"] == "Unannotated"
        assert out.loc[0, "kegg_function"] == "RTX toxin RtxA"
        assert out.loc[1, "kegg_function"] == "filamentous hemagglutinin; putative surface-exposed virulence protein"
        assert out.loc[2, "kegg_function"] == ""

    def test_no_op_without_cog_or_kegg(self):
        df = pd.DataFrame({"locus_tag": ["a"]})
        out = _add_functional_names(df)
        assert "cog_category_name" not in out.columns
        assert "kegg_function" not in out.columns
