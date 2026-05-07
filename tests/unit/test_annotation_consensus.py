"""Unit tests for annotation_consensus.py.

Exercises the pure-regex `classify_description()` categoriser and the
`compute_consensus()` voting function. Covers the 17 broad functional
categories the pipeline votes over, plus the agreeing-count math,
concordance ratio, confidence tiers, and evidence-keywords output.
"""

from annotation_consensus import (  # noqa: E402, F401
    CATEGORY_PATTERNS,
    classify_description,
    compute_consensus,
)


class TestClassifyDescription:
    """Spot-check every broad category by feeding a canonical description."""

    def test_empty_string_returns_empty_list(self):
        assert classify_description("") == []
        assert classify_description("   ") == []

    def test_none_returns_empty_list(self):
        assert classify_description(None) == []

    def test_adhesin_keyword(self):
        assert "Adhesin" in classify_description("fimbrial adhesin subunit")
        assert "Adhesin" in classify_description("type IV pilin")

    def test_autotransporter_keyword(self):
        assert "Autotransporter" in classify_description("Autotransporter domain-containing protein")

    def test_protease_keyword(self):
        assert "Protease" in classify_description("serine protease")
        assert "Protease" in classify_description("metallopeptidase M10")

    def test_lipase_keyword(self):
        assert "Lipase/Esterase" in classify_description("extracellular lipase")
        assert "Lipase/Esterase" in classify_description("phospholipase A1")

    def test_nuclease_keyword(self):
        assert "Nuclease" in classify_description("DNase I family nuclease")

    def test_glycoside_hydrolase_keyword(self):
        assert "Glycoside hydrolase" in classify_description("glycoside hydrolase family 18")
        assert "Glycoside hydrolase" in classify_description("lysozyme")

    def test_toxin_keyword(self):
        assert "Toxin" in classify_description("alpha-hemolysin HlyA")
        assert "Toxin" in classify_description("Rtx toxin")

    def test_transporter_keyword(self):
        assert "Transporter" in classify_description("ABC transporter permease")
        assert "Transporter" in classify_description("MFS efflux pump")

    def test_secretion_system_keyword(self):
        assert "Secretion system" in classify_description("type I secretion system ABC transporter")
        assert "Secretion system" in classify_description("VirB4 homologue")

    def test_oxidoreductase_keyword(self):
        assert "Oxidoreductase" in classify_description("NADH dehydrogenase subunit")

    def test_chaperone_keyword(self):
        assert "Chaperone" in classify_description("usher chaperone")

    def test_regulatory_keyword(self):
        assert "Regulatory" in classify_description("LysR-family transcriptional regulator")
        assert "Regulatory" in classify_description("two-component response regulator")

    def test_hypothetical_keyword(self):
        assert "Hypothetical" in classify_description("hypothetical protein")
        assert "Hypothetical" in classify_description("DUF1234 domain protein")

    def test_multiple_categories_possible(self):
        """A description can match more than one category."""
        cats = classify_description("hemolysin-like transporter")
        assert "Toxin" in cats
        assert "Transporter" in cats

    def test_fallback_for_unmatched_description(self):
        """Descriptions with no pattern match fall back to a title-cased stub."""
        cats = classify_description("WeirdlyNamed Unknown Thingy")
        assert cats  # must not be empty
        assert cats[0] != ""

    def test_fallback_ignores_nan_like_strings(self):
        assert classify_description("nan") == []
        assert classify_description("NaN") == []
        assert classify_description("None") == []

    def test_17_categories_accounted_for(self):
        """Sanity check: CATEGORY_PATTERNS hasn't silently lost a category."""
        category_names = {cat for cat, _pat in CATEGORY_PATTERNS}
        assert len(category_names) == 17


class TestComputeConsensusEmpty:
    def test_empty_input_returns_none_tier(self):
        result = compute_consensus({})
        assert result["confidence_tier"] == "None"
        assert result["n_tools_agreeing"] == 0
        assert result["n_tools_with_hits"] == 0
        assert result["broad_annotation"] == ""
        assert result["concordance_ratio"] == 0.0


class TestComputeConsensusSingleTool:
    def test_one_tool_one_category_low_tier(self):
        result = compute_consensus({"BLASTp": "serine protease"})
        assert result["broad_annotation"] == "Protease"
        assert result["n_tools_agreeing"] == 1
        assert result["n_tools_with_hits"] == 1
        assert result["concordance_ratio"] == 1.0
        assert result["confidence_tier"] == "Low"
        assert "BLASTp" in result["broad_consensus_annotation"]


class TestComputeConsensusAgreement:
    def test_two_tools_agree_medium_tier(self):
        result = compute_consensus({"BLASTp": "serine protease", "InterProScan": "peptidase domain"})
        assert result["broad_annotation"] == "Protease"
        assert result["n_tools_agreeing"] == 2
        assert result["confidence_tier"] == "Medium"

    def test_three_tools_agree_high_tier(self):
        result = compute_consensus(
            {
                "BLASTp": "serine protease",
                "InterProScan": "peptidase domain",
                "Bakta": "endopeptidase",
            }
        )
        assert result["broad_annotation"] == "Protease"
        assert result["n_tools_agreeing"] == 3
        assert result["confidence_tier"] == "High"

    def test_tools_disagree_most_common_wins(self):
        """When tools vote different categories, the most-supported wins."""
        result = compute_consensus(
            {
                "BLASTp": "serine protease",
                "InterProScan": "peptidase",
                "Bakta": "DUF1234 domain",  # Hypothetical
            }
        )
        assert result["broad_annotation"] == "Protease"
        assert result["n_tools_agreeing"] == 2
        assert result["n_tools_with_hits"] == 3

    def test_concordance_ratio_math(self):
        """concordance_ratio = n_agreeing / n_with_hits."""
        result = compute_consensus(
            {
                "BLASTp": "serine protease",
                "InterProScan": "peptidase",
                "Bakta": "DUF1234 domain",
                "EggNOG": "uncharacterized protein",
            }
        )
        # 2 tools vote Protease out of 4 total = 0.5
        assert result["n_tools_agreeing"] == 2
        assert result["n_tools_with_hits"] == 4
        assert result["concordance_ratio"] == 0.5


class TestComputeConsensusEvidence:
    def test_evidence_keywords_lists_all_categories_with_tools(self):
        result = compute_consensus(
            {
                "BLASTp": "serine protease",
                "InterProScan": "peptidase",
                "Bakta": "hemolysin",
            }
        )
        evidence = result["evidence_keywords"]
        assert "Protease" in evidence
        assert "Toxin" in evidence
        assert "BLASTp" in evidence
        assert "Bakta" in evidence

    def test_consensus_annotation_names_supporting_tools(self):
        result = compute_consensus(
            {
                "BLASTp": "hemolysin",
                "Bakta": "cytolysin",
            }
        )
        assert result["broad_annotation"] == "Toxin"
        # Both tools support Toxin; both should appear in the consensus
        assert "BLASTp" in result["broad_consensus_annotation"]
        assert "Bakta" in result["broad_consensus_annotation"]

    def test_detailed_annotation_caps_at_15(self):
        """If many distinct terms are present, only up to 15 appear."""
        many = {f"Tool{i}": f"unique term {i} something" for i in range(20)}
        result = compute_consensus(many)
        parts = result["detailed_annotation"].split(" | ")
        assert len(parts) <= 15


class TestComputeConsensusWithNewTools:
    """Regression: Bakta and EggNOG integrate cleanly into voting."""

    def test_bakta_and_eggnog_vote_like_any_other_tool(self):
        result = compute_consensus(
            {
                "Bakta": "ABC transporter permease",
                "EggNOG": "efflux pump",
                "BLASTp": "MFS transporter",
            }
        )
        assert result["broad_annotation"] == "Transporter"
        assert result["n_tools_agreeing"] == 3
        assert result["confidence_tier"] == "High"

    def test_mixed_new_and_old_tools(self):
        result = compute_consensus(
            {
                "Bakta": "alpha-hemolysin",
                "EggNOG": "Rtx toxin",
                "pLM-BLAST": "hemolysin-like protein",
                "InterProScan": "hypothetical",
            }
        )
        assert result["broad_annotation"] == "Toxin"
        assert result["n_tools_agreeing"] == 3
