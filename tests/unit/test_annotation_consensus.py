"""Unit tests for annotation_consensus.py (v2, tool-weighted).

Exercises the single-category `classify_description()` matcher and the
tool-weighted `compute_consensus()` voter: credibility weighting, the
structural-fold down-weight, machinery-by-component-identity (no `t[1-9]ss`
hijack), RTX-by-activity, the Hypothetical/Unclassified floors, and the
preserved output-dict schema.
"""

import pandas as pd
from annotation_consensus import (
    APPARATUS,
    CATEGORY_NAMES,
    CATEGORY_PATTERNS,
    UNCLASSIFIED,
    classify_description,
    compute_consensus,
)


class TestClassifyDescription:
    """classify_description returns a SINGLE category (most-specific first) or None."""

    def test_empty_and_none_return_none(self):
        assert classify_description("") is None
        assert classify_description("   ") is None
        assert classify_description(None) is None

    def test_no_match_returns_none_not_a_minted_label(self):
        # v1 minted a title-case stub from the first 3 words; v2 must not.
        assert classify_description("WeirdlyNamed Unknown Thingy") is None

    def test_returns_single_string_not_a_list(self):
        result = classify_description("serine protease")
        assert isinstance(result, str)
        assert result == "Protease/Peptidase"

    def test_specific_beats_generic_first_match_wins(self):
        # "protease" (specific) is ordered before "transporter/channel" (generic);
        # a serine-protease autotransporter classifies by its passenger FUNCTION.
        assert classify_description("serine protease autotransporter EspP") == "Protease/Peptidase"

    def test_effector_categories(self):
        assert classify_description("T6SS amidase effector Tae1") == "Peptidoglycan hydrolase"
        assert classify_description("mono-ADP-ribosyltransferase toxin") == "ADP-ribosyltransferase"
        assert classify_description("OspF phosphothreonine lyase") == "Phosphothreonine lyase"
        assert classify_description("NleB glycosyltransferase") == "Glycosyltransferase"

    def test_pore_forming_vs_generic_toxin(self):
        assert classify_description("alpha-hemolysin HlyA") == "Pore-forming toxin"
        assert classify_description("adenylate cyclase toxin CyaA") == "Pore-forming toxin"
        assert classify_description("colicin E1") == "Toxin (other)"

    def test_autotransporter_passenger_only_when_unspecific(self):
        assert classify_description("Autotransporter domain-containing protein") == "Autotransporter passenger"

    def test_hypothetical(self):
        assert classify_description("hypothetical protein") == "Hypothetical"
        assert classify_description("DUF1234 domain protein") == "Hypothetical"


class TestMachineryByIdentity:
    """Machinery is detected by component identity, never by naming an SS type."""

    def test_component_identity_routes_to_apparatus(self):
        assert classify_description("VgrG protein") == APPARATUS
        assert classify_description("Hcp secreted protein") == APPARATUS
        assert classify_description("PAAR-repeat protein") == APPARATUS
        assert classify_description("flagellar hook protein FlgE") == APPARATUS

    def test_translocator_context_routes_to_apparatus(self):
        assert classify_description("haemolysin secretion protein ShlB") == APPARATUS
        assert classify_description("two-partner secretion protein") == APPARATUS

    def test_effector_named_after_system_is_not_hijacked(self):
        # v1's `t[1-9]ss` keyword swallowed effectors into "Secretion system".
        assert classify_description("T6SS amidase effector (Tae)") == "Peptidoglycan hydrolase"
        assert classify_description("Type VI secretion lipase effector Tle1") == "Lipase/Phospholipase"


class TestRTXByActivity:
    """No standalone RTX bucket — RTX proteins classify by their real activity."""

    def test_rtx_metalloprotease_is_protease(self):
        assert classify_description("RTX serralysin metalloprotease") == "Protease/Peptidase"

    def test_rtx_cytolysin_is_pore_forming(self):
        assert classify_description("RTX cytolysin") == "Pore-forming toxin"


class TestWeightedVote:
    """Tool credibility weights the vote; structure never outvotes function."""

    def test_functional_tool_outweighs_structural(self):
        # InterProScan lipase (w2) vs HHpred_PDB pore-forming fold (w1) -> Lipase.
        result = compute_consensus({"InterProScan": "triacylglycerol lipase", "HHpred_PDB": "aerolysin pore-forming"})
        assert result["broad_annotation"] == "Lipase/Phospholipase"

    def test_two_tier2_beat_one_tier3(self):
        result = compute_consensus(
            {"HHpred_Pfam": "peptidase M4", "InterProScan": "endopeptidase", "pLM-BLAST": "beta-barrel"}
        )
        assert result["broad_annotation"] == "Protease/Peptidase"

    def test_fold_name_is_capped_and_cannot_win(self):
        # Bakta is Tier-1 (w3) but reports a structural SUPERFAMILY name (capped to
        # w1); a Tier-2 functional call (w2) must beat it.
        result = compute_consensus({"Bakta": "Peroxidase-like superfamily", "InterProScan": "serine protease"})
        assert result["broad_annotation"] == "Protease/Peptidase"

    def test_duplicate_descriptions_counted_once(self):
        # Two Tier-1 tools echo the SAME hemolysin string -> one deduped vote (w3),
        # so two independent Tier-2 lipase calls (w2+w2=4) win. Without dedup the
        # hemolysin would double to w6 and wrongly win.
        result = compute_consensus(
            {
                "BLASTp": "alpha-hemolysin",
                "Bakta": "alpha-hemolysin",
                "InterProScan": "lipase",
                "HHpred_Pfam": "esterase",
            }
        )
        assert result["broad_annotation"] == "Lipase/Phospholipase"

    def test_translocator_agreed_by_good_tools_wins(self):
        result = compute_consensus({"EggNOG": "haemolysin secretion protein", "HHpred_PDB": "beta-helix adhesin fold"})
        assert result["broad_annotation"] == APPARATUS


class TestFloors:
    """Hypothetical and Unclassified never outrank a functional call."""

    def test_single_function_beats_multiple_hypothetical(self):
        result = compute_consensus(
            {"BLASTp": "hypothetical protein", "Bakta": "uncharacterized protein", "EggNOG": "serine protease"}
        )
        assert result["broad_annotation"] == "Protease/Peptidase"

    def test_only_hypothetical_yields_hypothetical(self):
        # (DUF4123 is deliberately NOT generic — it is the Eag chaperone family —
        # so use a generic DUF here.)
        result = compute_consensus({"Bakta": "hypothetical protein", "GBFF": "DUF1234 domain protein"})
        assert result["broad_annotation"] == "Hypothetical"

    def test_no_functional_or_hypothetical_call_yields_unclassified(self):
        result = compute_consensus({"pLM-BLAST": "zztop random string"})
        assert result["broad_annotation"] == UNCLASSIFIED
        assert result["broad_consensus_annotation"] == UNCLASSIFIED  # no empty "()"
        assert result["n_tools_agreeing"] == 0


class TestConfidenceTierFromWeightedSupport:
    def test_two_tier1_agree_high(self):
        result = compute_consensus({"BLASTp": "serine protease", "Bakta": "metallopeptidase"})
        assert result["confidence_tier"] == "High"  # 3 + 3 = 6
        assert result["n_tools_agreeing"] == 2

    def test_one_tier1_medium(self):
        result = compute_consensus({"Bakta": "serine protease"})
        assert result["confidence_tier"] == "Medium"  # weight 3

    def test_one_tier2_low(self):
        result = compute_consensus({"InterProScan": "peptidase domain"})
        assert result["confidence_tier"] == "Low"  # weight 2

    def test_one_tier3_low(self):
        result = compute_consensus({"HHpred_PDB": "aerolysin"})
        assert result["confidence_tier"] == "Low"  # weight 1


class TestOutputSchema:
    """The output dict keys are an interface integrate_annotations depends on."""

    EXPECTED_KEYS = {
        "broad_annotation",
        "broad_consensus_annotation",
        "detailed_annotation",
        "detailed_consensus_annotation",
        "evidence_keywords",
        "n_tools_agreeing",
        "n_tools_with_hits",
        "concordance_ratio",
        "confidence_tier",
    }

    def test_keys_present_for_populated_input(self):
        result = compute_consensus({"BLASTp": "serine protease"})
        assert set(result.keys()) == self.EXPECTED_KEYS

    def test_keys_present_for_empty_input(self):
        assert set(compute_consensus({}).keys()) == self.EXPECTED_KEYS

    def test_empty_input_is_no_hits_sentinel(self):
        result = compute_consensus({})
        assert result["confidence_tier"] == "no_hits"
        assert result["broad_annotation"] == ""
        assert result["n_tools_with_hits"] == 0
        assert result["concordance_ratio"] == 0.0

    def test_no_hits_survives_csv_roundtrip(self, tmp_path):
        # Regression: the sentinel must be "no_hits", not "None" (pandas coerces
        # bare "None" to NaN on read, blanking the column in the master CSV).
        df = pd.DataFrame([compute_consensus({})])
        out = tmp_path / "consensus.csv"
        df.to_csv(out, index=False)
        roundtrip = pd.read_csv(out)
        assert roundtrip["confidence_tier"].iloc[0] == "no_hits"
        assert not pd.isna(roundtrip["confidence_tier"].iloc[0])

    def test_concordance_ratio_math(self):
        # 2 tools classify Protease, 1 Hypothetical, 1 no-match -> 2/4.
        result = compute_consensus(
            {
                "BLASTp": "serine protease",
                "InterProScan": "peptidase",
                "Bakta": "DUF1234 domain",
                "EggNOG": "zztop nonsense",
            }
        )
        assert result["broad_annotation"] == "Protease/Peptidase"
        assert result["n_tools_agreeing"] == 2
        assert result["n_tools_with_hits"] == 4
        assert result["concordance_ratio"] == 0.5

    def test_evidence_keywords_lists_categories_with_tools(self):
        result = compute_consensus(
            {"BLASTp": "serine protease", "InterProScan": "peptidase", "Bakta": "alpha-hemolysin"}
        )
        evidence = result["evidence_keywords"]
        assert "Protease/Peptidase" in evidence
        assert "Pore-forming toxin" in evidence
        assert "BLASTp" in evidence and "Bakta" in evidence

    def test_consensus_annotation_names_supporting_tools(self):
        result = compute_consensus({"BLASTp": "alpha-hemolysin", "Bakta": "cytolysin"})
        assert result["broad_annotation"] == "Pore-forming toxin"
        assert "BLASTp" in result["broad_consensus_annotation"]
        assert "Bakta" in result["broad_consensus_annotation"]

    def test_detailed_annotation_caps_at_15(self):
        many = {f"Tool{i}": f"unique term {i} something" for i in range(20)}
        result = compute_consensus(many)
        assert len(result["detailed_annotation"].split(" | ")) <= 15


class TestTaxonomy:
    def test_no_standalone_rtx_bucket(self):
        names = {cat for cat, _ in CATEGORY_PATTERNS}
        assert not any("RTX" in n for n in names)

    def test_new_effector_categories_present(self):
        names = {cat for cat, _ in CATEGORY_PATTERNS}
        for expected in (
            "Peptidoglycan hydrolase",
            "ADP-ribosyltransferase",
            "Glycosyltransferase",
            "Phosphothreonine lyase",
            "Ubiquitin-pathway",
            "Beta-lactamase",
            "Pore-forming toxin",
            "Phage/mobile element",
        ):
            assert expected in names

    def test_category_names_export_resolves_apparatus_sentinel(self):
        assert "__APPARATUS__" not in CATEGORY_NAMES
        assert APPARATUS in CATEGORY_NAMES
        assert UNCLASSIFIED not in CATEGORY_NAMES  # a floor, not a classify output
