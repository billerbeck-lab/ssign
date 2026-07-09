"""Unit tests for the controlled functional vocabularies (COG / KEGG / consensus)."""

from ssign_app.scripts.ssign_lib.functional_vocab import (
    APPARATUS_BUCKET,
    COG_CATEGORY_NAMES,
    cog_category_names,
    consensus_bucket,
    kegg_descriptions,
    kegg_ko_name,
    kegg_kos,
)


def test_cog_single_letter_maps_to_name():
    assert cog_category_names("M") == ["Cell wall/membrane/envelope biogenesis"]
    assert cog_category_names("U") == ["Intracellular trafficking, secretion, and vesicular transport"]


def test_cog_multi_letter_splits_each_letter():
    # EggNOG packs several categories into one string with no separator.
    assert cog_category_names("MU") == [COG_CATEGORY_NAMES["M"], COG_CATEGORY_NAMES["U"]]


def test_cog_blank_and_missing_are_unannotated():
    for v in ("", "-", None, "nan"):
        assert cog_category_names(v) == ["Unannotated"]


def test_cog_unknown_letter_falls_back_to_function_unknown():
    assert cog_category_names("@") == [COG_CATEGORY_NAMES["S"]]


def test_kegg_strips_ko_prefix_and_splits():
    assert kegg_kos("ko:K10953;ko:K12516;ko:K15125") == ["K10953", "K12516", "K15125"]


def test_kegg_blank_returns_empty():
    for v in ("", "-", None):
        assert kegg_kos(v) == []


def test_consensus_known_category_passes_through():
    assert consensus_bucket("Nuclease (InterProScan)") == "Nuclease"  # suffix stripped
    assert consensus_bucket("Autotransporter passenger") == "Autotransporter passenger"
    assert consensus_bucket("Protease/Peptidase (Bakta, EggNOG)") == "Protease/Peptidase"


def test_consensus_apparatus_passes_through():
    # The v2 voter emits "Apparatus-associated" directly for machinery.
    assert consensus_bucket("Apparatus-associated") == APPARATUS_BUCKET
    assert consensus_bucket("Apparatus-associated (EggNOG)") == APPARATUS_BUCKET


def test_kegg_ko_name_from_bundled_table():
    # Bundled KEGG list/ko table maps IDs to readable function names.
    assert kegg_ko_name("K10953") == "RTX toxin RtxA"
    assert kegg_ko_name("K00001") == "alcohol dehydrogenase"  # gene-symbol prefix + EC bracket dropped


def test_kegg_ko_name_unknown_falls_back_to_id():
    assert kegg_ko_name("K99999") == "K99999"


def test_kegg_descriptions_maps_and_dedupes():
    out = kegg_descriptions("ko:K10953;ko:K10953;ko:K15125")
    assert out == ["RTX toxin RtxA", "filamentous hemagglutinin"]  # deduped, order-preserving
    assert kegg_descriptions("-") == []


def test_consensus_noise_blank_and_unclassified_become_other():
    # Unrecognised strings, blanks, and the voter's honest "Unclassified" floor
    # all collapse to a single "Other" catch-all in the functional chart.
    assert consensus_bucket("Low Calcium Response") == "Other"
    assert consensus_bucket("Unclassified") == "Other"
    assert consensus_bucket("") == "Other"
    assert consensus_bucket(None) == "Other"
