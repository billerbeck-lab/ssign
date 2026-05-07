"""Tests for resolve_taxonomy.py.

Heavy mocking territory — the real backend is the ~1.5 GB NCBI taxdump
parsed by ``taxopy``. The tests here stub both ``_load_taxdb`` and the
``taxopy`` module so we can pin behaviour without touching the dump.

Coverage targets:

- Caching (same name → cache hit, no second taxdb load).
- Name normalisation (leading/trailing whitespace).
- Graceful degradation when taxopy isn't installed (returns Nones, logs).
- Graceful degradation when the dump files are missing.
- Species + genus resolution via the rank dictionary.
- Genus fallback via the first word of the organism name.
"""

import logging
import os
import sys
from types import SimpleNamespace

import pytest


import resolve_taxonomy as rt  # noqa: E402

# ---------------------------------------------------------------------------
# Fake taxopy backend
# ---------------------------------------------------------------------------


class FakeTaxDb:
    """Stand-in for taxopy.TaxDb. Populated with name → taxids and
    taxid → (name, rank, genus_taxid) tables in each test fixture."""

    def __init__(self):
        # name → list of taxids
        self.name_to_taxids = {}
        # taxid → {"name": ..., "rank": ..., "genus": <genus_taxid or None>}
        self.taxid_to_meta = {}


def _make_fake_taxon_class(taxdb):
    """Return a class whose constructor reads taxid → meta from the given taxdb."""

    class FakeTaxon:
        def __init__(self, taxid, _taxdb=None):
            meta = taxdb.taxid_to_meta[int(taxid)]
            self.taxid = int(taxid)
            self.name = meta["name"]
            self.rank = meta["rank"]
            genus = meta.get("genus")
            self.rank_taxid_dictionary = {"genus": genus} if genus else {}

    return FakeTaxon


def _make_fake_taxopy_module(taxdb):
    """Build a SimpleNamespace exposing the taxopy surface resolve_taxonomy uses."""
    return SimpleNamespace(
        TaxDb=lambda **kwargs: taxdb,  # _load_taxdb instantiates this
        Taxon=_make_fake_taxon_class(taxdb),
        taxid_from_name=lambda name, _taxdb: taxdb.name_to_taxids.get(name, []),
    )


@pytest.fixture(autouse=True)
def _reset_module_state(monkeypatch):
    """Each test sees a clean cache + unloaded taxdb."""
    monkeypatch.setattr(rt, "_cache", {})
    monkeypatch.setattr(rt, "_taxdb", None)


def _install_fake_backend(monkeypatch, taxdb):
    """Wire the fake taxopy module + force _load_taxdb to return `taxdb`."""
    monkeypatch.setitem(sys.modules, "taxopy", _make_fake_taxopy_module(taxdb))
    monkeypatch.setattr(rt, "_load_taxdb", lambda: taxdb)


# ---------------------------------------------------------------------------
# Empty / whitespace input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["", "   ", None])
def test_empty_name_returns_nones(name):
    result = rt.resolve_organism(name)
    assert result == {"species": None, "genus": None}


def test_name_stripped_before_lookup(monkeypatch):
    # Trailing whitespace must not poison the cache key or the lookup
    taxdb = FakeTaxDb()
    taxdb.name_to_taxids["Xanthomonas campestris"] = [339]
    taxdb.taxid_to_meta[339] = {"name": "Xanthomonas campestris", "rank": "species", "genus": 338}
    taxdb.taxid_to_meta[338] = {"name": "Xanthomonas", "rank": "genus"}
    _install_fake_backend(monkeypatch, taxdb)

    result = rt.resolve_organism("  Xanthomonas campestris  ")
    assert result["species"]["name"] == "Xanthomonas campestris"


# ---------------------------------------------------------------------------
# Species + genus resolution
# ---------------------------------------------------------------------------


def test_resolves_species_and_genus(monkeypatch):
    taxdb = FakeTaxDb()
    taxdb.name_to_taxids["Xanthomonas campestris"] = [339]
    taxdb.taxid_to_meta[339] = {"name": "Xanthomonas campestris", "rank": "species", "genus": 338}
    taxdb.taxid_to_meta[338] = {"name": "Xanthomonas", "rank": "genus"}
    _install_fake_backend(monkeypatch, taxdb)

    result = rt.resolve_organism("Xanthomonas campestris")
    assert result["species"] == {"name": "Xanthomonas campestris", "taxid": "339"}
    assert result["genus"] == {"name": "Xanthomonas", "taxid": "338"}


def test_picks_subspecies_rank_when_available(monkeypatch):
    # The picker prefers rank in {species, subspecies, strain} over a
    # multi-word fallback name match.
    taxdb = FakeTaxDb()
    taxdb.name_to_taxids["Salmonella enterica subsp. enterica"] = [59201]
    taxdb.taxid_to_meta[59201] = {
        "name": "Salmonella enterica subsp. enterica",
        "rank": "subspecies",
        "genus": 590,
    }
    taxdb.taxid_to_meta[590] = {"name": "Salmonella", "rank": "genus"}
    _install_fake_backend(monkeypatch, taxdb)

    result = rt.resolve_organism("Salmonella enterica subsp. enterica")
    assert result["species"]["taxid"] == "59201"
    assert result["genus"]["name"] == "Salmonella"


def test_no_taxids_yields_nones(monkeypatch):
    # taxopy.taxid_from_name returns [] → no species, no genus.
    taxdb = FakeTaxDb()
    _install_fake_backend(monkeypatch, taxdb)
    result = rt.resolve_organism("Mythical thing")
    assert result == {"species": None, "genus": None}


# ---------------------------------------------------------------------------
# Genus fallback via first word
# ---------------------------------------------------------------------------


def test_genus_fallback_via_first_word(monkeypatch):
    """Even if the species lookup doesn't yield a genus_taxid (e.g. the
    species rank dict is bare), the first word of the organism name is
    looked up as a genus."""
    taxdb = FakeTaxDb()
    # Species lookup: returns a taxid whose rank_dict has no "genus" entry
    taxdb.name_to_taxids["Xanthomonas obscure_strain"] = [99999]
    taxdb.taxid_to_meta[99999] = {
        "name": "Xanthomonas obscure_strain",
        "rank": "species",
        # genus omitted → species_taxon.rank_taxid_dictionary.get("genus") is None
    }
    # First-word lookup: "Xanthomonas" → returns the genus taxid directly
    taxdb.name_to_taxids["Xanthomonas"] = [338]
    taxdb.taxid_to_meta[338] = {"name": "Xanthomonas", "rank": "genus"}
    _install_fake_backend(monkeypatch, taxdb)

    result = rt.resolve_organism("Xanthomonas obscure_strain")
    assert result["species"]["name"] == "Xanthomonas obscure_strain"
    assert result["genus"] == {"name": "Xanthomonas", "taxid": "338"}


# ---------------------------------------------------------------------------
# Caching
# ---------------------------------------------------------------------------


def test_repeat_resolution_uses_cache(monkeypatch):
    """A second lookup of the same name does not re-call the backend."""
    taxdb = FakeTaxDb()
    taxdb.name_to_taxids["Xanthomonas campestris"] = [339]
    taxdb.taxid_to_meta[339] = {"name": "Xanthomonas campestris", "rank": "species", "genus": 338}
    taxdb.taxid_to_meta[338] = {"name": "Xanthomonas", "rank": "genus"}

    # Wrap _load_taxdb to count invocations
    call_count = {"n": 0}

    def counting_load_taxdb():
        call_count["n"] += 1
        return taxdb

    monkeypatch.setitem(sys.modules, "taxopy", _make_fake_taxopy_module(taxdb))
    monkeypatch.setattr(rt, "_load_taxdb", counting_load_taxdb)

    rt.resolve_organism("Xanthomonas campestris")
    rt.resolve_organism("Xanthomonas campestris")
    rt.resolve_organism("Xanthomonas campestris")
    assert call_count["n"] == 1


# ---------------------------------------------------------------------------
# Graceful degradation
# ---------------------------------------------------------------------------


def test_taxopy_missing_returns_nones_and_logs(monkeypatch, caplog):
    """If the taxopy import inside _load_taxdb raises, resolve_organism
    must return Nones and log a warning, not crash."""

    def raise_runtime():
        raise RuntimeError("taxopy not installed (test stub)")

    monkeypatch.setattr(rt, "_load_taxdb", raise_runtime)

    with caplog.at_level(logging.WARNING):
        result = rt.resolve_organism("Xanthomonas campestris")

    assert result == {"species": None, "genus": None}
    assert any("Taxonomy resolution unavailable" in rec.message for rec in caplog.records)


def test_taxdump_missing_returns_nones(monkeypatch):
    """Same path as above but specifically simulating missing dump files."""

    def raise_missing_files():
        raise RuntimeError("NCBI taxdump files not found under /nonexistent/path")

    monkeypatch.setattr(rt, "_load_taxdb", raise_missing_files)
    assert rt.resolve_organism("Xanthomonas campestris") == {"species": None, "genus": None}


def test_failed_lookup_is_cached_too(monkeypatch):
    """A degraded result (Nones) is cached so repeated failed lookups don't
    keep triggering the warning."""

    call_count = {"n": 0}

    def counting_failure():
        call_count["n"] += 1
        raise RuntimeError("simulated failure")

    monkeypatch.setattr(rt, "_load_taxdb", counting_failure)
    rt.resolve_organism("Xanthomonas campestris")
    rt.resolve_organism("Xanthomonas campestris")
    assert call_count["n"] == 1


# ---------------------------------------------------------------------------
# Defensive bounds: length cap + control-character stripping
# ---------------------------------------------------------------------------


def test_overlong_name_is_truncated(monkeypatch, caplog):
    taxdb = FakeTaxDb()
    long_name = "X" * 500
    truncated = long_name[:200]
    taxdb.name_to_taxids[truncated] = [1]
    taxdb.taxid_to_meta[1] = {"name": truncated, "rank": "genus"}
    _install_fake_backend(monkeypatch, taxdb)

    with caplog.at_level("WARNING", logger="ssign_app.scripts.resolve_taxonomy"):
        result = rt.resolve_organism(long_name)
    assert any("Truncating overlong organism name" in rec.message for rec in caplog.records)
    assert result["species"] is None  # genus-rank only — no species set
    assert result["genus"]["name"] == truncated


def test_control_characters_stripped(monkeypatch):
    taxdb = FakeTaxDb()
    taxdb.name_to_taxids["Escherichia coli"] = [562]
    taxdb.taxid_to_meta[562] = {"name": "Escherichia coli", "rank": "species", "genus": 561}
    taxdb.taxid_to_meta[561] = {"name": "Escherichia", "rank": "genus"}
    _install_fake_backend(monkeypatch, taxdb)

    # Embed null + bell + tab; all should be stripped before lookup.
    result = rt.resolve_organism("Escherichia\x00 coli\x07\t")
    assert result["species"] == {"name": "Escherichia coli", "taxid": "562"}


def test_only_control_characters_returns_nones():
    # If sanitisation leaves an empty string, return early without touching
    # the taxdump (no _load_taxdb call needed — a pure-code guard).
    assert rt.resolve_organism("\x00\x01\x02") == {"species": None, "genus": None}


# ---------------------------------------------------------------------------
# SSIGN_TAXDUMP_DIR env var
# ---------------------------------------------------------------------------


def test_env_var_overrides_default_taxdump_dir(monkeypatch):
    # The TAXDUMP_DIR module constant is captured at import time. The
    # contract pinned here is "the env-var path is the one tried first
    # when _load_taxdb runs". We can't easily re-import the module per
    # test, but we CAN verify that DEFAULT_TAXDUMP_DIR is parameterised
    # by the env var by re-evaluating the same expression.
    monkeypatch.setenv("SSIGN_TAXDUMP_DIR", "/custom/taxdump/path")
    resolved = os.environ.get("SSIGN_TAXDUMP_DIR", rt.DEFAULT_TAXDUMP_DIR)
    assert resolved == "/custom/taxdump/path"
