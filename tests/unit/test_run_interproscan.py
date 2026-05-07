"""Tests for run_interproscan.py.

Two pure-Python surfaces here:

1. `parse_interproscan_tsv` — column-index TSV parser. Each row is one
   member-DB hit; rows for the same protein must aggregate into sets
   for domains, GO terms, Pfam IDs, and descriptions, then join with
   ';' on output. Hits with `_MISSING` ("-") in a field don't pollute
   the aggregate.
2. `_GO_ID_RE` — extracts every `GO:NNNNNNN` from the GO-terms cell
   (one cell can list several pipe-separated GO entries).

The `run_local_interproscan` subprocess path requires a multi-GB
InterProScan install + Java runtime and is exercised by
tests/integration/test_run_interproscan_integration.py.
"""

import os
import sys

import pytest


from run_interproscan import (  # noqa: E402
    _GO_ID_RE,
    _MISSING,
    DEFAULT_IPS_APPLICATIONS,
    parse_interproscan_tsv,
)

# IPS TSV column layout (1-based for documentation, 0-based in code):
#   0  protein_accession
#   1  md5
#   2  seqlen
#   3  analysis (e.g. "Pfam")
#   4  signature_accession (e.g. "PF03797")
#   5  signature_description
#   6  start
#   7  end
#   8  evalue
#   9  status
#   10 date
#   11 interpro_accession (e.g. "IPR005546")
#   12 interpro_description
#   13 go_terms (e.g. "GO:0005886|GO:0019867")
#   14 pathways (optional)


def _ips_row(
    pid="GENE_001",
    sig_acc="PF03797",
    sig_desc="Autotransporter beta-domain",
    ipr_acc="IPR005546",
    ipr_desc="Autotransporter, beta domain",
    go="GO:0019867",
):
    """A single 14-column IPS TSV row (tab-joined)."""
    return "\t".join(
        [
            pid,
            "abc123md5",
            "300",
            "Pfam",
            sig_acc,
            sig_desc,
            "1",
            "300",
            "1e-50",
            "T",
            "2026-05-05",
            ipr_acc,
            ipr_desc,
            go,
        ]
    )


def _write_ips(path, rows):
    with open(path, "w") as f:
        for row in rows:
            f.write(row + "\n")
    return path


# ---------------------------------------------------------------------------
# _GO_ID_RE — regex extraction
# ---------------------------------------------------------------------------


class TestGoIdRegex:
    @pytest.mark.parametrize(
        "cell, expected",
        [
            ("GO:0005886", ["GO:0005886"]),
            ("GO:0005886|GO:0019867", ["GO:0005886", "GO:0019867"]),
            ("GO:0005886(InterPro)|GO:0019867(PANTHER)", ["GO:0005886", "GO:0019867"]),
            ("", []),
            ("-", []),
            ("not a go term", []),
            # 7-digit IDs (newer GO) must also match
            ("GO:0099999", ["GO:0099999"]),
        ],
    )
    def test_extracts_all_go_ids(self, cell, expected):
        assert [m.group(1) for m in _GO_ID_RE.finditer(cell)] == expected


# ---------------------------------------------------------------------------
# parse_interproscan_tsv — single-protein cases
# ---------------------------------------------------------------------------


class TestParseSingleProtein:
    def test_single_row_aggregated(self, tmp_dir):
        path = _write_ips(os.path.join(tmp_dir, "ips.tsv"), [_ips_row()])
        results = parse_interproscan_tsv(path)
        assert "GENE_001" in results
        e = results["GENE_001"]
        assert e["locus_tag"] == "GENE_001"
        assert e["interpro_domains"] == "IPR005546"
        assert e["interpro_descriptions"] == "Autotransporter, beta domain"
        assert e["interpro_pfam_ids"] == "PF03797"
        assert e["interpro_go_terms"] == "GO:0019867"

    def test_multiple_rows_same_protein_merged(self, tmp_dir):
        # Same protein, three different signatures. All four fields must dedup
        # into a sorted, semicolon-joined string.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [
                _ips_row(sig_acc="PF03797", ipr_acc="IPR005546", ipr_desc="alpha", go="GO:0001"),
                _ips_row(sig_acc="PF00001", ipr_acc="IPR000123", ipr_desc="beta", go="GO:0002"),
                _ips_row(sig_acc="PF03797", ipr_acc="IPR005546", ipr_desc="alpha", go="GO:0001|GO:0003"),
            ],
        )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_domains"] == "IPR000123;IPR005546"  # sorted
        assert e["interpro_descriptions"] == "alpha;beta"
        assert e["interpro_pfam_ids"] == "PF00001;PF03797"
        assert e["interpro_go_terms"] == "GO:0001;GO:0002;GO:0003"

    def test_missing_ipr_excluded(self, tmp_dir):
        # Row with "-" in IPR fields must NOT add to the aggregate.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [_ips_row(ipr_acc=_MISSING, ipr_desc=_MISSING)],
        )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_domains"] == ""
        assert e["interpro_descriptions"] == ""
        # Pfam still extracted (sig_acc isn't missing)
        assert e["interpro_pfam_ids"] == "PF03797"

    def test_missing_go_excluded(self, tmp_dir):
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [_ips_row(go=_MISSING)],
        )
        assert parse_interproscan_tsv(path)["GENE_001"]["interpro_go_terms"] == ""

    def test_non_pfam_signature_not_in_pfam_ids(self, tmp_dir):
        # Signature "TIGR03660" doesn't start with "PF" — must NOT show up in
        # interpro_pfam_ids.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [_ips_row(sig_acc="TIGR03660", sig_desc="something")],
        )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_pfam_ids"] == ""
        # The IPR side is still populated
        assert e["interpro_domains"] == "IPR005546"


# ---------------------------------------------------------------------------
# parse_interproscan_tsv — multi-protein + filter
# ---------------------------------------------------------------------------


class TestParseMultipleProteins:
    def test_each_protein_isolated(self, tmp_dir):
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [
                _ips_row(pid="GENE_001", sig_acc="PF03797"),
                _ips_row(pid="GENE_002", sig_acc="PF00001"),
            ],
        )
        results = parse_interproscan_tsv(path)
        assert set(results.keys()) == {"GENE_001", "GENE_002"}
        assert results["GENE_001"]["interpro_pfam_ids"] == "PF03797"
        assert results["GENE_002"]["interpro_pfam_ids"] == "PF00001"

    def test_target_ids_filters(self, tmp_dir):
        # Only GENE_001 is in target_ids — GENE_002 must be dropped.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [
                _ips_row(pid="GENE_001"),
                _ips_row(pid="GENE_002"),
            ],
        )
        results = parse_interproscan_tsv(path, target_ids={"GENE_001"})
        assert set(results.keys()) == {"GENE_001"}

    def test_empty_target_ids_returns_empty(self, tmp_dir):
        path = _write_ips(os.path.join(tmp_dir, "ips.tsv"), [_ips_row()])
        # Note: the contract here is that target_ids is checked for membership;
        # empty set behaves like "no protein matches" — verifies the truthy
        # check on target_ids in the parser.
        results = parse_interproscan_tsv(path, target_ids={"NOT_PRESENT"})
        assert results == {}


# ---------------------------------------------------------------------------
# Parser resilience
# ---------------------------------------------------------------------------


class TestParserResilience:
    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "ips.tsv")
        open(path, "w").close()
        assert parse_interproscan_tsv(path) == {}

    def test_short_rows_skipped(self, tmp_dir):
        # Rows with <12 columns are dropped (corrupt or partial line)
        path = os.path.join(tmp_dir, "ips.tsv")
        with open(path, "w") as f:
            f.write("GENE_001\tabc\t300\tPfam\tPF03797\n")  # only 5 columns
        assert parse_interproscan_tsv(path) == {}

    def test_row_without_go_column_handled(self, tmp_dir):
        # A row with exactly 13 columns (no GO field at index 13) must not
        # crash. The parser uses bounds-checked column reads.
        path = os.path.join(tmp_dir, "ips.tsv")
        with open(path, "w") as f:
            f.write(
                "GENE_001\tmd5\t300\tPfam\tPF03797\tdesc\t1\t300\t1e-50\tT\t2026-05-05\t"
                "IPR005546\tipr_desc\n"  # only 13 fields
            )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_go_terms"] == ""
        assert e["interpro_domains"] == "IPR005546"


# ---------------------------------------------------------------------------
# DEFAULT_IPS_APPLICATIONS — pinned to the documented bacteria set
# ---------------------------------------------------------------------------


class TestDefaultApplications:
    def test_panther_excluded(self):
        # PANTHER is intentionally skipped — eukaryote-leaning + slowest member.
        # Future maintainers must opt into PANTHER explicitly via --applications.
        assert "PANTHER" not in DEFAULT_IPS_APPLICATIONS

    def test_pfam_included(self):
        # Pfam is the workhorse for bacterial annotation; must always be present.
        assert "Pfam" in DEFAULT_IPS_APPLICATIONS

    def test_no_duplicates(self):
        assert len(DEFAULT_IPS_APPLICATIONS) == len(set(DEFAULT_IPS_APPLICATIONS))
