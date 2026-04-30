"""Integration tests for proximity_analysis.py.

The per-component proximity logic is load-bearing — the previous
version (system-boundary proximity, finding proteins within the full
span of a secretion system) produced 26 false positives. This test
suite exercises the real proximity_analysis.py via subprocess against
synthetic-but-realistic input files modelled on the T5aSS fixture
(BIMENO_04457 + neighbors) so we can construct controlled scenarios
and assert specific substrate-or-not outcomes for each.

Tests cover:
    - basic per-component proximity (neighbor within ±window appears)
    - boundary cases (exactly at window, just outside, far away)
    - SS components themselves are excluded as substrates
    - a single neighbor flanking two distinct SS components shows up once
      with both ss_types aggregated
    - DSE cross-genome leakage filter (DSE-only with type not in genome)
    - DSE type-match flagging (dse_type_match column)
    - same-contig invariant (proximity NEVER spans contigs)
    - default ±3 window vs --window override

Each test writes the three input TSVs (gene_order, ss_components,
predictions) into a tempdir, invokes the script, parses the output
TSV, and asserts on the rows.

Run with:
    pytest -m integration tests/integration/test_proximity_analysis_integration.py
"""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = PROJECT_ROOT / "src" / "ssign_app" / "scripts" / "proximity_analysis.py"


# ── helpers ─────────────────────────────────────────────────────────────


def _write_gene_order(path: Path, rows: list[dict]) -> None:
    """Write a gene_order.tsv. Each row needs gene_index + locus_tag + contig.

    Other columns (start, end, strand) are not required by
    proximity_analysis.py but we include them for realism — the file
    matches what extract_gene_order.py produces.
    """
    fieldnames = ["gene_index", "locus_tag", "contig", "start", "end", "strand"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _write_ss_components(path: Path, rows: list[dict]) -> None:
    """Write a ss_components.tsv. Required cols: locus_tag, ss_type, excluded."""
    fieldnames = ["locus_tag", "ss_type", "excluded", "system_id", "component"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "False" if k == "excluded" else "") for k in fieldnames})


def _write_predictions(path: Path, rows: list[dict]) -> None:
    """Write a predictions.tsv with the columns proximity_analysis reads."""
    fieldnames = [
        "locus_tag", "dlp_extracellular_prob", "predicted_localization",
        "dlp_max_localization", "dlp_max_probability",
        "dse_ss_type", "dse_max_prob",
        "signalp_prediction", "signalp_probability", "signalp_cs_position",
        "product",
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def _run(tmp_path: Path, gene_order_rows, ss_rows, pred_rows,
         window: int = 3, conf: float = 0.8) -> list[dict]:
    """Run proximity_analysis.py and return the parsed output rows."""
    go = tmp_path / "gene_order.tsv"
    sc = tmp_path / "ss_components.tsv"
    pr = tmp_path / "predictions.tsv"
    out = tmp_path / "substrates.tsv"

    _write_gene_order(go, gene_order_rows)
    _write_ss_components(sc, ss_rows)
    _write_predictions(pr, pred_rows)

    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--gene-order", str(go),
            "--ss-components", str(sc),
            "--predictions", str(pr),
            "--sample", "test",
            "--window", str(window),
            "--conf-threshold", str(conf),
            "--output", str(out),
        ],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, (
        f"proximity_analysis.py exit {result.returncode}\n"
        f"stderr: {result.stderr}"
    )

    with open(out) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _make_neighborhood(n: int = 11, contig: str = "c1") -> list[dict]:
    """Return n CDS on the same contig, indexed 0..n-1, locus G0..G{n-1}.

    Indexing matches what proximity_analysis expects (0-based gene_index
    on the same contig, sorted by position).
    """
    return [
        {
            "gene_index": str(i),
            "locus_tag": f"G{i:02d}",
            "contig": contig,
            "start": str(i * 1000),
            "end": str(i * 1000 + 800),
            "strand": "+",
        }
        for i in range(n)
    ]


# ── tests ────────────────────────────────────────────────────────────────


class TestPerComponentProximity:
    """The headline guarantee: a neighbor within ±window of an SS component
    on the same contig gets flagged; a neighbor outside the window does not.
    """

    def test_neighbor_within_window_is_found(self, tmp_dir):
        tmp = Path(tmp_dir)
        # G05 is the SS component. G03 and G07 are within ±2 (window=2).
        # Both have DLP extracellular hit. Both should appear.
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [
            {"locus_tag": "G03", "dlp_extracellular_prob": "0.95"},
            {"locus_tag": "G07", "dlp_extracellular_prob": "0.95"},
        ]
        rows = _run(tmp, genes, ss, preds, window=2)
        loci = {r["locus_tag"] for r in rows}
        assert loci == {"G03", "G07"}

    def test_neighbor_just_outside_window_is_excluded(self, tmp_dir):
        tmp = Path(tmp_dir)
        # G05 is the component. G02 (offset -3) is just outside ±2 window.
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [{"locus_tag": "G02", "dlp_extracellular_prob": "0.95"}]
        rows = _run(tmp, genes, ss, preds, window=2)
        assert rows == []

    def test_neighbor_far_away_is_excluded(self, tmp_dir):
        """The whole point of per-component proximity: distant proteins
        are not picked up just because something else far away is also
        an SS component."""
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(20)
        # Component at G05; far-away protein G18 has DLP hit.
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [{"locus_tag": "G18", "dlp_extracellular_prob": "0.95"}]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert rows == []

    def test_per_component_not_per_system_span(self, tmp_dir):
        """REGRESSION: the old span-based proximity counted everything
        between MIN(component_idx) and MAX(component_idx) as a substrate
        candidate. With two components at G02 and G18 and a target at
        G10 (10 genes away from both), per-component proximity correctly
        excludes G10. Span-based would (incorrectly) include it.
        """
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(20)
        ss = [
            {"locus_tag": "G02", "ss_type": "T5aSS", "excluded": "False"},
            {"locus_tag": "G18", "ss_type": "T5aSS", "excluded": "False"},
        ]
        # G10 is 8 away from G02 and 8 away from G18. Outside ±3 window
        # of each. Span-based would include it (between the two
        # components); per-component correctly excludes.
        preds = [{"locus_tag": "G10", "dlp_extracellular_prob": "0.95"}]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert rows == [], (
            "G10 is 8 genes from each component — per-component proximity "
            "MUST exclude it. If this test fails, the fix that removed 26 "
            "false positives has regressed."
        )


class TestSsComponentSelfExclusion:
    """SS components themselves are not substrates (they ARE the system).

    Note: T5SS substrates are added separately by t5ss_handler.py, not
    by proximity_analysis.py. proximity_analysis only finds NEIGHBORS.
    """

    def test_component_itself_is_not_in_output(self, tmp_dir):
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        # Component G05 itself has DLP+ — it should NOT be in proximity
        # output (proximity_analysis only finds neighbors).
        preds = [
            {"locus_tag": "G05", "dlp_extracellular_prob": "0.95"},
            {"locus_tag": "G06", "dlp_extracellular_prob": "0.95"},
        ]
        rows = _run(tmp, genes, ss, preds, window=2)
        loci = {r["locus_tag"] for r in rows}
        assert loci == {"G06"}, (
            "G05 is the SS component itself — proximity_analysis must NOT "
            "list components as substrates. T5SS-substrate-is-component "
            "logic lives in t5ss_handler.py."
        )

    def test_excluded_component_is_ignored(self, tmp_dir):
        """Components marked excluded=True (e.g. flagellum, T3SS) should
        not pull in substrates from their neighborhood."""
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [
            {"locus_tag": "G05", "ss_type": "Flagellum", "excluded": "True"},
            {"locus_tag": "G02", "ss_type": "T5aSS", "excluded": "False"},
        ]
        # G06 is in the flagellum component's window (±3) but NOT in T5aSS's.
        # Should be excluded because the flagellum is excluded.
        preds = [{"locus_tag": "G06", "dlp_extracellular_prob": "0.95"}]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert rows == []


class TestNeighborSharedAcrossComponents:
    def test_one_neighbor_two_components_aggregates_ss_types(self, tmp_dir):
        """If the same neighbor is in proximity to two components of
        different SS types, the substrate row should include BOTH types
        in nearby_ss_types (comma-joined)."""
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        # G03 is 1 away from G02 (T5aSS) and 2 away from G05 (T2SS).
        ss = [
            {"locus_tag": "G02", "ss_type": "T5aSS", "excluded": "False"},
            {"locus_tag": "G05", "ss_type": "T2SS", "excluded": "False"},
        ]
        preds = [{"locus_tag": "G03", "dlp_extracellular_prob": "0.95"}]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert len(rows) == 1
        nearby_types = set(rows[0]["nearby_ss_types"].split(","))
        assert nearby_types == {"T5aSS", "T2SS"}, (
            f"Expected both SS types aggregated, got {nearby_types}"
        )


class TestSameContigInvariant:
    """Proximity NEVER spans contigs — gene_index is per-contig."""

    def test_neighbor_on_different_contig_excluded(self, tmp_dir):
        tmp = Path(tmp_dir)
        # G00-G04 on contig c1, X00-X04 on contig c2. G02 is the SS
        # component on c1. X00 has gene_index=0 and DLP+ — same numerical
        # offset (-2) from G02 but on c2. Must be excluded.
        genes = (
            _make_neighborhood(5, contig="c1")
            + [
                {"gene_index": str(i), "locus_tag": f"X{i:02d}",
                 "contig": "c2", "start": str(i * 1000),
                 "end": str(i * 1000 + 800), "strand": "+"}
                for i in range(5)
            ]
        )
        ss = [{"locus_tag": "G02", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [
            {"locus_tag": "X00", "dlp_extracellular_prob": "0.95"},
            {"locus_tag": "G04", "dlp_extracellular_prob": "0.95"},
        ]
        rows = _run(tmp, genes, ss, preds, window=3)
        loci = {r["locus_tag"] for r in rows}
        # G04 is +2 from G02 on c1 — included.
        # X00 is on a different contig — excluded even though the gene
        # index difference would otherwise be in the window.
        assert loci == {"G04"}


class TestDseCrossGenomeLeakageFilter:
    """DSE-only substrates with a SS type that doesn't exist in the
    genome's MacSyFinder calls get filtered out (Session 11 fix)."""

    def test_dse_only_with_type_in_genome_kept(self, tmp_dir):
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T2SS", "excluded": "False"}]
        # G06: DSE T2SE-positive (matches genome's T2SS), no DLP
        preds = [{
            "locus_tag": "G06",
            "dlp_extracellular_prob": "0.0",
            "dse_ss_type": "T2SS",
            "dse_max_prob": "0.95",
        }]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert {r["locus_tag"] for r in rows} == {"G06"}

    def test_dse_only_with_type_not_in_genome_filtered(self, tmp_dir):
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        # Genome has T5aSS only; DSE on neighbor predicts T6SE.
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [{
            "locus_tag": "G06",
            "dlp_extracellular_prob": "0.0",
            "dse_ss_type": "T6SS",
            "dse_max_prob": "0.95",
        }]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert rows == [], (
            "DSE T6SE on a neighbor of T5aSS-only genome should be "
            "filtered (DSE cross-genome leakage fix per CLAUDE.md "
            "Critical Bug Fix #2)."
        )

    def test_dlp_plus_dse_mismatch_kept_because_dlp_carries(self, tmp_dir):
        """If DLP is positive AND DSE is type-mismatched, the substrate
        is still kept — DLP alone is enough. Only DSE-only mismatches
        get filtered."""
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [{
            "locus_tag": "G06",
            "dlp_extracellular_prob": "0.95",
            "dse_ss_type": "T6SS",
            "dse_max_prob": "0.95",
        }]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert {r["locus_tag"] for r in rows} == {"G06"}


class TestDseTypeMatchFlag:
    """The dse_type_match column flags whether DSE's predicted SS type
    matches the nearby MacSyFinder system type. Useful downstream for
    ranking substrate confidence."""

    def test_dse_type_match_true_when_aligned(self, tmp_dir):
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T2SS", "excluded": "False"}]
        preds = [{
            "locus_tag": "G06",
            "dlp_extracellular_prob": "0.95",
            "dse_ss_type": "T2SS",
            "dse_max_prob": "0.95",
        }]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert len(rows) == 1
        assert rows[0]["dse_type_match"] == "True"

    def test_dse_type_match_false_when_misaligned(self, tmp_dir):
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T2SS", "excluded": "False"}]
        # DLP rescues the row; DSE predicts T6SE which doesn't match T2SS.
        preds = [{
            "locus_tag": "G06",
            "dlp_extracellular_prob": "0.95",
            "dse_ss_type": "T6SS",
            "dse_max_prob": "0.95",
        }]
        rows = _run(tmp, genes, ss, preds, window=3)
        assert len(rows) == 1
        assert rows[0]["dse_type_match"] == "False"


class TestOutputSchema:
    """Confirm the output TSV has every column downstream consumers expect."""

    def test_output_columns_present(self, tmp_dir):
        tmp = Path(tmp_dir)
        genes = _make_neighborhood(11)
        ss = [{"locus_tag": "G05", "ss_type": "T5aSS", "excluded": "False"}]
        preds = [{"locus_tag": "G06", "dlp_extracellular_prob": "0.95"}]
        rows = _run(tmp, genes, ss, preds)
        assert len(rows) == 1

        required = {
            "locus_tag", "sample_id", "tool", "nearby_ss_types",
            "dlp_extracellular_prob", "predicted_localization",
            "dlp_max_localization", "dlp_max_probability",
            "dse_ss_type", "dse_max_prob",
            "signalp_prediction", "signalp_probability", "signalp_cs_position",
            "dse_type_match", "product",
        }
        assert required <= set(rows[0].keys())
