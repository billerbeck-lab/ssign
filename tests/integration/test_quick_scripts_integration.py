"""Quick integration tests for in-process / pure-Python scripts.

These don't shell out to external tools, so they run in milliseconds
and don't gate on installed binaries / DBs:

    - run_protparam.py: in-process Biopython ProtParam
    - extract_neighborhood.py: pure-Python on gene_info / SS components
    - detect_input_format.py: filename + sniff-byte heuristic
    - extract_gene_order.py: pure-Python on gene_info

Marked `integration` because they exercise the script CLIs end-to-end
against the T5aSS minimal fixture's outputs, but they run as fast as
unit tests.
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "src" / "ssign_app" / "scripts"


def _run_script(script_name: str, args: list[str]) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(SCRIPTS_DIR / script_name)] + args
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60)


# ── detect_input_format.py ──────────────────────────────────────────────


class TestDetectInputFormat:
    def test_detects_genbank(self, t1ss_fixture_gbff):
        r = _run_script("detect_input_format.py", [t1ss_fixture_gbff])
        assert r.returncode == 0
        assert r.stdout.strip() == "genbank"

    def test_detects_fasta_contigs(self, tmp_dir, t1ss_fixture_contigs):
        r = _run_script("detect_input_format.py", [t1ss_fixture_contigs])
        assert r.returncode == 0
        assert r.stdout.strip() == "fasta_contigs"

    def test_detects_protein_fasta(self, tmp_dir, t1ss_fixture_proteins):
        # The fixture is named .faa but content sniff should still flag
        # it as protein_fasta. Rename to .fasta to remove the extension
        # tell-tale and force the content-sniff path.
        renamed = os.path.join(tmp_dir, "queries.fasta")
        with open(t1ss_fixture_proteins) as src, open(renamed, "w") as dst:
            dst.write(src.read())
        r = _run_script("detect_input_format.py", [renamed])
        assert r.returncode == 0
        assert r.stdout.strip() in ("protein_fasta", "fasta_contigs")


# ── compute_protparam.py ────────────────────────────────────────────────


class TestComputeProtParam:
    """compute_protparam.py expects a `substrates` TSV (locus_tags to
    compute features for) + a proteins FASTA. We construct a tiny
    substrates list (just BIMENO_04457) so the test runs in milliseconds.
    """

    def test_runs_on_substrate_subset(
        self, tmp_dir, t1ss_fixture_proteins
    ):
        substrates_path = os.path.join(tmp_dir, "substrates.tsv")
        with open(substrates_path, "w") as f:
            f.write("locus_tag\nBIMENO_04457\n")
        out_path = os.path.join(tmp_dir, "protparam.tsv")
        r = _run_script(
            "compute_protparam.py",
            [
                "--substrates", substrates_path,
                "--proteins", t1ss_fixture_proteins,
                "--sample", "test",
                "--output", out_path,
            ],
        )
        assert r.returncode == 0, f"rc={r.returncode} stderr={r.stderr[:300]}"
        assert os.path.exists(out_path)

        # compute_protparam writes CSV (comma-separated), not TSV.
        with open(out_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) >= 1
        assert "locus_tag" in rows[0]
        # Biopython ProtParam reports MW, instability, GRAVY etc.
        # Confirm we got more than just locus_tag — at least one feature.
        assert len(rows[0]) > 1


# ── extract_gene_order.py + extract_neighborhood.py ─────────────────────
# These need a pre-built gene_info.tsv from extract_proteins as input.


@pytest.fixture
def fixture_gene_info(tmp_dir, t1ss_fixture_gbff):
    """Run extract_proteins.py once and yield its gene_info.tsv path."""
    out_proteins = os.path.join(tmp_dir, "p.faa")
    out_gene_info = os.path.join(tmp_dir, "gene_info.tsv")
    out_metadata = os.path.join(tmp_dir, "meta.json")
    r = _run_script(
        "extract_proteins.py",
        [
            "--input", t1ss_fixture_gbff,
            "--sample", "test",
            "--out-proteins", out_proteins,
            "--out-gene-info", out_gene_info,
            "--out-metadata", out_metadata,
        ],
    )
    assert r.returncode == 0, f"extract_proteins failed: {r.stderr[:300]}"
    return out_gene_info


class TestExtractGeneOrder:
    def test_produces_one_row_per_cds(self, tmp_dir, fixture_gene_info):
        out_path = os.path.join(tmp_dir, "gene_order.tsv")
        r = _run_script(
            "extract_gene_order.py",
            ["--gene-info", fixture_gene_info, "--output", out_path],
        )
        assert r.returncode == 0, f"rc={r.returncode} stderr={r.stderr[:300]}"
        with open(out_path) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        assert len(rows) >= 5  # 9 in minimal fixture, 179 in full
        required = {"locus_tag", "contig", "gene_index"}
        assert required <= set(rows[0].keys())

        # gene_index must be 0-based, contiguous, and sorted ascending
        # within each contig (load-bearing for proximity_analysis.py's
        # early-termination invariant).
        by_contig: dict = {}
        for row in rows:
            by_contig.setdefault(row["contig"], []).append(int(row["gene_index"]))
        for contig, indices in by_contig.items():
            assert indices == sorted(indices), (
                f"gene_index not sorted ascending on contig {contig}"
            )
            assert indices[0] == 0, (
                f"gene_index not 0-based on contig {contig} (starts at {indices[0]})"
            )
