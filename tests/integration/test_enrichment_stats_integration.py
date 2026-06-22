"""End-to-end integration test for the --enrichment-stats path.

Drives the per-genome circular-shift test (enrichment_testing.py) and the
cross-genome pooler (pool_enrichment_stats) as subprocesses / calls against two
synthetic two-contig genomes. Predictions (DLP/DSE) are hand-synthesised because
the stats engine just consumes their TSV outputs -- the actual tool binaries
aren't part of what's being tested here.

Mirrors the test_quick_scripts_integration.py style: small fixtures, no external
binaries / databases, runs in well under a second.
"""

import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.integration


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "src" / "ssign_app" / "scripts"


def _run_script(script_name: str, args: list[str]) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(SCRIPTS_DIR / script_name)] + args
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60)


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict]):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _build_genome(root: Path, prefix: str, positives_neigh: set):
    """Build gene_order + ss_components + whole-genome DLP/DSE for one genome.

    Layout: 15 proteins, contig_A genes 0-9 + contig_B genes 0-4. Two T2SS
    components on contig_A at gene_index 5, 6. Whole-genome predictions (every
    protein has a row) as enrichment now forces. Returns dict of file paths.
    """
    root.mkdir(parents=True, exist_ok=True)
    loci = [f"{prefix}_GA_{i:04d}" for i in range(10)] + [f"{prefix}_GB_{i:04d}" for i in range(5)]

    gene_order = root / "gene_order.tsv"
    rows = [
        {"contig": f"{prefix}_contig_A", "gene_index": str(i), "locus_tag": f"{prefix}_GA_{i:04d}"} for i in range(10)
    ]
    rows += [
        {"contig": f"{prefix}_contig_B", "gene_index": str(i), "locus_tag": f"{prefix}_GB_{i:04d}"} for i in range(5)
    ]
    _write_tsv(gene_order, ["contig", "gene_index", "locus_tag"], rows)

    ss_components = root / "ss_components.tsv"
    _write_tsv(
        ss_components,
        ["locus_tag", "ss_type", "sys_id", "excluded"],
        [
            {"locus_tag": f"{prefix}_GA_0005", "ss_type": "T2SS", "sys_id": f"{prefix}_T2SS_1", "excluded": "False"},
            {"locus_tag": f"{prefix}_GA_0006", "ss_type": "T2SS", "sys_id": f"{prefix}_T2SS_1", "excluded": "False"},
        ],
    )

    dlp_rows, dse_rows = [], []
    for sid in loci:
        pos = sid in positives_neigh
        dlp_rows.append({"locus_tag": sid, "dlp_extracellular_prob": "0.95" if pos else "0.10"})
        dse_rows.append(
            {
                "locus_tag": sid,
                "dse_ss_type": "T2SS" if pos else "Non-secreted",
                "dse_max_prob": "0.95" if pos else "0.10",
            }
        )
    dlp_tsv, dse_tsv = root / "dlp.tsv", root / "dse.tsv"
    _write_tsv(dlp_tsv, ["locus_tag", "dlp_extracellular_prob"], dlp_rows)
    _write_tsv(dse_tsv, ["locus_tag", "dse_ss_type", "dse_max_prob"], dse_rows)

    return {"gene_order": gene_order, "ss_components": ss_components, "dlp": dlp_tsv, "dse": dse_tsv}


class TestStatsEndToEnd:
    def test_two_genome_run_emits_per_genome_and_pooled(self, tmp_path):
        # Both genomes: GA_0003 + GA_0004 are positive, inside the +/-3 window of
        # the T2SS components at indices 5, 6 -> observed = 2 per genome per tool.
        ga = _build_genome(tmp_path / "ga", "A", {"A_GA_0003", "A_GA_0004"})
        gb = _build_genome(tmp_path / "gb", "B", {"B_GA_0003", "B_GA_0004"})

        # ── per-genome circular-shift test (no null sampling needed) ──
        per_genome_tsvs = []
        for prefix, g in (("A", ga), ("B", gb)):
            out = tmp_path / f"{prefix}_enrichment_stats.tsv"
            nulls = tmp_path / f"{prefix}_enrichment_nulls.npz"
            r = _run_script(
                "enrichment_testing.py",
                [
                    "--ss-components", str(g["ss_components"]),
                    "--gene-order", str(g["gene_order"]),
                    "--dlp", str(g["dlp"]),
                    "--dse", str(g["dse"]),
                    "--window", "3",
                    "--conf-threshold", "0.8",
                    "--sample", f"genome_{prefix}",
                    "--out", str(out),
                    "--nulls-out", str(nulls),
                ],
            )  # fmt: skip
            assert r.returncode == 0, r.stderr
            assert nulls.exists()
            per_genome_tsvs.append(out)

            rows = list(csv.DictReader(open(out), delimiter="\t"))
            assert len(rows) == 2  # one T2SS type x 2 tools
            assert {row["tool"] for row in rows} == {"DLP", "DSE"}
            assert {row["ss_type"] for row in rows} == {"T2SS"}
            assert {row["mode"] for row in rows} == {"window"}
            assert all(int(row["observed"]) == 2 for row in rows)  # GA_0003 + GA_0004
            # npz carries a null array per (type, tool)
            with np.load(nulls) as npz:
                assert {"T2SS__DLP", "T2SS__DSE"} <= set(npz.files)

        # ── pool across the two genomes ──
        sys.path.insert(0, str(PROJECT_ROOT / "src"))
        from ssign_app.core.runner import pool_enrichment_stats

        pooled_path = tmp_path / "pooled_enrichment_stats.tsv"
        pooled_nulls = tmp_path / "pooled_enrichment_nulls.npz"
        n_pooled = pool_enrichment_stats(
            [str(p) for p in per_genome_tsvs], str(pooled_path), nulls_output=str(pooled_nulls)
        )
        assert n_pooled == 2  # (T2SS, DLP) and (T2SS, DSE)
        assert pooled_nulls.exists()

        pooled = list(csv.DictReader(open(pooled_path), delimiter="\t"))
        assert {r["ss_type"] for r in pooled} == {"T2SS"}
        for r in pooled:
            assert int(r["observed"]) == 4  # 2 + 2 across the two genomes
            assert r["mode"] == "window"
            assert float(r["fold"]) > 1.0  # both genomes planted enrichment

        # ── combined pooled figure renders from the pooled TSV + nulls ──
        pooled_png = tmp_path / "pooled_enrichment_null_distributions.png"
        r = _run_script(
            "run_enrichment_figure.py",
            ["--stats", str(pooled_path), "--nulls", str(pooled_nulls), "--out", str(pooled_png)],
        )
        assert r.returncode == 0, r.stderr
        assert pooled_png.exists() and pooled_png.stat().st_size > 0

        # ── shared helper (used by the CLI multi-genome path + GUI) does pool+figure ──
        from ssign_app.core.runner import pool_and_plot_enrichment

        helper_dir = tmp_path / "helper_out"
        helper_dir.mkdir()
        # per-genome stats+npz must be co-located for the helper; reuse the genome files
        import shutil as _shutil

        for prefix in ("A", "B"):
            _shutil.copy(tmp_path / f"{prefix}_enrichment_stats.tsv", helper_dir)
            _shutil.copy(tmp_path / f"{prefix}_enrichment_nulls.npz", helper_dir)
        n = pool_and_plot_enrichment(
            [str(helper_dir / f"{p}_enrichment_stats.tsv") for p in ("A", "B")], str(helper_dir)
        )
        assert n == 2
        assert (helper_dir / "pooled_enrichment_stats.tsv").exists()
        assert (helper_dir / "figures" / "pooled_enrichment_null_distributions.png").exists()
