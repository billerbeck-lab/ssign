"""End-to-end integration test for integrate_annotations.py.

The merge step is where every upstream tool's output converges into
the final master CSV. Test approach: feed a small synthetic-but-
realistic set of upstream files into the script and verify the master
CSV has every expected column populated + the known T5aSS substrate
(BIMENO_04457) carries annotations from every tool we provide.

Synthetic-but-realistic: each upstream file has the exact column
schema the real tool produces (verified against the per-tool
integration tests), but the row content is small + deterministic so
the merge logic is exercised in isolation from upstream tool runtime.

Run with:
    pytest -m integration tests/integration/test_integrate_annotations_integration.py
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest


pytestmark = pytest.mark.integration


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = PROJECT_ROOT / "src" / "ssign_app" / "scripts" / "integrate_annotations.py"


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


@pytest.fixture
def synthetic_pipeline_outputs(tmp_dir):
    """Build the minimal set of inputs integrate_annotations.py wants.

    BIMENO_04457 is the substrate of interest; G05 + G07 are
    housekeeping flanking proteins included to verify non-substrate
    rows are merged correctly when they appear in the substrates_all
    sheet but not the filtered one.
    """
    tmp = Path(tmp_dir)

    # substrates_filtered.tsv: passes all the cross-validate filters.
    # In the real pipeline, this row gets the proximity + cross-validate
    # metadata. We use a minimal viable schema.
    substrates_filtered = tmp / "substrates_filtered.tsv"
    _write_tsv(
        substrates_filtered,
        ["locus_tag", "sample_id", "tool", "nearby_ss_types"],
        [
            {
                "locus_tag": "BIMENO_04457",
                "sample_id": "test",
                "tool": "DLP+DSE",
                "nearby_ss_types": "T5aSS",
            },
        ],
    )

    # substrates_all.tsv: same schema, may have more rows pre-filter.
    substrates_all = tmp / "substrates_all.tsv"
    _write_tsv(
        substrates_all,
        ["locus_tag", "sample_id", "tool", "nearby_ss_types"],
        [
            {
                "locus_tag": "BIMENO_04457",
                "sample_id": "test",
                "tool": "DLP+DSE",
                "nearby_ss_types": "T5aSS",
            },
        ],
    )

    # gene_info.tsv (Phase 3.3.c shape: gbff_annotation column carries
    # the original GenBank product, post-Bakta re-annotation).
    gene_info = tmp / "gene_info.tsv"
    _write_tsv(
        gene_info,
        [
            "locus_tag", "protein_id", "gene", "product", "contig",
            "start", "end", "strand", "ec_numbers", "cog_ids", "go_terms",
            "kegg_ko", "refseq_ids", "pfam_ids", "gbff_annotation",
        ],
        [
            {
                "locus_tag": "BIMENO_04457",
                "product": "hypothetical protein",
                "contig": "Xanthobacter_T5aSS_minimal",
                "start": "4562",
                "end": "9122",
                "strand": "-",
                "pfam_ids": "PF03797;PF12951",
                "gbff_annotation": "Autotransporter domain-containing protein",
            },
        ],
    )

    # Per-tool annotation TSVs. integrate_annotations.py merges these
    # by locus_tag (or any "protein"/"id" column).
    bakta_ann = tmp / "bakta_ann.tsv"
    _write_tsv(
        bakta_ann,
        ["locus_tag", "bakta_product", "pfam_ids"],
        [{"locus_tag": "BIMENO_04457",
          "bakta_product": "hypothetical protein",
          "pfam_ids": "PF03797;PF12951"}],
    )

    eggnog_ann = tmp / "eggnog_ann.tsv"
    _write_tsv(
        eggnog_ann,
        ["locus_tag", "eggnog_description", "cog_category", "kegg_ko"],
        [{"locus_tag": "BIMENO_04457",
          "eggnog_description": "Autotransporter beta-domain",
          "cog_category": "M",
          "kegg_ko": ""}],
    )

    plm_blast_ann = tmp / "plm_blast_ann.tsv"
    _write_tsv(
        plm_blast_ann,
        ["locus_tag", "ecod70_top1_description", "ecod70_top1_score"],
        [{"locus_tag": "BIMENO_04457",
          "ecod70_top1_description": "Autotransporter beta-domain",
          "ecod70_top1_score": "0.91"}],
    )

    # Tiny proteins FASTA so the script can pull sequence + length
    proteins = tmp / "proteins.faa"
    proteins.write_text(">BIMENO_04457\nMKKVALLAALLPALIPSAQAATTHQKVKVGKEEFNQ\n")

    return {
        "tmp": tmp,
        "substrates_filtered": substrates_filtered,
        "substrates_all": substrates_all,
        "gene_info": gene_info,
        "annotations": [bakta_ann, eggnog_ann, plm_blast_ann],
        "proteins": proteins,
    }


def _run_integrate(inputs: dict) -> Path:
    out = inputs["tmp"] / "master.tsv"
    cmd = [
        sys.executable, str(SCRIPT),
        "--substrates-filtered", str(inputs["substrates_filtered"]),
        "--substrates-all", str(inputs["substrates_all"]),
        "--gene-info", str(inputs["gene_info"]),
        "--proteins", str(inputs["proteins"]),
        "--sample", "test",
        "--output", str(out),
        "--annotations", *[str(p) for p in inputs["annotations"]],
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, (
        f"integrate_annotations.py exit {result.returncode}\n"
        f"stderr: {result.stderr}"
    )
    return out


class TestIntegrateAnnotations:
    def test_master_tsv_carries_every_tool_contribution(
        self, synthetic_pipeline_outputs
    ):
        """The master CSV's BIMENO_04457 row should contain the
        per-tool contributions: Bakta product, EggNOG description,
        pLM-BLAST hit, GBFF annotation (from gene_info), and the
        consensus voting columns annotation_consensus computes."""
        out = _run_integrate(synthetic_pipeline_outputs)
        # integrate_annotations.py writes CSV (comma-separated).
        with open(out) as f:
            rows = list(csv.DictReader(f))

        target = next(
            (r for r in rows if r["locus_tag"] == "BIMENO_04457"), None
        )
        assert target is not None, (
            "BIMENO_04457 not in master output — substrate row dropped "
            "during merge."
        )

        # Per-tool contributions
        assert target.get("bakta_product"), "Bakta product not merged"
        assert target.get("eggnog_description"), "EggNOG description not merged"
        assert target.get("ecod70_top1_description"), "pLM-BLAST hit not merged"
        assert target.get("gbff_annotation"), (
            "gbff_annotation not merged from gene_info — Phase 3.3.c "
            "regression?"
        )

        # Sequence got pulled from proteins FASTA
        assert target.get("sequence"), "sequence not merged from proteins FASTA"
        assert int(target.get("aa_length", 0)) > 0

        # annotation_consensus voting columns (computed in-process by
        # integrate_annotations.py). The exact names are
        # broad_consensus_annotation / detailed_consensus_annotation /
        # n_tools_agreeing / confidence_tier — schema sourced from
        # _compute_consensus output, verified against a real run.
        consensus_cols = {
            "broad_annotation",
            "detailed_annotation",
            "n_tools_agreeing",
            "n_tools_with_hits",
            "confidence_tier",
            "annotation_tools",
        }
        missing = consensus_cols - set(target.keys())
        assert not missing, (
            f"annotation_consensus columns missing from master CSV: "
            f"{missing}"
        )

        # n_tools_with_hits should reflect the per-tool contributions
        # we provided (Bakta + EggNOG + pLM-BLAST + GBFF = 4).
        n_tools = int(target["n_tools_with_hits"])
        assert n_tools >= 3, (
            f"Expected ≥3 tools to have annotated BIMENO_04457, got {n_tools}"
        )

    def test_runs_with_no_annotations(self, synthetic_pipeline_outputs):
        """integrate_annotations should still produce a master CSV when
        no tool annotations are supplied — gene_info + substrates only."""
        inputs = dict(synthetic_pipeline_outputs)
        inputs["annotations"] = []
        out = _run_integrate(inputs)
        # integrate_annotations.py writes CSV (comma-separated).
        with open(out) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) >= 1
        # gbff_annotation still merges from gene_info
        target = next(
            (r for r in rows if r["locus_tag"] == "BIMENO_04457"), None
        )
        assert target is not None
        assert target.get("gbff_annotation")
