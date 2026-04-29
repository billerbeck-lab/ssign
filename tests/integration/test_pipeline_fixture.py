"""Whole-pipeline integration test on the T1SS fixture genome.

Runs the full `PipelineRunner` end-to-end against the bundled
Xanthobacter tagetidis contig_87 fixture (213 kb, 179 CDS, one known
T1SS system with BIMENO_04457 as the expected substrate) and asserts
the high-level invariants ssign's downstream consumers care about:

    - Output directory is populated
    - Input parsing produced a protein FASTA
    - MacSyFinder detected at least one secretion system
    - BIMENO_04457 appears somewhere in the downstream outputs

Prerequisites (auto-skip when missing):
    - MacSyFinder + TXSScan installed (pip install macsyfinder + pyhmmer
      shim + `macsydata install TXSScan`)
    - DeepLocPro reachable — either via the DTU webface (default,
      requires network) or a local install (Phase 2.3, blocked on DTU
      redistribution approval)
    - Explicit opt-in: set SSIGN_RUN_FULL_PIPELINE=1. The test is
      expensive enough (~2-5 min) to skip by default even when all
      prerequisites are installed.

Run explicitly with:
    SSIGN_RUN_FULL_PIPELINE=1 pytest -m integration \\
        tests/integration/test_pipeline_fixture.py
"""

import importlib
import os
import shutil
import sys

import pytest


pytestmark = pytest.mark.integration


SRC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from conftest import dlp_pipeline_kwargs  # noqa: E402


def _skip_unless_pipeline_prereqs_ready():
    """Skip the current test unless every whole-pipeline prerequisite is met."""
    if os.environ.get("SSIGN_RUN_FULL_PIPELINE") != "1":
        pytest.skip(
            "SSIGN_RUN_FULL_PIPELINE=1 not set — whole-pipeline test is "
            "opt-in because it takes several minutes; by default it hits "
            "the DTU webface for DeepLocPro, or runs DLP locally if "
            "SSIGN_DEEPLOCPRO_PATH is set."
        )

    # MacSyFinder 2.x renamed the Python package: `macsypy` (legacy) →
    # `macsylib` (current). Accept either to keep this test working
    # across MacSyFinder versions ssign deps support.
    if importlib.util.find_spec("macsylib") is None and importlib.util.find_spec("macsypy") is None:
        pytest.skip("MacSyFinder (macsylib/macsypy) not installed")

    if shutil.which("macsyfinder") is None:
        pytest.skip("macsyfinder binary not on PATH")


def _make_t1ss_config(tmp_dir, t1ss_fixture_gbff):
    """Shared PipelineConfig for the two test methods.

    use_input_annotations=True skips Bakta re-annotation (Phase 3.3.c
    default) — keeps the test fast and self-contained.
    """
    from ssign_app.core.runner import PipelineConfig
    return PipelineConfig(
        input_path=t1ss_fixture_gbff,
        sample_id="t1ss_fixture",
        outdir=tmp_dir,
        use_input_annotations=True,
        wholeness_threshold=0.8,
        excluded_systems=["Flagellum", "Tad", "T3SS"],
        conf_threshold=0.8,
        proximity_window=3,
        # Keep the run cheap — skip every heavy annotation tool
        **dlp_pipeline_kwargs(),
        skip_deepsece=True,
        skip_signalp=True,
        skip_blastp=True,
        skip_hhsuite=True,
        skip_interproscan=True,
        skip_plmblast=True,
        skip_protparam=False,
        skip_structure=True,
    )


class TestPipelineOnT1SSFixture:
    def test_pipeline_runs_and_detects_t1ss(self, tmp_dir, t1ss_fixture_gbff):
        _skip_unless_pipeline_prereqs_ready()

        from ssign_app.core.runner import PipelineRunner

        runner = PipelineRunner(_make_t1ss_config(tmp_dir, t1ss_fixture_gbff))
        results = runner.run(resume=False)

        failed = [r for r in results if not r.success]
        assert not failed, "Pipeline steps failed:\n" + "\n".join(
            f"  - {r.name}: {r.message}" for r in failed
        )

        # Output directory should be populated
        assert os.listdir(tmp_dir), "Output directory is empty after pipeline run"

        # Core artefacts expected by downstream consumers. The runner
        # writes intermediate workfiles into a tempdir (`runner.work_dir`)
        # and copies final products into `outdir` (== tmp_dir for this
        # test). The proteins FASTA is the intermediate one — fetched
        # via the runner's file registry, not by guessing the path.
        proteins_path = runner.files.get("proteins")
        assert proteins_path and os.path.exists(proteins_path), (
            f"Expected protein FASTA missing from runner.files. "
            f"runner.files keys: {list(runner.files.keys())}"
        )
        assert os.path.getsize(proteins_path) > 0

    def test_bimeno_04457_appears_in_outputs(self, tmp_dir, t1ss_fixture_gbff):
        """Known T1SS substrate BIMENO_04457 should appear in at least one
        ssign output file. This is the headline biological assertion — if
        it fails, the pipeline has regressed on real secretion biology."""
        _skip_unless_pipeline_prereqs_ready()

        from ssign_app.core.runner import PipelineRunner

        PipelineRunner(_make_t1ss_config(tmp_dir, t1ss_fixture_gbff)).run(resume=False)

        found_in = []
        for root, _dirs, files in os.walk(tmp_dir):
            for name in files:
                path = os.path.join(root, name)
                try:
                    with open(path) as f:
                        if "BIMENO_04457" in f.read():
                            found_in.append(os.path.relpath(path, tmp_dir))
                except (UnicodeDecodeError, PermissionError):
                    continue

        assert found_in, (
            "BIMENO_04457 (the known T1SS substrate on the fixture contig) "
            "did not appear in any output file. Either the fixture has "
            "been corrupted or the pipeline dropped the protein somewhere "
            "between extraction and reporting."
        )
