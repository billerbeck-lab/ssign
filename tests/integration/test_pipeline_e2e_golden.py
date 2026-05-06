"""End-to-end golden-output regression test.

Runs the full ssign pipeline on the minimal T5aSS fixture and asserts every
produced TSV/CSV/text output is byte-identical to the frozen reference in
`tests/fixtures/golden/t5ass_minimal/`. This is the canonical scientific-
output regression check: a code change that quietly alters a downstream
column, drops a row, or shifts a probability fails this test.

The test runs offline. DeepLocPro must be installed locally and pointed at
via `SSIGN_DEEPLOCPRO_PATH` (DTU academic license). Every other prediction
and annotation tool is disabled so the pipeline finishes in ~2 minutes on
CPU.

Outputs covered:
  - User-facing (outdir): {sid}_results.csv, _results_raw.csv, _summary.txt
  - Intermediate (work_dir): gene_info, gene_order, valid_systems,
    ss_components, deeplocpro, predictions, substrates, substrates_filtered,
    t5ss_substrates, integrated, enrichment_fisher

Figures are checked for existence only — matplotlib output is not
byte-stable across OS / font versions.

When the diff is intentional, regenerate references per
`tests/fixtures/golden/REGENERATE.md`.
"""

import filecmp
import os
import sys

import pytest

pytestmark = pytest.mark.integration


SRC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from conftest import (  # noqa: E402
    T5ASS_MINIMAL_FIXTURE_GBFF,
    skip_unless_pipeline_prereqs_ready,
)

_GOLDEN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "fixtures", "golden", "t5ass_minimal"))

_EXPECTED_FIGURES = [
    "fig1_ss_type_distribution.png",
    "fig2_tool_coverage.png",
    "fig3_protein_lengths.png",
    "fig5_category_distribution.png",
    "fig7_functional_summary.png",
]


def _golden_populated() -> bool:
    if not os.path.isdir(_GOLDEN_DIR):
        return False
    return any(not name.startswith(".") for name in os.listdir(_GOLDEN_DIR))


def _make_config(tmp_outdir: str, dlp_path: str):
    """Build the offline-only PipelineConfig used to regenerate refs."""
    from ssign_app.core.runner import PipelineConfig

    return PipelineConfig(
        input_path=T5ASS_MINIMAL_FIXTURE_GBFF,
        sample_id="t5ass_minimal",
        outdir=tmp_outdir,
        use_input_annotations=True,
        wholeness_threshold=0.8,
        excluded_systems=["Flagellum", "Tad", "T3SS"],
        conf_threshold=0.8,
        proximity_window=3,
        deeplocpro_mode="local",
        deeplocpro_path=dlp_path,
        skip_deepsece=True,
        skip_signalp=True,
        skip_plm_effector=True,
        skip_blastp=True,
        skip_hhsuite=True,
        skip_interproscan=True,
        skip_plmblast=True,
        skip_protparam=True,
        skip_structure=True,
    )


def _resolve_produced_path(name: str, outdir: str, work_dir: str) -> str | None:
    """Find a produced file by name — outdir first, then work_dir."""
    for parent in (outdir, work_dir):
        candidate = os.path.join(parent, name)
        if os.path.isfile(candidate):
            return candidate
    return None


def test_pipeline_e2e_matches_golden_outputs(tmp_dir):
    """Run the pipeline on the minimal T5aSS fixture and diff every output
    against the frozen reference. Byte-identical match required."""
    if not _golden_populated():
        pytest.skip("tests/fixtures/golden/t5ass_minimal/ is empty — populate via tests/fixtures/golden/REGENERATE.md.")

    dlp_path = skip_unless_pipeline_prereqs_ready(require_dlp_local=True)

    from ssign_app.core.runner import PipelineRunner

    # Regen mode: SSIGN_GOLDEN_REGEN_DIR pins outdir so produced files
    # survive after the test exits and can be copied into the reference dir.
    outdir = os.environ.get("SSIGN_GOLDEN_REGEN_DIR", tmp_dir)
    if outdir != tmp_dir:
        os.makedirs(outdir, exist_ok=True)

    runner = PipelineRunner(_make_config(outdir, dlp_path))
    results = runner.run(resume=False)

    failed = [r for r in results if not r.success]
    assert not failed, "Pipeline steps failed:\n" + "\n".join(f"  - {r.name}: {r.message}" for r in failed)

    diffs: list[str] = []
    for ref_name in sorted(os.listdir(_GOLDEN_DIR)):
        if ref_name.startswith("."):
            continue
        ref_path = os.path.join(_GOLDEN_DIR, ref_name)
        produced = _resolve_produced_path(ref_name, outdir, runner.work_dir)
        if produced is None:
            diffs.append(f"  MISSING  {ref_name} (not produced in outdir or work_dir)")
            continue
        if not filecmp.cmp(ref_path, produced, shallow=False):
            diffs.append(f"  DIFFERS  {ref_name}\n           ref:      {ref_path}\n           produced: {produced}")

    # Figures: existence check only — matplotlib output is not byte-stable.
    figures_dir = os.path.join(outdir, "figures", "t5ass_minimal")
    if os.path.isdir(figures_dir):
        present = set(os.listdir(figures_dir))
        for fig_name in _EXPECTED_FIGURES:
            if fig_name not in present:
                diffs.append(f"  MISSING  figures/t5ass_minimal/{fig_name}")
    else:
        diffs.append("  MISSING  figures/t5ass_minimal/ directory")

    if diffs:
        diff_block = "\n".join(diffs)
        pytest.fail(
            f"Golden-output regression detected:\n{diff_block}\n\n"
            f"If the change is intentional, regenerate references per "
            f"tests/fixtures/golden/REGENERATE.md.\n"
            f"runner.work_dir={runner.work_dir}\n"
            f"outdir={outdir}\n"
        )
