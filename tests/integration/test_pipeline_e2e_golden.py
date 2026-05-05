"""End-to-end golden-output regression test (skeleton).

The real test is populated in Task 4.i — this file scaffolds the structure
so the reference directory layout, normalisation rules, and skip logic are
all in place when 4.i comes around. Until then, the test skips with a
pointer to `tests/fixtures/golden/REGENERATE.md`.

When implemented, the flow is:

    1. Run the full ssign pipeline on the minimal T5aSS fixture.
    2. Walk every file under tests/fixtures/golden/t5ass_minimal/.
    3. For each, normalise non-deterministic fields (timestamps, paths,
       run UUIDs, ProtParam float precision artefacts) on both sides.
    4. Assert the normalised contents match byte-for-byte.

CI failure on diff is the strongest signal that a code change altered
scientific output. Updating the references is gated by REGENERATE.md.
"""

import os

import pytest

pytestmark = pytest.mark.integration


_GOLDEN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "fixtures", "golden", "t5ass_minimal"))


def _golden_populated():
    """True once 4.i has copied reference outputs into tests/fixtures/golden/."""
    if not os.path.isdir(_GOLDEN_DIR):
        return False
    # Empty directory or directory with only dotfiles → not populated
    return any(not name.startswith(".") for name in os.listdir(_GOLDEN_DIR))


def test_pipeline_e2e_matches_golden_outputs():
    """Run the pipeline on the minimal T5aSS fixture and diff against frozen refs.

    Skipped until Task 4.i populates `tests/fixtures/golden/t5ass_minimal/`.
    Once populated, this is the canonical scientific-output regression check
    for v1.0.0+.
    """
    if not _golden_populated():
        pytest.skip(
            "tests/fixtures/golden/t5ass_minimal/ is empty — populate via "
            "tests/fixtures/golden/REGENERATE.md (Task 4.i)."
        )

    if os.environ.get("SSIGN_RUN_FULL_PIPELINE") != "1":
        pytest.skip(
            "SSIGN_RUN_FULL_PIPELINE=1 not set — golden-output test is "
            "opt-in because the full pipeline takes several minutes."
        )

    # 4.i implementation lands here:
    #   1. Build PipelineConfig pointing at the minimal T5aSS fixture
    #   2. Run PipelineRunner end-to-end into a tmp_dir
    #   3. For each golden ref file, normalise + assert match
    pytest.skip(
        "Golden-output diff body lands in Task 4.i — see tests/fixtures/golden/REGENERATE.md for the planned flow."
    )
