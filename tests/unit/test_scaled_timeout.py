"""Unit tests for the size-aware tool timeout helper.

Uses the real bundled coefficients.json (no mocks) so the tests exercise the
actual fits the runner will use.
"""

import math
import pytest

from ssign_app.runtime.effort_model import effort, load_coefficients, resolve_regime
from ssign_app.runtime.timeouts import (
    DEFAULT_MARGIN,
    LOW_CONFIDENCE_MARGIN,
    scaled_timeout,
)
from ssign_app.scripts.ssign_lib.constants import TOOL_TIMEOUT_S

COEFFS = load_coefficients()


def test_large_whole_genome_scales_above_floor():
    # 160,831-protein pool (the failing benchmark case): predicted ~5.4h,
    # 2x margin ~11h, well above the 4h floor.
    n = 160_831
    e = effort("deeplocpro", n, "whole_genome", COEFFS)
    assert e is not None and not e.low_confidence
    out = scaled_timeout("deeplocpro", n, "whole_genome")
    assert out == math.ceil(DEFAULT_MARGIN * e.seconds)
    assert out > TOOL_TIMEOUT_S


def test_small_input_stays_at_floor():
    # A few hundred neighborhood proteins: 2x predicted is far below the floor.
    out = scaled_timeout("deeplocpro", 200, "neighborhood")
    assert out == TOOL_TIMEOUT_S


def test_unmodelled_tool_returns_floor():
    # A tool with no coefficient block -> effort is None -> floor unchanged.
    assert effort("made_up_tool", 100_000, "substrates", COEFFS) is None
    assert scaled_timeout("made_up_tool", 100_000, "substrates") == TOOL_TIMEOUT_S


def test_unmodelled_tool_respects_custom_floor():
    assert scaled_timeout("made_up_tool", 100_000, "substrates", floor=999) == 999


def test_confident_fit_uses_default_margin():
    # low_confidence=False fit -> 2x. Force a tiny floor so scaling shows through.
    n = 160_831
    e = effort("deeplocpro", n, "whole_genome", COEFFS)
    assert not e.low_confidence
    out = scaled_timeout("deeplocpro", n, "whole_genome", floor=1)
    assert out == math.ceil(DEFAULT_MARGIN * e.seconds)


def test_low_confidence_fit_uses_wider_margin():
    # eggnog substrates is low_confidence=True -> wider margin than the default.
    n = 5_000
    e = effort("eggnog", n, "substrates", COEFFS)
    assert e is not None and e.low_confidence
    out = scaled_timeout("eggnog", n, "substrates", floor=1)
    assert out == math.ceil(LOW_CONFIDENCE_MARGIN * e.seconds)
    assert LOW_CONFIDENCE_MARGIN > DEFAULT_MARGIN


def test_hhsuite_large_substrate_set_scales_above_floor():
    # Pooled full-tier reruns (~593/768 substrates) overran the flat 4h floor.
    # hh_suite.substrates is a low_confidence prior (a=35, b=600) -> 3x margin
    # lifts a big pool well above the floor so it isn't killed mid-step.
    n = 768
    e = effort("hh_suite", n, "substrates", COEFFS)
    assert e is not None and e.low_confidence
    out = scaled_timeout("hh_suite", n, "substrates")
    assert out == math.ceil(LOW_CONFIDENCE_MARGIN * e.seconds)
    assert out > TOOL_TIMEOUT_S


def test_hhsuite_small_substrate_set_stays_at_floor():
    # A single genome's handful of substrates: 3x predicted is below the 4h
    # floor, so the floor is unchanged (no regression for small runs).
    out = scaled_timeout("hh_suite", 30, "substrates")
    assert out == TOOL_TIMEOUT_S


def test_blastp_is_substrate_scoped_and_scales_above_floor():
    # blastp BLASTs the substrate set (regime "substrates", not "fixed"). Its
    # prior (a=20, b=300, low_confidence) is sized for the slow opt-in NR path,
    # so a big substrate set gets headroom above the 4h floor.
    assert resolve_regime("blastp") == "substrates"
    n = 768
    e = effort("blastp", n, "substrates", COEFFS)
    assert e is not None and e.low_confidence
    out = scaled_timeout("blastp", n, "substrates")
    assert out == math.ceil(LOW_CONFIDENCE_MARGIN * e.seconds)
    assert out > TOOL_TIMEOUT_S


def test_blastp_small_substrate_set_stays_at_floor():
    # Swiss-Prot default finishes in minutes; a handful of substrates stays at
    # the floor (the prior only lifts large NR runs above it).
    out = scaled_timeout("blastp", 30, "substrates")
    assert out == TOOL_TIMEOUT_S


def test_unmodeled_tool_protected_by_floor():
    # A tool absent from the effort model (e.g. plm_effector, removed
    # 2026-07-15) has no fit, so scaled_timeout falls back to the floor
    # rather than under-sizing the timeout. (SignalP whole_genome used to
    # be a degenerate a=0 fit with the same outcome; it was proxied to the
    # DeepLocPro rate on 2026-07-03 after the pooled 160k run showed its
    # a=0 fit floored SignalP at 4h and killed it.)
    out = scaled_timeout("plm_effector", 160_831, "whole_genome")
    assert out == TOOL_TIMEOUT_S


def test_returns_int():
    assert isinstance(scaled_timeout("deeplocpro", 160_831, "whole_genome"), int)


def test_machine_rate_slow_widens_timeout():
    # A machine at half the reference GPU rate (0.5) should roughly double the
    # cap above what the reference estimate alone would give, when both clear the floor.
    n = 160_831
    ref = scaled_timeout("deeplocpro", n, "whole_genome")
    slow = scaled_timeout("deeplocpro", n, "whole_genome", machine_rate=0.5)
    assert slow > ref
    assert slow == pytest.approx(2 * ref, rel=0.01)


def test_machine_rate_fast_never_below_floor():
    # A fast machine (rate 4) shrinks the raw prediction, but the floor still
    # backstops it so a still-running tool can't be false-killed.
    out = scaled_timeout("deeplocpro", 160_831, "whole_genome", machine_rate=4.0)
    assert out >= TOOL_TIMEOUT_S


def test_machine_rate_none_is_unchanged():
    n = 160_831
    assert scaled_timeout("deeplocpro", n, "whole_genome") == scaled_timeout(
        "deeplocpro", n, "whole_genome", machine_rate=None
    )


def test_machine_rate_ignored_for_unmodeled_tool():
    # No fit -> floor, regardless of any machine rate.
    assert scaled_timeout("plm_effector", 160_831, "whole_genome", machine_rate=0.1) == TOOL_TIMEOUT_S
