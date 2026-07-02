"""Unit tests for the size-aware tool timeout helper.

Uses the real bundled coefficients.json (no mocks) so the tests exercise the
actual fits the runner will use.
"""

import math

from ssign_app.runtime.effort_model import effort, load_coefficients
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
    # blastp has no coefficient block -> effort is None -> floor unchanged.
    assert effort("blastp", 100_000, "fixed", COEFFS) is None
    assert scaled_timeout("blastp", 100_000, "fixed") == TOOL_TIMEOUT_S


def test_unmodelled_tool_respects_custom_floor():
    assert scaled_timeout("blastp", 100_000, "fixed", floor=999) == 999


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


def test_degenerate_low_fit_protected_by_floor():
    # SignalP whole_genome is a garbage flat fit (~13 min) but SignalP is fast,
    # so the floor keeps its cap safe rather than the tiny prediction.
    out = scaled_timeout("signalp", 160_831, "whole_genome")
    assert out == TOOL_TIMEOUT_S


def test_returns_int():
    assert isinstance(scaled_timeout("deeplocpro", 160_831, "whole_genome"), int)
