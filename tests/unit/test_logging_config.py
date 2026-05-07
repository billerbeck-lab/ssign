"""Tests for ssign_lib.logging_config."""

import logging

from ssign_lib.logging_config import LOG_FORMAT, init_script_logger


def test_log_format_matches_existing_scripts():
    # Pin: every run_*.py module currently calls
    # logging.basicConfig(level=INFO, format="%(levelname)s: %(message)s").
    # The constant must stay in lockstep so a single-edit migration of
    # the 27 pipeline scripts to init_script_logger() keeps stderr output
    # byte-identical to the pre-migration state.
    assert LOG_FORMAT == "%(levelname)s: %(message)s"


def test_init_script_logger_returns_named_logger():
    lg = init_script_logger("ssign_test_logger_module_a")
    assert isinstance(lg, logging.Logger)
    assert lg.name == "ssign_test_logger_module_a"


def test_init_script_logger_is_idempotent():
    # Two calls from different module names share the root handler set
    # rather than stacking duplicate handlers — basicConfig short-circuits
    # once the root logger already has at least one handler.
    init_script_logger("ssign_test_logger_module_b")
    handlers_before = list(logging.getLogger().handlers)
    init_script_logger("ssign_test_logger_module_c")
    handlers_after = list(logging.getLogger().handlers)
    assert len(handlers_after) == len(handlers_before)
