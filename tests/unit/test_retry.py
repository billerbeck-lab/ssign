"""Tests for ssign_lib/retry.py — `retry_once` wrapper.

Two-attempt contract: if the first call raises, sleep `delay` seconds
and try once more. Return `(result, "success")` if either attempt
succeeds; `(None, "failed: <reason>")` after two consecutive failures.

Used by per-protein web-API steps (DTU BioLib) where a single
retransmit is worth the wait but further retries waste time on a
genuinely broken request.
"""

import logging
import os
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from ssign_lib.retry import retry_once  # noqa: E402

# ---------------------------------------------------------------------------
# Success on first attempt
# ---------------------------------------------------------------------------


class TestFirstAttemptSucceeds:
    def test_returns_result_and_success_status(self):
        def func(item_id):
            return f"value_{item_id}"

        result, status = retry_once(func, "GENE_001")
        assert result == "value_GENE_001"
        assert status == "success"

    def test_no_sleep_when_first_attempt_succeeds(self, monkeypatch):
        # `time.sleep` only runs in the retry branch; first-attempt success
        # should never sleep.
        sleeps = []
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: sleeps.append(s))

        def func(item_id):
            return "ok"

        retry_once(func, "GENE_001")
        assert sleeps == []

    def test_kwargs_forwarded(self):
        # Extra kwargs are passed through to the wrapped callable.
        def func(item_id, *, multiplier):
            return item_id * multiplier

        result, status = retry_once(func, "X", multiplier=3)
        assert result == "XXX"
        assert status == "success"


# ---------------------------------------------------------------------------
# Retry path
# ---------------------------------------------------------------------------


class TestRetryAfterFailure:
    def test_succeeds_on_second_attempt(self, monkeypatch):
        # First call raises, second returns. Final result is the second-
        # attempt return; status is "success" (no per-attempt distinction).
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)
        attempts = []

        def func(item_id):
            attempts.append(item_id)
            if len(attempts) == 1:
                raise RuntimeError("transient")
            return "ok_on_retry"

        result, status = retry_once(func, "GENE_001")
        assert result == "ok_on_retry"
        assert status == "success"
        assert len(attempts) == 2

    def test_sleep_invoked_with_delay(self, monkeypatch):
        # The delay= argument controls the inter-attempt sleep duration.
        sleeps = []
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: sleeps.append(s))

        def func(item_id):
            if not sleeps:
                raise RuntimeError("first")
            return "ok"

        retry_once(func, "X", delay=2.5)
        assert sleeps == [2.5]


# ---------------------------------------------------------------------------
# Both attempts fail
# ---------------------------------------------------------------------------


class TestBothAttemptsFail:
    def test_returns_none_with_failure_status(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)

        def func(item_id):
            raise RuntimeError("boom")

        result, status = retry_once(func, "GENE_001")
        assert result is None
        assert status.startswith("failed:")
        assert "RuntimeError" in status
        assert "boom" in status

    def test_status_includes_exception_type_and_message(self, monkeypatch):
        # Different exception types produce identifiable statuses for the
        # caller's manifest log.
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)

        def func(item_id):
            raise ValueError("bad input")

        _, status = retry_once(func, "GENE_001")
        assert "ValueError" in status
        assert "bad input" in status

    def test_long_error_message_truncated_to_200(self, monkeypatch):
        # Defensive: an exception with a 5 KB message must not flood the log.
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)

        def func(item_id):
            raise RuntimeError("X" * 5000)

        _, status = retry_once(func, "GENE_001")
        # status = "failed: RuntimeError: <truncated msg>"
        # truncation cap is 200 on the message portion.
        msg_portion = status.split(": ", 2)[-1]
        assert len(msg_portion) <= 200

    def test_warning_then_error_logged(self, monkeypatch, caplog):
        # First failure → warning ("retrying"); second failure → error
        # ("skipping"). Pin the log-level contract.
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)

        def func(item_id):
            raise RuntimeError("boom")

        with caplog.at_level(logging.DEBUG, logger="ssign_lib.retry"):
            retry_once(func, "GENE_001")

        levels = [r.levelno for r in caplog.records]
        assert logging.WARNING in levels
        assert logging.ERROR in levels


# ---------------------------------------------------------------------------
# Exception types
# ---------------------------------------------------------------------------


class TestExceptionTypes:
    @pytest.mark.parametrize(
        "exc_type",
        [RuntimeError, ValueError, ConnectionError, TimeoutError, OSError],
    )
    def test_any_exception_caught_and_retried(self, monkeypatch, exc_type):
        # Catches `Exception` — every common per-item failure should retry.
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)
        attempts = []

        def func(item_id):
            attempts.append(item_id)
            if len(attempts) == 1:
                raise exc_type("first")
            return "ok"

        result, status = retry_once(func, "X")
        assert result == "ok"
        assert status == "success"

    def test_keyboard_interrupt_not_caught(self, monkeypatch):
        # KeyboardInterrupt inherits from BaseException, not Exception —
        # must propagate so Ctrl-C still works mid-retry.
        monkeypatch.setattr("ssign_lib.retry.time.sleep", lambda s: None)

        def func(item_id):
            raise KeyboardInterrupt()

        with pytest.raises(KeyboardInterrupt):
            retry_once(func, "X")
