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


from ssign_lib import retry as retry_with_backoff_module  # noqa: E402
from ssign_lib.retry import retry_once, retry_with_backoff  # noqa: E402

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


# ---------------------------------------------------------------------------
# retry_with_backoff
# ---------------------------------------------------------------------------


class TestRetryWithBackoff:
    """retry_with_backoff: N-attempt linear backoff for batch operations.

    Used by DTU SignalP + DLP remote wrappers — replaces the duplicated
    `for attempt in range(...): ... time.sleep(30 * attempt)` blocks.
    """

    def test_first_attempt_succeeds_returns_result(self):
        calls = []

        def func():
            calls.append(1)
            return "ok"

        assert retry_with_backoff(func, max_attempts=3, label="t") == "ok"
        assert len(calls) == 1

    def test_succeeds_on_third_attempt(self, monkeypatch):
        sleeps: list[float] = []
        monkeypatch.setattr(retry_with_backoff_module.time, "sleep", lambda s: sleeps.append(s))
        attempts = {"n": 0}

        def func():
            attempts["n"] += 1
            if attempts["n"] < 3:
                raise RuntimeError(f"flaky {attempts['n']}")
            return "finally"

        assert retry_with_backoff(func, max_attempts=3, initial_delay=10.0) == "finally"
        # Sleeps: 10s after attempt 1, 20s after attempt 2; nothing after 3.
        assert sleeps == [10.0, 20.0]

    def test_reraises_after_max_attempts(self, monkeypatch):
        monkeypatch.setattr(retry_with_backoff_module.time, "sleep", lambda s: None)

        def func():
            raise RuntimeError("doomed")

        with pytest.raises(RuntimeError, match="doomed"):
            retry_with_backoff(func, max_attempts=2, initial_delay=1.0)

    def test_only_retries_on_specified_exceptions(self, monkeypatch):
        # ValueError is NOT in retry_on=RuntimeError → no retry, immediate raise
        sleeps: list[float] = []
        monkeypatch.setattr(retry_with_backoff_module.time, "sleep", lambda s: sleeps.append(s))

        def func():
            raise ValueError("not retryable")

        with pytest.raises(ValueError):
            retry_with_backoff(func, max_attempts=3, retry_on=RuntimeError)
        assert sleeps == []

    def test_retry_on_tuple_of_exceptions(self, monkeypatch):
        monkeypatch.setattr(retry_with_backoff_module.time, "sleep", lambda s: None)
        attempts = {"n": 0}

        def func():
            attempts["n"] += 1
            if attempts["n"] == 1:
                raise RuntimeError("net")
            if attempts["n"] == 2:
                raise OSError("disk")
            return "ok"

        result = retry_with_backoff(func, max_attempts=3, retry_on=(RuntimeError, OSError))
        assert result == "ok"
        assert attempts["n"] == 3

    def test_invalid_max_attempts_rejected(self):
        with pytest.raises(ValueError, match="max_attempts"):
            retry_with_backoff(lambda: None, max_attempts=0)

    def test_max_attempts_one_means_no_retry(self):
        attempts = {"n": 0}

        def func():
            attempts["n"] += 1
            raise RuntimeError("nope")

        with pytest.raises(RuntimeError):
            retry_with_backoff(func, max_attempts=1)
        assert attempts["n"] == 1

    def test_label_appears_in_warning_log(self, monkeypatch, caplog):
        monkeypatch.setattr(retry_with_backoff_module.time, "sleep", lambda s: None)
        attempts = {"n": 0}

        def func():
            attempts["n"] += 1
            if attempts["n"] < 2:
                raise RuntimeError("flaky")
            return "ok"

        with caplog.at_level(logging.WARNING):
            retry_with_backoff(func, max_attempts=2, label="DTU SignalP", initial_delay=1.0)

        assert any("DTU SignalP" in rec.message for rec in caplog.records)
