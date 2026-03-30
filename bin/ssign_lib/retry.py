"""Retry-once-then-skip logic for recoverable per-item failures.

Provides a simple wrapper that attempts a callable twice before giving up.
Used for network requests, file operations, and other per-protein steps
where a single retry is reasonable but further retries waste time.
"""
from __future__ import annotations

import logging
import time
from typing import Any, Callable

logger = logging.getLogger(__name__)


def retry_once(
    func: Callable, item_id: str, delay: float = 1.0, **kwargs
) -> tuple[Any, str]:
    """Execute *func*, retry once on failure, return ``(result, status)``.

    Args:
        func: Callable that takes *item_id* as its first argument plus
            any additional keyword arguments.
        item_id: Identifier for logging (e.g. a protein ID).
        delay: Seconds to wait between attempts.
        **kwargs: Extra keyword arguments forwarded to *func*.

    Returns:
        ``(result, "success")`` on success.
        ``(None, "failed: <reason>")`` after two consecutive failures.
    """
    for attempt in range(2):
        try:
            result = func(item_id, **kwargs)
            return result, "success"
        except Exception as e:
            error_msg = f"{type(e).__name__}: {str(e)[:200]}"
            if attempt == 0:
                logger.warning(
                    "%s: attempt 1 failed (%s), retrying in %ss...",
                    item_id,
                    error_msg,
                    delay,
                )
                time.sleep(delay)
            else:
                logger.error(
                    "%s: attempt 2 failed (%s), skipping", item_id, error_msg
                )
                return None, f"failed: {error_msg}"

    # Should never reach here, but guard against it.
    return None, "failed: unknown"
