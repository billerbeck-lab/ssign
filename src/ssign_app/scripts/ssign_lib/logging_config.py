"""Single source of truth for the format string used by every pipeline
script's stderr logging.

Every `run_<tool>.py` and helper currently sets up logging with the same
two-line preamble:

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logger = logging.getLogger(__name__)

The format choice (``LEVEL: message``) is deliberate — it matches the
output of pyhmmer / MacSyFinder / BLAST so users see uniform-looking
log lines regardless of which tool produced them. This module exists
so a future maintainer who wants to change that decision (e.g. add a
timestamp prefix once the webserver lands) updates one place.

Usage
-----

In a script's module body::

    from ssign_lib.logging_config import init_script_logger
    logger = init_script_logger(__name__)

That replaces the two-line preamble. Existing call sites can adopt this
incrementally; the format string stays in sync as long as new scripts
use the helper.
"""

import logging

LOG_FORMAT = "%(levelname)s: %(message)s"


def init_script_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Apply the standard ssign script-level logging config and return
    a per-module logger. Idempotent — calling from multiple scripts in
    the same process won't duplicate handlers because basicConfig() is
    a no-op once the root logger has handlers attached."""
    logging.basicConfig(level=level, format=LOG_FORMAT)
    return logging.getLogger(name)
