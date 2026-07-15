"""Tests for runner.run_script's stream_stderr option.

The default (stream_stderr=False) buffers stderr until the subprocess exits.
stream_stderr=True forwards each line to the runner logger live so long
steps like pLM-BLAST embedding can surface incremental progress to the run log.
"""

import logging

import pytest

from ssign_app.core import runner


@pytest.fixture
def planted_script(tmp_path, monkeypatch):
    """Drop a tiny script in a fake BIN_DIR and patch the runner to find it."""

    def _plant(filename, source):
        bin_dir = tmp_path / "scripts"
        bin_dir.mkdir(exist_ok=True)
        script = bin_dir / filename
        script.write_text(source)
        monkeypatch.setattr(runner, "BIN_DIR", bin_dir)
        return filename

    return _plant


def test_default_does_not_log_stderr(planted_script, caplog):
    name = planted_script(
        "quiet_child.py",
        "import sys\nprint('chatter', file=sys.stderr)\n",
    )

    with caplog.at_level(logging.INFO, logger=runner.__name__):
        rc, stdout, stderr = runner.run_script(name, [])

    assert rc == 0
    assert "chatter" in stderr
    forwarded = [r for r in caplog.records if "chatter" in r.getMessage()]
    assert forwarded == [], "default path must not forward stderr to the logger"


def test_stream_stderr_forwards_each_line(planted_script, caplog):
    name = planted_script(
        "noisy_child.py",
        (
            "import sys, time\n"
            "for i in range(3):\n"
            "    print(f'line {i}', file=sys.stderr, flush=True)\n"
            "print('out-only', flush=True)\n"
        ),
    )

    with caplog.at_level(logging.INFO, logger=runner.__name__):
        rc, stdout, stderr = runner.run_script(name, [], stream_stderr=True)

    assert rc == 0
    assert "line 0" in stderr and "line 1" in stderr and "line 2" in stderr
    assert "out-only" in stdout
    forwarded = [r.getMessage() for r in caplog.records if "noisy_child.py" in r.getMessage()]
    assert any("line 0" in m for m in forwarded)
    assert any("line 2" in m for m in forwarded)
    assert all("out-only" not in m for m in forwarded), "stdout should not be forwarded"


def test_stream_stderr_logs_lines_on_nonzero_exit(planted_script, caplog):
    name = planted_script(
        "failing_child.py",
        ("import sys\nprint('about to fail', file=sys.stderr, flush=True)\nsys.exit(7)\n"),
    )

    with caplog.at_level(logging.INFO, logger=runner.__name__):
        rc, stdout, stderr = runner.run_script(name, [], stream_stderr=True)

    assert rc == 7
    assert "about to fail" in stderr
    forwarded = [r.getMessage() for r in caplog.records if "about to fail" in r.getMessage()]
    assert forwarded, "stderr lines should be logged even when the child fails"


def test_stream_stderr_timeout_kills_child(planted_script):
    name = planted_script(
        "sleeper_child.py",
        "import time\ntime.sleep(30)\n",
    )

    rc, stdout, stderr = runner.run_script(name, [], timeout=1, stream_stderr=True)

    assert rc == -1
    assert "Timeout" in stderr
