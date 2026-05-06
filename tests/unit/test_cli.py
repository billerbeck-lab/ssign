"""Tests for ssign_app/cli.py — the `ssign` entry point.

Mostly thin glue around `streamlit run Home.py`. Testable surfaces:

- `--version` short-circuit (prints version, never spawns Streamlit).
- argparse defaults for `--port` and `--no-browser`.
- App-file missing → exits 1.
- Streamlit-not-found (FileNotFoundError from `subprocess.Popen`) → exits 1.
- KeyboardInterrupt → graceful "stopped" message.

The full Streamlit subprocess path is exercised when the user actually
runs `ssign`, not in unit tests.
"""

import sys

import pytest

# cli imports as a sub-module of the installed package; no SCRIPTS_DIR shim.
from ssign_app import cli

# ---------------------------------------------------------------------------
# --version
# ---------------------------------------------------------------------------


class TestVersion:
    def test_version_flag_prints_and_returns(self, monkeypatch, capsys):
        # --version must not spawn a subprocess.
        spawned = []
        monkeypatch.setattr(cli.subprocess, "Popen", lambda *a, **k: spawned.append(a))
        monkeypatch.setattr(sys, "argv", ["ssign", "--version"])
        cli.main()
        out = capsys.readouterr().out
        assert "ssign" in out
        # No subprocess invocation.
        assert spawned == []


# ---------------------------------------------------------------------------
# argparse defaults
# ---------------------------------------------------------------------------


class TestArgparseDefaults:
    def test_help_exits_cleanly(self, monkeypatch):
        # `--help` triggers SystemExit(0) — argparse standard behaviour.
        monkeypatch.setattr(sys, "argv", ["ssign", "--help"])
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 0

    def test_unknown_flag_exits_with_error(self, monkeypatch):
        # argparse default: unknown flag → SystemExit(2)
        monkeypatch.setattr(sys, "argv", ["ssign", "--bogus-flag"])
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 2


# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------


class TestBanner:
    def test_banner_contains_tagline(self):
        # The "ssign" in the banner is ASCII art; the human-readable tagline
        # is the only literal substring we can reliably pin against.
        assert "Identification" in cli.BANNER
        assert "Gram Negatives" in cli.BANNER


# ---------------------------------------------------------------------------
# Subprocess and error paths
# ---------------------------------------------------------------------------


class _FakeSocket:
    def __init__(self, free=True):
        self._free = free

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def bind(self, addr):
        if not self._free:
            raise OSError("port in use")


class _FakePopen:
    """Stand-in for subprocess.Popen — empty stdout/stderr, wait()→0."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.stderr = iter([])
        self.stdout = iter([])

    def wait(self):
        return 0


@pytest.fixture
def streamlit_env(monkeypatch, tmp_path):
    """Stage the cli's environment for a fake-Streamlit run.

    Provides:
      - A fake `Home.py` next to the cli module (so the app-file existence
        check passes).
      - A free socket so the port-finder doesn't loop.
      - A `captured_cmd` list and a Popen that appends argv into it.
    Tests then monkeypatch `cli.subprocess.Popen` themselves if they need
    a custom mock (e.g., to raise an exception).
    """
    fake_home = tmp_path / "Home.py"
    fake_home.write_text("# stub")
    monkeypatch.setattr(cli.os.path, "abspath", lambda _: str(tmp_path / "cli.py"))
    monkeypatch.setattr(cli.socket, "socket", lambda *a, **k: _FakeSocket(free=True))

    captured_cmd: list[str] = []

    class _RecordingPopen(_FakePopen):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            captured_cmd.extend(args[0])

    monkeypatch.setattr(cli.subprocess, "Popen", _RecordingPopen)
    return captured_cmd


class TestStreamlitPath:
    def test_streamlit_not_found_exits_1(self, monkeypatch, streamlit_env):
        # Override Popen with one that raises FileNotFoundError.
        def _raise(*a, **k):
            raise FileNotFoundError("streamlit")

        monkeypatch.setattr(cli.subprocess, "Popen", _raise)
        monkeypatch.setattr(sys, "argv", ["ssign", "--no-browser"])
        with pytest.raises(SystemExit) as exc:
            cli.main()
        assert exc.value.code == 1

    def test_keyboard_interrupt_handled_cleanly(self, monkeypatch, streamlit_env, capsys):
        def _raise(*a, **k):
            raise KeyboardInterrupt()

        monkeypatch.setattr(cli.subprocess, "Popen", _raise)
        monkeypatch.setattr(sys, "argv", ["ssign", "--no-browser"])
        cli.main()  # KeyboardInterrupt must be swallowed
        assert "stopped" in capsys.readouterr().out.lower()


# ---------------------------------------------------------------------------
# Streamlit command construction
# ---------------------------------------------------------------------------


class TestStreamlitCmd:
    def test_streamlit_invoked_with_app_file_and_port(self, monkeypatch, streamlit_env):
        monkeypatch.setattr(sys, "argv", ["ssign", "--no-browser", "--port", "9999"])
        cli.main()
        assert "streamlit" in streamlit_env
        assert "run" in streamlit_env
        port_idx = streamlit_env.index("--server.port")
        assert streamlit_env[port_idx + 1] == "9999"
        headless_idx = streamlit_env.index("--server.headless")
        assert streamlit_env[headless_idx + 1] == "true"

    def test_streamlit_uses_python_executable(self, monkeypatch, streamlit_env):
        # `sys.executable -m streamlit ...` for venv-safety.
        monkeypatch.setattr(sys, "argv", ["ssign", "--no-browser"])
        cli.main()
        assert streamlit_env[0] == sys.executable
        assert streamlit_env[1] == "-m"
        assert streamlit_env[2] == "streamlit"

    def test_browser_mode_passes_headless_false(self, monkeypatch, streamlit_env):
        # Without --no-browser, headless must be "false" so Streamlit
        # auto-opens the browser.
        monkeypatch.setattr(sys, "argv", ["ssign"])
        cli.main()
        headless_idx = streamlit_env.index("--server.headless")
        assert streamlit_env[headless_idx + 1] == "false"
