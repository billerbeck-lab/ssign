#!/usr/bin/env python3
"""CLI entry point for ssign.

Usage:
    ssign              # Launch Streamlit GUI
    ssign --version    # Print version
    ssign --help       # Show help
"""

import argparse
import os
import re
import subprocess
import sys
import threading

BANNER = r"""
  ┌─────────────────────────────────────────┐
  │                                         │
  │    ___  ___  _  __ _ _ __               │
  │   / __|/ __|| |/ _` | '_ \              │
  │   \__ \\__ \| | (_| | | | |             │
  │   |___/|___/|_|\__, |_| |_|             │
  │                |___/                    │
  │                                         │
  │   Secretion-System Identification       │
  │   for Gram Negatives                    │
  │                                         │
  └─────────────────────────────────────────┘
"""


def main():
    parser = argparse.ArgumentParser(
        prog="ssign",
        description="Secretion-system Identification for Gram Negatives",
    )
    parser.add_argument(
        "--version", action="store_true",
        help="Print version and exit",
    )
    parser.add_argument(
        "--no-browser", action="store_true",
        help="Start the GUI server without opening a browser",
    )
    parser.add_argument(
        "--port", type=int, default=8501,
        help="Port for the Streamlit server (default: 8501)",
    )

    args = parser.parse_args()

    if args.version:
        from ssign_app import __version__
        print(f"ssign {__version__}")
        return

    print(BANNER, flush=True)

    # Find the Streamlit app file and config
    app_dir = os.path.dirname(os.path.abspath(__file__))
    app_file = os.path.join(app_dir, "Home.py")
    config_dir = os.path.join(app_dir, ".streamlit")

    if not os.path.exists(app_file):
        print(f"Error: Could not find {app_file}", file=sys.stderr)
        sys.exit(1)

    # Set STREAMLIT_CONFIG_DIR so it finds our bundled config.toml
    # This sets maxUploadSize=500MB and other defaults
    env = os.environ.copy()
    if os.path.isdir(config_dir):
        env["STREAMLIT_CONFIG_DIR"] = config_dir

    # Also set via env vars as a fallback
    env["STREAMLIT_SERVER_MAX_UPLOAD_SIZE"] = "500"
    env["STREAMLIT_SERVER_MAX_MESSAGE_SIZE"] = "500"
    env["STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION"] = "false"
    env["STREAMLIT_SERVER_ENABLE_CORS"] = "false"
    env["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"

    # Build Streamlit command
    cmd = [
        sys.executable, "-m", "streamlit", "run", app_file,
        "--server.port", str(args.port),
        "--server.headless", "true" if args.no_browser else "false",
        "--server.maxUploadSize", "500",
        "--server.maxMessageSize", "500",
        "--server.enableXsrfProtection", "false",
        "--server.enableCORS", "false",
    ]

    port = args.port
    if args.no_browser:
        print(f"  Open in browser: http://localhost:{port}", flush=True)
    else:
        print(f"  Opening http://localhost:{port} ...", flush=True)
    print(flush=True)

    try:
        # Filter Streamlit's startup banner from stderr while still
        # showing real errors (import failures, config problems, etc.)
        _banner_re = re.compile(
            r"You can now view|Network URL|External URL|Local URL"
        )

        proc = subprocess.Popen(
            cmd, env=env, stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE, text=True,
        )

        def _filter_stderr():
            for line in proc.stderr:
                if _banner_re.search(line) or not line.strip():
                    continue
                sys.stderr.write(line)
                sys.stderr.flush()

        t = threading.Thread(target=_filter_stderr, daemon=True)
        t.start()
        proc.wait()
        t.join(timeout=1)
    except KeyboardInterrupt:
        print("\nssign stopped.")
    except FileNotFoundError:
        print(
            "Error: Streamlit not found. Install with: pip install ssign",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
