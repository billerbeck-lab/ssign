"""Compatibility shims for external tool binaries.

These scripts emulate external command-line tools using pure-Python
libraries, so that `pip install ssign` works without system packages.

If a shim breaks due to upstream changes, the fallback is to install
the real tool (e.g., `sudo apt install hmmer`) which will take
precedence on PATH if installed to a system location.
"""
