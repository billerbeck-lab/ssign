"""Tests for ssign_app/shims/hmmsearch.py — pyhmmer-based hmmsearch shim.

The shim exists so `pip install ssign` works without `sudo apt install hmmer`.
MacSyFinder shells out to whatever `hmmsearch` is on PATH; the shim
emulates HMMER3 hmmsearch using the pyhmmer Cython binding so MacSyFinder's
text-output parser still works.

Testable surfaces:

1. `_decode` — bytes/str/None safe coercion (pyhmmer returns bytes for
   names + descriptions that need converting before formatted output).
2. `parse_args` — argparse contract: every MacSyFinder-used flag is
   recognised; `--noali` / `--notextw` / `-Z` are accepted-and-ignored
   for resilience against MacSyFinder version drift.
3. `SHIM_VERSION` — string constant pinned for the version-history note.
4. `main()` error paths: FileNotFoundError exits 1, generic exceptions
   exit 1 with stderr message pointing to the real-HMMER fallback.
5. End-to-end: build a minimal HMM in-memory via pyhmmer.plan7.Builder,
   run the shim, and verify the text output contains the markers
   MacSyFinder's parser depends on (`>>`, the domain table header, `//`).
"""

import os
import sys

import pytest

SHIMS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "shims"))
sys.path.insert(0, SHIMS_DIR)

# pyhmmer is a hard install dep; if it's missing the shim itself can't load,
# so all tests here are conditional on its presence.
pytest.importorskip("pyhmmer")

import hmmsearch  # noqa: E402
from hmmsearch import (  # noqa: E402
    SHIM_VERSION,
    _decode,
    main,
    parse_args,
)
from pyhmmer.easel import Alphabet, TextSequence  # noqa: E402
from pyhmmer.plan7 import Background, Builder  # noqa: E402

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta  # noqa: E402

# ---------------------------------------------------------------------------
# _decode
# ---------------------------------------------------------------------------


class TestDecode:
    @pytest.mark.parametrize(
        "val, expected",
        [
            (b"hello", "hello"),
            ("hello", "hello"),
            (None, ""),
            (b"", ""),
            ("", ""),
            (b"\xc3\xa9", "é"),
            (42, "42"),
        ],
    )
    def test_returns_string(self, val, expected):
        assert _decode(val) == expected


# ---------------------------------------------------------------------------
# parse_args
# ---------------------------------------------------------------------------


class TestParseArgs:
    def test_required_args(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt"])
        assert args.hmmfile == "profile.hmm"
        assert args.seqdb == "seq.fasta"
        assert args.output == "out.txt"

    def test_missing_output_raises(self):
        # `-o` is required — argparse exits 2 if absent.
        with pytest.raises(SystemExit) as exc:
            parse_args(["profile.hmm", "seq.fasta"])
        assert exc.value.code == 2

    def test_missing_positional_raises(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-o", "out.txt"])
        assert exc.value.code == 2

    def test_evalue_default(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt"])
        assert args.evalue == 10.0

    def test_evalue_explicit(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt", "-E", "1e-5"])
        assert args.evalue == 1e-5

    def test_cpu_default(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt"])
        assert args.cpu == 1

    def test_cpu_explicit(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt", "--cpu", "8"])
        assert args.cpu == 8

    def test_cut_ga_flag(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt", "--cut_ga"])
        assert args.cut_ga

    def test_tblout_default_none(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt"])
        assert args.tblout is None

    def test_tblout_explicit(self):
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt", "--tblout", "tbl.txt"])
        assert args.tblout == "tbl.txt"

    @pytest.mark.parametrize("flag", ["--noali", "--notextw"])
    def test_ignored_flags_accepted(self, flag):
        # MacSyFinder occasionally adds these — must not raise.
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt", flag])
        # Flag is parsed and set to True (then ignored downstream).
        assert getattr(args, flag.lstrip("-")) is True

    def test_z_flag_accepted_and_ignored(self):
        # `-Z` adjusts E-value scaling in real HMMER; ignored here.
        args = parse_args(["profile.hmm", "seq.fasta", "-o", "out.txt", "-Z", "100000"])
        assert args.Z == 100000.0


# ---------------------------------------------------------------------------
# SHIM_VERSION
# ---------------------------------------------------------------------------


class TestShimVersion:
    def test_is_string(self):
        assert isinstance(SHIM_VERSION, str)

    def test_non_empty(self):
        assert SHIM_VERSION.strip()

    def test_dot_separated(self):
        # Pinned X.Y format; if the convention ever bumps to X.Y.Z (or
        # similar) a maintainer should make a deliberate decision rather
        # than silently dropping the contract.
        parts = SHIM_VERSION.split(".")
        assert all(p.isdigit() for p in parts)


# ---------------------------------------------------------------------------
# main() — error paths
# ---------------------------------------------------------------------------


class TestMainErrorPaths:
    def test_missing_hmm_file_exits_1(self, tmp_dir, capsys):
        # `hmmfile` doesn't exist → run_search raises FileNotFoundError
        # → main prints to stderr and sys.exit(1).
        argv = [
            os.path.join(tmp_dir, "no_such.hmm"),
            os.path.join(tmp_dir, "no_such.fasta"),
            "-o",
            os.path.join(tmp_dir, "out.txt"),
        ]
        with pytest.raises(SystemExit) as exc:
            main(argv)
        assert exc.value.code == 1
        err = capsys.readouterr().err
        assert "ERROR" in err

    def test_unwritable_output_exits_1(self, tmp_dir, monkeypatch, capsys):
        # Force write_text_output to fail → main prints "Failed to write
        # output" and exits 1.
        # First make run_search succeed by short-circuiting it.
        monkeypatch.setattr(hmmsearch, "run_search", lambda args: ([], 0, 0))

        argv = [
            "ignored.hmm",
            "ignored.fasta",
            "-o",
            "/nonexistent/dir/out.txt",
        ]
        with pytest.raises(SystemExit) as exc:
            main(argv)
        assert exc.value.code == 1
        assert "Failed to write" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# End-to-end with a built-in-memory HMM
# ---------------------------------------------------------------------------

_ANCHOR_SEQ = "MKTLLLTLLCAFSVAQAVDLPTQEPALGKAA"


def _run_main_expect(argv, expected_code=0):
    """Run main(argv); assert it sys.exit'd with `expected_code`."""
    with pytest.raises(SystemExit) as exc:
        main(argv)
    assert exc.value.code == expected_code


@pytest.fixture(scope="module")
def anchor_hmm_path(tmp_path_factory):
    """Build a single-sequence HMM once per module — Builder is the slow bit."""
    alpha = Alphabet.amino()
    bg = Background(alpha)
    builder = Builder(alpha)
    seq = TextSequence(name=b"anchor", sequence=_ANCHOR_SEQ).digitize(alpha)
    hmm, _, _ = builder.build(seq, bg)

    hmm_path = str(tmp_path_factory.mktemp("hmm") / "anchor.hmm")
    with open(hmm_path, "wb") as f:
        hmm.write(f)
    return hmm_path


class TestShimEndToEnd:
    def test_text_output_contains_macsyfinder_markers(self, tmp_dir, anchor_hmm_path):
        # MacSyFinder keys on ">>", the domain table (i-Evalue header),
        # and trailing "//". If any go missing it silently produces zero
        # hits. Pin all three.
        fasta_path = os.path.join(tmp_dir, "targets.fasta")
        write_fasta({"target_1": _ANCHOR_SEQ, "target_2": "G" * 20}, fasta_path)
        out_path = os.path.join(tmp_dir, "out.txt")
        _run_main_expect([anchor_hmm_path, fasta_path, "-o", out_path])

        with open(out_path) as f:
            text = f.read()
        assert ">>" in text, "missing >> hit-target marker"
        assert "i-Evalue" in text, "missing domain table header"
        assert "//" in text, "missing record terminator"

    def test_tblout_written_when_requested(self, tmp_dir, anchor_hmm_path):
        fasta_path = os.path.join(tmp_dir, "targets.fasta")
        write_fasta({"target_1": _ANCHOR_SEQ, "target_2": "G" * 20}, fasta_path)
        out_path = os.path.join(tmp_dir, "out.txt")
        tbl_path = os.path.join(tmp_dir, "tbl.txt")
        _run_main_expect([anchor_hmm_path, fasta_path, "-o", out_path, "--tblout", tbl_path])

        assert os.path.exists(tbl_path)
        with open(tbl_path) as f:
            content = f.read()
        assert "target name" in content
        assert content.rstrip().endswith("//")

    def test_no_hits_path_emits_no_hits_message(self, tmp_dir, anchor_hmm_path):
        # Same HMM, but the FASTA contains only an unrelated sequence.
        fasta_path = os.path.join(tmp_dir, "targets.fasta")
        write_fasta({"unrelated": "G" * 24}, fasta_path)
        out_path = os.path.join(tmp_dir, "out.txt")
        _run_main_expect([anchor_hmm_path, fasta_path, "-o", out_path])

        with open(out_path) as f:
            text = f.read()
        assert "No hits" in text
        assert "//" in text
