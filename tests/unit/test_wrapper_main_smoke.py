"""Smoke tests: each external-tool wrapper's main() runs far enough to
reach its first external call without raising a missing-symbol error.

What this catches: NameError / ImportError / AttributeError inside
main() before the subprocess (or model) is invoked. The class of bug
that bit DSE + SignalP in commit 6db791c: the wrapper builds its cmd
list (or calls an env-setup helper), references an undefined name,
and crashes in <1 second on the cluster — yet the unit suite passes
because no test ever ran the cmd-construction path.

What this does NOT cover:
- Whether the external binary exists or works (integration tests).
- Whether the cmd flags are semantically valid for the tool.
- Whether the output file is parsed correctly (per-wrapper unit tests).

The shape:
- sys.argv is set per wrapper to a minimal valid invocation.
- subprocess.run is stubbed to return success.
- The wrapper's main() is imported and called.
- Any exception that is NOT NameError / ImportError / AttributeError
  is treated as expected (e.g. FileNotFoundError parsing the fake
  empty output, RuntimeError from "results file missing"). Those
  failures live downstream of the bug class we're guarding against.
"""

from __future__ import annotations

import importlib
import sys
from types import SimpleNamespace

import pytest

# The bug class we want this test to surface. Any other exception means
# main() got far enough to invoke the external entry point — which is
# the whole point.
GUARDED_BUGS = (NameError, ImportError, AttributeError)


def _fake_subprocess_run(*args, **kwargs):
    """Stand-in for subprocess.run that returns a successful result."""
    return SimpleNamespace(returncode=0, stdout="", stderr="", args=args[0] if args else [])


def _fake_completed_process_with_output(out_file_path: str):
    """Some wrappers check that the subprocess produced an output file
    before parsing it. Patching subprocess.run to also touch the
    expected output path lets main() reach the parse stage."""

    def _runner(*args, **kwargs):
        # The cmd list is args[0]; many tools write to a dir given via
        # a flag. Easier: just touch the expected output path and let
        # the parse stage choke on empty content (caught as non-bug-class).
        try:
            with open(out_file_path, "w") as f:
                f.write("")
        except OSError:
            pass
        return SimpleNamespace(returncode=0, stdout="", stderr="", args=args[0] if args else [])

    return _runner


@pytest.fixture
def tiny_fasta(tmp_path):
    """One protein, valid FASTA."""
    p = tmp_path / "in.faa"
    p.write_text(">test_protein\nMAGSKLVAVLFLLSLLLSPGSDA\n")
    return p


@pytest.fixture
def tiny_substrates_tsv(tmp_path):
    """Substrates TSV with locus_tag column (the only column most
    annotation wrappers need to find a substrate to annotate)."""
    p = tmp_path / "subs.tsv"
    p.write_text("locus_tag\tsequence\ntest_protein\tMAGSKLVAVL\n")
    return p


def _run_wrapper_main(module_name: str, argv: list[str], monkeypatch, extra_patches=None):
    """Common harness. Imports wrapper module fresh, sets argv, patches
    subprocess.run + any extras, runs main(). Re-raises only on the
    guarded bug class."""
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr("subprocess.run", _fake_subprocess_run)
    # Bring in or refresh the module — wrappers do sys.path mangling at
    # import time, so a fresh import is safest.
    if module_name in sys.modules:
        module = importlib.reload(sys.modules[module_name])
    else:
        module = importlib.import_module(module_name)

    if extra_patches:
        for patch in extra_patches:
            patch(monkeypatch, module)

    try:
        module.main()
    except GUARDED_BUGS:
        raise
    except SystemExit:
        # argparse exits cleanly on --help, or main() called sys.exit on
        # a known error — not a bug-class failure.
        pass
    except Exception:
        # Any other exception (FileNotFoundError, RuntimeError parsing
        # fake output, ConnectionError on remote mode, etc.) is expected.
        # The bug class is what matters here.
        pass


class TestPredictionWrappers:
    """The three secretion-prediction wrappers (DLP, DSE, SignalP).
    These are the parallel group whose CPU-sharing rework introduced the
    NameErrors we're guarding against."""

    def test_signalp_local_mode(self, tmp_path, tiny_fasta, monkeypatch):
        argv = [
            "run_signalp.py",
            "--input",
            str(tiny_fasta),
            "--sample",
            "test",
            "--output",
            str(tmp_path / "out.tsv"),
            "--mode",
            "local",
            "--signalp-path",
            "/nonexistent",
        ]
        _run_wrapper_main("ssign_app.scripts.run_signalp", argv, monkeypatch)

    def test_deepsece(self, tmp_path, tiny_fasta, monkeypatch):
        # Stub the in-module run_deepsece() function so we don't actually
        # load DeepSecE model weights. The bug class we're guarding lives
        # BEFORE this call: in main() at line 533 (torch.set_num_threads
        # without `import torch`).
        def _fake_run(*args, **kwargs):
            return str(tmp_path / "fake_dse_out.tsv")

        argv = [
            "run_deepsece.py",
            "--input",
            str(tiny_fasta),
            "--sample",
            "test",
            "--output",
            str(tmp_path / "out.tsv"),
        ]
        _run_wrapper_main(
            "ssign_app.scripts.run_deepsece",
            argv,
            monkeypatch,
            extra_patches=[lambda mp, mod: mp.setattr(mod, "run_deepsece", _fake_run)],
        )

    def test_deeplocpro(self, tmp_path, tiny_fasta, monkeypatch):
        argv = [
            "run_deeplocpro.py",
            "--input",
            str(tiny_fasta),
            "--sample",
            "test",
            "--output",
            str(tmp_path / "out.tsv"),
            "--mode",
            "local",
            "--deeplocpro-path",
            "/nonexistent",
        ]
        _run_wrapper_main("ssign_app.scripts.run_deeplocpro", argv, monkeypatch)


class TestAnnotationWrappers:
    """The annotation parallel group (IPS, EggNOG, BLASTp, pLM-BLAST,
    HH-suite). These all had CPU-sharing changes in commit 82cece9."""

    def test_interproscan(self, tmp_path, tiny_fasta, tiny_substrates_tsv, monkeypatch):
        argv = [
            "run_interproscan.py",
            "--substrates",
            str(tiny_substrates_tsv),
            "--proteins",
            str(tiny_fasta),
            "--sample",
            "test",
            "--output",
            str(tmp_path / "ips.tsv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_interproscan", argv, monkeypatch)

    def test_eggnog(self, tmp_path, tiny_fasta, tiny_substrates_tsv, monkeypatch):
        db = tmp_path / "eggnog_db"
        db.mkdir()
        argv = [
            "run_eggnog.py",
            "--substrates",
            str(tiny_substrates_tsv),
            "--proteins",
            str(tiny_fasta),
            "--db",
            str(db),
            "--sample",
            "test",
            "--out",
            str(tmp_path / "eggnog.tsv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_eggnog", argv, monkeypatch)

    def test_blastp(self, tmp_path, tiny_fasta, tiny_substrates_tsv, monkeypatch):
        argv = [
            "run_blastp.py",
            "--substrates",
            str(tiny_substrates_tsv),
            "--proteins",
            str(tiny_fasta),
            "--sample",
            "test",
            "--db",
            "/nonexistent/nr",
            "--output",
            str(tmp_path / "blast.tsv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_blastp", argv, monkeypatch)

    def test_plm_blast(self, tmp_path, tiny_fasta, tiny_substrates_tsv, monkeypatch):
        ecod = tmp_path / "ecod"
        ecod.mkdir()
        argv = [
            "run_plm_blast.py",
            "--substrates",
            str(tiny_substrates_tsv),
            "--proteins",
            str(tiny_fasta),
            "--ecod-db",
            str(ecod),
            "--out",
            str(tmp_path / "plmblast.tsv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_plm_blast", argv, monkeypatch)

    def test_hhsuite(self, tmp_path, tiny_fasta, tiny_substrates_tsv, monkeypatch):
        argv = [
            "run_hhsuite.py",
            "--substrates",
            str(tiny_substrates_tsv),
            "--proteins",
            str(tiny_fasta),
            "--sample",
            "test",
            "--output",
            str(tmp_path / "hh.tsv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_hhsuite", argv, monkeypatch)


class TestWholeGenomeWrappers:
    """Wrappers that operate on the whole genome (Bakta) or the combined
    substrate set across genomes (ortholog grouping)."""

    def test_bakta(self, tmp_path, tiny_fasta, monkeypatch):
        db = tmp_path / "bakta_db"
        db.mkdir()
        argv = [
            "run_bakta.py",
            "--input",
            str(tiny_fasta),
            "--db",
            str(db),
            "--sample",
            "test",
            "--out-proteins",
            str(tmp_path / "proteins.faa"),
            "--out-gene-info",
            str(tmp_path / "gene_info.tsv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_bakta", argv, monkeypatch)

    def test_ortholog_grouping(self, tmp_path, tiny_fasta, monkeypatch):
        argv = [
            "run_ortholog_grouping.py",
            "--substrates-fasta",
            str(tiny_fasta),
            "--output",
            str(tmp_path / "groups.csv"),
        ]
        _run_wrapper_main("ssign_app.scripts.run_ortholog_grouping", argv, monkeypatch)
