"""Tests for src/ssign_app/core/runner.py.

Surface tested here:

- PipelineConfig dataclass — defaults, mutable-default safety, factories
- StepResult dataclass — basic shape
- run_script() — subprocess wrapper: missing script, exit codes, timeout
- PipelineRunner.__init__ — callback default, semaphore wiring, state
- check_dependencies() — warning detection via shutil.which monkeypatch
- _elapsed_str() — pre-start vs post-start formatting
- _save_progress() + load_progress() — JSON round-trip
- _try_resume() — completed-step set computation, missing/stale work_dir handling

Deep orchestration (the run() main loop, parallel groups, per-step _step_*
methods) is out of unit-test scope; tests/integration/test_pipeline_fixture.py
exercises that on the real T1SS fixture.
"""

import json
import os
import shutil
import sys
import threading
import time
from pathlib import Path

import pytest

from ssign_app.core import runner

# runner.py is the orchestrator under src/ssign_app/core/. Import via the
# installed package — editable install puts src/ on sys.path.
from ssign_app.core.runner import (  # noqa: E402
    PipelineConfig,
    PipelineRunner,
    StepResult,
    run_script,
)

# ---------------------------------------------------------------------------
# PipelineConfig
# ---------------------------------------------------------------------------


class TestPipelineConfig:
    def test_defaults_are_sane(self):
        c = PipelineConfig()
        assert c.wholeness_threshold == 0.8
        assert c.conf_threshold == 0.8
        assert c.proximity_window == 3
        assert c.required_fraction_correct == 0.8
        assert c.deeplocpro_mode == "remote"
        assert c.signalp_mode == "remote"

    def test_excluded_systems_default(self):
        c = PipelineConfig()
        assert c.excluded_systems == ["Flagellum", "Tad", "T3SS"]

    def test_excluded_systems_isolated_per_instance(self):
        # Mutable defaults via field(default_factory=...) — each instance
        # must get its own list. Pin this so a refactor to plain `= [...]`
        # fails the test.
        a = PipelineConfig()
        b = PipelineConfig()
        a.excluded_systems.append("MUTATED")
        assert "MUTATED" not in b.excluded_systems

    def test_plm_effector_types_default(self):
        c = PipelineConfig()
        assert c.plm_effector_types == ["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"]

    def test_plm_effector_types_isolated_per_instance(self):
        a = PipelineConfig()
        b = PipelineConfig()
        a.plm_effector_types.append("EXTRA")
        assert "EXTRA" not in b.plm_effector_types

    def test_skip_flags_align_with_install_tier(self):
        c = PipelineConfig()
        # base tier defaults: HH-suite, pLM-BLAST, EggNOG, PLM-Effector skipped
        assert c.skip_hhsuite is True
        assert c.skip_plmblast is True
        assert c.skip_eggnog is True
        assert c.skip_plm_effector is True
        # base tier defaults: BLASTp, IPS, ProtParam, DeepSecE, SignalP active
        assert c.skip_blastp is False
        assert c.skip_interproscan is False
        assert c.skip_protparam is False
        assert c.skip_deepsece is False
        assert c.skip_signalp is False

    def test_cpu_per_genome_factory_yields_positive_int(self):
        c = PipelineConfig()
        assert isinstance(c.cpu_per_genome, int)
        assert c.cpu_per_genome >= 1


class TestPipelineConfigHHsuiteEnvFallback:
    """SSIGN_HHSUITE_* env vars fill in empty fields; explicit paths win.

    Documented in --hhsuite-*-db CLI help and docs/optional_tools.md.
    """

    @pytest.mark.parametrize(
        "env_name, attr",
        [
            ("SSIGN_HHSUITE_PFAM", "hhsuite_pfam_db"),
            ("SSIGN_HHSUITE_PDB70", "hhsuite_pdb70_db"),
            ("SSIGN_HHSUITE_UNICLUST", "hhsuite_uniclust_db"),
        ],
    )
    def test_empty_field_picks_up_env_var(self, monkeypatch, env_name, attr):
        monkeypatch.setenv(env_name, "/tmp/fake_db")
        c = PipelineConfig()
        assert getattr(c, attr) == "/tmp/fake_db"

    @pytest.mark.parametrize(
        "env_name, attr",
        [
            ("SSIGN_HHSUITE_PFAM", "hhsuite_pfam_db"),
            ("SSIGN_HHSUITE_PDB70", "hhsuite_pdb70_db"),
            ("SSIGN_HHSUITE_UNICLUST", "hhsuite_uniclust_db"),
        ],
    )
    def test_explicit_path_overrides_env_var(self, monkeypatch, env_name, attr):
        monkeypatch.setenv(env_name, "/tmp/from_env")
        c = PipelineConfig(**{attr: "/explicit/path"})
        assert getattr(c, attr) == "/explicit/path"

    def test_no_env_var_set_keeps_empty_default(self, monkeypatch):
        for env in ("SSIGN_HHSUITE_PFAM", "SSIGN_HHSUITE_PDB70", "SSIGN_HHSUITE_UNICLUST"):
            monkeypatch.delenv(env, raising=False)
        c = PipelineConfig()
        assert c.hhsuite_pfam_db == ""
        assert c.hhsuite_pdb70_db == ""
        assert c.hhsuite_uniclust_db == ""

    def test_env_var_fallback_logs_audit_trail(self, monkeypatch, caplog):
        # Stale env vars in a dev shell should be visible — log INFO when the
        # fallback fires so the user can spot they're picking up a leftover.
        monkeypatch.setenv("SSIGN_HHSUITE_PFAM", "/tmp/fake_pfam")
        with caplog.at_level("INFO", logger="ssign_app.core.runner"):
            PipelineConfig()
        assert any("hhsuite_pfam_db" in rec.message and "/tmp/fake_pfam" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# StepResult
# ---------------------------------------------------------------------------


class TestStepResult:
    def test_default_output_files_is_empty_dict(self):
        r = StepResult(name="step", success=True, message="ok")
        assert r.output_files == {}

    def test_output_files_isolated_per_instance(self):
        a = StepResult(name="a", success=True, message="ok")
        b = StepResult(name="b", success=True, message="ok")
        a.output_files["k"] = "v"
        assert b.output_files == {}


# ---------------------------------------------------------------------------
# run_script — subprocess wrapper
# ---------------------------------------------------------------------------


class TestRunScript:
    def test_missing_script_returns_minus_one(self):
        rc, stdout, stderr = run_script("does_not_exist.py", [])
        assert rc == -1
        assert "not found" in stderr.lower()

    def test_real_exit_code_returned(self, tmp_dir, monkeypatch):
        # Drop a tiny Python script into BIN_DIR so run_script picks it up

        script_path = os.path.join(tmp_dir, "exit_code.py")
        with open(script_path, "w") as f:
            f.write("import sys; sys.exit(int(sys.argv[1]))\n")
        monkeypatch.setattr(runner, "BIN_DIR", Path(tmp_dir))
        rc, _, _ = run_script("exit_code.py", ["7"])
        assert rc == 7

    def test_stdout_captured(self, tmp_dir, monkeypatch):

        script_path = os.path.join(tmp_dir, "echo.py")
        with open(script_path, "w") as f:
            f.write("print('hello from script')\n")
        monkeypatch.setattr(runner, "BIN_DIR", Path(tmp_dir))
        rc, stdout, _ = run_script("echo.py", [])
        assert rc == 0
        assert "hello from script" in stdout

    def test_timeout_yields_minus_one(self, tmp_dir, monkeypatch):

        script_path = os.path.join(tmp_dir, "slow.py")
        with open(script_path, "w") as f:
            f.write("import time; time.sleep(10)\n")
        monkeypatch.setattr(runner, "BIN_DIR", Path(tmp_dir))
        rc, _, stderr = run_script("slow.py", [], timeout=1)
        assert rc == -1
        assert "timeout" in stderr.lower()


# ---------------------------------------------------------------------------
# PipelineRunner __init__ + small accessors
# ---------------------------------------------------------------------------


class TestRunnerInit:
    def test_default_progress_callback_is_noop(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        r = PipelineRunner(c)
        # Default callback must be safely invokable without args
        r.progress("step", 50, "msg")  # no exception

    def test_results_files_start_empty(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        r = PipelineRunner(c)
        assert r.results == []
        assert r.files == {}

    def test_api_semaphores_default_to_empty_dict(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        r = PipelineRunner(c)
        assert r.api_sem == {}

    def test_api_semaphores_passed_through(self, tmp_dir):

        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        sem = threading.Semaphore(2)
        r = PipelineRunner(c, api_semaphores={"dtu": sem})
        assert r.api_sem == {"dtu": sem}


class TestElapsedStr:
    def test_returns_empty_before_start(self, tmp_dir):
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir))
        assert r._elapsed_str() == ""

    def test_returns_seconds_only_under_a_minute(self, tmp_dir):
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir))
        r.start_time = time.monotonic() - 5
        assert r._elapsed_str() == "5s"

    def test_returns_minutes_and_seconds(self, tmp_dir):
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir))
        r.start_time = time.monotonic() - 125  # 2m 5s
        assert r._elapsed_str() == "2m 5s"


# ---------------------------------------------------------------------------
# check_dependencies
# ---------------------------------------------------------------------------


class TestCheckDependencies:
    @pytest.mark.parametrize("tool", ["hmmsearch", "macsyfinder"])
    def test_warns_when_core_tool_missing(self, tmp_dir, monkeypatch, tool):
        c = PipelineConfig(outdir=tmp_dir, skip_deepsece=True)
        r = PipelineRunner(c)
        monkeypatch.setattr(shutil, "which", lambda name: None)
        warnings = r.check_dependencies()
        assert any(tool in w for w in warnings)

    def test_no_warnings_when_all_present(self, tmp_dir, monkeypatch):
        c = PipelineConfig(outdir=tmp_dir, skip_deepsece=True)
        r = PipelineRunner(c)
        monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
        warnings = r.check_dependencies()
        assert warnings == []

    def test_skip_deepsece_suppresses_torch_check(self, tmp_dir, monkeypatch):
        # Even with torch absent, skip_deepsece=True must produce no torch warning
        c = PipelineConfig(outdir=tmp_dir, skip_deepsece=True)
        r = PipelineRunner(c)
        monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
        # Force the torch import to fail by clobbering it in sys.modules
        monkeypatch.setitem(sys.modules, "torch", None)
        warnings = r.check_dependencies()
        # No torch-related warning should appear when DSE is skipped
        assert not any("torch" in w for w in warnings)


# ---------------------------------------------------------------------------
# _save_progress + load_progress + _try_resume
# ---------------------------------------------------------------------------


class TestProgressRoundTrip:
    def test_save_then_load_recovers_steps(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r.results = [
            StepResult(name="extract_proteins", success=True, message="ok"),
            StepResult(name="macsyfinder", success=True, message="ok"),
            StepResult(name="deeplocpro", success=False, message="failed"),
        ]
        r.files = {"proteins": os.path.join(tmp_dir, "proteins.faa")}
        r._save_progress()

        loaded_results, loaded_files, loaded_work, _ = PipelineRunner.load_progress(
            tmp_dir,
            sample_id="my_sample",
        )
        assert [s.name for s in loaded_results] == [
            "extract_proteins",
            "macsyfinder",
            "deeplocpro",
        ]
        assert {s.name: s.success for s in loaded_results} == {
            "extract_proteins": True,
            "macsyfinder": True,
            "deeplocpro": False,
        }
        assert loaded_files["proteins"] == r.files["proteins"]
        assert loaded_work == tmp_dir

    def test_load_returns_none_when_no_progress_file(self, tmp_dir):
        results, files, work, cfg = PipelineRunner.load_progress(
            tmp_dir,
            sample_id="nope",
        )
        assert results is None
        assert files == {}
        assert work == ""
        assert cfg == {}

    def test_save_creates_hidden_ssign_dir(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r._save_progress()
        progress_file = os.path.join(tmp_dir, ".ssign", "my_sample_progress.json")
        assert os.path.exists(progress_file)
        with open(progress_file) as f:
            data = json.load(f)
        # Includes sample_id, work_dir, steps, files, config
        assert data["sample_id"] == "my_sample"
        assert "config" in data


class TestTryResume:
    def test_no_prior_progress_returns_empty_set(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="fresh")
        r = PipelineRunner(c)
        assert r._try_resume() == set()

    def test_resumes_only_successful_steps(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        r1.results = [
            StepResult(name="extract_proteins", success=True, message="ok"),
            StepResult(name="macsyfinder", success=True, message="ok"),
            StepResult(name="deeplocpro", success=False, message="failed"),
        ]
        r1._save_progress()

        r2 = PipelineRunner(c)
        completed = r2._try_resume()
        assert "extract_proteins" in completed
        assert "macsyfinder" in completed
        assert "deeplocpro" not in completed

    def test_resume_restores_work_dir_and_files(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        # Use a real file path so the resume validator keeps it
        proteins = os.path.join(tmp_dir, "proteins.faa")
        with open(proteins, "w") as f:
            f.write(">x\nMKT\n")
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1.files = {"proteins": proteins}
        r1._save_progress()

        r2 = PipelineRunner(c)
        r2._try_resume()
        assert r2.work_dir == tmp_dir
        assert r2.files["proteins"] == proteins

    def test_resume_skipped_if_work_dir_gone(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        # Reference a work_dir that doesn't exist
        r1.work_dir = os.path.join(tmp_dir, "deleted_work")
        os.makedirs(r1.work_dir)
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1._save_progress()
        # Wipe the work_dir before resume
        os.rmdir(r1.work_dir)

        r2 = PipelineRunner(c)
        completed = r2._try_resume()
        assert completed == set()


class TestReportFiguresUpstreamGuards:
    """`_step_report` and `_step_figures` must skip cleanly when their upstream
    `integrated` CSV is missing — otherwise they invoke generate_report.py /
    generate_figures.py with an empty path and produce a cryptic subprocess
    error."""

    @pytest.mark.parametrize("step_name", ["report", "figures"])
    def test_skip_when_integrated_missing(self, tmp_dir, step_name):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r.files = {}
        result = r._step_report() if step_name == "report" else r._step_figures()
        assert result.success is False
        assert result.name == step_name
        assert "skipping" in result.message.lower()

    @pytest.mark.parametrize("step_name", ["report", "figures"])
    def test_skip_when_integrated_path_does_not_exist(self, tmp_dir, step_name):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r.files = {"integrated": os.path.join(tmp_dir, "never_written.csv")}
        result = r._step_report() if step_name == "report" else r._step_figures()
        assert result.success is False
        assert "skipping" in result.message.lower()

    def test_resume_aborted_if_persisted_output_missing(self, tmp_dir):
        # Manifest says step succeeded with file X; X was deleted between
        # runs. Resume must refuse — better fresh than corrupt.
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        proteins = os.path.join(tmp_dir, "proteins.faa")
        with open(proteins, "w") as f:
            f.write(">x\nMKT\n")
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1.files = {"proteins": proteins}
        r1._save_progress()

        os.remove(proteins)

        r2 = PipelineRunner(c)
        assert r2._try_resume() == set()

    def test_resume_aborted_if_persisted_output_empty(self, tmp_dir):
        # Zero-byte TSV looks like exists==True but breaks downstream
        # parsers. Must not be silently accepted.
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        empty_path = os.path.join(tmp_dir, "ss_components.tsv")
        open(empty_path, "w").close()
        r1.results = [StepResult(name="macsyfinder", success=True, message="ok")]
        r1.files = {"ss_components": empty_path}
        r1._save_progress()

        r2 = PipelineRunner(c)
        assert r2._try_resume() == set()

    def test_resume_ignores_intentionally_blank_paths(self, tmp_dir):
        # Some steps record an empty-string path to mean "no output produced
        # this run" (e.g., skipped optional tool). Validator must skip those.
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        proteins = os.path.join(tmp_dir, "proteins.faa")
        with open(proteins, "w") as f:
            f.write(">x\nMKT\n")
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1.files = {"proteins": proteins, "blastp_results": ""}
        r1._save_progress()

        r2 = PipelineRunner(c)
        completed = r2._try_resume()
        assert "extract_proteins" in completed
        assert r2.files["blastp_results"] == ""


# ---------------------------------------------------------------------------
# Constants integration
# ---------------------------------------------------------------------------


def test_hhsuite_min_prob_default_matches_constants():
    """PipelineConfig.hhsuite_min_prob must mirror the canonical value in
    ssign_lib.constants — single source of truth."""
    from ssign_app.scripts.ssign_lib.constants import HHSUITE_MIN_PROB

    assert PipelineConfig().hhsuite_min_prob == HHSUITE_MIN_PROB
