"""Tests for the `ssign run` subcommand CLI parsing.

Coverage targets:
  - Every PipelineConfig field reachable through a CLI flag (argparse
    roundtrip).
  - Defaults match PipelineConfig defaults.
  - argparse.BooleanOptionalAction emits a `--no-<flag>` inverse.
  - List-typed flags accept multiple values and reach the config as a list.
  - Missing input file returns exit code 2 without invoking the runner.
  - Sample-id auto-derives from input filename stem when not provided.

We don't run the pipeline here — PipelineRunner is monkeypatched so each
test captures the config that would have been used.
"""

import sys

import pytest

from ssign_app import cli
from ssign_app.core.runner import PipelineConfig


@pytest.fixture
def fake_input(tmp_path):
    """Create an empty file standing in for an input genome."""
    p = tmp_path / "GENOME.gbff"
    p.write_text("LOCUS stub")
    return str(p)


@pytest.fixture
def captured_config(monkeypatch):
    """Replace PipelineRunner with a stub that captures the config and
    returns an empty results list. Returns a list[PipelineConfig] that
    gets populated on each main() invocation."""
    bag: list[PipelineConfig] = []

    class _StubRunner:
        def __init__(self, config, **kwargs):
            bag.append(config)

        def run(self, resume=False):
            return []

    import ssign_app.core.runner as runner_module

    monkeypatch.setattr(runner_module, "PipelineRunner", _StubRunner)
    return bag


def _invoke(monkeypatch, argv):
    """Run cli.main() with given argv, returning its exit code."""
    monkeypatch.setattr(sys, "argv", argv)
    return cli.main()


class TestDefaults:
    def test_minimal_invocation_uses_pipeline_defaults(self, monkeypatch, captured_config, fake_input):
        rc = _invoke(monkeypatch, ["ssign", "run", fake_input])
        assert rc == 0
        assert len(captured_config) == 1
        cfg = captured_config[0]
        defaults = PipelineConfig()
        # Spot-check every flag-controlled field against PipelineConfig defaults.
        assert cfg.outdir == defaults.outdir
        assert cfg.wholeness_threshold == defaults.wholeness_threshold
        assert cfg.excluded_systems == defaults.excluded_systems
        assert cfg.macsyfinder_db_type == defaults.macsyfinder_db_type
        assert cfg.conf_threshold == defaults.conf_threshold
        assert cfg.proximity_window == defaults.proximity_window
        assert cfg.required_fraction_correct == defaults.required_fraction_correct
        assert cfg.deeplocpro_mode == defaults.deeplocpro_mode
        assert cfg.signalp_mode == defaults.signalp_mode
        assert cfg.skip_signalp == defaults.skip_signalp
        assert cfg.skip_deepsece == defaults.skip_deepsece
        assert cfg.skip_blastp == defaults.skip_blastp
        assert cfg.skip_hhsuite == defaults.skip_hhsuite
        assert cfg.skip_interproscan == defaults.skip_interproscan
        assert cfg.skip_plmblast == defaults.skip_plmblast
        assert cfg.skip_eggnog == defaults.skip_eggnog
        assert cfg.skip_plm_effector == defaults.skip_plm_effector
        assert cfg.skip_protparam == defaults.skip_protparam
        assert cfg.dpi == defaults.dpi

    def test_sample_id_derives_from_input_filename_stem(self, monkeypatch, captured_config, fake_input):
        # fake_input is .../GENOME.gbff → sample_id = "GENOME"
        _invoke(monkeypatch, ["ssign", "run", fake_input])
        assert captured_config[0].sample_id == "GENOME"

    def test_explicit_sample_id_overrides_filename_derivation(self, monkeypatch, captured_config, fake_input):
        _invoke(monkeypatch, ["ssign", "run", fake_input, "--sample-id", "custom"])
        assert captured_config[0].sample_id == "custom"


class TestNumericFlags:
    @pytest.mark.parametrize(
        "flag, attr, value",
        [
            ("--wholeness-threshold", "wholeness_threshold", 0.5),
            ("--conf-threshold", "conf_threshold", 0.95),
            ("--proximity-window", "proximity_window", 7),
            ("--required-fraction-correct", "required_fraction_correct", 0.6),
            ("--deepsece-min-prob", "deepsece_min_prob", 0.7),
            ("--signalp-min-prob", "signalp_min_prob", 0.4),
            ("--bakta-threads", "bakta_threads", 8),
            ("--blastp-min-pident", "blastp_min_pident", 70.0),
            ("--blastp-min-qcov", "blastp_min_qcov", 75.0),
            ("--blastp-evalue", "blastp_evalue", 1e-10),
            ("--interproscan-min-evalue", "interproscan_min_evalue", 1e-3),
            ("--ortholog-min-pident", "ortholog_min_pident", 50.0),
            ("--ortholog-min-qcov", "ortholog_min_qcov", 80.0),
            ("--cpu-per-genome", "cpu_per_genome", 16),
            ("--dpi", "dpi", 200),
        ],
    )
    def test_numeric_flag_roundtrips_to_config(self, monkeypatch, captured_config, fake_input, flag, attr, value):
        _invoke(monkeypatch, ["ssign", "run", fake_input, flag, str(value)])
        assert getattr(captured_config[0], attr) == value


class TestBooleanFlags:
    @pytest.mark.parametrize(
        "attr, default",
        [
            ("skip_signalp", False),
            ("skip_deepsece", False),
            ("skip_blastp", False),
            ("skip_hhsuite", True),
            ("skip_interproscan", False),
            ("skip_plmblast", True),
            ("skip_eggnog", True),
            ("skip_plm_effector", True),
            ("skip_protparam", False),
            ("skip_structure", True),
            ("filter_dse_type_mismatch", True),
            ("use_input_annotations", False),
            ("run_bakta", False),
            ("dlp_whole_genome", False),
            ("dse_whole_genome", False),
            ("sp_whole_genome", False),
        ],
    )
    def test_boolean_flag_can_flip_either_way(self, monkeypatch, captured_config, fake_input, attr, default):
        # Affirmative — `--<flag>` always sets True.
        flag = "--" + attr.replace("_", "-")
        _invoke(monkeypatch, ["ssign", "run", fake_input, flag])
        assert getattr(captured_config[-1], attr) is True

        # Inverse — `--no-<flag>` always sets False.
        no_flag = "--no-" + attr.replace("_", "-")
        _invoke(monkeypatch, ["ssign", "run", fake_input, no_flag])
        assert getattr(captured_config[-1], attr) is False

    def test_default_when_neither_form_passed(self, monkeypatch, captured_config, fake_input):
        # No flag → PipelineConfig default applies.
        _invoke(monkeypatch, ["ssign", "run", fake_input])
        defaults = PipelineConfig()
        for attr in ("skip_blastp", "skip_hhsuite", "use_input_annotations"):
            assert getattr(captured_config[0], attr) == getattr(defaults, attr)


class TestListFlags:
    def test_excluded_systems_accepts_multiple_values(self, monkeypatch, captured_config, fake_input):
        _invoke(
            monkeypatch,
            ["ssign", "run", fake_input, "--excluded-systems", "Tad", "Flagellum"],
        )
        assert captured_config[0].excluded_systems == ["Tad", "Flagellum"]

    def test_plm_effector_types_accepts_multiple_values(self, monkeypatch, captured_config, fake_input):
        _invoke(
            monkeypatch,
            ["ssign", "run", fake_input, "--plm-effector-types", "T1SE", "T6SE"],
        )
        assert captured_config[0].plm_effector_types == ["T1SE", "T6SE"]


class TestChoiceFlags:
    @pytest.mark.parametrize("mode", ["local", "remote"])
    def test_deeplocpro_mode_accepts_local_or_remote(self, monkeypatch, captured_config, fake_input, mode):
        _invoke(
            monkeypatch,
            ["ssign", "run", fake_input, "--deeplocpro-mode", mode],
        )
        assert captured_config[0].deeplocpro_mode == mode

    @pytest.mark.parametrize("db_type", ["ordered_replicon", "unordered"])
    def test_macsyfinder_db_type_accepts_documented_values(self, monkeypatch, captured_config, fake_input, db_type):
        _invoke(
            monkeypatch,
            ["ssign", "run", fake_input, "--macsyfinder-db-type", db_type],
        )
        assert captured_config[0].macsyfinder_db_type == db_type

    def test_invalid_choice_rejected(self, monkeypatch, fake_input):
        # argparse should emit SystemExit(2) on an invalid --choices value.
        monkeypatch.setattr(
            sys,
            "argv",
            ["ssign", "run", fake_input, "--deeplocpro-mode", "garbage"],
        )
        with pytest.raises(SystemExit):
            cli.main()


class TestEarlyExit:
    def test_missing_input_file_returns_2(self, monkeypatch, captured_config):
        rc = _invoke(monkeypatch, ["ssign", "run", "/does/not/exist.gbff"])
        assert rc == 2
        # Runner must NOT have been invoked.
        assert captured_config == []
