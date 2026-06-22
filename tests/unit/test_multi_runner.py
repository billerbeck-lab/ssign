"""Tests for ssign_app.core.multi_runner.MultiGenomeRunner.

Covers:
- N=1 dispatches straight to PipelineRunner.run() (no segment-splitting code path)
- Constructor validation (empty, duplicate sample_ids, reserved name,
  separator-in-sample_id, heterogeneous skip flags)
- slice_stages_by_segment correctly partitions a stages list
- N>1 flow calls pool/split helpers and per-segment runners in the
  right order (verified via mocked PipelineRunner methods)
"""

from __future__ import annotations

from dataclasses import dataclass
from unittest.mock import MagicMock, patch

import pytest

from ssign_app.core._pool_utils import SEPARATOR
from ssign_app.core.multi_runner import (
    POOL_SAMPLE_ID,
    MultiGenomeRunner,
    slice_stages_by_segment,
)

# ---------------------------------------------------------------------------
# Lightweight config + stages stand-ins
# ---------------------------------------------------------------------------


@dataclass
class _StubConfig:
    """A drop-in PipelineConfig stand-in. Includes every field that
    MultiGenomeRunner reaches into (validation flags + the fields
    overridden in ``_make_pool_runner`` via dataclasses.replace())."""

    sample_id: str
    outdir: str = "/tmp/out"
    input_path: str = ""
    monitor_resources: bool = True
    skip_deeplocpro: bool = False
    skip_deepsece: bool = False
    skip_signalp: bool = False
    skip_plm_effector: bool = False
    skip_blastp: bool = False
    skip_hhsuite: bool = True  # match runner default for tests
    skip_interproscan: bool = False
    skip_eggnog: bool = False
    skip_plmblast: bool = False
    skip_protparam: bool = False
    enrichment_stats: bool = False
    t5ass_annotate_whole: bool = False


def _stub_step(step_id: str):
    """Return a callable whose __name__ mirrors a real _step_* method."""
    fn = MagicMock()
    fn.__name__ = f"_step_{step_id}"
    return fn


# ---------------------------------------------------------------------------
# slice_stages_by_segment
# ---------------------------------------------------------------------------


class TestSliceStagesBySegment:
    def test_canonical_layout(self):
        # Mimics the shape returned by PipelineRunner._build_stages():
        # 5 segment-A tuples, 1 parallel group (B prediction), 1 segment-B
        # tuple (plm_effector), 4 segment-C tuples, 1 parallel group (D
        # annotation), 5 segment-E tuples.
        stages = [
            ("Detecting input format", _stub_step("detect_format")),
            ("Extracting proteins", _stub_step("extract_proteins")),
            ("Running MacSyFinder", _stub_step("macsyfinder")),
            ("Validating secretion systems", _stub_step("validate_systems")),
            ("Extracting SS neighborhood", _stub_step("extract_neighborhood")),
            [
                ("DLP", _stub_step("deeplocpro")),
                ("DSE", _stub_step("deepsece")),
                ("SignalP", _stub_step("signalp")),
            ],
            ("PLM-Effector", _stub_step("plm_effector")),
            ("Cross-validating", _stub_step("cross_validate")),
            ("Proximity", _stub_step("proximity")),
            ("T5SS", _stub_step("t5ss")),
            ("Filtering", _stub_step("filtering")),
            [
                ("BLASTp", _stub_step("blastp")),
                ("EggNOG", _stub_step("eggnog")),
                ("IPS", _stub_step("interproscan")),
                ("pLM-BLAST", _stub_step("plm_blast")),
                ("ProtParam", _stub_step("protparam")),
            ],
            ("Integrate", _stub_step("integrate")),
            ("Orthologs", _stub_step("orthologs")),
            ("Enrichment", _stub_step("enrichment")),
            ("Report", _stub_step("report")),
            ("Figures", _stub_step("figures")),
        ]
        segments = slice_stages_by_segment(stages)
        # A is the 5 tuples before the first parallel group
        assert len(segments["A"]) == 5
        # B is the prediction parallel group + plm_effector
        assert len(segments["B"]) == 2
        assert isinstance(segments["B"][0], list)
        # C is the 4 sequential tuples between the parallel groups
        assert len(segments["C"]) == 4
        # D is the annotation parallel group
        assert len(segments["D"]) == 1
        assert isinstance(segments["D"][0], list)
        # E is the 5 trailing tuples
        assert len(segments["E"]) == 5

    def test_pre_prediction_steps_route_to_A(self):
        # Every pre-prediction step (before the first parallel group) lands in A.
        stages = [
            ("Detecting input format", _stub_step("detect_format")),
            ("Extracting SS neighborhood", _stub_step("extract_neighborhood")),
            [("DLP", _stub_step("deeplocpro"))],
            ("Cross-validating", _stub_step("cross_validate")),
            [("EggNOG", _stub_step("eggnog"))],
            ("Integrate", _stub_step("integrate")),
        ]
        segments = slice_stages_by_segment(stages)
        assert len(segments["A"]) == 2  # detect_format + extract_neighborhood

    def test_empty_parallel_group_skipped(self):
        # All predictions skipped → empty prediction parallel list.
        stages = [
            ("Detecting input format", _stub_step("detect_format")),
            [],  # empty parallel group (all predictions skipped)
            ("Cross-validating", _stub_step("cross_validate")),
        ]
        segments = slice_stages_by_segment(stages)
        assert len(segments["A"]) == 1
        assert segments["B"] == []
        assert len(segments["C"]) == 1


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------


class TestConstructorValidation:
    def test_empty_configs_raises(self):
        with pytest.raises(ValueError, match="at least one"):
            MultiGenomeRunner([])

    def test_duplicate_sample_ids_raises(self):
        configs = [_StubConfig(sample_id="ecoli"), _StubConfig(sample_id="ecoli")]
        with pytest.raises(ValueError, match="unique"):
            MultiGenomeRunner(configs)

    def test_pool_sample_id_reserved(self):
        configs = [
            _StubConfig(sample_id="ecoli"),
            _StubConfig(sample_id=POOL_SAMPLE_ID),
        ]
        with pytest.raises(ValueError, match="reserved"):
            MultiGenomeRunner(configs)

    def test_separator_in_sample_id_raises(self):
        configs = [_StubConfig(sample_id=f"bad{SEPARATOR}id")]
        with pytest.raises(ValueError, match="separator"):
            MultiGenomeRunner(configs)

    def test_heterogeneous_skip_flags_raises(self):
        configs = [
            _StubConfig(sample_id="ecoli", skip_eggnog=False),
            _StubConfig(sample_id="pao1", skip_eggnog=True),
        ]
        with pytest.raises(ValueError, match="skip_eggnog"):
            MultiGenomeRunner(configs)

    def test_uniform_configs_pass(self):
        configs = [
            _StubConfig(sample_id="ecoli"),
            _StubConfig(sample_id="pao1"),
            _StubConfig(sample_id="vc"),
        ]
        # Should not raise.
        runner = MultiGenomeRunner(configs)
        assert [c.sample_id for c in runner.configs] == ["ecoli", "pao1", "vc"]


# ---------------------------------------------------------------------------
# N=1 dispatch — delegates straight to PipelineRunner
# ---------------------------------------------------------------------------


class TestSingleGenomeDispatch:
    def test_n1_delegates_to_pipeline_runner(self, monkeypatch):
        # Patch PipelineRunner so we can capture how MultiGenomeRunner uses it.
        fake_instance = MagicMock()
        fake_instance.run.return_value = ["fake_result_1", "fake_result_2"]
        fake_cls = MagicMock(return_value=fake_instance)
        monkeypatch.setattr("ssign_app.core.multi_runner.PipelineRunner", fake_cls)

        config = _StubConfig(sample_id="ecoli")
        mgr = MultiGenomeRunner([config], progress_callback=None)
        out = mgr.run(resume=True)

        # PipelineRunner constructed exactly once, with the single config.
        fake_cls.assert_called_once()
        args, kwargs = fake_cls.call_args
        assert args[0] is config
        # run() called with resume=True
        fake_instance.run.assert_called_once_with(resume=True)
        # Output keyed by sample_id
        assert out == {"ecoli": ["fake_result_1", "fake_result_2"]}

    def test_n1_does_not_invoke_multi_path(self, monkeypatch):
        # If N=1 accidentally fell into _run_multi, it would try to pool —
        # which would touch the pool_utils helpers. Patch them so we'd notice.
        fake_runner = MagicMock()
        fake_runner.run.return_value = []
        monkeypatch.setattr(
            "ssign_app.core.multi_runner.PipelineRunner",
            MagicMock(return_value=fake_runner),
        )
        pool_fastas_spy = MagicMock()
        monkeypatch.setattr("ssign_app.core.multi_runner.pool_fastas", pool_fastas_spy)

        MultiGenomeRunner([_StubConfig(sample_id="solo")]).run()
        pool_fastas_spy.assert_not_called()


# ---------------------------------------------------------------------------
# N>1 flow — exercise the segment dance with mocked PipelineRunner methods
# ---------------------------------------------------------------------------


@pytest.fixture
def n2_runner(monkeypatch, tmp_path):
    """Build a MultiGenomeRunner with N=2 stub configs and patched PipelineRunner.

    Returns ``(mgr, fake_runners_by_sid, fake_pool_runner, call_log)``:
    - mgr: the MultiGenomeRunner under test
    - fake_runners_by_sid: dict mapping sample_id -> mocked per-genome runner
    - fake_pool_runner: the mocked pool runner (the third PipelineRunner instance)
    - call_log: list of (instance, method, args, kwargs) for ordering assertions
    """
    configs = [
        _StubConfig(sample_id="g1", outdir=str(tmp_path / "g1")),
        _StubConfig(sample_id="g2", outdir=str(tmp_path / "g2")),
    ]
    for c in configs:
        # The constructor for the *real* PipelineRunner needs outdir to exist
        # only when run() is called; our mock doesn't care, but make the dirs
        # so _make_pool_runner's commonpath logic produces tmp_path.
        import os as _os

        _os.makedirs(c.outdir, exist_ok=True)

    call_log: list[tuple] = []

    def _make_fake_runner(sample_id: str) -> MagicMock:
        fake = MagicMock()
        fake.config = MagicMock()
        fake.config.sample_id = sample_id
        fake.config.outdir = str(tmp_path / sample_id) if sample_id != POOL_SAMPLE_ID else str(tmp_path / "_pool")
        fake.files = {}
        fake.work_dir = str(tmp_path / f"work_{sample_id}")
        import os as _os

        _os.makedirs(fake.work_dir, exist_ok=True)
        fake.results = []

        # Stub _build_stages to return a minimal stages list with one tuple
        # per segment letter so all five segments are populated.
        stages = [
            ("A1", _stub_step("detect_format")),
            [("B1", _stub_step("deeplocpro"))],
            ("B2", _stub_step("plm_effector")),
            ("C1", _stub_step("cross_validate")),
            [("D1", _stub_step("eggnog"))],
            ("E1", _stub_step("integrate")),
        ]
        fake._build_stages = MagicMock(return_value=stages)

        def _exec(stages_arg, skip_steps):
            call_log.append((sample_id, "_execute_stages", len(stages_arg)))
            return False  # core_failed=False

        fake._execute_stages = MagicMock(side_effect=_exec)
        return fake

    fake_g1 = _make_fake_runner("g1")
    fake_g2 = _make_fake_runner("g2")
    fake_pool = _make_fake_runner(POOL_SAMPLE_ID)

    # Sequence the PipelineRunner constructor: 2 per-genome, then 1 pool.
    instances = iter([fake_g1, fake_g2, fake_pool])
    fake_cls = MagicMock(side_effect=lambda *a, **kw: next(instances))
    monkeypatch.setattr("ssign_app.core.multi_runner.PipelineRunner", fake_cls)

    # Patch pool/split helpers so we can verify they're called at the
    # right boundaries without needing real FASTA/TSV inputs.
    pool_fastas_spy = MagicMock(return_value=0)
    pool_tsvs_spy = MagicMock(return_value=0)
    split_tsv_spy = MagicMock(return_value={})
    monkeypatch.setattr("ssign_app.core.multi_runner.pool_fastas", pool_fastas_spy)
    monkeypatch.setattr("ssign_app.core.multi_runner.pool_tsvs", pool_tsvs_spy)
    monkeypatch.setattr("ssign_app.core.multi_runner.split_tsv_by_source", split_tsv_spy)

    mgr = MultiGenomeRunner(configs, write_combined_summary=False)

    return (
        mgr,
        {"g1": fake_g1, "g2": fake_g2},
        fake_pool,
        call_log,
        {
            "pool_fastas": pool_fastas_spy,
            "pool_tsvs": pool_tsvs_spy,
            "split_tsv_by_source": split_tsv_spy,
        },
    )


def _seed_two_genome_files(per_genome: dict) -> None:
    """Seed neighborhood / substrates / proteins paths so the pooling
    boundaries in _run_multi find something to gather. Used by every
    test in TestMultiGenomeFlow."""
    for sid in ("g1", "g2"):
        per_genome[sid].files["neighborhood_proteins"] = f"/tmp/{sid}_nb.faa"
        per_genome[sid].files["substrates_filtered"] = f"/tmp/{sid}_sub.tsv"
        per_genome[sid].files["proteins"] = f"/tmp/{sid}_p.faa"


class TestMultiGenomeFlow:
    def test_per_genome_runners_execute_segments_a_c_e(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        # Each per-genome runner ran _execute_stages 3 times (segments A, C, E).
        per_sid_calls = {sid: 0 for sid in per_genome}
        for sid, method, _n in log:
            if method == "_execute_stages" and sid in per_genome:
                per_sid_calls[sid] += 1
        assert per_sid_calls == {"g1": 3, "g2": 3}

    def test_pool_runner_executes_segments_b_d(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        pool_calls = [entry for entry in log if entry[0] == POOL_SAMPLE_ID]
        # Pool runner runs once for segment B and once for segment D.
        assert len(pool_calls) == 2

    def test_pool_helpers_called(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        # Pool boundary 1: pool_fastas called for neighborhood
        assert spies["pool_fastas"].called
        # Pool boundary 2: pool_tsvs called for substrates
        assert spies["pool_tsvs"].called

    def test_result_keyed_by_sample_id(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        per_genome["g1"].results = ["g1_result"]
        per_genome["g2"].results = ["g2_result"]

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            out = mgr.run(resume=True)

        assert set(out.keys()) == {"g1", "g2"}
        assert out["g1"] == ["g1_result"]
        assert out["g2"] == ["g2_result"]

    def test_copy_outputs_called_per_genome(self, n2_runner):
        """Regression for #107: multi-genome runs were producing
        combined_results.csv at top level but no per-genome
        ``<sid>_results.csv`` because ``_copy_outputs`` only fires inside
        ``PipelineRunner.run()`` — and multi-genome bypasses that and
        goes straight to ``_execute_stages``. The fix is to call
        ``runner._copy_outputs()`` per genome after segment E.
        """
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        for sid, runner in per_genome.items():
            runner._copy_outputs.assert_called_once()
        # Pool runner has no user-facing outputs — it must NOT be touched.
        pool._copy_outputs.assert_not_called()

    def test_copy_outputs_failure_isolated_per_genome(self, n2_runner):
        """One genome's _copy_outputs raising must not stop the others
        from getting theirs written. Mirrors the
        ``_write_combined_summary`` skip-with-warning posture."""
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        per_genome["g1"]._copy_outputs.side_effect = RuntimeError("boom")

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        per_genome["g1"]._copy_outputs.assert_called_once()
        per_genome["g2"]._copy_outputs.assert_called_once()

    def test_each_runner_builds_its_own_stages(self, n2_runner):
        """Regression for the 2026-06-05 bound-method bug.

        ``_build_stages`` tuples hold bound methods; reusing one
        runner's slice across others silently runs the first runner's
        methods (manifested as 4 K-12 copies in a 4-genome batched CX3
        run). Each runner must build its own — assert _build_stages
        was hit on every runner instance individually.
        """
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        for sid, runner in per_genome.items():
            assert runner._build_stages.call_count >= 3, (
                f"runner {sid} _build_stages called "
                f"{runner._build_stages.call_count}x; expected >= 3 (segments A, C, E)"
            )
        assert pool._build_stages.call_count >= 2, (
            f"pool runner _build_stages called {pool._build_stages.call_count}x; expected >= 2 (segments B, D)"
        )


class TestT5assWholePoolPass:
    """The pool t5ass_whole pass between segments D and E runs ONLY when
    at least one config has --t5ass-annotate-whole on. When it runs, the
    per-genome t5ass_whole step is skipped via skip_steps so it doesn't
    overwrite the split paths.
    """

    def test_skipped_when_no_config_enables_whole(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)
        # Default configs all have t5ass_annotate_whole=False.
        mgr._pool_t5ass_whole_inputs = MagicMock()
        mgr._run_pool_t5ass_whole_segment = MagicMock()
        mgr._split_t5ass_whole_outputs = MagicMock()

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        mgr._pool_t5ass_whole_inputs.assert_not_called()
        mgr._run_pool_t5ass_whole_segment.assert_not_called()
        mgr._split_t5ass_whole_outputs.assert_not_called()

    def test_inputs_called_when_flag_on(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)
        # Toggle the flag on one config — that's enough to trigger the pool path.
        mgr.configs[0].t5ass_annotate_whole = True

        # Simulate "no genome has clean T5aSS substrates" via the inputs
        # method returning False; the downstream methods then stay quiet.
        mgr._pool_t5ass_whole_inputs = MagicMock(return_value=False)
        mgr._run_pool_t5ass_whole_segment = MagicMock()
        mgr._split_t5ass_whole_outputs = MagicMock()

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        mgr._pool_t5ass_whole_inputs.assert_called_once()
        mgr._run_pool_t5ass_whole_segment.assert_not_called()
        mgr._split_t5ass_whole_outputs.assert_not_called()

    def test_full_pool_pass_fires_when_inputs_pool_succeeds(self, n2_runner):
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)
        mgr.configs[0].t5ass_annotate_whole = True
        mgr.configs[1].t5ass_annotate_whole = True

        mgr._pool_t5ass_whole_inputs = MagicMock(return_value=True)
        mgr._run_pool_t5ass_whole_segment = MagicMock()
        mgr._split_t5ass_whole_outputs = MagicMock()

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        mgr._pool_t5ass_whole_inputs.assert_called_once()
        mgr._run_pool_t5ass_whole_segment.assert_called_once()
        mgr._split_t5ass_whole_outputs.assert_called_once()

    def test_per_genome_segment_e_skips_t5ass_whole_when_pool_ran(self, n2_runner):
        """When the pool ran the second pass, per-genome _execute_stages
        for segment E must be invoked with the t5ass_whole_annotations
        step in its skip_steps set so it doesn't re-do (and clobber) the
        split paths.
        """
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)
        mgr.configs[0].t5ass_annotate_whole = True
        mgr.configs[1].t5ass_annotate_whole = True

        mgr._pool_t5ass_whole_inputs = MagicMock(return_value=True)
        mgr._run_pool_t5ass_whole_segment = MagicMock()
        mgr._split_t5ass_whole_outputs = MagicMock()

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        # Look at the per-genome _execute_stages calls. The third call
        # per runner is segment E; its skip_steps kwarg must include
        # t5ass_whole_annotations.
        for sid, runner in per_genome.items():
            e_calls = runner._execute_stages.call_args_list[-1]
            skip = e_calls.kwargs.get("skip_steps") or (e_calls.args[1] if len(e_calls.args) > 1 else None)
            assert skip is not None, f"runner {sid} segment-E call had no skip_steps"
            assert "t5ass_whole_annotations" in skip, (
                f"runner {sid} segment-E skip_steps missing t5ass_whole_annotations: {skip}"
            )

    def test_per_genome_segment_e_no_skip_when_flag_off(self, n2_runner):
        """When the flag is off, segment E runs unrestricted (empty skip)."""
        mgr, per_genome, pool, log, spies = n2_runner
        _seed_two_genome_files(per_genome)

        with patch("ssign_app.scripts.ssign_lib.substrates.load_substrate_ids", return_value=set()):
            mgr.run(resume=True)

        for sid, runner in per_genome.items():
            e_calls = runner._execute_stages.call_args_list[-1]
            skip = e_calls.kwargs.get("skip_steps") or (e_calls.args[1] if len(e_calls.args) > 1 else None)
            # skip can be an empty set or None — either way, t5ass_whole shouldn't be in it.
            if skip:
                assert "t5ass_whole_annotations" not in skip
