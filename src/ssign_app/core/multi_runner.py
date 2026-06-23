"""MultiGenomeRunner: orchestrate ssign for N>1 genomes in one invocation.

Pools predictions over neighborhoods and annotations over substrates so
the heavy startup tax (InterProScan JVM, EggNOG DIAMOND, pLM-BLAST
embeddings, PLM-Effector models) is paid once per batch instead of once
per genome.

Five-segment design: A (per-genome detect → neighborhood) → B (pooled
predictions) → C (per-genome cross-validate → filtering) → D (pooled
annotations) → E (per-genome integrate → report).

For N=1, ``run()`` delegates straight to ``PipelineRunner`` so the
single-genome path stays bit-identical to today.
"""

from __future__ import annotations

import logging
import os
import tempfile
from dataclasses import replace
from pathlib import Path
from typing import Callable, Literal, Optional

SegmentLetter = Literal["A", "B", "C", "D", "E"]

from ssign_app.core._pool_utils import (
    pool_fastas,
    pool_tsvs,
    split_tsv_by_source,
    validate_sample_id,
)
from ssign_app.core.runner import PipelineConfig, PipelineRunner, StepResult

logger = logging.getLogger(__name__)

POOL_SAMPLE_ID = "_pool"

_SEGMENT_BY_STEP = {
    "detect_format": "A",
    "extract_proteins": "A",
    "macsyfinder": "A",
    "validate_systems": "A",
    "extract_neighborhood": "A",
    "deeplocpro": "B",
    "deepsece": "B",
    "signalp": "B",
    "plm_effector": "B",
    "cross_validate": "C",
    "proximity": "C",
    "t5ss": "C",
    "filtering": "C",
    "build_passenger_fasta": "C",
    "blastp": "D",
    "hhsuite": "D",
    "interproscan": "D",
    "eggnog": "D",
    "plm_blast": "D",
    "protparam": "D",
    "t5ass_whole_annotations": "E",
    "integrate": "E",
    "orthologs": "E",
    "enrichment": "E",
    "report": "E",
    "figures": "E",
    "enrichment_figure": "E",
}

# Output keys (in PipelineRunner.files) for each pooled segment's tools.
# Used at split time to wire pool TSVs back into per-genome runners.
_SEGMENT_B_OUTPUT_KEYS = ("deeplocpro", "deepsece", "signalp", "plm_effector")
_SEGMENT_D_OUTPUT_KEYS = (
    "blastp",
    "hhsuite",
    "interproscan",
    "eggnog",
    "plm_blast",
    "protparam",
)

# Skip-flag fields whose values must agree across all genome configs in
# a multi-genome run. Heterogeneous flags would require per-genome pool
# routing, which is deferred.
_REQUIRED_UNIFORM_FLAGS = (
    "skip_deeplocpro",
    "skip_deepsece",
    "skip_signalp",
    "skip_plm_effector",
    "skip_blastp",
    "skip_hhsuite",
    "skip_interproscan",
    "skip_eggnog",
    "skip_plmblast",
    "skip_protparam",
    "enrichment_stats",
)


def _step_id(func) -> str:
    return func.__name__.replace("_step_", "")


def slice_stages_by_segment(stages: list) -> dict[str, list]:
    """Partition a ``PipelineRunner._build_stages()`` list into A..E segments.

    Each stage is either a sequential ``(name, callable)`` tuple or a
    parallel-group list of such tuples. The segment a stage belongs to
    is determined by the step_id of its callable (single step) or its
    first member (parallel group).
    """
    segments: dict[str, list] = {seg: [] for seg in "ABCDE"}
    for stage in stages:
        if isinstance(stage, list):
            if not stage:
                continue
            first_step_id = _step_id(stage[0][1])
            segments[_SEGMENT_BY_STEP[first_step_id]].append(stage)
            continue
        _, func = stage
        segments[_SEGMENT_BY_STEP[_step_id(func)]].append(stage)
    return segments


class MultiGenomeRunner:
    """Orchestrate ssign for one or more genomes in a single process.

    For N=1, ``run()`` is functionally identical to
    ``PipelineRunner(configs[0]).run()``. For N>1, prediction tools
    (DLP/DSE/SignalP/PLM-E) run once on a pooled neighborhood FASTA
    and annotation tools (BLASTp/EggNOG/IPS/pLM-BLAST/HHsuite/ProtParam)
    run once on a pooled substrate set, while context-dependent steps
    (cross_validate, proximity, T5SS handling, filtering, integrate,
    report, figures) still run per-genome.
    """

    def __init__(
        self,
        configs: list[PipelineConfig],
        progress_callback: Optional[Callable] = None,
        api_semaphores: Optional[dict] = None,
        write_combined_summary: bool = True,
    ):
        if not configs:
            raise ValueError("MultiGenomeRunner needs at least one PipelineConfig")
        ids = [c.sample_id for c in configs]
        if len(set(ids)) != len(ids):
            raise ValueError(f"sample_ids must be unique across configs; got {ids}")
        if POOL_SAMPLE_ID in ids:
            raise ValueError(f"sample_id {POOL_SAMPLE_ID!r} is reserved for the pool runner")
        for c in configs:
            validate_sample_id(c.sample_id)
        if len(configs) > 1:
            ref = configs[0]
            for c in configs[1:]:
                for flag in _REQUIRED_UNIFORM_FLAGS:
                    if getattr(c, flag) != getattr(ref, flag):
                        raise ValueError(
                            f"Configs disagree on {flag}: multi-genome v1 requires "
                            f"uniform skip flags across all genomes."
                        )

        self.configs = list(configs)
        self._progress_callback = progress_callback
        self.api_semaphores = api_semaphores
        self.write_combined_summary = write_combined_summary
        self.results: dict[str, list[StepResult]] = {}

    def run(self, resume: bool = True) -> dict[str, list[StepResult]]:
        if len(self.configs) == 1:
            runner = PipelineRunner(
                self.configs[0],
                progress_callback=self._progress_callback,
                api_semaphores=self.api_semaphores,
            )
            self.results = {self.configs[0].sample_id: runner.run(resume=resume)}
            return self.results
        return self._run_multi()

    def _run_multi(self) -> dict[str, list[StepResult]]:
        runners: dict[str, PipelineRunner] = {
            c.sample_id: PipelineRunner(
                c,
                progress_callback=self._progress_callback,
                api_semaphores=self.api_semaphores,
            )
            for c in self.configs
        }

        # Top-level outdir is the shared parent of all per-genome outdirs.
        # Pool runner artefacts land in <top_outdir>/_pool/.
        top_outdir = Path(os.path.commonpath([c.outdir for c in self.configs])).resolve()
        pool_outdir = top_outdir / "_pool"
        pool_outdir.mkdir(parents=True, exist_ok=True)

        # === Segment A (per-genome) ===
        self._run_per_genome_segment(runners, "A")

        # === Pool boundary 1: pool neighborhoods ===
        pool_runner = self._make_pool_runner(top_outdir)
        self._pool_segment_b_inputs(runners, pool_runner, pool_outdir)

        # === Segment B (pooled) ===
        # split is unconditional: if segment B was empty (all predictions
        # skipped), pool_runner.files won't have any segment-B keys, and
        # _split_pooled_outputs no-ops over the empty key set.
        self._run_pool_segment(pool_runner, "B")
        self._split_pooled_outputs(runners, pool_runner, pool_outdir, _SEGMENT_B_OUTPUT_KEYS)

        # === Segment C (per-genome) ===
        self._run_per_genome_segment(runners, "C")

        # === Pool boundary 2: pool substrates + substrate-only proteomes ===
        self._pool_segment_d_inputs(runners, pool_runner, pool_outdir)

        # === Segment D (pooled) ===
        # Same no-op-on-empty rule as the segment B split above.
        self._run_pool_segment(pool_runner, "D")
        self._split_pooled_outputs(runners, pool_runner, pool_outdir, _SEGMENT_D_OUTPUT_KEYS)

        # === Optional pool t5ass_whole pass (between D and E) ===
        # When any config has --t5ass-annotate-whole on and at least one
        # genome has clean T5aSS substrates, the second annotation pass
        # runs ONCE on the pool instead of N times per-genome. The
        # per-genome t5ass_whole step in segment E is then skipped via
        # the skip_steps argument below — it would otherwise re-do the
        # work for each genome AND overwrite the split paths we just
        # wired into runners[sid].files["t5ass_whole_<tool>"].
        skip_e_steps: set[str] = set()
        if any(c.t5ass_annotate_whole for c in self.configs):
            if self._pool_t5ass_whole_inputs(runners, pool_runner, pool_outdir):
                self._run_pool_t5ass_whole_segment(pool_runner, pool_outdir)
                self._split_t5ass_whole_outputs(runners, pool_runner, pool_outdir)
                skip_e_steps.add("t5ass_whole_annotations")

        # === Segment E (per-genome) ===
        self._run_per_genome_segment(runners, "E", skip_steps=skip_e_steps)
        self._copy_per_genome_outputs(runners)

        if self.write_combined_summary:
            self._write_combined_summary(runners, top_outdir)

        # Combined cross-genome enrichment view (pooled table + figure) when the
        # enrichment test ran. Each genome's _enrichment_stats.tsv + sibling .npz
        # come straight from its runner.files.
        if any(c.enrichment_stats for c in self.configs):
            from ssign_app.core.runner import pool_and_plot_enrichment

            tsvs = [r.files.get("enrichment_stats", "") for r in runners.values()]
            try:
                n = pool_and_plot_enrichment(tsvs, str(top_outdir))
                if n:
                    logger.info("pooled enrichment: %d cross-genome rows + figure -> %s", n, top_outdir)
            except Exception as e:
                logger.warning("pooled enrichment failed: %s", str(e)[:160])

        # Pooled cross-genome regular figures (P01-P03) from every genome's
        # integrated CSV. generate_figures.py's pooled mode no-ops for <2 genomes;
        # gate here too to skip a needless subprocess.
        integrated_csvs = [r.files.get("integrated", "") for r in runners.values()]
        integrated_csvs = [p for p in integrated_csvs if p and os.path.exists(p)]
        if len(integrated_csvs) >= 2:
            from ssign_app.core.runner import run_script

            pooled_fig_dir = os.path.join(str(top_outdir), "figures")
            os.makedirs(pooled_fig_dir, exist_ok=True)
            rc, _out, err = run_script(
                "generate_figures.py",
                [
                    "--master-csvs",
                    *integrated_csvs,
                    "--outdir",
                    pooled_fig_dir,
                    "--mode",
                    "pooled",
                    "--dpi",
                    str(self.configs[0].dpi),
                ],
            )
            if rc == 0:
                logger.info("pooled figures -> %s", pooled_fig_dir)
            else:
                logger.warning("pooled figures failed: %s", err[:160])

        self.results = {sid: r.results for sid, r in runners.items()}
        return self.results

    def _make_pool_runner(self, top_outdir: Path) -> PipelineRunner:
        """Build a PipelineRunner whose work_dir holds pooled artefacts."""
        pool_config = replace(
            self.configs[0],
            sample_id=POOL_SAMPLE_ID,
            outdir=str(top_outdir / "_pool"),
            # No input genome; we seed self.files manually before each segment.
            input_path="",
            # Avoid a second concurrent sampler thread — per-genome runners
            # could each have their own, multiplying the I/O cost.
            monitor_resources=False,
        )
        pool_runner = PipelineRunner(
            pool_config,
            progress_callback=self._progress_callback,
            api_semaphores=self.api_semaphores,
        )
        pool_runner.work_dir = tempfile.mkdtemp(prefix="ssign_pool_")
        os.makedirs(pool_runner.config.outdir, exist_ok=True)
        return pool_runner

    def _stages_for_segment(self, runner: PipelineRunner, segment_letter: SegmentLetter) -> list:
        """Slice the named segment from this runner's own ``_build_stages``.

        Stages must come from the executing runner because the
        ``(name, callable)`` tuples hold bound methods — reusing one
        runner's slice across others silently executes the first
        runner's methods (caught 2026-06-05: a 4-genome batched CX3 run
        produced four copies of K-12's results because every runner
        re-ran ref._step_extract_proteins against ref.config).
        """
        return slice_stages_by_segment(runner._build_stages())[segment_letter]

    def _run_per_genome_segment(
        self,
        runners: dict[str, PipelineRunner],
        segment_letter: SegmentLetter,
        skip_steps: Optional[set[str]] = None,
    ) -> None:
        """Run the named segment for each per-genome runner.

        ``skip_steps`` is forwarded to ``_execute_stages`` so callers can
        suppress steps that have already run in a multi-genome pool
        pass — e.g. ``t5ass_whole_annotations`` runs once on the pool
        between segments D and E and should not re-fire per-genome.
        """
        skip = set(skip_steps) if skip_steps else set()
        for sid, runner in runners.items():
            if not runner.work_dir:
                runner.work_dir = tempfile.mkdtemp(prefix=f"ssign_{sid}_")
                os.makedirs(runner.config.outdir, exist_ok=True)
            stages = self._stages_for_segment(runner, segment_letter)
            if not stages:
                continue
            runner._execute_stages(stages, skip_steps=skip)

    def _run_pool_segment(self, pool_runner: PipelineRunner, segment_letter: SegmentLetter) -> None:
        stages = self._stages_for_segment(pool_runner, segment_letter)
        if not stages:
            return
        pool_runner._execute_stages(stages, skip_steps=set())

    def _pool_segment_b_inputs(
        self,
        runners: dict[str, PipelineRunner],
        pool_runner: PipelineRunner,
        pool_outdir: Path,
    ) -> None:
        sources = [
            (sid, Path(r.files["neighborhood_proteins"]))
            for sid, r in runners.items()
            if "neighborhood_proteins" in r.files
        ]
        if not sources:
            return
        pooled = pool_outdir / "pooled_neighborhood.faa"
        n = pool_fastas(sources, pooled)
        logger.info(
            "Segment B input: pooled %d neighborhood proteins from %d genomes -> %s",
            n,
            len(sources),
            pooled,
        )
        pool_runner.files["neighborhood_proteins"] = str(pooled)
        # Default `proteins` to the neighborhood pool; the whole-genome block below
        # overrides it for DLP/DSE when enrichment needs the full proteome.
        pool_runner.files["proteins"] = str(pooled)

        # --enrichment-stats forces whole-genome DLP/DSE (Config.__post_init__) so
        # the circular-shift null's background is genome-wide. In batched mode the
        # pooled DLP/DSE run must therefore see the FULL proteome, not just the
        # neighborhood: _resolve_step_input_fasta(whole_genome=True) returns
        # `proteins`, so we pool the whole proteomes and stage them under that key.
        # SignalP/PLM-E keep `neighborhood_proteins` (whole_genome=False), matching
        # single-genome scope. Without this the background is computed from
        # neighborhood-only positivity and the null is ~7x too low (fold/p inflated).
        if self.configs[0].dlp_whole_genome or self.configs[0].dse_whole_genome:
            wg_sources = [(sid, Path(r.files["proteins"])) for sid, r in runners.items() if "proteins" in r.files]
            if wg_sources:
                pooled_wg = pool_outdir / "pooled_whole_genome.faa"
                n_wg = pool_fastas(wg_sources, pooled_wg)
                pool_runner.files["proteins"] = str(pooled_wg)
                logger.info(
                    "Segment B: pooled %d whole-genome proteins for DLP/DSE (enrichment) -> %s",
                    n_wg,
                    pooled_wg,
                )

    def _split_pooled_outputs(
        self,
        runners: dict[str, PipelineRunner],
        pool_runner: PipelineRunner,
        pool_outdir: Path,
        output_keys: tuple[str, ...],
    ) -> None:
        """Demux each pooled tool TSV into per-genome split files and rewire
        ``runners[sid].files[key]`` to point at them."""
        for key in output_keys:
            pooled = pool_runner.files.get(key)
            if not pooled or not os.path.exists(pooled):
                continue
            id_col = "seq_id" if key == "plm_effector" else "locus_tag"
            split_dir = pool_outdir / f"split_{key}"
            paths = split_tsv_by_source(Path(pooled), split_dir, id_column=id_col)
            for sid, path in paths.items():
                if sid in runners:
                    runners[sid].files[key] = str(path)

    def _pool_segment_d_inputs(
        self,
        runners: dict[str, PipelineRunner],
        pool_runner: PipelineRunner,
        pool_outdir: Path,
    ) -> None:
        from ssign_app.scripts.ssign_lib.substrates import (
            load_substrate_ids,
            write_substrates_only_fasta,
        )

        substrates_sources = [
            (sid, Path(r.files["substrates_filtered"]))
            for sid, r in runners.items()
            if "substrates_filtered" in r.files
        ]
        if not substrates_sources:
            return

        pooled_substrates_tsv = pool_outdir / "pooled_substrates_filtered.tsv"
        n_rows = pool_tsvs(substrates_sources, pooled_substrates_tsv)
        pool_runner.files["substrates_filtered"] = str(pooled_substrates_tsv)
        logger.info(
            "Segment D input: pooled %d substrate rows from %d genomes -> %s",
            n_rows,
            len(substrates_sources),
            pooled_substrates_tsv,
        )

        # Per-genome substrate-only FASTAs (pre-filtered so the pooled
        # FASTA is small), then pool. Build TWO pools:
        #   - "proteins"                       full-protein sequences (IPS reads this)
        #   - "proteins_for_passenger_tools"   T5aSS substrates carry passenger seqs;
        #                                      everything else carries the full seq
        # The same routing fires inside the pool runner via
        # PipelineRunner._annotation_input_proteins.
        substrate_fasta_sources = []
        passenger_fasta_sources = []
        for sid, r in runners.items():
            sf = r.files.get("substrates_filtered")
            proteins = r.files.get("proteins")
            if not sf or not proteins:
                continue
            sub_ids = load_substrate_ids(sf)
            if not sub_ids:
                continue
            per_genome_fasta = Path(r.work_dir) / f"{sid}_substrate_proteins.faa"
            write_substrates_only_fasta(proteins, sub_ids, str(per_genome_fasta))
            substrate_fasta_sources.append((sid, per_genome_fasta))

            # Passenger-substituted substrate-only FASTA. Falls back to the
            # full substrate FASTA when build_passenger_fasta wasn't staged.
            passenger_proteins = r.files.get("proteins_for_passenger_tools") or proteins
            per_genome_pfasta = Path(r.work_dir) / f"{sid}_substrate_proteins_passenger.faa"
            write_substrates_only_fasta(passenger_proteins, sub_ids, str(per_genome_pfasta))
            passenger_fasta_sources.append((sid, per_genome_pfasta))

        if substrate_fasta_sources:
            pooled_proteins = pool_outdir / "pooled_substrate_proteins.faa"
            pool_fastas(substrate_fasta_sources, pooled_proteins)
            pool_runner.files["proteins"] = str(pooled_proteins)

        if passenger_fasta_sources:
            pooled_passenger = pool_outdir / "pooled_substrate_proteins_passenger.faa"
            pool_fastas(passenger_fasta_sources, pooled_passenger)
            pool_runner.files["proteins_for_passenger_tools"] = str(pooled_passenger)

    def _pool_t5ass_whole_inputs(
        self,
        runners: dict[str, PipelineRunner],
        pool_runner: PipelineRunner,
        pool_outdir: Path,
    ) -> bool:
        """Pool T5aSS-only-clean substrates + full proteins across genomes.

        Returns True if at least one genome contributed substrates (and
        the pool t5ass_whole pass should fire), False otherwise.
        Mirrors ``_pool_segment_d_inputs`` but on the T5aSS subset and
        with the FULL protein sequence (not passenger-substituted) so
        the second-pass tools see the whole autotransporter.
        """
        from ssign_app.scripts.ssign_lib.substrates import (
            load_substrate_ids,
            write_substrates_only_fasta,
        )

        t5_substrates_sources: list[tuple[str, Path]] = []
        t5_full_fasta_sources: list[tuple[str, Path]] = []
        for sid, r in runners.items():
            t5_sub_path = r._build_t5ass_only_substrates_tsv()
            if t5_sub_path is None:
                continue
            t5_substrates_sources.append((sid, Path(t5_sub_path)))

            proteins = r.files.get("proteins")
            if not proteins:
                continue
            t5_ids = load_substrate_ids(t5_sub_path)
            if not t5_ids:
                continue
            per_genome_fasta = Path(r.work_dir) / f"{sid}_t5ass_substrate_proteins_full.faa"
            write_substrates_only_fasta(proteins, t5_ids, str(per_genome_fasta))
            t5_full_fasta_sources.append((sid, per_genome_fasta))

        if not t5_substrates_sources:
            return False

        pooled_t5_subs = pool_outdir / "pooled_t5ass_substrates.tsv"
        pool_tsvs(t5_substrates_sources, pooled_t5_subs)

        pooled_t5_fasta = pool_outdir / "pooled_t5ass_substrate_proteins_full.faa"
        if t5_full_fasta_sources:
            pool_fastas(t5_full_fasta_sources, pooled_t5_fasta)

        # Stash on pool_runner under disambiguated keys; the actual
        # swap into substrates_filtered + proteins happens inside
        # _run_pool_t5ass_whole_segment so the segment-D outputs aren't
        # clobbered before split.
        pool_runner.files["t5ass_only_substrates"] = str(pooled_t5_subs)
        if t5_full_fasta_sources:
            pool_runner.files["t5ass_full_proteins"] = str(pooled_t5_fasta)
        return True

    def _run_pool_t5ass_whole_segment(
        self,
        pool_runner: PipelineRunner,
        pool_outdir: Path,
    ) -> None:
        """Run the 5 routed annotation tools once on the pooled T5aSS inputs.

        Same hijack pattern as the single-genome
        ``_step_t5ass_whole_annotations``: swap work_dir + substrates
        + drop the passenger FASTA so the step methods see the full
        pooled FASTA, run each tool, column-prefix the output, store
        under ``pool_runner.files["t5ass_whole_<tool>"]``. Restores
        the segment-D state in ``finally``.
        """
        import shutil

        from ssign_app.core.runner import _rename_csv_columns_with_prefix

        t5_sub_path = pool_runner.files.get("t5ass_only_substrates", "")
        t5_full_fasta = pool_runner.files.get("t5ass_full_proteins", "")
        if not t5_sub_path or not t5_full_fasta:
            return

        tools = PipelineRunner._T5ASS_WHOLE_TOOLS
        second_work_dir = tempfile.mkdtemp(prefix="ssign_pool_t5ass_whole_")
        saved_work_dir = pool_runner.work_dir
        saved_substrates = pool_runner.files.get("substrates_filtered")
        saved_proteins = pool_runner.files.get("proteins")
        saved_passenger = pool_runner.files.get("proteins_for_passenger_tools")
        saved_outputs: dict[str, object] = {tool: pool_runner.files.get(tool) for tool in tools}

        pool_runner.work_dir = second_work_dir
        pool_runner.files["substrates_filtered"] = t5_sub_path
        pool_runner.files["proteins"] = t5_full_fasta
        pool_runner.files.pop("proteins_for_passenger_tools", None)

        try:
            for tool in tools:
                if saved_outputs[tool] is None:
                    continue  # main pool pass didn't run this tool
                try:
                    step_method = getattr(pool_runner, f"_step_{tool}")
                    result = step_method()
                    if not result.success:
                        logger.warning("pool t5ass_whole %s failed: %s", tool, result.message[:120])
                        continue
                    tool_output = pool_runner.files.get(tool, "")
                    if not tool_output or not os.path.exists(tool_output):
                        continue
                    ext = os.path.splitext(tool_output)[1]
                    final_path = pool_outdir / f"_pool_t5ass_whole_{tool}{ext}"
                    _rename_csv_columns_with_prefix(
                        tool_output, str(final_path), prefix="t5ass_whole_", keep_columns={"locus_tag"}
                    )
                    saved_outputs[tool] = (saved_outputs[tool], str(final_path))
                except Exception as e:
                    logger.warning("pool t5ass_whole %s raised: %s", tool, str(e)[:120])
        finally:
            pool_runner.work_dir = saved_work_dir
            if saved_substrates is not None:
                pool_runner.files["substrates_filtered"] = saved_substrates
            if saved_proteins is not None:
                pool_runner.files["proteins"] = saved_proteins
            if saved_passenger is not None:
                pool_runner.files["proteins_for_passenger_tools"] = saved_passenger
            for tool, info in saved_outputs.items():
                if isinstance(info, tuple):
                    main_path, t5_path = info
                    pool_runner.files[tool] = main_path
                    pool_runner.files[f"t5ass_whole_{tool}"] = t5_path
                elif info is not None:
                    pool_runner.files[tool] = info
            try:
                shutil.rmtree(second_work_dir)
            except OSError:
                pass

    def _split_t5ass_whole_outputs(
        self,
        runners: dict[str, PipelineRunner],
        pool_runner: PipelineRunner,
        pool_outdir: Path,
    ) -> None:
        """Demux each ``t5ass_whole_<tool>`` pool file into per-genome splits."""
        for tool in PipelineRunner._T5ASS_WHOLE_TOOLS:
            key = f"t5ass_whole_{tool}"
            pooled = pool_runner.files.get(key)
            if not pooled or not os.path.exists(pooled):
                continue
            split_dir = pool_outdir / f"split_{key}"
            paths = split_tsv_by_source(Path(pooled), split_dir, id_column="locus_tag")
            for sid, path in paths.items():
                if sid in runners:
                    runners[sid].files[key] = str(path)

    def _copy_per_genome_outputs(self, runners: dict[str, PipelineRunner]) -> None:
        """Materialise per-genome user-facing files.

        Single-genome runs get this via ``PipelineRunner.run()``;
        multi-genome bypasses ``run()`` and drives ``_execute_stages``
        directly, so the finaliser has to be invoked explicitly. One
        genome failing shouldn't block the others.
        """
        for sid, runner in runners.items():
            try:
                runner._copy_outputs()
            except Exception as e:
                logger.warning("copy_outputs failed for genome %s: %s", sid, str(e)[:120])

    def _write_combined_summary(self, runners: dict[str, PipelineRunner], top_outdir: Path) -> None:
        """Concatenate per-genome master CSVs into ``combined_results.csv``.

        Each row gets a ``source_genome`` column tagging which genome the
        substrate came from. Genomes whose ``integrated`` CSV is missing
        (e.g. segment-E failed for that genome) are skipped with a warning.

        Matches the single-genome ``<sample>_results.csv`` schema +
        delimiter so downstream consumers can use the same parser. Method
        name and ``write_combined_summary`` flag still say "summary" for
        backward compat with callers; the file on disk is the results table.
        """
        import csv

        out_path = top_outdir / "combined_results.csv"
        all_fields: list[str] = ["source_genome"]
        seen: set[str] = {"source_genome"}
        rows: list[dict[str, str]] = []

        for sid, runner in runners.items():
            integrated = runner.files.get("integrated")
            if not integrated or not os.path.exists(integrated):
                logger.warning("combined_results: no integrated CSV for genome %s", sid)
                continue
            with open(integrated) as f:
                # integrated_annotations.py writes a comma-separated CSV.
                reader = csv.DictReader(f)
                for col in reader.fieldnames or []:
                    if col not in seen:
                        all_fields.append(col)
                        seen.add(col)
                for row in reader:
                    row["source_genome"] = sid
                    rows.append(row)

        if not rows:
            logger.warning("combined_results: no rows to write")
            return

        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_fields)
            writer.writeheader()
            writer.writerows(rows)
        logger.info(
            "combined_results: wrote %d rows across %d genomes -> %s",
            len(rows),
            len(runners),
            out_path,
        )
