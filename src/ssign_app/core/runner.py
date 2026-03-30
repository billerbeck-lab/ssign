"""Pipeline runner — orchestrates ssign steps without Nextflow.

For ssign-lite, this runs the Python scripts directly in sequence.
Each step is a function that returns (success, message) and updates
a callback with progress info.
"""

import csv
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# The bin/ directory containing all pipeline scripts
# Scripts are packaged inside ssign_app/scripts/ (installed with pip)
# Fallback to bin/ for development mode
_PACKAGE_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
_DEV_BIN = Path(__file__).resolve().parents[3] / "bin"
BIN_DIR = _PACKAGE_SCRIPTS if _PACKAGE_SCRIPTS.exists() else _DEV_BIN


@dataclass
class PipelineConfig:
    """All configuration for a pipeline run."""
    # Input
    input_path: str = ""
    original_filename: str = ""  # Original filename when input is a temp upload
    sample_id: str = ""
    outdir: str = "./results"

    # Phase 2: SS detection
    wholeness_threshold: float = 0.8
    excluded_systems: list = field(default_factory=lambda: ["Flagellum", "Tad", "T3SS"])

    # Phase 3: Prediction
    conf_threshold: float = 0.8
    proximity_window: int = 3
    required_fraction_correct: float = 0.8

    # Phase 1: ORF prediction options (for FASTA contigs input)
    run_bakta: bool = False     # Use Bakta instead of Prodigal for ORF prediction
    bakta_db: str = ""          # Path to Bakta database (required if run_bakta=True)
    bakta_threads: int = 4

    # Phase 3: Tool paths (DTU licensed)
    deeplocpro_mode: str = "remote"  # "local" or "remote"
    deeplocpro_path: str = ""
    signalp_mode: str = "remote"
    signalp_path: str = ""
    skip_signalp: bool = False
    skip_deepsece: bool = False

    # Phase 5: Annotation tools
    skip_blastp: bool = False
    blastp_mode: str = "remote"
    blastp_db: str = ""
    blastp_exclude_taxid: str = ""
    blastp_min_pident: float = 80.0
    blastp_min_qcov: float = 80.0
    blastp_evalue: float = 1e-5

    skip_hhsuite: bool = True
    hhsuite_mode: str = "remote"
    hhsuite_pfam_db: str = ""
    hhsuite_pdb70_db: str = ""
    hhsuite_uniclust_db: str = ""
    hhpred_min_probability: float = 40.0  # HHpred/Pfam min probability to keep hit

    skip_interproscan: bool = False
    interproscan_mode: str = "remote"
    interproscan_db: str = ""
    interproscan_min_evalue: float = 1e-5

    skip_foldseek: bool = True
    foldseek_db: str = ""

    skip_plmblast: bool = True
    plmblast_db: str = ""

    skip_protparam: bool = False
    skip_structure: bool = True

    # DSE type-match filter: remove DSE-only substrates where predicted
    # SS type doesn't match nearby MacSyFinder system
    filter_dse_type_mismatch: bool = True

    # Figures
    dpi: int = 300
    fig_category: bool = True
    fig_ss_comp: bool = True
    fig_tool_heatmap: bool = True
    fig_substrate_count: bool = True
    fig_func_summary: bool = True

    # InterProScan threshold
    interproscan_evalue: float = 1e-5

    # Foldseek thresholds
    foldseek_evalue: float = 1e-3
    foldseek_min_tmscore: float = 0.5

    # DeepSecE threshold
    deepsece_min_prob: float = 0.8

    # SignalP threshold
    signalp_min_prob: float = 0.5

    # Ortholog group thresholds
    ortholog_min_pident: float = 40.0
    ortholog_min_qcov: float = 70.0


@dataclass
class StepResult:
    """Result of a pipeline step."""
    name: str
    success: bool
    message: str
    output_files: dict = field(default_factory=dict)


def run_script(script_name: str, args: list, timeout: int = 7200) -> tuple:
    """Run a bin/ script with arguments. Returns (returncode, stdout, stderr)."""
    script_path = BIN_DIR / script_name
    if not script_path.exists():
        return (-1, "", f"Script not found: {script_path}")

    cmd = [sys.executable, str(script_path)] + args
    logger.info(f"Running: {' '.join(cmd[:4])}...")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
        )
        return (result.returncode, result.stdout, result.stderr)
    except subprocess.TimeoutExpired:
        return (-1, "", f"Timeout after {timeout}s")
    except Exception as e:
        return (-1, "", str(e))


class PipelineRunner:
    """Runs the ssign pipeline step by step with progress callbacks."""

    def __init__(self, config: PipelineConfig,
                 progress_callback: Optional[Callable] = None):
        self.config = config
        self.progress = progress_callback or (lambda step, pct, msg: None)
        self.results: list[StepResult] = []
        self.work_dir = ""
        self.files = {}  # Track intermediate file paths
        self.start_time: float | None = None

    def check_dependencies(self) -> list[str]:
        """Pre-flight check for required and optional dependencies.

        Returns a list of warning strings. Does NOT block the pipeline —
        individual steps produce definitive errors if a tool is missing.

        FRAGILE: This function checks for external tools and packages that
        can change with updates. Each check can fail independently. These
        warnings are informational — the pipeline will still attempt to run.
        """
        warnings = []

        # ── Core dependencies (pipeline cannot run without these) ──

        # hmmsearch: needed by MacSyFinder. Should be provided by our
        # pyhmmer shim (installed as console_script via pip install ssign).
        # FRAGILE: if pip install didn't create the console_script (e.g.,
        # --no-deps, broken PATH, or venv not activated), this will be missing.
        # Fallback: sudo apt install hmmer
        if not shutil.which("hmmsearch"):
            warnings.append(
                "hmmsearch not found on PATH. MacSyFinder needs it.\n"
                "  This should have been installed automatically with ssign "
                "(pyhmmer shim).\n"
                "  How to fix:\n"
                "    - Reinstall: pip install --force-reinstall ssign\n"
                "    - Or install real HMMER: sudo apt install hmmer\n"
                "    - Or: conda install -c bioconda hmmer"
            )

        # macsyfinder: pip-installed, should be on PATH
        if not shutil.which("macsyfinder"):
            warnings.append(
                "macsyfinder not found on PATH.\n"
                "  Install: pip install macsyfinder\n"
                "  It should have been installed as a dependency of ssign."
            )

        # macsydata: needed to install TXSScan models on first run
        if not shutil.which("macsydata"):
            warnings.append(
                "macsydata not found on PATH (needed to install TXSScan "
                "models).\n"
                "  Install: pip install macsyfinder"
            )

        # ── Optional dependencies (checked based on config) ──

        # DeepSecE: large optional dependency (torch + ESM + DeepSecE)
        if not self.config.skip_deepsece:
            for pkg, install in [
                ("torch", "pip install 'ssign[full]'  # or: pip install torch"),
                ("esm", "pip install fair-esm"),
            ]:
                try:
                    __import__(pkg)
                except ImportError:
                    warnings.append(
                        f"Python package '{pkg}' not found (needed for "
                        f"DeepSecE).\n"
                        f"  Install: {install}\n"
                        f"  Or skip DeepSecE in the GUI settings."
                    )

        return warnings

    def _elapsed_str(self) -> str:
        """Format elapsed time since pipeline start."""
        if self.start_time is None:
            return ""
        elapsed = int(time.monotonic() - self.start_time)
        m, s = divmod(elapsed, 60)
        return f"{m}m {s}s" if m else f"{s}s"

    def run(self, resume: bool = True) -> list[StepResult]:
        """Run the full pipeline. Returns list of step results.

        If resume=True (default), attempts to skip steps that completed
        successfully in a previous run (detected via ssign_progress.json
        in the output directory).
        """
        self.start_time = time.monotonic()

        # Pre-flight dependency check
        dep_warnings = self.check_dependencies()
        if dep_warnings:
            for w in dep_warnings:
                logger.warning(f"Dependency check: {w}")
            self.progress(
                "Dependency check", 0,
                f"{len(dep_warnings)} warning(s) — see log for details"
            )

        # Create output directory
        os.makedirs(self.config.outdir, exist_ok=True)

        # Try to resume from previous progress
        skip_steps = set()
        if resume:
            skip_steps = self._try_resume()

        # Only create new work_dir if not resuming
        if not self.work_dir:
            self.work_dir = tempfile.mkdtemp(prefix="ssign_")

        # ── Build pipeline stages ──
        # Stages are either a single step (sequential) or a list of steps
        # that can run in parallel. Parallel steps must write to different
        # keys in self.files and not depend on each other.
        #
        # Dependency graph:
        #   detect_format → extract_proteins → macsyfinder → validate →
        #   extract_neighborhood →
        #     [deeplocpro || deepsece || signalp]  ← PARALLEL GROUP 1
        #   → cross_validate → proximity → t5ss → filtering →
        #     [blastp || hhsuite || interproscan || protparam]  ← PARALLEL GROUP 2
        #   → integrate → orthologs → enrichment → report → figures

        # Prediction tools (parallel group 1)
        prediction_steps = [
            ("Predicting localization (DeepLocPro)", self._step_deeplocpro),
        ]
        if not self.config.skip_deepsece:
            prediction_steps.append(
                ("Predicting secretion type (DeepSecE)", self._step_deepsece))
        if not self.config.skip_signalp:
            prediction_steps.append(
                ("Predicting signal peptides (SignalP)", self._step_signalp))

        # Annotation tools (parallel group 2)
        annotation_steps = []
        if not self.config.skip_blastp:
            annotation_steps.append(("Running BLASTp", self._step_blastp))
        if not self.config.skip_hhsuite:
            annotation_steps.append(("Running HH-suite", self._step_hhsuite))
        if not self.config.skip_interproscan:
            annotation_steps.append(
                ("Running InterProScan", self._step_interproscan))
        if not self.config.skip_protparam:
            annotation_steps.append(
                ("Computing physicochemical properties", self._step_protparam))

        # Full pipeline as stages: each stage is either:
        #   ("name", func)         — single sequential step
        #   [("name", func), ...]  — parallel group
        stages = [
            ("Detecting input format", self._step_detect_format),
            ("Extracting proteins", self._step_extract_proteins),
            ("Running MacSyFinder", self._step_macsyfinder),
            ("Validating secretion systems", self._step_validate_systems),
            ("Extracting SS neighborhood", self._step_extract_neighborhood),
            prediction_steps,  # PARALLEL: deeplocpro + deepsece + signalp
            ("Cross-validating predictions", self._step_cross_validate),
            ("Running proximity analysis", self._step_proximity),
            ("Handling T5SS autotransporters", self._step_t5ss),
            ("Filtering systems", self._step_filtering),
            annotation_steps,  # PARALLEL: blastp + hhsuite + interproscan + protparam
            ("Integrating annotations", self._step_integrate),
            ("Assigning ortholog groups", self._step_orthologs),
            ("Running enrichment analysis", self._step_enrichment),
            ("Generating report", self._step_report),
            ("Generating figures", self._step_figures),
        ]

        # Flatten for counting and core-step tracking
        all_steps = []
        for stage in stages:
            if isinstance(stage, list):
                all_steps.extend(stage)
            else:
                all_steps.append(stage)

        CORE_STEPS = {
            "detect_format", "extract_proteins", "macsyfinder",
            "validate_systems", "extract_neighborhood",
            "deeplocpro", "cross_validate",
            "proximity", "t5ss", "filtering",
            "integrate",
        }

        total = len(all_steps)
        core_failed = False
        n_skipped = 0
        step_counter = 0

        from concurrent.futures import ThreadPoolExecutor, as_completed

        for stage in stages:
            # Normalize: single step becomes a one-element list
            if isinstance(stage, tuple):
                stage_steps = [stage]
                parallel = False
            else:
                stage_steps = stage
                parallel = len(stage_steps) > 1

            if not stage_steps:
                continue

            # Check resume / core-failed for all steps in this stage
            steps_to_run = []
            for name, func in stage_steps:
                step_counter += 1
                step_id = func.__name__.replace("_step_", "")

                if step_id in skip_steps:
                    n_skipped += 1
                    pct = int(100 * step_counter / total)
                    self.progress(name, pct, f"Skipped (already done) | {self._elapsed_str()} elapsed")
                    self.results.append(StepResult(
                        step_id, True, "Resumed (already completed)"))
                    continue

                if core_failed and step_id in CORE_STEPS:
                    self.results.append(StepResult(
                        step_id, False, "Skipped (earlier step failed)"))
                    continue

                steps_to_run.append((name, func, step_id, step_counter))

            if not steps_to_run:
                continue

            # Run steps — in parallel if multiple, sequential if single
            if parallel and len(steps_to_run) > 1:
                step_names = ", ".join(n for n, _, _, _ in steps_to_run)
                pct = int(100 * steps_to_run[0][3] / total)
                self.progress(
                    f"Running in parallel: {step_names}",
                    pct,
                    f"{len(steps_to_run)} tools running simultaneously | {self._elapsed_str()} elapsed",
                )
                print(
                    f"[ssign] Starting parallel group: {step_names}",
                    flush=True,
                )

                with ThreadPoolExecutor(
                    max_workers=len(steps_to_run)
                ) as executor:
                    futures = {}
                    for name, func, step_id, sc in steps_to_run:
                        futures[executor.submit(func)] = (
                            name, step_id, sc)

                    for future in as_completed(futures):
                        name, step_id, sc = futures[future]
                        try:
                            result = future.result()
                            self.results.append(result)
                            print(
                                f"[ssign] Finished (parallel): {name} -> "
                                f"{'OK' if result.success else 'FAILED: ' + result.message[:100]}",
                                flush=True,
                            )
                            if not result.success:
                                logger.error(
                                    f"Step '{name}' failed: "
                                    f"{result.message}"
                                )
                                if step_id in CORE_STEPS:
                                    core_failed = True
                        except Exception as e:
                            print(
                                f"[ssign] EXCEPTION (parallel): "
                                f"{name} -> {e}",
                                flush=True,
                            )
                            self.results.append(
                                StepResult(step_id, False, str(e)))
                            logger.exception(
                                f"Step '{name}' raised exception")
                            if step_id in CORE_STEPS:
                                core_failed = True

                self._save_progress()
            else:
                # Sequential execution (single step or single-element group)
                for name, func, step_id, sc in steps_to_run:
                    pct = int(100 * sc / total)
                    self.progress(name, pct, f"Step {sc}/{total} | {self._elapsed_str()} elapsed")
                    print(
                        f"[ssign] Starting step {sc}/{total}: "
                        f"{name} ({step_id})",
                        flush=True,
                    )

                    try:
                        result = func()
                        self.results.append(result)
                        self._save_progress()
                        print(
                            f"[ssign] Finished step {sc}/{total}: "
                            f"{name} -> "
                            f"{'OK' if result.success else 'FAILED: ' + result.message[:100]}",
                            flush=True,
                        )
                        if not result.success:
                            self.progress(
                                name, pct,
                                f"Failed: {result.message}")
                            logger.error(
                                f"Step '{name}' failed: "
                                f"{result.message}")
                            if step_id in CORE_STEPS:
                                core_failed = True
                    except Exception as e:
                        print(
                            f"[ssign] EXCEPTION in step {sc}/{total}: "
                            f"{name} -> {e}",
                            flush=True,
                        )
                        self.results.append(
                            StepResult(step_id, False, str(e)))
                        self._save_progress()
                        logger.exception(
                            f"Step '{name}' raised exception")
                        if step_id in CORE_STEPS:
                            core_failed = True

        if n_skipped:
            logger.info(f"Resumed: skipped {n_skipped} previously completed steps")

        self.progress("Complete", 100, f"Pipeline finished in {self._elapsed_str()}")

        # Copy final outputs to outdir
        self._copy_outputs()
        self._save_progress()

        return self.results

    def _wf(self, name):
        """Get work file path."""
        return os.path.join(self.work_dir, name)

    # ── Phase 1: Input Processing ──

    def _step_detect_format(self) -> StepResult:
        rc, stdout, stderr = run_script("detect_input_format.py", [
            self.config.input_path,
        ])
        if rc == 0:
            fmt = stdout.strip()
            self.files['format'] = fmt
            self.files['input'] = self.config.input_path
            return StepResult("detect_format", True, f"Format: {fmt}")
        return StepResult("detect_format", False, stderr[:500])

    def _step_extract_proteins(self) -> StepResult:
        proteins_path = self._wf(f"{self.config.sample_id}_proteins.faa")
        gene_info_path = self._wf(f"{self.config.sample_id}_gene_info.tsv")
        fmt = self.files.get('format', '')

        if fmt == 'fasta_contigs' and self.config.run_bakta:
            # Use Bakta for richer annotation of raw contig FASTA
            if not self.config.bakta_db:
                return StepResult("extract_proteins", False,
                                  "Bakta mode requires --bakta-db path. "
                                  "Download with: bakta_db download --output /path/to/db --type light")
            rc, stdout, stderr = run_script("run_bakta.py", [
                "--input", self.config.input_path,
                "--db", self.config.bakta_db,
                "--sample", self.config.sample_id,
                "--threads", str(self.config.bakta_threads),
                "--out-proteins", proteins_path,
                "--out-gene-info", gene_info_path,
            ], timeout=14400)
            tool_name = "Bakta"
        else:
            # Handles GenBank, GFF3, and FASTA contigs (via Prodigal)
            metadata_path = self._wf(f"{self.config.sample_id}_metadata.json")
            extract_args = [
                "--input", self.config.input_path,
                "--sample", self.config.sample_id,
                "--out-proteins", proteins_path,
                "--out-gene-info", gene_info_path,
                "--out-metadata", metadata_path,
            ]
            if self.config.original_filename:
                extract_args.extend(["--original-filename", self.config.original_filename])
            rc, stdout, stderr = run_script("extract_proteins.py", extract_args)
            tool_name = "Prodigal" if fmt == 'fasta_contigs' else "GenBank parser"

        if rc == 0:
            self.files['proteins'] = proteins_path
            self.files['gene_info'] = gene_info_path

            # Read organism name from metadata if available
            metadata_path = self._wf(f"{self.config.sample_id}_metadata.json")
            if os.path.exists(metadata_path):
                try:
                    import json
                    with open(metadata_path) as mf:
                        meta = json.load(mf)
                    self.files['organism'] = meta.get('organism', '')
                except Exception:
                    pass

            # Also extract gene order
            gene_order_path = self._wf(f"{self.config.sample_id}_gene_order.tsv")
            rc2, _, stderr2 = run_script("extract_gene_order.py", [
                "--gene-info", gene_info_path,
                "--output", gene_order_path,
            ])
            if rc2 == 0:
                self.files['gene_order'] = gene_order_path

            n = sum(1 for line in open(proteins_path) if line.startswith('>'))
            return StepResult("extract_proteins", True, f"{tool_name}: extracted {n} proteins")
        return StepResult("extract_proteins", False, stderr[:500])

    # ── Phase 2: SS Detection ──

    def _step_macsyfinder(self) -> StepResult:
        """Run MacSyFinder v2 to detect secretion systems."""
        msf_out = self._wf("macsyfinder_out")
        # Do NOT pre-create: MacSyFinder v2 requires --out-dir to not exist
        if os.path.exists(msf_out):
            shutil.rmtree(msf_out)

        proteins = self.files.get('proteins', '')
        if not proteins or not os.path.exists(proteins):
            return StepResult("macsyfinder", False, "No proteins file from previous step")

        # FRAGILE: macsydata install — downloads TXSScan HMM models from the
        # macsy-models GitHub repository. Can fail behind corporate firewalls
        # or if GitHub is unreachable. Models persist in ~/.macsyfinder/ after
        # first successful install, so this usually only runs once.
        # If this breaks: run manually: macsydata install --user TXSScan
        install_cmd = ["macsydata", "install", "--user", "TXSScan"]
        try:
            install_result = subprocess.run(
                install_cmd, capture_output=True, text=True, timeout=120
            )
            if install_result.returncode != 0:
                logger.warning(
                    f"macsydata install returned exit code "
                    f"{install_result.returncode}. TXSScan models may "
                    f"already be installed — continuing.\n"
                    f"  If MacSyFinder fails, run manually:\n"
                    f"    macsydata install --user TXSScan"
                )
        except FileNotFoundError:
            return StepResult(
                "macsyfinder", False,
                "macsydata command not found.\n"
                "  Install: pip install macsyfinder\n"
                "  It should have been installed as a dependency of ssign."
            )
        except subprocess.TimeoutExpired:
            logger.warning(
                "macsydata install timed out — TXSScan models may already "
                "exist. Continuing."
            )
        except Exception:
            pass  # May already be installed

        # FRAGILE: macsyfinder CLI — internally calls hmmsearch via subprocess.
        # Our pyhmmer shim (installed as console_script 'hmmsearch') should be
        # found on PATH. If not, MacSyFinder will fail with a message about
        # hmmsearch not being found.
        # If this breaks: sudo apt install hmmer (or: conda install -c bioconda hmmer)
        cmd = [
            "macsyfinder",
            "--sequence-db", proteins,
            "--db-type", "ordered_replicon",
            "--models", "TXSScan", "all",
            "--out-dir", msf_out,
            "--mute",
        ]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600,
            )
            if result.returncode == 0:
                self.files['macsyfinder_out'] = msf_out
                return StepResult("macsyfinder", True, "MacSyFinder v2 complete")

            error_msg = result.stderr[:500]
            # Detect common MacSyFinder failures and provide specific advice
            if "hmmsearch" in error_msg.lower() or "hmmer" in error_msg.lower():
                return StepResult(
                    "macsyfinder", False,
                    f"MacSyFinder failed — hmmsearch not available.\n"
                    f"  The pyhmmer shim should provide this automatically.\n"
                    f"  How to fix:\n"
                    f"    - pip install --force-reinstall ssign\n"
                    f"    - Or install real HMMER: sudo apt install hmmer"
                )
            if "txscan" in error_msg.lower() or "model" in error_msg.lower():
                return StepResult(
                    "macsyfinder", False,
                    f"MacSyFinder failed — TXSScan models may not be "
                    f"installed.\n"
                    f"  Run: macsydata install --user TXSScan\n"
                    f"  Original error: {error_msg}"
                )
            return StepResult("macsyfinder", False, error_msg)
        except FileNotFoundError:
            return StepResult(
                "macsyfinder", False,
                "macsyfinder command not found.\n"
                "  Install: pip install macsyfinder\n"
                "  It should have been installed as a dependency of ssign."
            )
        except subprocess.TimeoutExpired:
            return StepResult("macsyfinder", False, "Timeout after 60 min")
        except Exception as e:
            return StepResult("macsyfinder", False, str(e))

    def _step_validate_systems(self) -> StepResult:
        """Validate MacSyFinder results and extract SS components."""
        msf_out = self.files.get('macsyfinder_out', '')
        gene_info = self.files.get('gene_info', '')
        components_path = self._wf(f"{self.config.sample_id}_ss_components.tsv")
        systems_path = self._wf(f"{self.config.sample_id}_valid_systems.tsv")

        if not msf_out or not os.path.isdir(msf_out):
            return StepResult("validate_systems", False, "No MacSyFinder output from previous step")

        rc, stdout, stderr = run_script("validate_macsyfinder_systems.py", [
            "--msf-dir", msf_out,
            "--gene-info", gene_info,
            "--sample", self.config.sample_id,
            "--wholeness-threshold", str(self.config.wholeness_threshold),
            "--excluded-systems", ",".join(self.config.excluded_systems),
            "--out-components", components_path,
            "--out-systems", systems_path,
        ])

        if rc == 0:
            self.files['ss_components'] = components_path
            self.files['valid_systems'] = systems_path
            return StepResult("validate_systems", True, "Systems validated")
        return StepResult("validate_systems", False, stderr[:500])

    def _step_extract_neighborhood(self) -> StepResult:
        """Extract proteins near SS components for focused prediction."""
        gene_order = self.files.get('gene_order', '')
        ss_components = self.files.get('ss_components', '')
        proteins = self.files.get('proteins', '')

        if not ss_components or not os.path.exists(ss_components):
            return StepResult("extract_neighborhood", False,
                              "No SS components from previous step")

        neighborhood_fasta = self._wf(f"{self.config.sample_id}_neighborhood.faa")

        rc, stdout, stderr = run_script("extract_neighborhood.py", [
            "--gene-order", gene_order,
            "--ss-components", ss_components,
            "--proteins", proteins,
            "--window", str(self.config.proximity_window),
            "--output", neighborhood_fasta,
        ])

        if rc == 0:
            self.files['neighborhood_proteins'] = neighborhood_fasta
            return StepResult("extract_neighborhood", True,
                              "Neighborhood proteins extracted")
        return StepResult("extract_neighborhood", False, stderr[:500])

    # ── Phase 3: Prediction ──

    def _step_deeplocpro(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_deeplocpro.tsv")

        # Use neighborhood proteins (focused) if available, else full proteome
        input_proteins = self.files.get('neighborhood_proteins',
                                         self.files.get('proteins', ''))

        args = [
            "--input", input_proteins,
            "--sample", self.config.sample_id,
            "--conf-threshold", str(self.config.conf_threshold),
            "--output", output,
        ]

        if self.config.deeplocpro_mode == "local" and self.config.deeplocpro_path:
            args.extend(["--mode", "local", "--deeplocpro-path", self.config.deeplocpro_path])
        else:
            args.extend(["--mode", "remote"])

        rc, stdout, stderr = run_script("run_deeplocpro.py", args, timeout=14400)
        if rc == 0:
            self.files['deeplocpro'] = output
            return StepResult("deeplocpro", True, "DeepLocPro complete")
        return StepResult("deeplocpro", False, stderr[:500])

    def _step_deepsece(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_deepsece.tsv")

        # Use neighborhood proteins if available
        input_proteins = self.files.get('neighborhood_proteins',
                                         self.files.get('proteins', ''))

        rc, stdout, stderr = run_script("run_deepsece.py", [
            "--input", input_proteins,
            "--sample", self.config.sample_id,
            "--output", output,
        ], timeout=14400)

        if rc == 0:
            self.files['deepsece'] = output
            return StepResult("deepsece", True, "DeepSecE complete")
        return StepResult("deepsece", False, stderr[:500])

    def _step_signalp(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_signalp.tsv")

        # Use neighborhood proteins if available
        input_proteins = self.files.get('neighborhood_proteins',
                                         self.files.get('proteins', ''))

        args = [
            "--input", input_proteins,
            "--sample", self.config.sample_id,
            "--output", output,
        ]
        if self.config.signalp_mode == "local" and self.config.signalp_path:
            args.extend(["--mode", "local", "--signalp-path", self.config.signalp_path])
        else:
            args.extend(["--mode", "remote"])

        rc, stdout, stderr = run_script("run_signalp.py", args, timeout=14400)
        if rc == 0:
            self.files['signalp'] = output
            return StepResult("signalp", True, "SignalP complete")
        return StepResult("signalp", False, stderr[:500])

    def _step_cross_validate(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_predictions.tsv")

        dlp = self.files.get('deeplocpro', '')
        if not dlp or not os.path.exists(dlp):
            return StepResult("cross_validate", False, "No DeepLocPro output from previous step")

        dse = self.files.get('deepsece', '')
        # DeepSecE is optional — cross_validate works with DLP alone

        valid_sys = self.files.get('valid_systems', '')
        if not valid_sys or not os.path.exists(valid_sys):
            return StepResult("cross_validate", False,
                              "No valid_systems file from previous step")

        args = [
            "--deeplocpro", dlp,
            "--valid-systems", valid_sys,
            "--sample", self.config.sample_id,
            "--conf-threshold", str(self.config.conf_threshold),
            "--output", output,
        ]
        if dse and os.path.exists(dse):
            args.extend(["--deepsece", dse])

        rc, stdout, stderr = run_script("cross_validate_predictions.py", args)
        if rc == 0:
            self.files['predictions'] = output
            return StepResult("cross_validate", True, "Predictions validated")
        return StepResult("cross_validate", False, stderr[:500])

    # ── Phase 4: Substrate ID ──

    def _step_proximity(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_substrates.tsv")

        rc, stdout, stderr = run_script("proximity_analysis.py", [
            "--gene-order", self.files.get('gene_order', ''),
            "--ss-components", self.files.get('ss_components', ''),
            "--predictions", self.files.get('predictions', ''),
            "--sample", self.config.sample_id,
            "--window", str(self.config.proximity_window),
            "--conf-threshold", str(self.config.conf_threshold),
            "--output", output,
        ])

        if rc == 0:
            self.files['substrates'] = output
            n = sum(1 for _ in open(output)) - 1  # minus header
            return StepResult("proximity", True, f"Found {n} substrate candidates")
        return StepResult("proximity", False, stderr[:500])

    def _step_t5ss(self) -> StepResult:
        out_sub = self._wf(f"{self.config.sample_id}_t5ss_substrates.tsv")
        out_dom = self._wf(f"{self.config.sample_id}_t5ss_domains.tsv")

        rc, stdout, stderr = run_script("t5ss_handler.py", [
            "--ss-components", self.files.get('ss_components', ''),
            "--predictions", self.files.get('predictions', ''),
            "--sample", self.config.sample_id,
            "--out-substrates", out_sub,
            "--out-domains", out_dom,
        ])

        if rc == 0:
            self.files['t5ss_substrates'] = out_sub
            return StepResult("t5ss", True, "T5SS handled")
        return StepResult("t5ss", False, stderr[:500])

    def _step_filtering(self) -> StepResult:
        out_filtered = self._wf(f"{self.config.sample_id}_substrates_filtered.tsv")
        out_all = self._wf(f"{self.config.sample_id}_substrates_all.tsv")

        substrates = self.files.get('substrates', '')
        t5ss = self.files.get('t5ss_substrates', '')
        if not substrates or not os.path.exists(substrates):
            return StepResult("filtering", False, "No substrates from previous step")

        filter_args = [
            "--proximity-substrates", substrates,
            "--t5ss-substrates", t5ss if t5ss and os.path.exists(t5ss) else substrates,
            "--valid-systems", self.files.get('valid_systems', ''),
            "--predictions", self.files.get('predictions', ''),
            "--sample", self.config.sample_id,
            "--excluded-systems", ",".join(self.config.excluded_systems),
            "--out-filtered", out_filtered,
            "--out-all", out_all,
        ]
        if self.config.filter_dse_type_mismatch:
            filter_args.append("--filter-dse-type-mismatch")

        rc, stdout, stderr = run_script("system_filtering.py", filter_args)

        if rc == 0:
            self.files['substrates_filtered'] = out_filtered
            self.files['substrates_all'] = out_all
            return StepResult("filtering", True, "Filtering complete")
        return StepResult("filtering", False, stderr[:500])

    # ── Phase 5: Annotation ──

    def _check_substrates_exist(self, step_name):
        """Check that upstream substrate files exist. Returns error StepResult or None."""
        sf = self.files.get('substrates_filtered', '')
        if not sf or not os.path.exists(sf):
            return StepResult(step_name, False, "Skipped — no substrates from upstream steps")
        return None

    def _step_blastp(self) -> StepResult:
        err = self._check_substrates_exist("blastp")
        if err:
            return err

        output = self._wf(f"{self.config.sample_id}_blastp.csv")

        args = [
            "--mode", self.config.blastp_mode,
            "--substrates", self.files.get('substrates_filtered', ''),
            "--proteins", self.files.get('proteins', ''),
            "--sample", self.config.sample_id,
            "--output", output,
            "--min-pident", str(self.config.blastp_min_pident),
            "--min-qcov", str(self.config.blastp_min_qcov),
            "--evalue", str(self.config.blastp_evalue),
        ]
        if self.config.blastp_mode == "local" and self.config.blastp_db:
            args.extend(["--db", self.config.blastp_db])
        if self.config.blastp_exclude_taxid:
            args.extend(["--exclude-taxid", self.config.blastp_exclude_taxid])

        rc, stdout, stderr = run_script("run_blastp.py", args, timeout=7200)
        if rc == 0:
            self.files['blastp'] = output
            return StepResult("blastp", True, "BLASTp complete")
        return StepResult("blastp", False, stderr[:500])

    def _step_hhsuite(self) -> StepResult:
        err = self._check_substrates_exist("hhsuite")
        if err:
            return err
        output = self._wf(f"{self.config.sample_id}_hhsuite.csv")

        args = [
            "--mode", self.config.hhsuite_mode,
            "--substrates", self.files.get('substrates_filtered', ''),
            "--proteins", self.files.get('proteins', ''),
            "--sample", self.config.sample_id,
            "--output", output,
        ]
        if self.config.hhsuite_mode == "local":
            if self.config.hhsuite_pfam_db:
                args.extend(["--pfam-db", self.config.hhsuite_pfam_db])
            if self.config.hhsuite_pdb70_db:
                args.extend(["--pdb70-db", self.config.hhsuite_pdb70_db])
            if self.config.hhsuite_uniclust_db:
                args.extend(["--uniclust-db", self.config.hhsuite_uniclust_db])

        rc, stdout, stderr = run_script("run_hhsuite.py", args, timeout=14400)
        if rc == 0:
            self.files['hhsuite'] = output
            return StepResult("hhsuite", True, "HH-suite complete")
        return StepResult("hhsuite", False, stderr[:500])

    def _step_interproscan(self) -> StepResult:
        err = self._check_substrates_exist("interproscan")
        if err:
            return err
        output = self._wf(f"{self.config.sample_id}_interproscan.csv")

        args = [
            "--mode", self.config.interproscan_mode,
            "--substrates", self.files.get('substrates_filtered', ''),
            "--proteins", self.files.get('proteins', ''),
            "--sample", self.config.sample_id,
            "--output", output,
        ]
        if self.config.interproscan_mode == "local" and self.config.interproscan_db:
            args.extend(["--db", self.config.interproscan_db])

        rc, stdout, stderr = run_script("run_interproscan.py", args, timeout=7200)
        if rc == 0:
            self.files['interproscan'] = output
            return StepResult("interproscan", True, "InterProScan complete")
        return StepResult("interproscan", False, stderr[:500])

    def _step_protparam(self) -> StepResult:
        err = self._check_substrates_exist("protparam")
        if err:
            return err
        output = self._wf(f"{self.config.sample_id}_protparam.csv")

        rc, stdout, stderr = run_script("compute_protparam.py", [
            "--substrates", self.files.get('substrates_filtered', ''),
            "--proteins", self.files.get('proteins', ''),
            "--sample", self.config.sample_id,
            "--output", output,
        ])

        if rc == 0:
            self.files['protparam'] = output
            return StepResult("protparam", True, "ProtParam complete")
        return StepResult("protparam", False, stderr[:500])

    # ── Phase 6: Integration ──

    def _step_integrate(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_integrated.csv")

        annotation_files = []
        for key in ['blastp', 'hhsuite', 'interproscan', 'protparam', 'signalp']:
            if key in self.files and os.path.exists(self.files[key]):
                annotation_files.append(self.files[key])

        args = [
            "--substrates-filtered", self.files.get('substrates_filtered', ''),
            "--substrates-all", self.files.get('substrates_all', ''),
            "--sample", self.config.sample_id,
            "--output", output,
        ]
        if annotation_files:
            args.extend(["--annotations"] + annotation_files)

        rc, stdout, stderr = run_script("integrate_annotations.py", args)
        if rc == 0:
            self.files['integrated'] = output
            return StepResult("integrate", True, "Annotations integrated")
        return StepResult("integrate", False, stderr[:500])

    def _step_orthologs(self) -> StepResult:
        """Assign ortholog groups via all-vs-all BLASTp + Union-Find clustering."""
        sf = self.files.get('substrates_filtered', '')
        proteins = self.files.get('proteins', '')

        if not sf or not os.path.exists(sf):
            return StepResult("orthologs", True, "No substrates — skipping orthologs")

        # Count substrates
        n_sub = sum(1 for _ in open(sf)) - 1
        if n_sub < 2:
            return StepResult("orthologs", True,
                              f"Only {n_sub} substrate(s) — skipping ortholog grouping")

        # Extract substrate sequences into a dedicated FASTA
        substrate_fasta = self._wf(f"{self.config.sample_id}_substrates_for_ortho.faa")
        substrate_ids = set()
        with open(sf) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                lid = row.get('locus_tag', row.get('protein_id', ''))
                if lid:
                    substrate_ids.add(lid)

        from Bio import SeqIO
        n_written = 0
        with open(substrate_fasta, 'w') as out:
            for rec in SeqIO.parse(proteins, 'fasta'):
                if rec.id in substrate_ids:
                    SeqIO.write(rec, out, 'fasta')
                    n_written += 1

        if n_written < 2:
            return StepResult("orthologs", True,
                              f"Only {n_written} substrate sequences — skipping")

        output = self._wf(f"{self.config.sample_id}_orthologs.csv")
        output_groups = self._wf(f"{self.config.sample_id}_ortholog_groups.csv")

        rc, stdout, stderr = run_script("run_ortholog_grouping.py", [
            "--substrates-fasta", substrate_fasta,
            "--min-pident", str(self.config.ortholog_min_pident),
            "--min-qcov", str(self.config.ortholog_min_qcov),
            "--output", output,
            "--output-groups", output_groups,
        ], timeout=3600)

        if rc == 0:
            self.files['orthologs'] = output
            self.files['ortholog_groups'] = output_groups

            # Merge ortholog assignments into the integrated CSV
            integrated = self.files.get('integrated', '')
            if integrated and os.path.exists(integrated) and os.path.exists(output):
                try:
                    import pandas as pd
                    df_int = pd.read_csv(integrated)
                    df_og = pd.read_csv(output)
                    # Merge on locus_tag
                    df_merged = df_int.merge(df_og, on='locus_tag', how='left')
                    df_merged.to_csv(integrated, index=False)
                    logger.info(f"Merged ortholog groups into {integrated}")
                except Exception as e:
                    logger.warning(f"Could not merge orthologs into integrated CSV: {e}")

            return StepResult("orthologs", True, "Ortholog groups assigned")
        return StepResult("orthologs", False, stderr[:500])

    def _step_enrichment(self) -> StepResult:
        """Run enrichment testing (Fisher's exact + permutation)."""
        integrated = self.files.get('integrated', '')
        if not integrated or not os.path.exists(integrated):
            return StepResult("enrichment", False, "No integrated CSV — skipping enrichment")

        out_fisher = self._wf(f"{self.config.sample_id}_enrichment_fisher.csv")
        out_perm = self._wf(f"{self.config.sample_id}_enrichment_permutation.csv")
        out_summary = self._wf(f"{self.config.sample_id}_enrichment_summary.txt")

        rc, stdout, stderr = run_script("enrichment_testing.py", [
            "--integrated-csv", integrated,
            "--n-permutations", "10000",
            "--out-fisher", out_fisher,
            "--out-permutation", out_perm,
            "--out-summary", out_summary,
        ])

        if rc == 0:
            self.files['enrichment_fisher'] = out_fisher
            self.files['enrichment_perm'] = out_perm
            self.files['enrichment_summary'] = out_summary
            return StepResult("enrichment", True, "Enrichment analysis complete")
        return StepResult("enrichment", False, stderr[:500])

    def _step_report(self) -> StepResult:
        html_out = self._wf(f"{self.config.sample_id}_report.html")
        txt_out = self._wf(f"{self.config.sample_id}_report.txt")

        rc, stdout, stderr = run_script("generate_report.py", [
            "--master-csvs", self.files.get('integrated', ''),
            "--out-html", html_out,
            "--out-txt", txt_out,
        ])

        if rc == 0:
            self.files['report_html'] = html_out
            self.files['report_txt'] = txt_out
            return StepResult("report", True, "Report generated")
        return StepResult("report", False, stderr[:500])

    def _step_figures(self) -> StepResult:
        fig_dir = self._wf("figures")
        os.makedirs(fig_dir, exist_ok=True)

        fig_args = [
            "--master-csvs", self.files.get('integrated', ''),
            "--outdir", fig_dir,
            "--dpi", str(self.config.dpi),
        ]
        # Pass figure toggles
        if not self.config.fig_category:
            fig_args.append("--no-category")
        if not self.config.fig_ss_comp:
            fig_args.append("--no-ss-comp")
        if not self.config.fig_tool_heatmap:
            fig_args.append("--no-tool-heatmap")
        if not self.config.fig_substrate_count:
            fig_args.append("--no-substrate-count")
        if not self.config.fig_func_summary:
            fig_args.append("--no-func-summary")

        rc, stdout, stderr = run_script("generate_figures.py", fig_args)

        if rc == 0:
            self.files['figures_dir'] = fig_dir
            return StepResult("figures", True, "Figures generated")
        return StepResult("figures", False, stderr[:500])

    def _save_progress(self):
        """Save pipeline progress to JSON for resume capability.

        Saves both step results and intermediate file paths so a resumed
        run can skip already-completed steps and find their output files.
        """
        try:
            outdir = Path(self.config.outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            progress = {
                "timestamp": datetime.now().isoformat(),
                "sample_id": self.config.sample_id,
                "work_dir": self.work_dir,
                "steps": [
                    {"name": r.name, "success": r.success, "message": r.message}
                    for r in self.results
                ],
                "files": self.files,
            }
            progress_path = outdir / "ssign_progress.json"
            with open(progress_path, 'w') as f:
                json.dump(progress, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save progress: {e}")

    @staticmethod
    def load_progress(outdir):
        """Load pipeline progress from a previous run.

        Returns (list[StepResult], dict_files, work_dir) or (None, {}, "").
        """
        progress_path = Path(outdir) / "ssign_progress.json"
        if not progress_path.exists():
            return None, {}, ""
        try:
            with open(progress_path) as f:
                data = json.load(f)
            results = [
                StepResult(s["name"], s["success"], s["message"])
                for s in data.get("steps", [])
            ]
            files = data.get("files", {})
            work_dir = data.get("work_dir", "")
            return results, files, work_dir
        except Exception:
            return None, {}, ""

    def _try_resume(self) -> set:
        """Attempt to resume from previous progress. Returns set of step
        names to skip (already completed successfully with valid outputs).
        """
        prev_results, prev_files, prev_work_dir = self.load_progress(
            self.config.outdir
        )
        if prev_results is None:
            return set()

        # Validate that the previous work_dir still exists with files
        if not prev_work_dir or not os.path.isdir(prev_work_dir):
            logger.info("Previous work directory gone — starting fresh")
            return set()

        # Check which steps completed successfully AND have valid output files
        completed = set()
        prev_succeeded = {r.name for r in prev_results if r.success}

        for step_name in prev_succeeded:
            completed.add(step_name)

        if completed:
            # Restore file paths and work directory from previous run
            self.work_dir = prev_work_dir
            # Only restore files that still exist
            for key, path in prev_files.items():
                if isinstance(path, str) and (
                    os.path.exists(path) or not path
                ):
                    self.files[key] = path
            logger.info(
                f"Resuming: {len(completed)} steps already done: "
                f"{', '.join(sorted(completed))}"
            )

        return completed

    def _copy_outputs(self):
        """Copy consolidated outputs to the user's output directory.

        Produces:
          {sample_id}_results.csv  - Master CSV (systems on top, substrates below)
          {sample_id}_summary.txt  - Report + enrichment stats
          figures/                 - Publication-quality plots
        """
        outdir = Path(self.config.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        sid = self.config.sample_id

        # 1. Master results CSV — systems + substrates in one file
        self._build_master_csv(outdir / f"{sid}_results.csv")

        # 2. Summary — report text + enrichment summary + Fisher table
        self._build_summary(outdir / f"{sid}_summary.txt")

        # 3. Figures directory
        if 'figures_dir' in self.files and os.path.exists(self.files['figures_dir']):
            dest = outdir / "figures"
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(self.files['figures_dir'], dest)

    def _build_master_csv(self, output_path: Path):
        """Build a single CSV: secretion systems on top, substrates below.

        Systems section: system-level rows then component rows, sorted by
        genome → ss_type. Substrates section: annotation rows sorted by
        genome → nearby_ss_types, with a nearby_sys_ids column linking
        back to the system rows. A blank row separates the two sections.
        """
        import pandas as pd

        # ── Systems section ──
        df_sys = None
        fpath = self.files.get('valid_systems', '')
        if fpath and os.path.exists(fpath):
            try:
                df_sys = pd.read_csv(fpath, sep='\t')
                df_sys.insert(0, 'record_type', 'system')
            except Exception:
                pass

        df_comp = None
        fpath = self.files.get('ss_components', '')
        if fpath and os.path.exists(fpath):
            try:
                df_comp = pd.read_csv(fpath, sep='\t')
                df_comp.insert(0, 'record_type', 'component')
            except Exception:
                pass

        # Build sys_id lookup: ss_type → list of sys_ids (for substrate linkage)
        sys_id_by_type: dict[str, list[str]] = {}
        if df_comp is not None and 'ss_type' in df_comp.columns and 'sys_id' in df_comp.columns:
            for ss_type, group in df_comp.groupby('ss_type'):
                sys_id_by_type[ss_type] = sorted(group['sys_id'].unique().tolist())

        # Combine and sort systems: by sample_id, ss_type, then record_type
        # (system rows before component rows for each SS type)
        sys_frames = [df for df in [df_sys, df_comp] if df is not None]
        df_systems = None
        if sys_frames:
            df_systems = pd.concat(sys_frames, ignore_index=True)
            type_order = {'system': 0, 'component': 1}
            df_systems['_sort'] = df_systems['record_type'].map(type_order)
            sort_cols = []
            if 'sample_id' in df_systems.columns:
                sort_cols.append('sample_id')
            if 'ss_type' in df_systems.columns:
                sort_cols.append('ss_type')
            sort_cols.append('_sort')
            df_systems = df_systems.sort_values(sort_cols).drop(columns=['_sort'])

        # ── Substrates section ──
        df_subs = None
        fpath = self.files.get('integrated', '')
        if fpath and os.path.exists(fpath):
            try:
                df_subs = pd.read_csv(fpath)
            except Exception:
                pass

        if df_subs is not None and sys_id_by_type:
            # Add nearby_sys_ids column: map each nearby_ss_type to its sys_ids
            def _resolve_sys_ids(nearby_types):
                if pd.isna(nearby_types) or not nearby_types:
                    return ''
                ids = []
                for t in str(nearby_types).split(','):
                    t = t.strip()
                    ids.extend(sys_id_by_type.get(t, []))
                return ', '.join(sorted(set(ids)))

            if 'nearby_ss_types' in df_subs.columns:
                # Insert nearby_sys_ids right after nearby_ss_types
                idx = df_subs.columns.get_loc('nearby_ss_types') + 1
                df_subs.insert(idx, 'nearby_sys_ids',
                               df_subs['nearby_ss_types'].apply(_resolve_sys_ids))

            # Sort by genome, then SS type
            sort_cols = []
            if 'sample_id' in df_subs.columns:
                sort_cols.append('sample_id')
            if 'nearby_ss_types' in df_subs.columns:
                sort_cols.append('nearby_ss_types')
            if sort_cols:
                df_subs = df_subs.sort_values(sort_cols)

        # ── Write combined CSV ──
        with open(output_path, 'w', newline='') as f:
            wrote_something = False
            if df_systems is not None and not df_systems.empty:
                df_systems.to_csv(f, index=False)
                wrote_something = True

            if df_subs is not None and not df_subs.empty:
                if wrote_something:
                    f.write('\n')  # blank separator row
                df_subs.to_csv(f, index=False)

    def _build_summary(self, output_path: Path):
        """Combine report text, enrichment summary, and Fisher results."""
        parts = []

        # Report text
        if 'report_txt' in self.files and os.path.exists(self.files['report_txt']):
            with open(self.files['report_txt']) as f:
                parts.append(f.read())

        # Enrichment summary
        if 'enrichment_summary' in self.files and os.path.exists(self.files['enrichment_summary']):
            with open(self.files['enrichment_summary']) as f:
                parts.append(
                    "\n\n" + "=" * 60 + "\n"
                    "ENRICHMENT ANALYSIS\n"
                    + "=" * 60 + "\n\n" + f.read()
                )

        # Fisher results table
        if 'enrichment_fisher' in self.files and os.path.exists(self.files['enrichment_fisher']):
            try:
                import pandas as pd
                df = pd.read_csv(self.files['enrichment_fisher'])
                if not df.empty:
                    parts.append(
                        "\n\n" + "-" * 60 + "\n"
                        "Fisher's Exact Test Results\n"
                        + "-" * 60 + "\n\n"
                        + df.to_string(index=False) + "\n"
                    )
            except Exception:
                pass

        if parts:
            with open(output_path, 'w') as f:
                f.write(''.join(parts))


def run_cross_genome_orthologs(
    genome_outdirs: list[str],
    output_dir: str,
    min_pident: float = 40.0,
    min_qcov: float = 70.0,
    progress_callback: Optional[Callable] = None,
) -> dict:
    """Run cross-genome ortholog grouping on substrates from multiple genomes.

    After each genome is processed individually, this function:
    1. Collects all substrate sequences from all genome output directories
    2. Runs all-vs-all BLASTp on the combined set
    3. Clusters into ortholog groups using Union-Find
    4. Merges cross-genome ortholog group IDs back into each genome's integrated CSV

    Returns dict with n_proteins, n_groups, genomes_updated, file paths.
    """
    from Bio import SeqIO
    import pandas as pd

    result = {
        'n_proteins': 0, 'n_groups': 0, 'genomes_updated': 0,
        'combined_fasta': '', 'orthologs_csv': '', 'ortholog_groups_csv': '',
    }

    if progress_callback:
        progress_callback("Cross-genome orthologs", 10,
                          "Collecting substrates from all genomes...")

    os.makedirs(output_dir, exist_ok=True)
    combined_fasta = os.path.join(output_dir, "all_substrates_combined.faa")

    # Step 1: Collect all substrate sequences
    all_substrate_ids = {}  # locus_tag -> genome_outdir
    n_written = 0

    with open(combined_fasta, 'w') as out_f:
        for gdir in genome_outdirs:
            gpath = Path(gdir)

            # Find integrated CSV to get substrate locus_tags
            integrated_csvs = list(gpath.glob("*_integrated.csv")) + \
                              list(gpath.glob("*integrated*.csv"))
            if not integrated_csvs:
                logger.warning(f"No integrated CSV in {gdir} -- skipping")
                continue

            integrated_csv = str(integrated_csvs[0])
            substrate_ids = set()
            try:
                df = pd.read_csv(integrated_csv)
                if 'locus_tag' in df.columns:
                    substrate_ids = set(df['locus_tag'].dropna().astype(str))
            except Exception as e:
                logger.warning(f"Could not read {integrated_csv}: {e}")
                continue

            if not substrate_ids:
                continue

            # Find protein FASTA — check progress.json first, then glob
            protein_fastas = []
            progress_path = gpath / "ssign_progress.json"
            if progress_path.exists():
                try:
                    with open(progress_path) as pf:
                        prog = json.load(pf)
                    proteins_path = prog.get("files", {}).get("proteins", "")
                    if proteins_path and os.path.exists(proteins_path):
                        protein_fastas.insert(0, Path(proteins_path))
                except Exception:
                    pass

            protein_fastas.extend(list(gpath.glob("*.faa")))

            for pf in protein_fastas:
                try:
                    for rec in SeqIO.parse(str(pf), 'fasta'):
                        if rec.id in substrate_ids and rec.id not in all_substrate_ids:
                            SeqIO.write(rec, out_f, 'fasta')
                            all_substrate_ids[rec.id] = gdir
                            n_written += 1
                except Exception:
                    continue

    result['combined_fasta'] = combined_fasta
    result['n_proteins'] = n_written

    if n_written < 2:
        logger.info(f"Only {n_written} substrate(s) across all genomes -- "
                     "skipping cross-genome ortholog grouping")
        return result

    if progress_callback:
        progress_callback("Cross-genome orthologs", 40,
                          f"Running all-vs-all BLASTp on {n_written} substrates...")

    # Step 2: Run ortholog grouping
    orthologs_csv = os.path.join(output_dir, "cross_genome_orthologs.csv")
    groups_csv = os.path.join(output_dir, "cross_genome_ortholog_groups.csv")

    rc, stdout, stderr = run_script("run_ortholog_grouping.py", [
        "--substrates-fasta", combined_fasta,
        "--min-pident", str(min_pident),
        "--min-qcov", str(min_qcov),
        "--output", orthologs_csv,
        "--output-groups", groups_csv,
    ], timeout=3600)

    if rc != 0:
        logger.error(f"Cross-genome ortholog grouping failed: {stderr[:500]}")
        return result

    result['orthologs_csv'] = orthologs_csv
    result['ortholog_groups_csv'] = groups_csv

    if progress_callback:
        progress_callback("Cross-genome orthologs", 70,
                          "Merging ortholog groups into per-genome results...")

    # Step 3: Read ortholog assignments and merge into each genome's integrated CSV
    try:
        df_og = pd.read_csv(orthologs_csv)
        result['n_groups'] = df_og['ortholog_group'].nunique()
    except Exception as e:
        logger.warning(f"Could not read ortholog results: {e}")
        return result

    # Prefix with "xg_" (cross-genome) to distinguish from per-genome groups
    df_og = df_og.rename(columns={
        'ortholog_group': 'xg_ortholog_group',
        'og_n_members': 'xg_og_n_members',
        'og_mean_pident': 'xg_og_mean_pident',
    })

    genomes_updated = 0
    for gdir in genome_outdirs:
        gpath = Path(gdir)
        integrated_csvs = list(gpath.glob("*_integrated.csv")) + \
                          list(gpath.glob("*integrated*.csv"))
        if not integrated_csvs:
            continue

        integrated_csv = str(integrated_csvs[0])
        try:
            df_int = pd.read_csv(integrated_csv)
            for col in ['xg_ortholog_group', 'xg_og_n_members', 'xg_og_mean_pident']:
                if col in df_int.columns:
                    df_int = df_int.drop(columns=[col])

            df_merged = df_int.merge(df_og, on='locus_tag', how='left')
            df_merged.to_csv(integrated_csv, index=False)
            genomes_updated += 1
            logger.info(f"Merged cross-genome orthologs into {integrated_csv}")
        except Exception as e:
            logger.warning(f"Could not merge into {integrated_csv}: {e}")

    result['genomes_updated'] = genomes_updated

    if progress_callback:
        progress_callback("Cross-genome orthologs", 100,
                          f"Done: {result['n_groups']} groups across "
                          f"{n_written} substrates from {genomes_updated} genomes")

    return result
