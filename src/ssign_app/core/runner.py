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
import threading
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Optional

from ssign_app.runtime import resolve_regime, scaled_timeout
from ssign_app.scripts.ssign_lib.constants import (
    DEFAULT_EXCLUDED_SYSTEMS,
    DEFAULT_TIER,
    HHSUITE_MIN_PROB,
    PLMBLAST_TIMEOUT_S,
    TIER_TOOL_DEFAULTS,
    TOOL_TIMEOUT_S,
)
from ssign_app.scripts.ssign_lib.dependency_manifest import find_db_by_env_var
from ssign_app.scripts.ssign_lib.fasta_io import count_sequences
from ssign_app.scripts.ssign_lib.resources import effective_cpu_count
from ssign_app.scripts.ssign_lib.substrates import load_substrate_ids

logger = logging.getLogger(__name__)

# The bin/ directory containing all pipeline scripts
# Scripts are packaged inside ssign_app/scripts/ (installed with pip)
# Fallback to bin/ for development mode
_PACKAGE_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
_DEV_BIN = Path(__file__).resolve().parents[3] / "bin"
BIN_DIR = _PACKAGE_SCRIPTS if _PACKAGE_SCRIPTS.exists() else _DEV_BIN

# Upper bound on how long we'll wait to enter the DTU rate-limit semaphore.
# 4h matches the longest-running per-step subprocess timeout, so a single
# stalled DTU job can't deadlock parallel-genome runs forever.
DTU_SEMAPHORE_TIMEOUT_S = 14400

# conda/mamba env roots scanned as a last-ditch binary lookup when PATH and
# configured paths both miss. Covers the conventional install locations; a
# user with a custom CONDA_ENVS_PATH still gets PATH-based discovery once
# they activate. Order matters only for log messages.
_CONDA_ENV_ROOTS = (
    "~/.conda/envs",
    "~/miniconda3/envs",
    "~/anaconda3/envs",
    "~/miniforge3/envs",
    "~/mambaforge/envs",
)


def _find_in_conda_envs(binary: str) -> Optional[str]:
    """Return the bin/ directory of a conda env holding ``binary``, or None.

    Scans the conventional env roots (see ``_CONDA_ENV_ROOTS``). Returns
    the directory rather than the full path so callers can hand it to
    ``signalp_path`` / ``deeplocpro_path`` unchanged — those fields already
    expect a directory, and the wrapper scripts append the binary name.
    """
    for root in _CONDA_ENV_ROOTS:
        envs_dir = os.path.expanduser(root)
        if not os.path.isdir(envs_dir):
            continue
        try:
            env_names = os.listdir(envs_dir)
        except OSError:
            continue
        for env_name in sorted(env_names):
            bin_dir = os.path.join(envs_dir, env_name, "bin")
            candidate = os.path.join(bin_dir, binary)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return bin_dir
    return None


def _detect_dtu_mode(binary: str, configured_path: str, display_name: str) -> tuple[str, Optional[str]]:
    """Pick local vs remote for a DTU-licensed tool.

    Tries (in order): explicit ``configured_path/binary`` (CLI flag or the
    matching ``SSIGN_*_PATH`` env var already filled in by __post_init__'s
    env loop), ``binary`` on PATH, then a scan of conventional conda env
    roots. Falls back to remote with a warning so an offline-only user
    doesn't silently get network-submitted jobs.

    Returns ``(mode, discovered_dir)``. ``discovered_dir`` is non-None only
    when the conda-env scan was the thing that succeeded, so the caller can
    persist it into the matching ``*_path`` config field; in every other
    case it's ``None``.
    """
    if configured_path:
        candidate = os.path.join(configured_path, binary)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return "local", None
    if shutil.which(binary):
        return "local", None
    conda_dir = _find_in_conda_envs(binary)
    if conda_dir:
        logger.info(
            "Auto-discovered %s in conda env: %s",
            display_name,
            os.path.join(conda_dir, binary),
        )
        return "local", conda_dir
    roots_for_msg = " or ".join(_CONDA_ENV_ROOTS)
    logger.warning(
        "No local %s binary found (%r not on PATH, no configured path holds it, "
        "and no conda env under %s has it). "
        "Falling back to the DTU webserver — pass --%s-mode local once you've "
        "installed it locally, or --%s-mode remote to silence this warning.",
        display_name,
        binary,
        roots_for_msg,
        display_name.lower(),
        display_name.lower(),
    )
    return "remote", None


def _resolve_db_root_for_runner() -> str:
    """Mirror of doctor's db-root resolver.

    Reads ``~/.ssign/db_root`` (written by ``scripts/fetch_databases.sh``)
    when present, falls back to ``~/.ssign/databases`` otherwise.

    Doctor has its own version because it takes ``data_root`` as a CLI flag;
    the runner doesn't expose that knob, so this helper is fixed to the
    standard home-dir location.
    """
    marker = os.path.expanduser("~/.ssign/db_root")
    if os.path.isfile(marker):
        try:
            with open(marker, encoding="utf-8") as f:
                recorded = f.read().strip()
            if recorded and os.path.isdir(recorded):
                return recorded
        except OSError:
            pass
    return os.path.expanduser("~/.ssign/databases")


def _read_tier_marker() -> Optional[str]:
    """Return the tier ``fetch_databases.sh`` recorded, or None.

    Looks at ``~/.ssign/tier`` (a one-line file with ``base`` /
    ``extended`` / ``full``). Returns None when the marker doesn't exist
    or holds an unrecognised value — caller falls back to ``DEFAULT_TIER``.
    """
    marker = os.path.expanduser("~/.ssign/tier")
    if not os.path.isfile(marker):
        return None
    try:
        with open(marker, encoding="utf-8") as f:
            value = f.read().strip().lower()
    except OSError:
        return None
    return value if value in TIER_TOOL_DEFAULTS else None


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
    # T3SS detected by default; DeepSecE is excluded from T3SS calls in
    # cross_validate_predictions (DSE unreliable on T3SS), so flagellar FPs stay out.
    excluded_systems: list = field(default_factory=lambda: list(DEFAULT_EXCLUDED_SYSTEMS))
    # MacSyFinder --db-type: "ordered_replicon" preserves gene-order signal
    # (more sensitive — recommended default). Switch to "unordered" for
    # highly fragmented assemblies where contig boundaries would otherwise
    # produce false-positive proximity calls.
    macsyfinder_db_type: str = "ordered_replicon"

    # CPU budget per genome for parallel sub-tools (e.g. macsyfinder -w).
    # When N genomes run concurrently, set this to cpu_per_genome // N to
    # avoid oversubscribing cores. Defaults to the cgroup-allocated CPU
    # count (not host total): os.cpu_count() reads /proc/cpuinfo and
    # ignores PBS/SLURM CPU pinning, so on a shared HPC node a 4-CPU job
    # would otherwise tell tools to spawn 24+ threads thrashing on 4 cores.
    cpu_per_genome: int = field(default_factory=effective_cpu_count)

    # Phase 3: Prediction
    conf_threshold: float = 0.8
    proximity_window: int = 3
    required_fraction_correct: float = 0.8
    # Localization-correctness gate (#58): per-component DLP max_probability
    # below this is treated as low-confidence and excluded from the gate's
    # numerator AND denominator. Separately tunable from conf_threshold,
    # which is the substrate-trigger DLP cutoff (different semantic).
    dlp_confidence_threshold: float = 0.8
    skip_localization_gate: bool = False

    # Phase 1: ORF prediction options
    # GenBank input is re-annotated through Bakta by default (Phase 3.3.c).
    # Original GenBank product strings are preserved as `gbff_annotation`
    # for annotation-consensus voting via map_gbff_to_bakta_cds.py.
    # Set use_input_annotations=True to skip Bakta and trust the input
    # annotations as-is — useful for hand-curated GenBank files.
    use_input_annotations: bool = False
    # FASTA contigs input always needs an ORF caller. run_bakta=True picks
    # Bakta (richer); otherwise Pyrodigal runs via extract_proteins.py.
    # Default True per plan addendum A.6 ("always run Bakta on the input").
    # Set --no-run-bakta to fall back to Pyrodigal (FASTA contigs only;
    # GenBank is governed by use_input_annotations).
    run_bakta: bool = True
    bakta_db: str = ""  # Required for any Bakta run (FASTA or GenBank re-annotation)
    # 0 => derive from cpu_per_genome in __post_init__. Bakta runs sequentially
    # in input processing, so it can take the full per-genome budget.
    bakta_threads: int = 0

    # --- Install tier (governs which optional tools default on) ------------
    # If unset (None), __post_init__ reads ~/.ssign/tier (written by
    # fetch_databases.sh) and falls back to DEFAULT_TIER. CLI --tier wins.
    # Every `skip_*` field below defaults to None ("use tier default");
    # __post_init__ resolves each one to a concrete bool via
    # TIER_TOOL_DEFAULTS. Passing --skip-X / --no-skip-X on the CLI bypasses
    # the tier-default lookup for that one tool.
    tier: Optional[str] = None

    # Phase 3: Tool paths (DTU licensed). ssign is offline-first: local
    # installs are the canonical path; remote submits to the DTU webserver
    # as a fallback. mode=None means "auto-detect" — __post_init__ picks
    # local if a binary is configured/on PATH, remote (with a warning)
    # otherwise. Explicit "local" / "remote" force the choice regardless.
    deeplocpro_mode: Optional[str] = None  # "local" / "remote" / None=auto
    deeplocpro_path: str = ""
    signalp_mode: Optional[str] = None
    signalp_path: str = ""
    skip_deeplocpro: Optional[bool] = None
    skip_signalp: Optional[bool] = None
    skip_deepsece: Optional[bool] = None
    dlp_whole_genome: bool = False  # Run on all proteins, not just neighborhood
    dse_whole_genome: bool = False
    sp_whole_genome: bool = False
    plme_whole_genome: bool = False

    # Resource sampling — writes outdir/runtime_data/resource_samples.csv during a run with
    # system CPU/RAM/GPU/disk samples tagged with the active step. Disable
    # if a tier-1 install is missing psutil or you don't want the
    # background thread.
    monitor_resources: bool = True
    monitor_interval_s: float = 5.0

    # --- Enrichment stats (opt-in) ------------------------------------------
    # When True, run the per-SS-type circular-shift permutation enrichment test
    # (enrichment_testing.py): fold + permutation p + BH q per system type, plus a
    # per-type null-distribution figure. This forces whole-genome DLP/DSE (see
    # __post_init__) because the rotation null needs every gene's positivity in gene
    # order. openspec: enrichment-circular-shift-per-run.
    enrichment_stats: bool = False

    # Phase 5: Annotation tools
    skip_blastp: Optional[bool] = None
    blastp_db: str = ""
    blastp_exclude_taxid: str = ""
    blastp_min_pident: float = 80.0
    blastp_min_qcov: float = 80.0
    blastp_evalue: float = 1e-5

    skip_hhsuite: Optional[bool] = None
    hhsuite_pfam_db: str = ""
    hhsuite_pdb70_db: str = ""
    hhsuite_uniclust_db: str = ""
    hhsuite_min_prob: float = HHSUITE_MIN_PROB

    skip_interproscan: Optional[bool] = None
    interproscan_db: str = ""
    interproscan_min_evalue: float = 1e-5

    skip_plmblast: Optional[bool] = None
    plmblast_db: str = ""
    # pLM-BLAST cosine percentile cutoff. 90 is the Kaminski 2023 paper
    # default — keeps top 10% of cosine pairs in pre-screening. Drop to
    # 70-80 for more permissive matching on short proteins (<200 aa where
    # cpc=90 often misses) at the cost of slower search.
    plmblast_cpc: int = 90

    # Phase 3.2.d: EggNOG-mapper (annotation-tier). Tier-driven default:
    # off at base (no DB shipped), on at extended/full.
    skip_eggnog: Optional[bool] = None
    eggnog_db: str = ""
    # None => let run_eggnog._autodetect_dbmem() decide based on host RAM
    # (44 GB resident, so only safe on >=50 GB hosts). Override with True/False.
    eggnog_dbmem: Optional[bool] = None
    # DIAMOND sensitivity. "sensitive" is ~10× DIAMOND default; "more-sensitive"
    # is ~2× that again and rarely rescues additional hits on the substrate
    # subset. See run_eggnog.run_emapper for the tradeoff.
    eggnog_sensmode: str = "sensitive"
    # Directory to stage the ~50 GB eggnog DB to before invoking emapper.
    # Mandatory on shared FS (gpfs/nfs/lustre): random-access mmap on the
    # 41 GB SQLite stalls for tens of minutes otherwise. Empty string ("")
    # resolves to PBS/SLURM $TMPDIR if set (local SSD on most HPCs).
    # Pass "off" (literal) to disable staging.
    eggnog_local_cache_dir: str = ""

    # Phase 3.2.d: PLM-Effector (prediction-tier, equal to DLP/DSE per
    # the cross-validate refactor in 3.2.b). Tier-driven default: on at
    # every tier (weights ship with the base bundle), but the user can
    # --skip-plm-effector on CPU-only nodes where the per-type ESM forward
    # is too slow.
    skip_plm_effector: Optional[bool] = None
    plm_effector_weights_dir: str = ""
    plm_effector_types: list = field(default_factory=lambda: ["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"])
    plm_chunk_size: int = 256

    skip_protparam: Optional[bool] = None

    # When True, run EggNOG/BLASTp/pLM-BLAST/HHsuite/ProtParam a SECOND time
    # on T5aSS substrates using the full protein FASTA (not the passenger-
    # substituted one) and emit the results under t5ass_whole_* columns in
    # the master CSV. Lets users compare passenger-only annotation (default,
    # often more functional) against full-AT annotation (often dominated by
    # the conserved β-barrel) side by side. IPS unchanged.
    t5ass_annotate_whole: bool = False

    # DSE type-match filter: remove DSE-only substrates where predicted
    # SS type doesn't match nearby MacSyFinder system
    filter_dse_type_mismatch: bool = True

    # Figures (per-genome curated set; pooled cross-genome figures are emitted
    # automatically by the multi-genome runner)
    dpi: int = 300
    fig_ss_comp: bool = True  # 01 secreted proteins by SS type (per genome for a group)
    fig_physicochemical: bool = True  # 02 size & physicochemical properties by type (incl. length)
    fig_func_summary: bool = True  # 03-06 functional categories by type (COG/KEGG/EggNOG/consensus)

    # DeepSecE threshold
    deepsece_min_prob: float = 0.8

    # SignalP threshold
    signalp_min_prob: float = 0.5

    # Ortholog group thresholds
    ortholog_min_pident: float = 40.0
    ortholog_min_qcov: float = 70.0

    def __post_init__(self) -> None:
        # Resolve the install tier first (every skip_* default depends on it):
        # explicit `tier` wins, then the marker at ~/.ssign/tier written by
        # fetch_databases.sh, then DEFAULT_TIER.
        if self.tier is None:
            self.tier = _read_tier_marker() or DEFAULT_TIER
        if self.tier not in TIER_TOOL_DEFAULTS:
            raise ValueError(f"Unknown tier {self.tier!r}; expected one of {sorted(TIER_TOOL_DEFAULTS)}")

        # Resolve each unset `skip_*` field via TIER_TOOL_DEFAULTS. The
        # dataclass field name is `skip_<tool>`, the tier-defaults key is
        # `<tool>`, and the tier table is keyed by "is enabled?" (True),
        # which inverts to skip=False. If the user explicitly set the flag
        # (True or False on the CLI), the dataclass arrives here non-None
        # and we leave it alone.
        tier_defaults = TIER_TOOL_DEFAULTS[self.tier]
        for tool, enabled_default in tier_defaults.items():
            field_name = f"skip_{tool}"
            if getattr(self, field_name) is None:
                setattr(self, field_name, not enabled_default)

        # bakta_threads=0 sentinel → take the full per-genome budget.
        if self.bakta_threads == 0:
            self.bakta_threads = self.cpu_per_genome

        # Circular-shift enrichment needs every gene's positivity in gene order, so
        # turning on --enrichment-stats forces whole-genome DLP/DSE/SignalP. That is
        # the runtime cost the flag warns about (~13 min/genome on a ~5k-gene
        # proteome). SignalP runs locally; no webserver. openspec:
        # enrichment-circular-shift-per-run, signalp-enrichment-track.
        if self.enrichment_stats and not (self.dlp_whole_genome and self.dse_whole_genome and self.sp_whole_genome):
            self.dlp_whole_genome = True
            self.dse_whole_genome = True
            self.sp_whole_genome = True
            logger.info(
                "Enrichment stats on: running DeepLocPro + DeepSecE + SignalP on the whole "
                "genome (needed for the circular-shift null; adds ~13 min/genome)."
            )

        # Step A — env-var verbatim. Trust whatever path the user (or HPC
        # session script) exported, even if the layout doesn't match
        # ssign's expectations yet. Tools that ship under their own
        # conventional env var (BAKTA_DB, EGGNOG_DATA_DIR) and our
        # SSIGN_*_DB family both get a chance here. Documented in the
        # matching CLI --*-db / --*-dir help and docs/how-to/install.md.
        for attr, env in (
            ("bakta_db", "BAKTA_DB"),
            ("hhsuite_pfam_db", "SSIGN_HHSUITE_PFAM"),
            ("hhsuite_pdb70_db", "SSIGN_HHSUITE_PDB70"),
            ("hhsuite_uniclust_db", "SSIGN_HHSUITE_UNICLUST"),
            ("interproscan_db", "SSIGN_INTERPROSCAN_PATH"),
            ("eggnog_db", "EGGNOG_DATA_DIR"),
            ("plmblast_db", "SSIGN_ECOD_DB"),
            ("plm_effector_weights_dir", "SSIGN_PLM_EFFECTOR_WEIGHTS"),
            ("signalp_path", "SSIGN_SIGNALP_PATH"),
            ("deeplocpro_path", "SSIGN_DEEPLOCPRO_PATH"),
        ):
            if not getattr(self, attr) and os.environ.get(env):
                value = os.environ[env]
                setattr(self, attr, value)
                logger.info("Using %s from env var %s: %s", attr, env, value)

        # Step B — marker fill-in + sentinel-driven descent. Uses
        # DatabasePath.resolve_path (the same resolver doctor consumes),
        # which globs the sentinel inside each candidate path and returns
        # the dir containing the first match. For each DB:
        #   - If the attr is unset, resolve_path returns the right path
        #     under the marker root.
        #   - If the attr is set (Step A) but points at the parent of the
        #     actual DB dir (e.g. BAKTA_DB=<dir>/bakta when the version-
        #     stamped subdir is at <dir>/bakta/db-light), resolve_path
        #     descends into the right inner dir and we update the attr.
        # This is what closes the "doctor green but Bakta crashes" gap.
        db_root = _resolve_db_root_for_runner()
        for attr, manifest_env_var in (
            ("bakta_db", "SSIGN_BAKTA_DB"),
            ("hhsuite_pfam_db", "SSIGN_HHSUITE_PFAM"),
            ("hhsuite_pdb70_db", "SSIGN_HHSUITE_PDB70"),
            ("hhsuite_uniclust_db", "SSIGN_HHSUITE_UNICLUST"),
            ("interproscan_db", "SSIGN_INTERPROSCAN_PATH"),
            ("eggnog_db", "SSIGN_EGGNOG_DB"),
            ("plmblast_db", "SSIGN_ECOD_DB"),
        ):
            entry = find_db_by_env_var(manifest_env_var)
            if entry is None:
                continue
            current = getattr(self, attr)
            extras = (current,) if current else ()
            resolved = entry.resolve_path(db_root, *extras)
            if resolved and resolved != current:
                setattr(self, attr, resolved)
                if current:
                    logger.info("Normalized %s: %s → %s", attr, current, resolved)
                else:
                    logger.info("Resolved %s → %s", attr, resolved)

        # BLASTp's DB is addressed by a name *prefix* (blastp -db
        # .../blast_swissprot/swissprot), not the directory resolve_path
        # returns. Resolve the Swiss-Prot dir (via SSIGN_BLAST_SWISSPROT or the
        # marker root) and append the `swissprot` basename. Swiss-Prot is the
        # full-tier default (curated + fast); NR — too slow on a real substrate
        # set — or any custom DB is opt-in via --blastp-db, which wins because
        # this only fills an unset blastp_db.
        if not self.blastp_db:
            sp_entry = find_db_by_env_var("SSIGN_BLAST_SWISSPROT")
            sp_dir = sp_entry.resolve_path(db_root) if sp_entry else None
            if sp_dir:
                self.blastp_db = os.path.join(sp_dir, "swissprot")
                logger.info("Resolved blastp_db → %s", self.blastp_db)

        # PLM-Effector weights also live under db_root when the user used
        # fetch_databases.sh's default layout. ModelWeights entry with
        # under_db_root=True; not part of DATABASE_PATHS so handled here.
        if not self.plm_effector_weights_dir:
            candidate = os.path.join(db_root, "plm_effector_weights")
            if os.path.isdir(candidate):
                self.plm_effector_weights_dir = candidate
                logger.info("Resolved plm_effector_weights_dir → %s", candidate)

        # Auto-detect DTU-tool mode when the user didn't pin it. Defaults
        # to local if a binary is discoverable (configured path or PATH);
        # falls back to remote with a warning otherwise so an offline-only
        # install doesn't silently start network-submitting jobs.
        if self.signalp_mode is None:
            self.signalp_mode, discovered = _detect_dtu_mode("signalp6", self.signalp_path, "SignalP")
            if discovered and not self.signalp_path:
                self.signalp_path = discovered
        if self.deeplocpro_mode is None:
            self.deeplocpro_mode, discovered = _detect_dtu_mode("deeplocpro", self.deeplocpro_path, "DeepLocPro")
            if discovered and not self.deeplocpro_path:
                self.deeplocpro_path = discovered


def _rename_csv_columns_with_prefix(src_path: str, dst_path: str, prefix: str, keep_columns: set[str]) -> None:
    """Copy a CSV/TSV from ``src_path`` to ``dst_path``, prefixing column
    headers except those in ``keep_columns``.

    Used by the t5ass_whole annotation second pass so each tool's output
    merges as side-by-side columns in the integrated CSV (e.g.
    ``eggnog_description`` vs ``t5ass_whole_eggnog_description``)
    instead of clobbering the passenger annotations on join.

    Delimiter for both files is inferred from the path extension via
    ``_pool_utils._delimiter_for`` so CSV/TSV mixed roundtrips keep
    each format. The ``keep_columns`` set typically holds the join
    column (``locus_tag``) that integrate_annotations merges on.
    """
    import csv as _csv

    from ssign_app.core._pool_utils import _delimiter_for

    in_delim = _delimiter_for(src_path)
    out_delim = _delimiter_for(dst_path)

    with open(src_path) as f:
        reader = _csv.DictReader(f, delimiter=in_delim)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    rename_map = {c: (c if c in keep_columns else f"{prefix}{c}") for c in fieldnames}
    new_fieldnames = [rename_map[c] for c in fieldnames]

    with open(dst_path, "w", newline="") as f:
        writer = _csv.DictWriter(f, fieldnames=new_fieldnames, delimiter=out_delim)
        writer.writeheader()
        for row in rows:
            writer.writerow({rename_map[k]: v for k, v in row.items()})


@dataclass
class StepResult:
    """Result of a pipeline step."""

    name: str
    success: bool
    message: str
    output_files: dict = field(default_factory=dict)
    duration_s: float = 0.0


def run_script(
    script_name: str,
    args: list,
    timeout: int = 7200,
    *,
    stream_stderr: bool = False,
) -> tuple:
    """Run a bin/ script with arguments. Returns (returncode, stdout, stderr).

    ``stream_stderr=True`` forwards each stderr line to the runner logger as it
    arrives, instead of buffering until the subprocess exits. Use it for long
    steps where intermediate progress (e.g. per-PLM-type completion in
    PLM-Effector) is what tells the user "still alive, on type N of 5". Default
    off — most tool wrappers are silent on success or chatty enough to flood
    the log.
    """
    script_path = BIN_DIR / script_name
    if not script_path.exists():
        return (-1, "", f"Script not found: {script_path}")

    cmd = [sys.executable, str(script_path)] + args
    logger.info(f"Running: {' '.join(cmd[:4])}...")

    if stream_stderr:
        return _run_script_streaming(cmd, script_name, timeout)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode != 0:
            _log_nonzero_exit(script_name, result.returncode, result.stdout, result.stderr)
        return (result.returncode, result.stdout, result.stderr)
    except subprocess.TimeoutExpired:
        logger.warning(
            "%s hit its %ss timeout and was killed — a hard cap, not a tool crash. If the input is "
            "legitimately large, the size-aware cap (openspec size-aware-tool-timeouts) needs a wider margin.",
            script_name,
            timeout,
        )
        return (-1, "", f"Timeout after {timeout}s (killed {script_name})")
    except Exception as e:
        # Log the traceback so the underlying cause (PermissionError, OSError,
        # encoding glitches in subprocess output) is recoverable from the log
        # even though the caller only sees the str(e) summary.
        logger.exception("run_script(%s) raised unexpectedly", script_name)
        return (-1, "", str(e))


def _log_nonzero_exit(script_name: str, rc: int, stdout: str, stderr: str) -> None:
    """Log a script's nonzero exit with the tail of both streams.

    StepResult callers truncate stderr to ~500 chars; benign HF / tokenizer
    warnings eat that budget before the real traceback. Logging the tail
    keeps the stack frame around without exploding the run log when tools
    (IPS, BLASTp) dump megabytes on failure.
    """
    logger.error(
        "%s exited with code %s\n--- stdout (tail) ---\n%s\n--- stderr (tail) ---\n%s",
        script_name,
        rc,
        stdout[-8000:],
        stderr[-8000:],
    )


def _run_script_streaming(cmd: list, script_name: str, timeout: int) -> tuple:
    """Popen variant of run_script that forwards stderr lines to the logger live.

    Both stdout and stderr must be drained in separate threads so a full pipe
    buffer on one stream doesn't deadlock the child. ``PYTHONUNBUFFERED=1`` in
    the child env defeats Python's default block-buffering of piped stderr so
    lines reach us as the child emits them, not in 4 KB chunks at exit.
    ``errors='replace'`` matters on HPC nodes whose locale is ``C``/``POSIX``;
    the default strict decoder would raise inside a daemon drain thread and
    silently kill it on the first non-ASCII byte.
    """
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
            env=env,
        )
    except (OSError, ValueError) as e:
        logger.exception("run_script(%s) failed to start subprocess", script_name)
        return (-1, "", str(e))

    stdout_chunks: list = []
    stderr_chunks: list = []

    def _drain(pipe, sink, log_prefix=None):
        for line in pipe:
            sink.append(line)
            if log_prefix is not None:
                msg = line.rstrip("\n")
                if msg:
                    logger.info("[%s] %s", log_prefix, msg)

    t_out = threading.Thread(target=_drain, args=(proc.stdout, stdout_chunks), daemon=True)
    t_err = threading.Thread(target=_drain, args=(proc.stderr, stderr_chunks, script_name), daemon=True)
    t_out.start()
    t_err.start()

    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        t_out.join(timeout=1)
        t_err.join(timeout=1)
        logger.warning(
            "%s hit its %ss timeout and was killed — a hard cap, not a tool crash. If the input is "
            "legitimately large, the size-aware cap (openspec size-aware-tool-timeouts) needs a wider margin.",
            script_name,
            timeout,
        )
        return (-1, "".join(stdout_chunks), f"Timeout after {timeout}s (killed {script_name})")

    t_out.join()
    t_err.join()

    rc = proc.returncode
    stdout_text = "".join(stdout_chunks)
    stderr_text = "".join(stderr_chunks)

    if rc != 0:
        _log_nonzero_exit(script_name, rc, stdout_text, stderr_text)
    return (rc, stdout_text, stderr_text)


class PipelineRunner:
    """Runs the ssign pipeline step by step with progress callbacks."""

    def __init__(
        self,
        config: PipelineConfig,
        progress_callback: Optional[Callable] = None,
        api_semaphores: Optional[dict] = None,
    ):
        self.config = config
        self._progress_callback = progress_callback or (lambda step, pct, msg: None)
        # Clamp percentages to monotonically non-decreasing — parallel-group
        # `as_completed` order is non-deterministic, so the raw `100*sc/total`
        # of a step that finishes after a higher-ordinal sibling could regress
        # the displayed bar. Hold the highest pct seen on this runner.
        self._max_progress_pct = 0
        self.results: list[StepResult] = []
        self.work_dir = ""
        self.files = {}  # Track intermediate file paths
        self.start_time: float | None = None
        # Lazily set in run() if monitor_resources is on. Initialised here
        # so callers that drive _execute_stages directly (MultiGenomeRunner)
        # don't AttributeError on the `if self._sampler is not None` checks
        # inside the step loop.
        self._sampler = None
        # Per-API semaphores for multi-genome concurrency control.
        # Keys: 'dtu', 'ncbi', 'mpi', 'ebi'. Values: threading.Semaphore
        self.api_sem = api_semaphores or {}
        # Guards self.results.append(...) + reads of self.files used for
        # progress JSON snapshots. Step functions write to self.files
        # under distinct keys (atomic under the GIL), but the result
        # list and progress-JSON serialiser need a coarse lock so a
        # parallel-group append doesn't race with another thread's
        # _save_progress() iteration.
        self._state_lock = threading.Lock()

    def progress(self, step: str, pct: int, msg: str) -> None:
        with self._state_lock:
            if pct < self._max_progress_pct:
                pct = self._max_progress_pct
            else:
                self._max_progress_pct = pct
        self._progress_callback(step, pct, msg)

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
                "macsydata not found on PATH (needed to install TXSScan models).\n  Install: pip install macsyfinder"
            )

        # ── Optional dependencies (checked based on config) ──

        # DeepSecE: large optional dependency (torch + ESM + DeepSecE)
        if not self.config.skip_deepsece:
            for pkg, install in [
                ("torch", "pip install 'ssign[deepsece]'  # or: pip install torch"),
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

    def check_required_executables(self) -> list[str]:
        """Pre-flight gate for optional tools that have no pip-extra.

        EggNOG-mapper and pLM-BLAST are deliberately outside ssign's pip
        extras (biopython hard-pin clash; alntools not on PyPI) so they
        are easy to forget to install. Their per-step crashes surface
        late: in the June 2026 K-12 PBS run, EggNOG crashed 56 minutes
        in. Catching the missing binary here aborts before any work.

        Returns a list of fatal error messages. An empty list means
        every enabled optional step has its binary on PATH.
        """
        errors = []

        # Match the conda-env auto-discovery already used for signalp6 /
        # deeplocpro — users who follow the docs' "conda install -n
        # bakta-deps -c bioconda eggnog-mapper" path get a binary that
        # isn't on PATH unless they activate the env first.
        if not self.config.skip_eggnog and not (shutil.which("emapper.py") or _find_in_conda_envs("emapper.py")):
            errors.append(
                "emapper.py not found on PATH or in any conda env (needed for EggNOG-mapper).\n"
                "  Install (recommended, avoids biopython pin clash):\n"
                "    pip install --no-deps eggnog-mapper\n"
                "  Or via conda:\n"
                "    conda install -c bioconda eggnog-mapper\n"
                "  Or pass --skip-eggnog to disable this step."
            )

        if not self.config.skip_plmblast:
            from ssign_app.scripts.run_plm_blast import find_plmblast_script

            if not find_plmblast_script():
                errors.append(
                    "plmblast.py not found on PATH and SSIGN_PLMBLAST_SCRIPT "
                    "is unset (needed for pLM-BLAST).\n"
                    "  Install: git clone https://github.com/labstructbioinf/pLM-BLAST.git\n"
                    "  Then: export SSIGN_PLMBLAST_SCRIPT=/path/to/pLM-BLAST/scripts/plmblast.py\n"
                    "  Or pass --skip-plmblast to disable this step."
                )

        return errors

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

        # Pre-flight: hard-fail on missing optional binaries before any
        # step runs. Catches the "forgot to install emapper.py / pLM-BLAST"
        # class of bug that previously wasted ~1h of HPC time per genome.
        dep_errors = self.check_required_executables()
        if dep_errors:
            for e in dep_errors:
                logger.error(f"Pre-flight: {e}")
            raise RuntimeError(
                "Required tool(s) for enabled optional steps are missing. "
                "Install them or pass the matching --skip-* flag. "
                f"{len(dep_errors)} issue(s) — see log for install instructions."
            )

        # Pre-flight dependency check
        dep_warnings = self.check_dependencies()
        if dep_warnings:
            for w in dep_warnings:
                logger.warning(f"Dependency check: {w}")
            self.progress(
                "Dependency check",
                0,
                f"{len(dep_warnings)} warning(s) — see log for details",
            )

        # Create output directory
        os.makedirs(self.config.outdir, exist_ok=True)

        # Spin up the in-process resource sampler. Writes
        # runtime_data/resource_samples.csv to outdir in real time
        # (system CPU/RAM/GPU/disk, tagged with the active step).
        # Daemon thread; sampling failures never crash the pipeline.
        # Disable with `monitor_resources=False`.
        self._sampler = None
        if self.config.monitor_resources:
            from ssign_app.core.resource_sampler import ResourceSampler

            self._sampler = ResourceSampler(
                out_path=os.path.join(self.config.outdir, "runtime_data", "resource_samples.csv"),
                interval=self.config.monitor_interval_s,
            )
            self._sampler.start()

        # Try to resume from previous progress
        skip_steps = set()
        if resume:
            skip_steps = self._try_resume()

        # Only create new work_dir if not resuming
        if not self.work_dir:
            self.work_dir = tempfile.mkdtemp(prefix="ssign_")

        stages = self._build_stages()

        core_failed = self._execute_stages(stages, skip_steps)

        # Copy final outputs to outdir BEFORE reporting 100% — otherwise the
        # progress bar hits 100 while files are still being written.
        # runtime_data/step_timings.csv is up-to-date already; _record_result
        # wrote it incrementally as each step finished.
        self._copy_outputs()
        self._save_progress()
        if self._sampler is not None:
            self._sampler.stop()
        self.progress("Complete", 100, f"Pipeline finished in {self._elapsed_str()}")

        # On clean success, drop the temp work_dir. On failure, retain it so
        # the user can `--resume` after fixing the underlying issue.
        if not core_failed and self.work_dir and os.path.isdir(self.work_dir):
            try:
                shutil.rmtree(self.work_dir)
                logger.info(f"Cleaned up work directory: {self.work_dir}")
            except OSError as e:
                logger.warning(f"Could not remove work_dir {self.work_dir}: {e}")

        return self.results

    def _build_stages(self) -> list:
        """Return the ordered pipeline stages list (memoised on first call).

        Each element is either a ``(name, callable)`` tuple (sequential
        step) or a list of such tuples (parallel group). Steps disabled
        in ``self.config`` (skip_*) are omitted.

        Dependency graph:
          detect_format → extract_proteins → macsyfinder → validate →
          extract_neighborhood →
            [deeplocpro || deepsece || signalp]  ← PARALLEL GROUP 1
          → plm_effector →
          → cross_validate → proximity → t5ss → filtering →
            [blastp || hhsuite || interproscan || plm_blast || protparam]  ← PARALLEL GROUP 2
          → integrate → orthologs → enrichment → report → figures

        PLM-Effector is deliberately serial after parallel group 1: it
        loads 3-4 PLMs (~6 GB resident each), and on a memory-constrained
        allocation the parallel group's shared memory ceiling kills its
        first model load while DLP/DSE/SignalP are still resident.

        Memoised because MultiGenomeRunner re-asks each per-genome
        runner for its stages once per per-genome segment (3 calls);
        recomputing each time is cheap (~ms) but pointless when config
        is frozen post-init.
        """
        cached = getattr(self, "_cached_stages", None)
        if cached is not None:
            return cached
        prediction_steps = []
        if not self.config.skip_deeplocpro:
            prediction_steps.append(("Predicting localization (DeepLocPro)", self._step_deeplocpro))
        if not self.config.skip_deepsece:
            prediction_steps.append(("Predicting secretion type (DeepSecE)", self._step_deepsece))
        if not self.config.skip_signalp:
            prediction_steps.append(("Predicting signal peptides (SignalP)", self._step_signalp))

        # HHpred listed first so it starts immediately — it's the bottleneck
        annotation_steps = []
        if not self.config.skip_hhsuite:
            annotation_steps.append(("Running HH-suite", self._step_hhsuite))
        if not self.config.skip_blastp:
            annotation_steps.append(("Running BLASTp", self._step_blastp))
        if not self.config.skip_interproscan:
            annotation_steps.append(("Running InterProScan", self._step_interproscan))
        if not self.config.skip_eggnog:
            annotation_steps.append(("Running EggNOG-mapper", self._step_eggnog))
        if not self.config.skip_plmblast:
            annotation_steps.append(("Running pLM-BLAST", self._step_plm_blast))
        if not self.config.skip_protparam:
            annotation_steps.append(("Computing physicochemical properties", self._step_protparam))

        stages = [
            ("Detecting input format", self._step_detect_format),
            ("Extracting proteins", self._step_extract_proteins),
            ("Running MacSyFinder", self._step_macsyfinder),
            ("Validating secretion systems", self._step_validate_systems),
            ("Extracting SS neighborhood", self._step_extract_neighborhood),
        ]
        # The circular-shift enrichment test needs no null sample: enrichment_stats
        # forces whole-genome DLP/DSE (see Config.__post_init__), so every gene
        # already has a prediction and the rotation null is the background.
        stages.append(prediction_steps)
        if not self.config.skip_plm_effector:
            stages.append(("Predicting effectors (PLM-Effector, 5 types)", self._step_plm_effector))
        stages.extend(
            [
                ("Cross-validating predictions", self._step_cross_validate),
                ("Running proximity analysis", self._step_proximity),
                ("Handling T5SS autotransporters", self._step_t5ss),
                ("Filtering systems", self._step_filtering),
                ("Building passenger-substituted FASTA", self._step_build_passenger_fasta),
                annotation_steps,
                ("Re-annotating T5aSS substrates whole-protein", self._step_t5ass_whole_annotations),
                ("Integrating annotations", self._step_integrate),
                ("Assigning ortholog groups", self._step_orthologs),
                ("Running enrichment analysis", self._step_enrichment),
                ("Generating report", self._step_report),
                ("Generating figures", self._step_figures),
            ]
        )
        if self.config.enrichment_stats:
            stages.append(("Generating enrichment figure", self._step_enrichment_figure))
        self._cached_stages = stages
        return stages

    def _execute_stages(self, stages: list, skip_steps: set) -> bool:
        """Execute the given stages list. Returns True iff a core step failed.

        Steps in ``skip_steps`` are recorded as resumed without re-running,
        EXCEPT downstream steps (integrate/orthologs/enrichment/report/
        figures) which are forced to re-run if any upstream step in this
        stages list ran (their outputs depend on upstream results).

        Progress percentages are relative to the passed stages list. For
        single-genome runs this is the full pipeline; for multi-genome
        callers a slice represents one phase.
        """
        all_steps = []
        for stage in stages:
            if isinstance(stage, list):
                all_steps.extend(stage)
            else:
                all_steps.append(stage)

        CORE_STEPS = {
            "detect_format",
            "extract_proteins",
            "macsyfinder",
            "validate_systems",
            "extract_neighborhood",
            "deeplocpro",
            "cross_validate",
            "proximity",
            "t5ss",
            "filtering",
            "integrate",
        }

        total = len(all_steps)
        core_failed = False
        n_skipped = 0
        step_counter = 0
        any_step_ran = False  # Track if any step ran (for forcing downstream re-runs)

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

                # Downstream steps (integrate, orthologs, enrichment, report,
                # figures) must re-run if ANY upstream step ran, since their
                # output depends on all upstream results
                _downstream = {
                    "integrate",
                    "orthologs",
                    "enrichment",
                    "report",
                    "figures",
                }
                _force_rerun = step_id in _downstream and any_step_ran

                if step_id in skip_steps and not _force_rerun:
                    n_skipped += 1
                    pct = int(100 * step_counter / total)
                    self.progress(
                        name,
                        pct,
                        f"Skipped (already done) | {self._elapsed_str()} elapsed",
                    )
                    self._record_result(StepResult(step_id, True, "Resumed (already completed)"))
                    continue

                if core_failed:
                    self._record_result(StepResult(step_id, False, "Skipped (earlier core step failed)"))
                    continue

                steps_to_run.append((name, func, step_id, step_counter))
                any_step_ran = True

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
                    f"[ssign] [{self.config.sample_id}] Starting parallel group: {step_names}",
                    flush=True,
                )
                if self._sampler is not None:
                    self._sampler.set_step(f"parallel:{','.join(sorted(n for n, _, _, _ in steps_to_run))}")

                # Tell every wrapper subprocess how many siblings it's
                # racing against. Every wrapper that calls resolve_threads()
                # (DLP/DSE/SignalP in the prediction group; IPS/EggNOG/
                # BLASTp/pLM-BLAST in the annotation group) reads
                # SSIGN_PARALLEL_GROUP_SIZE via parallel_share_cpus() and
                # divides the CPU budget evenly instead of each grabbing
                # the whole allocation. The context manager restores the
                # env var when the group exits (or raises).
                from ssign_app.scripts.ssign_lib.resources import parallel_group

                with parallel_group(len(steps_to_run)), ThreadPoolExecutor(max_workers=len(steps_to_run)) as executor:
                    futures = {}
                    for name, func, step_id, sc in steps_to_run:
                        futures[executor.submit(func)] = (name, step_id, sc, time.monotonic())

                    for future in as_completed(futures):
                        name, step_id, sc, t_step_start = futures[future]
                        try:
                            result = future.result()
                            result.duration_s = time.monotonic() - t_step_start
                            with self._state_lock:
                                self._record_result(result)
                            pct = int(100 * sc / total)
                            print(
                                f"[ssign] [{self.config.sample_id}] Finished (parallel): {name} -> "
                                f"{'OK' if result.success else 'FAILED: ' + result.message[:100]}",
                                flush=True,
                            )
                            if result.success:
                                self.progress(
                                    name,
                                    pct,
                                    f"Done: {result.message} | {self._elapsed_str()} elapsed",
                                )
                            else:
                                self.progress(name, pct, f"Failed: {result.message}")
                                logger.error(f"Step '{name}' failed: {result.message}")
                                if step_id in CORE_STEPS:
                                    core_failed = True
                        except Exception as e:
                            print(
                                f"[ssign] [{self.config.sample_id}] EXCEPTION (parallel): {name} -> {e}",
                                flush=True,
                            )
                            with self._state_lock:
                                self._record_result(
                                    StepResult(
                                        step_id,
                                        False,
                                        str(e),
                                        duration_s=time.monotonic() - t_step_start,
                                    )
                                )
                            logger.exception(f"Step '{name}' raised exception")
                            if step_id in CORE_STEPS:
                                core_failed = True

                self._save_progress()
            else:
                # Sequential execution (single step or single-element group)
                for name, func, step_id, sc in steps_to_run:
                    pct = int(100 * sc / total)
                    self.progress(name, pct, f"Step {sc}/{total} | {self._elapsed_str()} elapsed")
                    print(
                        f"[ssign] [{self.config.sample_id}] Starting step {sc}/{total}: {name} ({step_id})",
                        flush=True,
                    )
                    if self._sampler is not None:
                        self._sampler.set_step(f"{sc}/{total}:{step_id}")

                    t_step = time.monotonic()
                    try:
                        result = func()
                        result.duration_s = time.monotonic() - t_step
                        self._record_result(result)
                        self._save_progress()
                        print(
                            f"[ssign] [{self.config.sample_id}] Finished step {sc}/{total}: "
                            f"{name} -> "
                            f"{'OK' if result.success else 'FAILED: ' + result.message[:100]}",
                            flush=True,
                        )
                        if result.success:
                            self.progress(
                                name,
                                pct,
                                f"Done: {result.message} | {self._elapsed_str()} elapsed",
                            )
                        else:
                            self.progress(name, pct, f"Failed: {result.message}")
                            logger.error(f"Step '{name}' failed: {result.message}")
                            if step_id in CORE_STEPS:
                                core_failed = True
                    except Exception as e:
                        print(
                            f"[ssign] [{self.config.sample_id}] EXCEPTION in step {sc}/{total}: {name} -> {e}",
                            flush=True,
                        )
                        self._record_result(StepResult(step_id, False, str(e), duration_s=time.monotonic() - t_step))
                        self._save_progress()
                        logger.exception(f"Step '{name}' raised exception")
                        if step_id in CORE_STEPS:
                            core_failed = True

        if n_skipped:
            logger.info(f"Resumed: skipped {n_skipped} previously completed steps")

        return core_failed

    def _record_result(self, result: "StepResult") -> None:
        """Append a step result and flush runtime_data/step_timings.csv to disk.

        Callers in the parallel group hold ``self._state_lock`` already;
        callers in the sequential path are single-threaded. Either way
        the append + write pair must not be split across an OS-level
        kill, which is why every site uses this helper instead of the
        bare append.
        """
        self.results.append(result)
        self._write_step_timings()

    def _write_step_timings(self) -> None:
        """Write per-step wallclock data to outdir/runtime_data/step_timings.csv.

        Companion to runtime_data/resource_samples.csv: that one is
        sampled over time; step_timings is one row per step with
        start/duration. Together they let the operator slice the
        per-second sampler data by step without grep-tagging.

        Called after every step finishes (not just at pipeline end) so
        a mid-pipeline crash or PBS walltime kill still leaves timings
        for the steps that did complete. Writes via temp-file +
        os.replace so a kill mid-write can't leave a corrupt CSV.
        Caller holds self._state_lock when writing from inside a
        parallel-group, so two threads can't race here.
        """
        if not self.results:
            return
        import csv

        out_path = os.path.join(self.config.outdir, "runtime_data", "step_timings.csv")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        tmp_path = f"{out_path}.tmp.{os.getpid()}"
        try:
            with open(tmp_path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["step", "success", "duration_s", "message"])
                for r in self.results:
                    w.writerow([r.name, int(r.success), f"{r.duration_s:.2f}", r.message[:200]])
            os.replace(tmp_path, out_path)
        except OSError as e:
            logger.warning("could not write step_timings.csv: %s", e)

    def _wf(self, name):
        """Get work file path."""
        return os.path.join(self.work_dir, name)

    # ── Phase 1: Input Processing ──

    def _step_detect_format(self) -> StepResult:
        rc, stdout, stderr = run_script(
            "detect_input_format.py",
            [
                self.config.input_path,
            ],
        )
        if rc == 0:
            fmt = stdout.strip()
            self.files["format"] = fmt
            self.files["input"] = self.config.input_path
            return StepResult("detect_format", True, f"Format: {fmt}")
        return StepResult("detect_format", False, stderr[:500])

    def _genbank_to_contigs_fasta(self, gbff_path: str, out_fasta: str) -> None:
        """Write a contigs FASTA from a GenBank file (one record = one contig).

        Bakta only ingests nucleotide FASTA, so to re-annotate GenBank
        input we strip out the contig sequences first. SeqIO handles the
        format conversion natively.
        """
        from Bio import SeqIO

        with open(out_fasta, "w") as f:
            for record in SeqIO.parse(gbff_path, "genbank"):
                f.write(f">{record.id}\n{str(record.seq)}\n")

    def _find_companion_fasta(self, gff_path: str) -> Optional[str]:
        """Locate the genomic FASTA that pairs with a GFF3 input.

        NCBI ships GFF3 + FASTA together with the same stem (e.g.
        ``ecoli_k12.gff`` + ``ecoli_k12.fna``); extract_proteins.py
        translates CDS coordinates against the FASTA, so the runner
        has to hand both files over.
        """
        base, _ = os.path.splitext(gff_path)
        for ext in (".fna", ".fa", ".fasta"):
            candidate = base + ext
            if os.path.isfile(candidate):
                return candidate
        return None

    def _run_extract_proteins_script(self, proteins_out: str, gene_info_out: str, metadata_out: str) -> tuple:
        """Invoke extract_proteins.py with this run's input + sample id.

        Returns (rc, stderr). All three output paths are required by the
        script.
        """
        args = [
            "--input",
            self.config.input_path,
            "--sample",
            self.config.sample_id,
            "--out-proteins",
            proteins_out,
            "--out-gene-info",
            gene_info_out,
            "--out-metadata",
            metadata_out,
        ]
        if self.config.original_filename:
            args.extend(["--original-filename", self.config.original_filename])
        if self.files.get("format", "") == "gff3":
            fasta = self._find_companion_fasta(self.config.input_path)
            if not fasta:
                return 1, (
                    "GFF3 input requires a companion nucleotide FASTA "
                    "with the same stem (e.g. ecoli_k12.gff + ecoli_k12.fna). "
                    "Place one alongside the GFF3 file."
                )
            args.extend(["--fasta", fasta])
        rc, _, stderr = run_script("extract_proteins.py", args)
        return rc, stderr

    def _step_extract_proteins(self) -> StepResult:
        proteins_path = self._wf(f"{self.config.sample_id}_proteins.faa")
        gene_info_path = self._wf(f"{self.config.sample_id}_gene_info.tsv")
        metadata_path = self._wf(f"{self.config.sample_id}_metadata.json")
        fmt = self.files.get("format", "")

        # Decide the path:
        # - GenBank with use_input_annotations=False (default): re-annotate
        #   via Bakta, preserve original products as gbff_annotation.
        # - GenBank with use_input_annotations=True: parse GenBank only.
        # - FASTA contigs with run_bakta=True: Bakta directly (no overlap map).
        # - FASTA contigs with run_bakta=False: Pyrodigal via extract_proteins.
        # - protein_fasta or gff3: extract_proteins handles them; Bakta can't.
        reannotate_gbff = fmt == "genbank" and not self.config.use_input_annotations
        bakta_only = fmt == "fasta_contigs" and self.config.run_bakta

        if reannotate_gbff or bakta_only:
            if not self.config.bakta_db:
                opt_out = (
                    "--use-input-annotations (preserve the input's annotations)"
                    if reannotate_gbff
                    else "--no-run-bakta (fall back to Pyrodigal)"
                )
                return StepResult(
                    "extract_proteins",
                    False,
                    f"Bakta runs by default (per plan A.6) and requires --bakta-db "
                    f"(or BAKTA_DB env var). Pass {opt_out} to skip Bakta. "
                    f"Download a DB with: bakta_db download --output /path/to/db --type light",
                )

            # For GenBank: capture the original annotations into a side
            # gene_info.tsv first (to be mapped onto Bakta's CDS later).
            # The companion _gbff_proteins.faa is kept on disk for
            # provenance/debugging; downstream uses Bakta's proteins.
            gbff_gene_info_path = ""
            if reannotate_gbff:
                gbff_gene_info_path = self._wf(f"{self.config.sample_id}_gbff_gene_info.tsv")
                gbff_proteins_path = self._wf(f"{self.config.sample_id}_gbff_proteins.faa")
                rc, stderr = self._run_extract_proteins_script(gbff_proteins_path, gbff_gene_info_path, metadata_path)
                if rc != 0:
                    return StepResult(
                        "extract_proteins",
                        False,
                        f"GenBank parse failed: {stderr[:300]}",
                    )

                # GenBank → contigs FASTA for Bakta
                bakta_input = self._wf(f"{self.config.sample_id}_contigs.fna")
                try:
                    self._genbank_to_contigs_fasta(self.config.input_path, bakta_input)
                except Exception as e:
                    return StepResult(
                        "extract_proteins",
                        False,
                        f"GenBank → FASTA conversion failed: {e}",
                    )
            else:
                bakta_input = self.config.input_path

            # Run Bakta. For GenBank, write to a side file so the overlap
            # mapper can layer gbff_annotation in afterwards. For FASTA
            # contigs (no GenBank source), write straight to gene_info_path.
            bakta_gene_info_path = (
                self._wf(f"{self.config.sample_id}_bakta_gene_info.tsv") if reannotate_gbff else gene_info_path
            )
            bakta_args = [
                "--input",
                bakta_input,
                "--db",
                self.config.bakta_db,
                "--sample",
                self.config.sample_id,
                "--threads",
                str(self.config.bakta_threads),
                "--out-proteins",
                proteins_path,
                "--out-gene-info",
                bakta_gene_info_path,
            ]
            # Stage Bakta DB to local SSD if it lives on a network filesystem.
            # Bakta does many small reads against ~30 GB of indexes; on gpfs
            # with --threads 4 the per-read latency stacks linearly (4 min
            # → 27 min regression observed). Same pattern as EggNOG.
            cache = os.environ.get("TMPDIR", "")
            if cache:
                bakta_args.extend(["--local-cache-dir", cache])
            rc, stdout, stderr = run_script("run_bakta.py", bakta_args, timeout=14400)
            if rc != 0:
                return StepResult("extract_proteins", False, stderr[:500])

            # Layer the original GenBank product strings onto Bakta's CDS
            # via reciprocal coordinate overlap (Phase 3.3.b).
            if reannotate_gbff:
                rc, _, stderr = run_script(
                    "map_gbff_to_bakta_cds.py",
                    [
                        "--bakta-gene-info",
                        bakta_gene_info_path,
                        "--genbank-gene-info",
                        gbff_gene_info_path,
                        "--out",
                        gene_info_path,
                    ],
                )
                if rc != 0:
                    return StepResult(
                        "extract_proteins",
                        False,
                        f"GenBank → Bakta annotation mapping failed: {stderr[:300]}",
                    )

            tool_name = "Bakta (re-annotated)" if reannotate_gbff else "Bakta"
        else:
            # GenBank with --use-input-annotations, GFF3, FASTA contigs
            # without Bakta, or protein FASTA — extract_proteins.py handles all.
            rc, stderr = self._run_extract_proteins_script(proteins_path, gene_info_path, metadata_path)
            tool_name = {
                "fasta_contigs": "Prodigal",
                "protein_fasta": "Protein FASTA",
                "genbank": "GenBank parser (input annotations preserved)",
            }.get(fmt, "GenBank parser")

        if rc == 0:
            self.files["proteins"] = proteins_path
            self.files["gene_info"] = gene_info_path

            # Read organism name from metadata if available
            metadata_path = self._wf(f"{self.config.sample_id}_metadata.json")
            if os.path.exists(metadata_path):
                try:
                    import json

                    with open(metadata_path) as mf:
                        meta = json.load(mf)
                    self.files["organism"] = meta.get("organism", "")
                except Exception:
                    pass

            # Also extract gene order
            gene_order_path = self._wf(f"{self.config.sample_id}_gene_order.tsv")
            rc2, _, stderr2 = run_script(
                "extract_gene_order.py",
                [
                    "--gene-info",
                    gene_info_path,
                    "--output",
                    gene_order_path,
                ],
            )
            if rc2 == 0:
                self.files["gene_order"] = gene_order_path

            n = sum(1 for line in open(proteins_path) if line.startswith(">"))
            return StepResult("extract_proteins", True, f"{tool_name}: extracted {n} proteins")
        return StepResult("extract_proteins", False, stderr[:500])

    # ── Phase 2: SS Detection ──

    def _step_macsyfinder(self) -> StepResult:
        """Run MacSyFinder v2 to detect secretion systems."""
        msf_out = self._wf("macsyfinder_out")
        # Do NOT pre-create: MacSyFinder v2 requires --out-dir to not exist
        if os.path.exists(msf_out):
            shutil.rmtree(msf_out)

        proteins = self.files.get("proteins", "")
        if not proteins or not os.path.exists(proteins):
            return StepResult("macsyfinder", False, "No proteins file from previous step")

        # FRAGILE: macsydata install — downloads TXSScan HMM models from the
        # macsy-models GitHub repository. Can fail behind corporate firewalls
        # or if GitHub is unreachable. Models persist in ~/.macsyfinder/ after
        # first successful install, so this usually only runs once.
        # Pin to TXSScan==1.1.4 — the version ssign was validated against.
        # If this breaks: run manually: macsydata install --user TXSScan==1.1.4
        txsscan_meta = Path.home() / ".macsyfinder" / "models" / "TXSScan" / "metadata.yml"
        if txsscan_meta.exists() and "vers: 1.1.4" in txsscan_meta.read_text():
            logger.info("TXSScan 1.1.4 already installed — skipping macsydata install")
        else:
            install_cmd = ["macsydata", "install", "--user", "TXSScan==1.1.4"]
            try:
                install_result = subprocess.run(install_cmd, capture_output=True, text=True, timeout=120)
                if install_result.returncode != 0:
                    logger.warning(
                        f"macsydata install returned exit code "
                        f"{install_result.returncode}. TXSScan models may "
                        f"already be installed — continuing.\n"
                        f"  If MacSyFinder fails, run manually:\n"
                        f"    macsydata install --user TXSScan==1.1.4"
                    )
            except FileNotFoundError:
                return StepResult(
                    "macsyfinder",
                    False,
                    "macsydata command not found.\n"
                    "  Install: pip install macsyfinder\n"
                    "  It should have been installed as a dependency of ssign.",
                )
            except subprocess.TimeoutExpired:
                logger.warning("macsydata install timed out — TXSScan models may already exist. Continuing.")
            except Exception:
                pass  # May already be installed

        # FRAGILE: macsyfinder CLI — internally calls hmmsearch via subprocess.
        # Our pyhmmer shim (installed as console_script 'hmmsearch') should be
        # found on PATH. If not, MacSyFinder will fail with a message about
        # hmmsearch not being found.
        # If this breaks: sudo apt install hmmer (or: conda install -c bioconda hmmer)
        cmd = [
            "macsyfinder",
            "--sequence-db",
            proteins,
            "--db-type",
            self.config.macsyfinder_db_type,
            "--models",
            "TXSScan",
            "all",
            "--out-dir",
            msf_out,
            "-w",
            str(self.config.cpu_per_genome),
            "--mute",
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600,
            )
            if result.returncode == 0:
                self.files["macsyfinder_out"] = msf_out
                return StepResult("macsyfinder", True, "MacSyFinder v2 complete")

            error_msg = result.stderr[:500]
            # Detect common MacSyFinder failures and provide specific advice
            if "hmmsearch" in error_msg.lower() or "hmmer" in error_msg.lower():
                return StepResult(
                    "macsyfinder",
                    False,
                    "MacSyFinder failed — hmmsearch not available.\n"
                    "  The pyhmmer shim should provide this automatically.\n"
                    "  How to fix:\n"
                    "    - pip install --force-reinstall ssign\n"
                    "    - Or install real HMMER: sudo apt install hmmer",
                )
            if "txscan" in error_msg.lower() or "model" in error_msg.lower():
                return StepResult(
                    "macsyfinder",
                    False,
                    f"MacSyFinder failed — TXSScan models may not be "
                    f"installed.\n"
                    f"  Run: macsydata install --user TXSScan\n"
                    f"  Original error: {error_msg}",
                )
            return StepResult("macsyfinder", False, error_msg)
        except FileNotFoundError:
            return StepResult(
                "macsyfinder",
                False,
                "macsyfinder command not found.\n"
                "  Install: pip install macsyfinder\n"
                "  It should have been installed as a dependency of ssign.",
            )
        except subprocess.TimeoutExpired:
            return StepResult("macsyfinder", False, "Timeout after 60 min")
        except Exception as e:
            return StepResult("macsyfinder", False, str(e))

    def _step_validate_systems(self) -> StepResult:
        """Validate MacSyFinder results and extract SS components."""
        msf_out = self.files.get("macsyfinder_out", "")
        gene_info = self.files.get("gene_info", "")
        components_path = self._wf(f"{self.config.sample_id}_ss_components.tsv")
        systems_path = self._wf(f"{self.config.sample_id}_valid_systems.tsv")

        if not msf_out or not os.path.isdir(msf_out):
            return StepResult("validate_systems", False, "No MacSyFinder output from previous step")

        rc, stdout, stderr = run_script(
            "validate_macsyfinder_systems.py",
            [
                "--msf-dir",
                msf_out,
                "--gene-info",
                gene_info,
                "--sample",
                self.config.sample_id,
                "--wholeness-threshold",
                str(self.config.wholeness_threshold),
                "--excluded-systems",
                ",".join(self.config.excluded_systems),
                "--out-components",
                components_path,
                "--out-systems",
                systems_path,
            ],
        )

        if rc == 0:
            self.files["ss_components"] = components_path
            self.files["valid_systems"] = systems_path
            return StepResult("validate_systems", True, "Systems validated")
        return StepResult("validate_systems", False, stderr[:500])

    def _step_extract_neighborhood(self) -> StepResult:
        """Extract proteins near SS components for focused prediction."""
        gene_order = self.files.get("gene_order", "")
        ss_components = self.files.get("ss_components", "")
        proteins = self.files.get("proteins", "")

        if not ss_components or not os.path.exists(ss_components):
            return StepResult("extract_neighborhood", False, "No SS components from previous step")

        neighborhood_fasta = self._wf(f"{self.config.sample_id}_neighborhood.faa")
        neighborhood_ids = self._wf(f"{self.config.sample_id}_neighborhood_ids.tsv")

        rc, stdout, stderr = run_script(
            "extract_neighborhood.py",
            [
                "--gene-order",
                gene_order,
                "--ss-components",
                ss_components,
                "--proteins",
                proteins,
                "--window",
                str(self.config.proximity_window),
                "--output",
                neighborhood_fasta,
                "--output-ids",
                neighborhood_ids,
            ],
        )

        if rc == 0:
            self.files["neighborhood_proteins"] = neighborhood_fasta
            self.files["neighborhood_ids"] = neighborhood_ids
            n_neigh = sum(1 for line in open(neighborhood_fasta) if line.startswith(">"))
            return StepResult("extract_neighborhood", True, f"{n_neigh} neighborhood proteins")
        return StepResult("extract_neighborhood", False, stderr[:500])

    # ── Phase 3: Prediction ──

    def _resolve_step_input_fasta(self, whole_genome: bool) -> str:
        """Pick the input FASTA path for a prediction step.

        - whole_genome=True forces the full proteome.

        Falls back to full proteome only when no neighborhood file is
        staged (e.g. MacSyFinder found no systems). Raises RuntimeError
        if neither neighborhood nor proteins are staged — that means
        input-processing never ran and the caller has no business
        launching a subprocess with --input "".
        """
        proteins = self.files.get("proteins", "")
        if whole_genome:
            if not proteins:
                raise RuntimeError("No proteins FASTA staged; input-processing step must run first")
            return proteins
        chosen = self.files.get("neighborhood_proteins") or proteins
        if not chosen:
            raise RuntimeError("No neighborhood or proteins FASTA staged; input-processing step must run first")
        return chosen

    def _scaled_tool_timeout(
        self, tool: str, size: int, *, whole_genome: bool = False, floor: int = TOOL_TIMEOUT_S
    ) -> int:
        """Size-aware subprocess timeout (seconds) for `tool`, from the effort model.

        `size` is the sequence count of what the tool actually processes (the
        proteome for predictors, the substrate count for annotation). Scales the
        cap so large inputs (e.g. a pooled 160k-protein whole-genome prediction)
        aren't killed by the flat 4h floor; unmodelled tools / unknown size fall
        back to `floor`. openspec: size-aware-tool-timeouts.
        """
        regime = resolve_regime(tool, whole_genome=whole_genome)
        if regime is None or size is None or size < 0:
            return floor
        return scaled_timeout(tool, size, regime, floor=floor)

    def _substrate_count(self) -> int:
        """Number of filtered substrates staged for this runner (0 if none/unreadable).

        This is the effort-model `size` axis for the substrate-scoped annotation
        tools; in a batched run the pool runner stages the pooled substrate set
        here, so the cap scales with the whole panel's substrate load.
        """
        path = self.files.get("substrates_filtered", "")
        if not path:
            return 0
        try:
            return len(load_substrate_ids(path))
        except (OSError, ValueError):
            return 0

    def _step_deeplocpro(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_deeplocpro.tsv")

        input_proteins = self._resolve_step_input_fasta(self.config.dlp_whole_genome)

        args = [
            "--input",
            input_proteins,
            "--sample",
            self.config.sample_id,
            "--conf-threshold",
            str(self.config.conf_threshold),
            "--output",
            output,
        ]

        if self.config.deeplocpro_mode == "local":
            args.extend(["--mode", "local"])
            if self.config.deeplocpro_path:
                args.extend(["--deeplocpro-path", self.config.deeplocpro_path])
        else:
            args.extend(["--mode", "remote"])

        sem = self.api_sem.get("dtu")
        held = False
        if sem:
            held = sem.acquire(timeout=DTU_SEMAPHORE_TIMEOUT_S)
            if not held:
                logger.warning(
                    "DTU semaphore acquire timed out after %ss — proceeding without "
                    "rate-limit hold; expect possible API throttling.",
                    DTU_SEMAPHORE_TIMEOUT_S,
                )
        timeout_s = self._scaled_tool_timeout(
            "deeplocpro", count_sequences(input_proteins), whole_genome=self.config.dlp_whole_genome
        )
        args.extend(["--timeout", str(timeout_s)])
        try:
            rc, stdout, stderr = run_script("run_deeplocpro.py", args, timeout=timeout_s)
        finally:
            if held:
                sem.release()
        if rc == 0:
            self.files["deeplocpro"] = output
            return StepResult("deeplocpro", True, "DeepLocPro complete")
        return StepResult("deeplocpro", False, stderr[:500])

    def _step_deepsece(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_deepsece.tsv")

        input_proteins = self._resolve_step_input_fasta(self.config.dse_whole_genome)

        rc, stdout, stderr = run_script(
            "run_deepsece.py",
            [
                "--input",
                input_proteins,
                "--sample",
                self.config.sample_id,
                "--output",
                output,
            ],
            timeout=self._scaled_tool_timeout(
                "deepsece", count_sequences(input_proteins), whole_genome=self.config.dse_whole_genome
            ),
        )

        if rc == 0:
            self.files["deepsece"] = output
            return StepResult("deepsece", True, "DeepSecE complete")
        return StepResult("deepsece", False, stderr[:500])

    def _step_signalp(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_signalp.tsv")

        input_proteins = self._resolve_step_input_fasta(self.config.sp_whole_genome)

        args = [
            "--input",
            input_proteins,
            "--sample",
            self.config.sample_id,
            "--output",
            output,
        ]
        if self.config.signalp_mode == "local":
            args.extend(["--mode", "local"])
            if self.config.signalp_path:
                args.extend(["--signalp-path", self.config.signalp_path])
        else:
            args.extend(["--mode", "remote"])

        sem = self.api_sem.get("dtu")
        held = False
        if sem:
            held = sem.acquire(timeout=DTU_SEMAPHORE_TIMEOUT_S)
            if not held:
                logger.warning(
                    "DTU semaphore acquire timed out after %ss — proceeding without "
                    "rate-limit hold; expect possible API throttling.",
                    DTU_SEMAPHORE_TIMEOUT_S,
                )
        timeout_s = self._scaled_tool_timeout(
            "signalp", count_sequences(input_proteins), whole_genome=self.config.sp_whole_genome
        )
        args.extend(["--timeout", str(timeout_s)])
        try:
            rc, stdout, stderr = run_script("run_signalp.py", args, timeout=timeout_s)
        finally:
            if held:
                sem.release()
        if rc == 0:
            self.files["signalp"] = output
            return StepResult("signalp", True, "SignalP complete")
        return StepResult("signalp", False, stderr[:500])

    def _step_cross_validate(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_predictions.tsv")

        dlp = self.files.get("deeplocpro", "")
        if not dlp or not os.path.exists(dlp):
            return StepResult("cross_validate", False, "No DeepLocPro output from previous step")

        dse = self.files.get("deepsece", "")
        plm_e = self.files.get("plm_effector", "")
        sp = self.files.get("signalp", "")

        # cross_validate itself no longer needs valid_systems (DSE-T3SS is flagged
        # unconditionally), but its presence still gates that the SS-detection step ran.
        valid_sys = self.files.get("valid_systems", "")
        if not valid_sys or not os.path.exists(valid_sys):
            return StepResult("cross_validate", False, "No valid_systems file from previous step")

        args = [
            "--deeplocpro",
            dlp,
            "--sample",
            self.config.sample_id,
            "--conf-threshold",
            str(self.config.conf_threshold),
            "--output",
            output,
        ]
        if dse and os.path.exists(dse):
            args.extend(["--deepsece", dse])
        if plm_e and os.path.exists(plm_e):
            args.extend(["--plm-effector", plm_e])
        if sp and os.path.exists(sp):
            args.extend(["--signalp", sp])
        ss_components = self.files.get("ss_components", "")
        if ss_components and os.path.exists(ss_components):
            args.extend(["--ss-components", ss_components])

        rc, stdout, stderr = run_script("cross_validate_predictions.py", args)
        if rc == 0:
            self.files["predictions"] = output
            return StepResult("cross_validate", True, "Predictions validated")
        return StepResult("cross_validate", False, stderr[:500])

    # ── Phase 4: Substrate ID ──

    def _step_proximity(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_substrates.tsv")

        rc, stdout, stderr = run_script(
            "proximity_analysis.py",
            [
                "--gene-order",
                self.files.get("gene_order", ""),
                "--ss-components",
                self.files.get("ss_components", ""),
                "--predictions",
                self.files.get("predictions", ""),
                "--sample",
                self.config.sample_id,
                "--window",
                str(self.config.proximity_window),
                "--conf-threshold",
                str(self.config.conf_threshold),
                "--output",
                output,
            ],
        )

        if rc == 0:
            self.files["substrates"] = output
            n = sum(1 for _ in open(output)) - 1  # minus header
            return StepResult("proximity", True, f"Found {n} substrate candidates")
        return StepResult("proximity", False, stderr[:500])

    def _step_t5ss(self) -> StepResult:
        out_sub = self._wf(f"{self.config.sample_id}_t5ss_substrates.tsv")
        out_dom = self._wf(f"{self.config.sample_id}_t5ss_domains.tsv")

        rc, stdout, stderr = run_script(
            "t5ss_handler.py",
            [
                "--ss-components",
                self.files.get("ss_components", ""),
                "--predictions",
                self.files.get("predictions", ""),
                "--proteins",
                self.files.get("proteins", self._wf(f"{self.config.sample_id}_proteins.faa")),
                "--sample",
                self.config.sample_id,
                "--out-substrates",
                out_sub,
                "--out-domains",
                out_dom,
            ],
        )

        if rc == 0:
            self.files["t5ss_substrates"] = out_sub
            self.files["t5ss_domains"] = out_dom
            return StepResult("t5ss", True, "T5SS handled")
        return StepResult("t5ss", False, stderr[:500])

    def _step_filtering(self) -> StepResult:
        out_filtered = self._wf(f"{self.config.sample_id}_substrates_filtered.tsv")
        out_all = self._wf(f"{self.config.sample_id}_substrates_all.tsv")

        substrates = self.files.get("substrates", "")
        t5ss = self.files.get("t5ss_substrates", "")
        if not substrates or not os.path.exists(substrates):
            return StepResult("filtering", False, "No substrates from previous step")

        gate_verdicts = self._wf(f"{self.config.sample_id}_gate_verdicts.tsv")
        filter_args = [
            "--proximity-substrates",
            substrates,
            "--t5ss-substrates",
            t5ss if t5ss and os.path.exists(t5ss) else substrates,
            "--valid-systems",
            self.files.get("valid_systems", ""),
            "--ss-components",
            self.files.get("ss_components", ""),
            "--predictions",
            self.files.get("predictions", ""),
            "--sample",
            self.config.sample_id,
            "--excluded-systems",
            ",".join(self.config.excluded_systems),
            "--required-fraction-correct",
            str(self.config.required_fraction_correct),
            "--dlp-confidence-threshold",
            str(self.config.dlp_confidence_threshold),
            "--out-filtered",
            out_filtered,
            "--out-all",
            out_all,
            "--out-gate-verdicts",
            gate_verdicts,
        ]
        if self.config.filter_dse_type_mismatch:
            filter_args.append("--filter-dse-type-mismatch")
        if self.config.skip_localization_gate:
            filter_args.append("--skip-localization-gate")

        rc, stdout, stderr = run_script("system_filtering.py", filter_args)

        if rc == 0:
            self.files["substrates_filtered"] = out_filtered
            self.files["substrates_all"] = out_all
            if os.path.exists(gate_verdicts):
                self.files["gate_verdicts"] = gate_verdicts
            n_subs = sum(1 for line in open(out_filtered) if not line.startswith("locus_tag") and line.strip()) - 1
            n_subs = max(0, n_subs)
            return StepResult("filtering", True, f"{n_subs} secreted proteins")
        return StepResult("filtering", False, stderr[:500])

    def _step_build_passenger_fasta(self) -> StepResult:
        """Build a passenger-substituted FASTA for routed annotation tools.

        For each clean T5aSS substrate (non-empty t5_quality_flag absent,
        passenger_length >= MIN_PASSENGER_FOR_ANNOTATION), the entry's
        sequence is replaced with the passenger subsequence
        ``seq[sp_end : barrel_start - LINKER_LENGTH]``. Every other
        entry (non-T5aSS substrates, T5aSS with quality flags, T5aSS
        with passenger too short, non-substrate proteins) carries its
        full sequence unchanged. The 5 annotation tools that benefit
        from passenger-only annotation (EggNOG / BLASTp / pLM-BLAST /
        HHsuite / ProtParam) read this FASTA instead of the full
        proteins FASTA; IPS continues to consume the full FASTA
        because it is already domain-aware.

        Always runs (even when there are no T5aSS substrates) so the
        ``proteins_for_passenger_tools`` key is always staged. When the
        t5ss step found no T5aSS components, the output FASTA is
        functionally a copy of the full proteins FASTA.
        """
        from ssign_app.scripts.ssign_lib.t5_passenger import build_passenger_substituted_fasta

        proteins = self.files.get("proteins", "")
        if not proteins:
            return StepResult(
                "build_passenger_fasta",
                False,
                "No proteins FASTA staged; extract_proteins step must run first",
            )
        classifications = self.files.get("t5ss_domains", "")

        out_fasta = self._wf(f"{self.config.sample_id}_proteins_passenger.faa")
        out_source = self._wf(f"{self.config.sample_id}_t5_annotation_source.tsv")
        routing = build_passenger_substituted_fasta(proteins, classifications, out_fasta, out_source)

        self.files["proteins_for_passenger_tools"] = out_fasta
        self.files["t5_annotation_source"] = out_source

        n_pass = sum(1 for v in routing.values() if v == "passenger")
        n_full = sum(1 for v in routing.values() if v == "full")
        return StepResult(
            "build_passenger_fasta",
            True,
            f"{n_pass} T5aSS substrates routed to passenger, {n_full} fell back to full protein",
        )

    # ── Phase 5: Annotation ──

    def _annotation_cpu_budget(self) -> int:
        """Per-tool CPU budget when annotation_steps run in parallel.

        BLASTp, InterProScan, HH-suite, and EggNOG all saturate CPU during
        their DIAMOND/hmmsearch/Java phases. When they fire concurrently in
        the parallel annotation_steps group (runner.py L591), naively giving
        each one `cpu_per_genome` produces N× oversubscription. Divide the
        budget by the count of CPU-heavy annotators currently enabled.
        pLM-BLAST and ProtParam aren't counted: the former is GPU-bound for
        its long phase, the latter is microseconds of pure-Python work.
        """
        n_heavy = sum(
            [
                not self.config.skip_hhsuite,
                not self.config.skip_blastp,
                not self.config.skip_interproscan,
                not self.config.skip_eggnog,
            ]
        )
        return max(1, self.config.cpu_per_genome // max(1, n_heavy))

    def _check_substrates_exist(self, step_name):
        """Check that upstream substrate files exist. Returns error StepResult or None."""
        sf = self.files.get("substrates_filtered", "")
        if not sf or not os.path.exists(sf):
            return StepResult(step_name, False, "Skipped — no substrates from upstream steps")
        return None

    # Annotation tools that benefit from passenger-only annotation on
    # T5aSS substrates. IPS is excluded — it's already domain-aware and
    # would lose its barrel-Pfam confirmation if pointed at passenger seqs.
    _PASSENGER_ROUTED_TOOLS = frozenset({"blastp", "hhsuite", "eggnog", "plm_blast", "protparam"})

    def _annotation_input_proteins(self, tool: str) -> str:
        """FASTA path that ``--proteins`` should point to for an annotation tool.

        Routed tools get the passenger-substituted FASTA produced by
        ``_step_build_passenger_fasta`` when it is staged; IPS and any
        future non-routed tool get the full proteins FASTA. Falls back
        to the full FASTA whenever the passenger one is missing, so
        callers degrade gracefully if the build step was skipped.
        """
        if tool in self._PASSENGER_ROUTED_TOOLS:
            return self.files.get("proteins_for_passenger_tools", "") or self.files.get("proteins", "")
        return self.files.get("proteins", "")

    def _step_blastp(self) -> StepResult:
        err = self._check_substrates_exist("blastp")
        if err:
            return err

        output = self._wf(f"{self.config.sample_id}_blastp.csv")

        if not self.config.blastp_db:
            return StepResult(
                "blastp",
                False,
                "BLASTp requires a local database. Set `blastp_db` in the "
                "config to the path of your BLAST+ database (e.g. NR).",
            )

        args = [
            "--substrates",
            self.files.get("substrates_filtered", ""),
            "--proteins",
            self._annotation_input_proteins("blastp"),
            "--sample",
            self.config.sample_id,
            "--output",
            output,
            "--db",
            self.config.blastp_db,
            "--min-pident",
            str(self.config.blastp_min_pident),
            "--min-qcov",
            str(self.config.blastp_min_qcov),
            "--evalue",
            str(self.config.blastp_evalue),
        ]
        if self.config.blastp_exclude_taxid:
            args.extend(["--exclude-taxid", self.config.blastp_exclude_taxid])
        args.extend(["--threads", str(self._annotation_cpu_budget())])

        # Size-aware timeout on the substrate count (like the other annotation
        # tools). Swiss-Prot (the full-tier default) finishes in minutes and stays
        # well under any cap; the prior is sized for the slow opt-in --blastp-db NR
        # run (the old hardcoded 7200s killed NR mid-search — full-tier smoke test).
        # Small runs keep the 4h floor.
        timeout_s = self._scaled_tool_timeout("blastp", self._substrate_count())
        rc, stdout, stderr = run_script("run_blastp.py", args, timeout=timeout_s)
        if rc == 0:
            self.files["blastp"] = output
            return StepResult("blastp", True, "BLASTp complete")
        return StepResult("blastp", False, stderr[:500])

    def _step_hhsuite(self) -> StepResult:
        err = self._check_substrates_exist("hhsuite")
        if err:
            return err

        if not self.config.hhsuite_uniclust_db:
            return StepResult(
                "hhsuite",
                False,
                "HH-suite requires a local UniRef30/UniClust30 database for hhblits. "
                "Set `hhsuite_uniclust_db` in the config.",
            )
        if not self.config.hhsuite_pfam_db and not self.config.hhsuite_pdb70_db:
            return StepResult(
                "hhsuite",
                False,
                "HH-suite requires at least one of `hhsuite_pfam_db` or `hhsuite_pdb70_db` to be set.",
            )

        output = self._wf(f"{self.config.sample_id}_hhsuite.csv")
        args = [
            "--substrates",
            self.files.get("substrates_filtered", ""),
            "--proteins",
            self._annotation_input_proteins("hhsuite"),
            "--sample",
            self.config.sample_id,
            "--output",
            output,
            "--uniclust-db",
            self.config.hhsuite_uniclust_db,
        ]
        if self.config.hhsuite_pfam_db:
            args.extend(["--pfam-db", self.config.hhsuite_pfam_db])
        if self.config.hhsuite_pdb70_db:
            args.extend(["--pdb70-db", self.config.hhsuite_pdb70_db])
        args.extend(["--min-prob", str(self.config.hhsuite_min_prob)])

        # Two-layer parallelism: workers × cpu_per_job. hhblits/hhsearch scale
        # sub-linearly past 2-4 threads per process, so spend the budget on
        # workers and pin cpu_per_job=2.
        budget = self._annotation_cpu_budget()
        args.extend(
            [
                "--max-workers",
                str(max(1, budget // 2)),
                "--cpu-per-job",
                "2",
            ]
        )
        # Stage uniclust/pfam/pdb70 to local SSD on network filesystems.
        # hhblits does many small random reads per query — fast on local
        # SSD, pathologically slow on gpfs/nfs without thread fan-out.
        cache = os.environ.get("TMPDIR", "")
        if cache:
            args.extend(["--local-cache-dir", cache])

        # Size-aware timeout (like the other annotation tools): a big pooled
        # substrate set legitimately needs many hours of per-protein hhblits/
        # hhsearch, so the flat 4h floor killed 74/57-genome pooled runs mid-step.
        # HH-suite has no overall --timeout flag (per-protein caps live inside the
        # wrapper), so only the outer subprocess cap is scaled.
        timeout_s = self._scaled_tool_timeout("hh_suite", self._substrate_count())
        rc, stdout, stderr = run_script("run_hhsuite.py", args, timeout=timeout_s)
        if rc == 0:
            self.files["hhsuite"] = output
            return StepResult("hhsuite", True, "HH-suite complete")
        return StepResult("hhsuite", False, stderr[:500])

    def _step_interproscan(self) -> StepResult:
        err = self._check_substrates_exist("interproscan")
        if err:
            return err
        output = self._wf(f"{self.config.sample_id}_interproscan.csv")

        args = [
            "--substrates",
            self.files.get("substrates_filtered", ""),
            "--proteins",
            self.files.get("proteins", ""),
            "--sample",
            self.config.sample_id,
            "--output",
            output,
        ]
        if self.config.interproscan_db:
            args.extend(["--db", self.config.interproscan_db])
        args.extend(["--cpu", str(self._annotation_cpu_budget())])
        # Stage IPS install tree to local SSD on network filesystems.
        cache = os.environ.get("TMPDIR", "")
        if cache:
            args.extend(["--local-cache-dir", cache])

        timeout_s = self._scaled_tool_timeout("interproscan", self._substrate_count())
        args.extend(["--timeout", str(timeout_s)])
        rc, stdout, stderr = run_script("run_interproscan.py", args, timeout=timeout_s)
        if rc == 0:
            self.files["interproscan"] = output
            return StepResult("interproscan", True, "InterProScan complete")
        return StepResult("interproscan", False, stderr[:500])

    def _step_protparam(self) -> StepResult:
        err = self._check_substrates_exist("protparam")
        if err:
            return err
        output = self._wf(f"{self.config.sample_id}_protparam.csv")

        rc, stdout, stderr = run_script(
            "compute_protparam.py",
            [
                "--substrates",
                self.files.get("substrates_filtered", ""),
                "--proteins",
                self._annotation_input_proteins("protparam"),
                "--sample",
                self.config.sample_id,
                "--output",
                output,
            ],
            timeout=self._scaled_tool_timeout("protparam", self._substrate_count()),
        )

        if rc == 0:
            self.files["protparam"] = output
            return StepResult("protparam", True, "ProtParam complete")
        return StepResult("protparam", False, stderr[:500])

    def _step_eggnog(self) -> StepResult:
        err = self._check_substrates_exist("eggnog")
        if err:
            return err
        if not self.config.eggnog_db:
            return StepResult(
                "eggnog",
                False,
                "EggNOG-mapper requires a database. Set `eggnog_db` in the "
                "config (extended/full install tier fetches it automatically).",
            )

        output = self._wf(f"{self.config.sample_id}_eggnog.tsv")
        args = [
            "--substrates",
            self.files.get("substrates_filtered", ""),
            "--proteins",
            self._annotation_input_proteins("eggnog"),
            "--db",
            self.config.eggnog_db,
            "--sample",
            self.config.sample_id,
            "--out",
            output,
            "--threads",
            str(self._annotation_cpu_budget()),
            "--sensmode",
            self.config.eggnog_sensmode,
        ]
        # None => let the wrapper auto-decide based on host RAM.
        if self.config.eggnog_dbmem is True:
            args.append("--dbmem")
        elif self.config.eggnog_dbmem is False:
            args.append("--no-dbmem")
        # Resolve local-cache-dir: "" => use $TMPDIR if set, "off" => skip.
        cache_dir = self.config.eggnog_local_cache_dir
        if cache_dir == "":
            cache_dir = os.environ.get("TMPDIR", "")
        if cache_dir and cache_dir != "off":
            args.extend(["--local-cache-dir", cache_dir])
        timeout_s = self._scaled_tool_timeout("eggnog", self._substrate_count())
        args.extend(["--timeout", str(timeout_s)])
        rc, stdout, stderr = run_script("run_eggnog.py", args, timeout=timeout_s)
        if rc == 0:
            self.files["eggnog"] = output
            return StepResult("eggnog", True, "EggNOG-mapper complete")
        return StepResult("eggnog", False, stderr[:500])

    def _step_plm_blast(self) -> StepResult:
        err = self._check_substrates_exist("plm_blast")
        if err:
            return err
        if not self.config.plmblast_db:
            return StepResult(
                "plm_blast",
                False,
                "pLM-BLAST requires an ECOD database. Set `plmblast_db` "
                "in the config; fetch from ftp.tuebingen.mpg.de.",
            )

        output = self._wf(f"{self.config.sample_id}_plm_blast.tsv")
        # Pass the per-genome CPU budget through to plmblast.py's
        # ProcessPoolExecutor (-workers). Without this the wrapper falls
        # back to its `default=4` regardless of node size, so the search
        # step only ever used 4 cores even on a 24-core HPC node.
        # Verified 2026-05-21: 4-worker pLM-BLAST hit the 14400s timeout
        # on a single E. coli K-12 substrate set; full cpu_per_genome
        # should bring that to ~40-60 min.
        args = [
            "--substrates",
            self.files.get("substrates_filtered", ""),
            "--proteins",
            self._annotation_input_proteins("plm_blast"),
            "--ecod-db",
            self.config.plmblast_db,
            "--out",
            output,
            "--threads",
            str(self.config.cpu_per_genome),
            "--cpc",
            str(self.config.plmblast_cpc),
        ]
        # Stage ECOD DB tree to local SSD on network filesystems (same
        # pattern as IPS / HH-suite / Bakta — see _step_interproscan).
        cache = os.environ.get("TMPDIR", "")
        if cache:
            args.extend(["--local-cache-dir", cache])

        timeout_s = self._scaled_tool_timeout("plm_blast", self._substrate_count(), floor=PLMBLAST_TIMEOUT_S)
        args.extend(["--timeout", str(timeout_s)])
        rc, stdout, stderr = run_script("run_plm_blast.py", args, timeout=timeout_s)
        if rc == 0:
            self.files["plm_blast"] = output
            return StepResult("plm_blast", True, "pLM-BLAST complete")
        return StepResult("plm_blast", False, stderr[:500])

    def _step_plm_effector(self) -> StepResult:
        """Run PLM-Effector across all requested SS types in one subprocess,
        then merge the per-type TSVs into a single combined TSV.

        The merged TSV has one row per protein with `passes_threshold=1`
        iff the ensemble flagged it for at least one SS type — shape that
        cross_validate_predictions expects.
        """
        if not self.config.plm_effector_weights_dir:
            return StepResult(
                "plm_effector",
                False,
                "PLM-Effector requires a weights directory. Set "
                "`plm_effector_weights_dir` in the config; extended/full "
                "install tier fetches the weights automatically.",
            )

        # Stage the ~19 GB PLM-Effector weights tree to local SSD once.
        # Each PLM subprocess re-reads ~4-12 GB of model weights; on gpfs
        # that's pathologically slow. Local staging makes reads 5-10×
        # faster. No-op when weights are already on local FS.
        weights_dir = self.config.plm_effector_weights_dir
        cache = os.environ.get("TMPDIR", "")
        if cache:
            from ssign_app.scripts.ssign_lib.resources import stage_directory_tree_to_local_ssd_if_remote

            weights_dir = stage_directory_tree_to_local_ssd_if_remote(weights_dir, cache)

        input_proteins = self._resolve_step_input_fasta(self.config.plme_whole_genome)

        out_dir = self._wf(f"{self.config.sample_id}_plm_effector")
        os.makedirs(out_dir, exist_ok=True)
        args = [
            "--input",
            input_proteins,
            "--weights-dir",
            weights_dir,
            "--effector-types",
            *self.config.plm_effector_types,
            "--out-dir",
            out_dir,
            "--chunk-size",
            str(self.config.plm_chunk_size),
        ]
        rc, stdout, stderr = run_script(
            "run_plm_effector.py",
            args,
            timeout=self._scaled_tool_timeout(
                "plm_effector", count_sequences(input_proteins), whole_genome=self.config.plme_whole_genome
            ),
            stream_stderr=True,
        )
        if rc != 0:
            return StepResult(
                "plm_effector",
                False,
                f"PLM-Effector failed: {stderr[:400]}",
            )

        per_type_paths = [os.path.join(out_dir, f"{eff_type}.tsv") for eff_type in self.config.plm_effector_types]
        merged = self._wf(f"{self.config.sample_id}_plm_effector_merged.tsv")
        rc, stdout, stderr = run_script(
            "merge_plm_effector_outputs.py",
            ["--inputs"] + per_type_paths + ["--out", merged],
        )
        if rc == 0:
            self.files["plm_effector"] = merged
            return StepResult(
                "plm_effector",
                True,
                f"PLM-Effector complete ({len(per_type_paths)} types merged)",
            )
        return StepResult("plm_effector", False, stderr[:500])

    # ── Phase 6: Integration ──

    _T5ASS_WHOLE_TOOLS = ("eggnog", "blastp", "plm_blast", "hhsuite", "protparam")

    def _build_t5ass_only_substrates_tsv(self) -> str | None:
        """Filter substrates_filtered.tsv to clean T5aSS rows; write to a new TSV.

        Returns the path to the filtered TSV, or None if no clean T5aSS
        substrates exist. A T5aSS substrate is recognised by either
        ``tool == "T5SS-self"`` (from t5ss_handler) or "T5aSS" in
        ``nearby_ss_types``; rows with a non-empty ``t5_quality_flag``
        are excluded since we don't trust their domain geometry.
        """
        import csv as _csv

        substrates_filtered = self.files.get("substrates_filtered", "")
        if not substrates_filtered or not os.path.exists(substrates_filtered):
            return None

        with open(substrates_filtered) as f:
            reader = _csv.DictReader(f, delimiter="\t")
            fieldnames = list(reader.fieldnames or [])
            t5_rows = []
            for row in reader:
                ss_types = (row.get("nearby_ss_types") or "").strip()
                tool = (row.get("tool") or "").strip()
                quality = (row.get("t5_quality_flag") or "").strip()
                is_t5a = (tool == "T5SS-self") or ("T5aSS" in ss_types)
                if is_t5a and not quality:
                    t5_rows.append(row)

        if not t5_rows:
            return None

        out_path = self._wf(f"{self.config.sample_id}_t5ass_only_substrates.tsv")
        with open(out_path, "w", newline="") as f:
            writer = _csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(t5_rows)
        return out_path

    def _step_t5ass_whole_annotations(self) -> StepResult:
        """Second annotation pass on T5aSS substrates using the full protein FASTA.

        Re-runs each routed annotation tool (EggNOG / BLASTp / pLM-BLAST /
        HHsuite / ProtParam) that succeeded in the main pass — this time
        against the full protein sequence instead of the passenger
        subsequence. Output files have their column headers prefixed
        with ``t5ass_whole_`` so integrate_annotations merges them as
        side-by-side columns rather than overwriting the passenger
        annotations.

        Skipped when ``config.t5ass_annotate_whole`` is False (default).
        Implementation re-uses the existing step methods by swapping
        ``self.work_dir`` to a sibling directory and overriding
        ``substrates_filtered`` and ``proteins_for_passenger_tools`` in
        ``self.files`` for the duration of each call; original state is
        restored at the end so downstream integrate sees the original
        main-pass outputs untouched plus the new t5ass_whole_* keys.

        Sequential across tools (not the parallel-group treatment the
        main annotation pass gets). T5aSS substrate counts are typically
        small (<10 per genome), so the per-tool startup cost dominates
        wallclock and serial sequencing is a few minutes of overhead
        rather than tens — fine for v1.
        """
        if not self.config.t5ass_annotate_whole:
            return StepResult(
                "t5ass_whole_annotations",
                True,
                "Skipped (--t5ass-annotate-whole not set)",
            )

        t5_sub_path = self._build_t5ass_only_substrates_tsv()
        if t5_sub_path is None:
            return StepResult(
                "t5ass_whole_annotations",
                True,
                "Skipped (no clean T5aSS substrates)",
            )

        second_work_dir = tempfile.mkdtemp(prefix=f"ssign_t5ass_whole_{self.config.sample_id}_")
        saved_work_dir = self.work_dir
        saved_substrates = self.files.get("substrates_filtered")
        saved_passenger = self.files.get("proteins_for_passenger_tools")
        saved_outputs = {tool: self.files.get(tool) for tool in self._T5ASS_WHOLE_TOOLS}

        # Hijack inputs: swap work_dir so step methods don't clobber main-
        # pass outputs, point substrates_filtered at the T5aSS subset, and
        # drop proteins_for_passenger_tools so _annotation_input_proteins
        # falls back to the full proteins FASTA.
        self.work_dir = second_work_dir
        self.files["substrates_filtered"] = t5_sub_path
        self.files.pop("proteins_for_passenger_tools", None)

        n_run = 0
        failures: list[str] = []
        try:
            for tool in self._T5ASS_WHOLE_TOOLS:
                if saved_outputs[tool] is None:
                    continue  # main pass didn't run this tool; nothing to compare to
                try:
                    step_method = getattr(self, f"_step_{tool}")
                    result = step_method()
                    if not result.success:
                        failures.append(f"{tool}: {result.message[:80]}")
                        continue
                    tool_output = self.files.get(tool, "")
                    if not tool_output or not os.path.exists(tool_output):
                        failures.append(f"{tool}: no output file")
                        continue
                    ext = os.path.splitext(tool_output)[1]
                    final_path = os.path.join(
                        saved_work_dir,
                        f"{self.config.sample_id}_t5ass_whole_{tool}{ext}",
                    )
                    _rename_csv_columns_with_prefix(
                        tool_output, final_path, prefix="t5ass_whole_", keep_columns={"locus_tag"}
                    )
                    # Stash the renamed-column path under the t5ass_whole key
                    # AFTER the finally-block restoration would overwrite it.
                    saved_outputs[tool] = (saved_outputs[tool], final_path)
                    n_run += 1
                except Exception as e:
                    failures.append(f"{tool}: {str(e)[:80]}")
        finally:
            self.work_dir = saved_work_dir
            if saved_substrates is not None:
                self.files["substrates_filtered"] = saved_substrates
            if saved_passenger is not None:
                self.files["proteins_for_passenger_tools"] = saved_passenger
            for tool, info in saved_outputs.items():
                if isinstance(info, tuple):
                    main_path, t5_path = info
                    self.files[tool] = main_path
                    self.files[f"t5ass_whole_{tool}"] = t5_path
                elif info is not None:
                    self.files[tool] = info
            try:
                shutil.rmtree(second_work_dir)
            except OSError:
                pass

        msg = f"Ran {n_run}/{len(self._T5ASS_WHOLE_TOOLS)} tools on T5aSS substrates"
        if failures:
            msg += f"; {len(failures)} failure(s): " + "; ".join(failures[:3])
        return StepResult("t5ass_whole_annotations", True, msg)

    def _step_integrate(self) -> StepResult:
        output = self._wf(f"{self.config.sample_id}_integrated.csv")

        # Annotation tools only (NOT SignalP — that's a prediction tool,
        # already included via cross_validate → substrates_filtered).
        # EggNOG + pLM-BLAST added in 3.2.d; PLM-Effector stays prediction-
        # tier and is consumed by cross_validate, not integrate.
        annotation_files = []
        for key in [
            "blastp",
            "hhsuite",
            "interproscan",
            "eggnog",
            "plm_blast",
            "protparam",
            # t5ass_whole_* files from the optional second annotation pass
            # carry t5ass_whole_-prefixed column headers, so they merge as
            # side-by-side columns rather than overwriting the main-pass
            # passenger annotations. Empty when --t5ass-annotate-whole=False.
            "t5ass_whole_blastp",
            "t5ass_whole_hhsuite",
            "t5ass_whole_eggnog",
            "t5ass_whole_plm_blast",
            "t5ass_whole_protparam",
            # Side-car from _step_build_passenger_fasta: one row per T5aSS
            # substrate with t5_annotation_source ∈ {passenger, full}.
            "t5_annotation_source",
        ]:
            if key in self.files and os.path.exists(self.files[key]):
                annotation_files.append(self.files[key])

        # Also pass gene_info for GBFF annotations
        gene_info = self.files.get("gene_info", "")

        args = [
            "--substrates-filtered",
            self.files.get("substrates_filtered", ""),
            "--substrates-all",
            self.files.get("substrates_all", ""),
            "--sample",
            self.config.sample_id,
            "--output",
            output,
        ]
        proteins = self.files.get("proteins", "")
        if gene_info and os.path.exists(gene_info):
            args.extend(["--gene-info", gene_info])
        if proteins and os.path.exists(proteins):
            args.extend(["--proteins", proteins])
        if annotation_files:
            args.extend(["--annotations"] + annotation_files)

        rc, stdout, stderr = run_script("integrate_annotations.py", args)
        if rc == 0:
            self.files["integrated"] = output
            return StepResult("integrate", True, "Annotations integrated")
        return StepResult("integrate", False, stderr[:500])

    def _step_orthologs(self) -> StepResult:
        """Assign ortholog groups via all-vs-all BLASTp + Union-Find clustering."""
        sf = self.files.get("substrates_filtered", "")
        proteins = self.files.get("proteins", "")

        if not sf or not os.path.exists(sf):
            return StepResult("orthologs", True, "No substrates — skipping orthologs")

        # Count substrates
        n_sub = sum(1 for _ in open(sf)) - 1
        if n_sub < 2:
            return StepResult(
                "orthologs",
                True,
                f"Only {n_sub} substrate(s) — skipping ortholog grouping",
            )

        # Extract substrate sequences into a dedicated FASTA
        substrate_fasta = self._wf(f"{self.config.sample_id}_substrates_for_ortho.faa")
        substrate_ids = set()
        with open(sf) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                lid = row.get("locus_tag", row.get("protein_id", ""))
                if lid:
                    substrate_ids.add(lid)

        from Bio import SeqIO

        n_written = 0
        with open(substrate_fasta, "w") as out:
            for rec in SeqIO.parse(proteins, "fasta"):
                if rec.id in substrate_ids:
                    SeqIO.write(rec, out, "fasta")
                    n_written += 1

        if n_written < 2:
            return StepResult("orthologs", True, f"Only {n_written} substrate sequences — skipping")

        output = self._wf(f"{self.config.sample_id}_orthologs.csv")
        output_groups = self._wf(f"{self.config.sample_id}_ortholog_groups.csv")

        rc, stdout, stderr = run_script(
            "run_ortholog_grouping.py",
            [
                "--substrates-fasta",
                substrate_fasta,
                "--min-pident",
                str(self.config.ortholog_min_pident),
                "--min-qcov",
                str(self.config.ortholog_min_qcov),
                "--output",
                output,
                "--output-groups",
                output_groups,
                # Sequential step, gets the full per-genome budget.
                "--threads",
                str(self.config.cpu_per_genome),
            ],
            timeout=3600,
        )

        if rc == 0:
            self.files["orthologs"] = output
            self.files["ortholog_groups"] = output_groups

            # Merge ortholog assignments into the integrated CSV
            integrated = self.files.get("integrated", "")
            if integrated and os.path.exists(integrated) and os.path.exists(output):
                try:
                    import pandas as pd

                    df_int = pd.read_csv(integrated)
                    df_og = pd.read_csv(output)
                    # Merge on locus_tag
                    df_merged = df_int.merge(df_og, on="locus_tag", how="left")
                    df_merged.to_csv(integrated, index=False)
                    logger.info(f"Merged ortholog groups into {integrated}")
                except Exception as e:
                    logger.warning(f"Could not merge orthologs into integrated CSV: {e}")

            return StepResult("orthologs", True, "Ortholog groups assigned")
        return StepResult("orthologs", False, stderr[:500])

    def _step_enrichment(self) -> StepResult:
        """Per-SS-type circular-shift enrichment test.

        Opt-in: only runs when --enrichment-stats is on. That flag forces
        whole-genome DLP/DSE (Config.__post_init__), so every gene has a
        prediction and the rotation null is the background — no null sample
        needed. PLM-Effector is deliberately not an enrichment predictor (it
        over-predicts at genome scale; openspec enrichment-background-and-plme-default-off).
        """
        if not self.config.enrichment_stats:
            return StepResult("enrichment", True, "Skipped (--enrichment-stats not set)")

        ss_components = self.files.get("ss_components", "")
        gene_order = self.files.get("gene_order", "")
        dlp = self.files.get("deeplocpro", "")
        dse = self.files.get("deepsece", "")
        signalp = self.files.get("signalp", "")
        for name, path in (
            ("ss_components", ss_components),
            ("gene_order", gene_order),
            ("dlp", dlp),
            ("dse", dse),
            ("signalp", signalp),
        ):
            if not path or not os.path.exists(path):
                return StepResult("enrichment", False, f"Missing upstream file for enrichment: {name}")

        out = self._wf(f"{self.config.sample_id}_enrichment_stats.tsv")
        nulls_out = self._wf(f"{self.config.sample_id}_enrichment_nulls.npz")
        enr_args = [
            "--ss-components",
            ss_components,
            "--gene-order",
            gene_order,
            "--dlp",
            dlp,
            "--dse",
            dse,
            "--signalp",
            signalp,
            "--window",
            str(self.config.proximity_window),
            "--conf-threshold",
            str(self.config.conf_threshold),
            "--sample",
            self.config.sample_id,
            "--out",
            out,
            "--nulls-out",
            nulls_out,
        ]
        rc, _stdout, stderr = run_script("enrichment_testing.py", enr_args)

        if rc == 0:
            self.files["enrichment_stats"] = out
            self.files["enrichment_nulls"] = nulls_out
            return StepResult("enrichment", True, "Enrichment analysis complete")
        return StepResult("enrichment", False, stderr[:500])

    def _step_enrichment_figure(self) -> StepResult:
        """Render the circular-shift fold-enrichment bar chart.

        Optional step (gated on --enrichment-stats). Reads the per-type enrichment
        stats TSV from _step_enrichment (the null-array dump is not needed here; it
        is only used by the cross-genome pooling math).
        """
        stats = self.files.get("enrichment_stats", "")
        if not stats or not os.path.exists(stats):
            return StepResult("enrichment_figure", True, "Skipped (no enrichment output)")

        fig_dir = self._wf("figures")
        os.makedirs(fig_dir, exist_ok=True)
        # Two figures: per-tool (DLP + DSE bars) and the combined DLP-or-DSE bar.
        out_png = os.path.join(fig_dir, f"{self.config.sample_id}_enrichment_fold.png")
        combined_png = os.path.join(fig_dir, f"{self.config.sample_id}_enrichment_fold_combined.png")
        rc, _stdout, stderr = run_script(
            "run_enrichment_figure.py",
            ["--stats", stats, "--out", out_png, "--dpi", str(self.config.dpi)],
        )
        rc_c, _o, err_c = run_script(
            "run_enrichment_figure.py",
            ["--stats", stats, "--out", combined_png, "--combined", "--dpi", str(self.config.dpi)],
        )
        if rc_c != 0:
            logger.warning("combined enrichment figure skipped: %s", err_c[:160])
        if rc == 0:
            self.files["enrichment_figure"] = out_png
            return StepResult("enrichment_figure", True, "Enrichment figure generated")
        # Non-core: a plotting failure shouldn't fail the genome
        return StepResult("enrichment_figure", True, f"Enrichment figure skipped: {stderr[:200]}")

    def _step_report(self) -> StepResult:
        integrated = self.files.get("integrated", "")
        if not integrated or not os.path.exists(integrated):
            return StepResult("report", False, "No integrated CSV — skipping report")

        html_out = self._wf(f"{self.config.sample_id}_report.html")
        txt_out = self._wf(f"{self.config.sample_id}_report.txt")

        rc, stdout, stderr = run_script(
            "generate_report.py",
            [
                "--master-csvs",
                integrated,
                "--out-html",
                html_out,
                "--out-txt",
                txt_out,
            ],
        )

        if rc == 0:
            self.files["report_html"] = html_out
            self.files["report_txt"] = txt_out
            return StepResult("report", True, "Report generated")
        return StepResult("report", False, stderr[:500])

    def _step_figures(self) -> StepResult:
        integrated = self.files.get("integrated", "")
        if not integrated or not os.path.exists(integrated):
            return StepResult("figures", False, "No integrated CSV — skipping figures")

        fig_dir = self._wf("figures")
        os.makedirs(fig_dir, exist_ok=True)

        fig_args = [
            "--master-csvs",
            integrated,
            "--outdir",
            fig_dir,
            "--dpi",
            str(self.config.dpi),
        ]
        # Pass figure toggles (config field -> generate_figures --no-* flag)
        toggles = {
            "fig_ss_comp": "--no-ss-comp",
            "fig_physicochemical": "--no-physicochemical",
            "fig_func_summary": "--no-func-summary",
        }
        for field, flag in toggles.items():
            if not getattr(self.config, field):
                fig_args.append(flag)

        rc, stdout, stderr = run_script("generate_figures.py", fig_args)

        if rc == 0:
            self.files["figures_dir"] = fig_dir
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
            # Snapshot self.results / self.files under the lock so a parallel
            # step's append doesn't trip a "list changed during iteration".
            with self._state_lock:
                steps_snapshot = [{"name": r.name, "success": r.success, "message": r.message} for r in self.results]
                files_snapshot = dict(self.files)
            progress = {
                "sample_id": self.config.sample_id,
                "work_dir": self.work_dir,
                "steps": steps_snapshot,
                "files": files_snapshot,
                "config": asdict(self.config),
            }
            # Per-genome progress file in hidden .ssign/ subdirectory
            sid = self.config.sample_id
            progress_dir = outdir / ".ssign"
            progress_dir.mkdir(exist_ok=True)
            progress_path = progress_dir / f"{sid}_progress.json"
            with open(progress_path, "w") as f:
                json.dump(progress, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save progress: {e}")

    @staticmethod
    def load_progress(outdir, sample_id=None):
        """Load pipeline progress from a previous run.

        Returns (list[StepResult], dict_files, work_dir, config_dict) or
        (None, {}, "", {}).
        """
        # Try per-genome progress file in .ssign/ subdirectory, then legacy locations
        if sample_id:
            progress_path = Path(outdir) / ".ssign" / f"{sample_id}_progress.json"
            if not progress_path.exists():
                # Legacy: progress file in outdir root
                progress_path = Path(outdir) / f"{sample_id}_progress.json"
        else:
            progress_path = Path(outdir) / ".ssign" / "ssign_progress.json"
        if not progress_path.exists():
            progress_path = Path(outdir) / "ssign_progress.json"
        if not progress_path.exists():
            return None, {}, "", {}
        try:
            with open(progress_path) as f:
                data = json.load(f)
            results = [StepResult(s["name"], s["success"], s["message"]) for s in data.get("steps", [])]
            files = data.get("files", {})
            work_dir = data.get("work_dir", "")
            saved_config = data.get("config", {})
            return results, files, work_dir, saved_config
        except Exception:
            return None, {}, "", {}

    def _try_resume(self) -> set:
        """Attempt to resume from previous progress. Returns set of step
        names to skip (already completed successfully with valid outputs).

        Validates that every persisted output file is still present and
        non-empty before agreeing to resume — a manifest that survives
        but with corrupt / truncated / deleted outputs is worse than a
        fresh run, because downstream steps would silently read garbage.
        """
        prev_results, prev_files, prev_work_dir, _saved_config = self.load_progress(
            self.config.outdir, self.config.sample_id
        )
        if prev_results is None:
            return set()

        if not prev_work_dir or not os.path.isdir(prev_work_dir):
            logger.info("Previous work directory gone — starting fresh")
            return set()

        # Validate every recorded output: must exist and be non-empty. A
        # zero-byte TSV looks valid to os.path.exists but breaks downstream.
        # If any are bad, refuse to resume — better a clean re-run than a
        # corrupt resume that taints the report.
        bad: list[str] = []
        for key, path in prev_files.items():
            if not isinstance(path, str) or not path:
                continue
            if not os.path.exists(path):
                bad.append(f"{key}: missing ({path})")
                continue
            try:
                if os.path.isfile(path) and os.path.getsize(path) == 0:
                    bad.append(f"{key}: empty ({path})")
            except OSError as e:
                bad.append(f"{key}: stat failed ({e})")

        if bad:
            logger.warning(
                "Resume aborted: %d previous output(s) missing or empty — starting fresh.\n  %s",
                len(bad),
                "\n  ".join(bad[:10]),
            )
            return set()

        completed = {r.name for r in prev_results if r.success}
        if completed:
            self.work_dir = prev_work_dir
            for key, path in prev_files.items():
                if isinstance(path, str):
                    self.files[key] = path
            logger.info(f"Resuming: {len(completed)} steps already done: {', '.join(sorted(completed))}")

        return completed

    def _copy_outputs(self):
        """Copy consolidated outputs to the user's output directory.

        Produces:
          {sample_id}_results.csv     - Normal CSV (secreted proteins, their SS, other SS)
          {sample_id}_results_raw.csv - Raw CSV with all tool data
          {sample_id}_summary.txt     - Report + enrichment stats
          figures/                    - Publication-quality plots
        """
        outdir = Path(self.config.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        sid = self.config.sample_id

        # 1. Normal results CSV (chunked: secreted proteins → their SS → other SS)
        self._build_master_csv(outdir / f"{sid}_results.csv")

        # 2. Raw results CSV (all columns from all tools)
        self._build_raw_csv(outdir / f"{sid}_results_raw.csv")

        # 3. Summary — report text + enrichment summary + Fisher table
        self._build_summary(outdir / f"{sid}_summary.txt")

        # 4. Figures directory (per-sample subfolder for parallel safety)
        if "figures_dir" in self.files and os.path.exists(self.files["figures_dir"]):
            dest = outdir / "figures" / sid
            dest.mkdir(parents=True, exist_ok=True)
            for fig_file in Path(self.files["figures_dir"]).iterdir():
                if fig_file.is_file():
                    shutil.copy2(fig_file, dest / fig_file.name)

        # 5. Enrichment stats + null dump (when --enrichment-stats ran): a standalone
        # table for the user, and the inputs cross-genome pooling reads back.
        for key in ("enrichment_stats", "enrichment_nulls"):
            src = self.files.get(key, "")
            if src and os.path.exists(src):
                shutil.copy2(src, outdir / Path(src).name)

    def _load_systems(self):
        """Load and combine system/component DataFrames, filtered by excluded."""
        import pandas as pd

        excluded = set(self.config.excluded_systems or [])

        df_sys = None
        fpath = self.files.get("valid_systems", "")
        if fpath and os.path.exists(fpath):
            try:
                df_sys = pd.read_csv(fpath, sep="\t")
                df_sys.insert(0, "record_type", "system")
            except Exception:
                pass

        df_comp = None
        fpath = self.files.get("ss_components", "")
        if fpath and os.path.exists(fpath):
            try:
                df_comp = pd.read_csv(fpath, sep="\t")
                df_comp.insert(0, "record_type", "component")
            except Exception:
                pass

        # Combine
        frames = [df for df in [df_sys, df_comp] if df is not None]
        if not frames:
            return pd.DataFrame()

        df = pd.concat(frames, ignore_index=True)

        # Sort: sample_id → ss_type → record_type (system before component)
        type_order = {"system": 0, "component": 1}
        df["_sort"] = df["record_type"].map(type_order)
        sort_cols = []
        if "sample_id" in df.columns:
            sort_cols.append("sample_id")
        if "ss_type" in df.columns:
            sort_cols.append("ss_type")
        sort_cols.append("_sort")
        df = df.sort_values(sort_cols).drop(columns=["_sort"])

        return df, excluded

    def _build_master_csv(self, output_path: Path):
        """Build the normal results CSV with three chunks:

        1. Secreted proteins (full annotation columns from integrated CSV)
        2. Secretion systems that have associated secreted proteins
        3. Remaining (non-excluded) secretion systems

        Excluded systems and their associated proteins are omitted.
        Blank rows separate the three chunks.
        """
        import pandas as pd

        excluded = set(self.config.excluded_systems or [])

        # ── Load systems ──
        df_sys = None
        fpath = self.files.get("valid_systems", "")
        if fpath and os.path.exists(fpath):
            try:
                df_sys = pd.read_csv(fpath, sep="\t")
                df_sys.insert(0, "record_type", "system")
            except Exception:
                pass

        df_comp = None
        fpath = self.files.get("ss_components", "")
        if fpath and os.path.exists(fpath):
            try:
                df_comp = pd.read_csv(fpath, sep="\t")
                df_comp.insert(0, "record_type", "component")
            except Exception:
                pass

        # Combine and filter out excluded
        sys_frames = [df for df in [df_sys, df_comp] if df is not None]
        df_systems = pd.DataFrame()
        if sys_frames:
            df_systems = pd.concat(sys_frames, ignore_index=True)
            if "ss_type" in df_systems.columns and excluded:
                df_systems = df_systems[~df_systems["ss_type"].isin(excluded)]
            type_order = {"system": 0, "component": 1}
            df_systems["_sort"] = df_systems["record_type"].map(type_order)
            sort_cols = []
            if "sample_id" in df_systems.columns:
                sort_cols.append("sample_id")
            if "ss_type" in df_systems.columns:
                sort_cols.append("ss_type")
            sort_cols.append("_sort")
            df_systems = df_systems.sort_values(sort_cols).drop(columns=["_sort"])

        # ── Load secreted proteins (integrated annotations) ──
        df_subs = pd.DataFrame()
        fpath = self.files.get("integrated", "")
        if fpath and os.path.exists(fpath):
            try:
                df_subs = pd.read_csv(fpath)
            except Exception:
                pass

        # Filter proteins from excluded systems
        if not df_subs.empty and "nearby_ss_types" in df_subs.columns and excluded:

            def _has_non_excluded_ss(val):
                if pd.isna(val) or not val:
                    return False
                types = {t.strip() for t in str(val).split(",")}
                return bool(types - excluded)

            df_subs = df_subs[df_subs["nearby_ss_types"].apply(_has_non_excluded_ss)]

        # Sort secreted proteins
        sort_cols = []
        if not df_subs.empty:
            if "sample_id" in df_subs.columns:
                sort_cols.append("sample_id")
            if "nearby_ss_types" in df_subs.columns:
                sort_cols.append("nearby_ss_types")
            if sort_cols:
                df_subs = df_subs.sort_values(sort_cols)

        # ── Organize columns like reference CSV ──
        # Priority column order (present columns are moved to front)
        # Column order matching the reference CSV structure
        priority_cols = [
            # Identity
            "locus_tag",
            "sample_id",
            # Consensus annotations (Phase 3 — computed if available)
            "broad_consensus_annotation",
            "broad_annotation",
            "detailed_annotation",
            "detailed_consensus_annotation",
            "evidence_keywords",
            "n_tools_agreeing",
            "n_tools_with_hits",
            "concordance_ratio",
            "confidence_tier",
            # Physicochemical properties
            "aa_length",
            "gravy",
            "mw_da",
            "isoelectric_point",
            "charge_ph7",
            "instability_index",
            "aromaticity",
            # Secretion system context
            "nearby_ss_types",
            "secretion_evidence",
            "is_secreted",
            # Localization predictions (DeepLocPro)
            "predicted_localization",
            "dlp_extracellular_prob",
            "dlp_max_localization",
            "dlp_max_probability",
            "periplasmic_prob",
            "outer_membrane_prob",
            "cytoplasmic_prob",
            # Secretion type prediction (DeepSecE)
            "dse_ss_type",
            "dse_max_prob",
            # Signal peptide prediction (SignalP)
            "signalp_prediction",
            "signalp_probability",
            "signalp_cs_position",
            # GBFF original annotation
            "gbff_annotation",
            # BLASTp hits
            "blastp_hit_accession",
            "blastp_hit_description",
            "blastp_pident",
            "blastp_qcov",
            "blastp_evalue",
            # HHpred Pfam
            "pfam_top1_id",
            "pfam_top1_description",
            "pfam_top1_probability",
            "pfam_top1_evalue",
            "pfam_top1_score",
            # HHpred PDB
            "pdb_top1_id",
            "pdb_top1_description",
            "pdb_top1_probability",
            "pdb_top1_evalue",
            "pdb_top1_score",
            # InterProScan
            "interpro_domains",
            "interpro_go_terms",
            "interpro_pfam_ids",
            "interpro_descriptions",
            # Ortholog groups
            "ortholog_group",
            "og_n_members",
            "og_mean_pident",
            # Annotation tool count
            "annotation_tools",
        ]
        if not df_subs.empty:
            existing_priority = [c for c in priority_cols if c in df_subs.columns]
            remaining = [c for c in df_subs.columns if c not in existing_priority and c != "sequence"]
            # Put sequence at the end
            col_order = existing_priority + sorted(remaining)
            if "sequence" in df_subs.columns:
                col_order.append("sequence")
            df_subs = df_subs[col_order]

        # ── Determine which SS have associated secreted proteins ──
        ss_types_with_subs = set()
        if not df_subs.empty and "nearby_ss_types" in df_subs.columns:
            for val in df_subs["nearby_ss_types"].dropna():
                for t in str(val).split(","):
                    ss_types_with_subs.add(t.strip())

        df_systems_with_subs = pd.DataFrame()
        df_systems_other = pd.DataFrame()
        if not df_systems.empty and "ss_type" in df_systems.columns:
            mask = df_systems["ss_type"].isin(ss_types_with_subs)
            df_systems_with_subs = df_systems[mask]
            df_systems_other = df_systems[~mask]

        # ── Write chunked CSV ──
        with open(output_path, "w", newline="") as f:
            chunk_written = False

            # Chunk 1: Secreted proteins
            if not df_subs.empty:
                f.write("# Secreted Proteins\n")
                df_subs.to_csv(f, index=False)
                chunk_written = True

            # Chunk 2: SS with associated secreted proteins
            if not df_systems_with_subs.empty:
                if chunk_written:
                    f.write("\n")
                f.write("# Secretion Systems (with secreted proteins)\n")
                df_systems_with_subs.to_csv(f, index=False)
                chunk_written = True

            # Chunk 3: Other non-excluded SS
            if not df_systems_other.empty:
                if chunk_written:
                    f.write("\n")
                f.write("# Secretion Systems (other)\n")
                df_systems_other.to_csv(f, index=False)

    def _build_raw_csv(self, output_path: Path):
        """Build raw CSV with ALL data from ALL tools — every protein, every
        per-tool intermediate column, no filtering, no column pruning.

        Distinct from the master ``_results.csv`` (which is the
        substrate-filtered consolidated view): the raw file is a debug /
        post-mortem dump where every column from every tool intermediate
        is left-joined onto the gene_info base, keyed by locus_tag.

        Joined sources (each contributes its full per-tool column set;
        empty values where the tool didn't run or skipped a protein):
          - ``gene_info``      — all proteins extracted from the genome
          - ``predictions``    — cross_validate output (DLP/DSE/SignalP/PLM-E)
          - ``plm_effector``   — merged per-type PLM-E TSV
          - ``substrates_all`` — proximity + T5SS-self substrate flags
          - per-tool annotation outputs (blastp, hhsuite, interproscan,
            eggnog, plm_blast, protparam)
        """
        import pandas as pd

        def _read_tsv_or_csv(path: str):
            """Auto-detect TSV vs CSV and read; return empty DataFrame on failure."""
            if not path or not os.path.exists(path):
                return pd.DataFrame()
            try:
                return pd.read_csv(path, sep=None, engine="python")
            except Exception as e:
                logger.warning("Could not read %s for raw CSV: %s", path, e)
                return pd.DataFrame()

        def _left_join(base: pd.DataFrame, addition: pd.DataFrame, label: str) -> pd.DataFrame:
            """Left-join `addition` onto `base` on locus_tag. Prefix-disambiguate
            overlapping columns (other than the join key) with the source label
            so we never silently drop per-tool intermediates."""
            if addition.empty or "locus_tag" not in addition.columns:
                return base
            if base.empty:
                return addition
            overlap = (set(addition.columns) & set(base.columns)) - {"locus_tag"}
            if overlap:
                addition = addition.rename(columns={c: f"{label}__{c}" for c in overlap})
            addition = addition.drop_duplicates(subset="locus_tag", keep="first")
            return base.merge(addition, on="locus_tag", how="left")

        # ── Base: every protein in the genome ──
        df = _read_tsv_or_csv(self.files.get("gene_info", ""))
        if df.empty or "locus_tag" not in df.columns:
            # Fallback to integrated CSV path (preserves old behaviour
            # for runs that never ran extract_proteins, e.g. failure tests).
            df = _read_tsv_or_csv(self.files.get("integrated", ""))
            if not df.empty:
                df.to_csv(output_path, index=False)
            else:
                with open(output_path, "w") as f:
                    f.write("locus_tag,sample_id\n")
            return

        df["sample_id"] = self.config.sample_id

        # ── Layer in every intermediate ──
        # Order matters only for label-prefixing of overlapping columns.
        sources = [
            ("pred", self.files.get("predictions", "")),
            ("plme", self.files.get("plm_effector", "")),
            ("substrates_all", self.files.get("substrates_all", "")),
            ("blastp", self.files.get("blastp", "")),
            ("hhsuite", self.files.get("hhsuite", "")),
            ("ips", self.files.get("interproscan", "")),
            ("eggnog", self.files.get("eggnog", "")),
            ("plmblast", self.files.get("plm_blast", "")),
            ("protparam", self.files.get("protparam", "")),
        ]
        for label, path in sources:
            df = _left_join(df, _read_tsv_or_csv(path), label)

        df.to_csv(output_path, index=False)

    def _build_summary(self, output_path: Path):
        """Combine report text + enrichment-stats table (when present)."""
        parts = []

        if "report_txt" in self.files and os.path.exists(self.files["report_txt"]):
            with open(self.files["report_txt"]) as f:
                parts.append(f.read())

        if "enrichment_stats" in self.files and os.path.exists(self.files["enrichment_stats"]):
            try:
                import pandas as pd

                df = pd.read_csv(self.files["enrichment_stats"], sep="\t")
                if not df.empty:
                    parts.append(
                        "\n\n"
                        + "=" * 60
                        + "\nENRICHMENT ANALYSIS\n"
                        + "=" * 60
                        + "\n\n"
                        + df.to_string(index=False)
                        + "\n"
                    )
            except Exception:
                pass

        if parts:
            with open(output_path, "w") as f:
                f.write("".join(parts))


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
    import pandas as pd
    from Bio import SeqIO

    result = {
        "n_proteins": 0,
        "n_groups": 0,
        "genomes_updated": 0,
        "combined_fasta": "",
        "orthologs_csv": "",
        "ortholog_groups_csv": "",
    }

    if progress_callback:
        progress_callback("Cross-genome orthologs", 10, "Collecting substrates from all genomes...")

    os.makedirs(output_dir, exist_ok=True)
    combined_fasta = os.path.join(output_dir, "all_substrates_combined.faa")

    # Step 1: Collect all substrate sequences
    all_substrate_ids = {}  # locus_tag -> genome_outdir
    n_written = 0

    with open(combined_fasta, "w") as out_f:
        for gdir in genome_outdirs:
            gpath = Path(gdir)

            # Find integrated CSV to get substrate locus_tags
            integrated_csvs = list(gpath.glob("*_integrated.csv")) + list(gpath.glob("*integrated*.csv"))
            if not integrated_csvs:
                logger.warning(f"No integrated CSV in {gdir} -- skipping")
                continue

            integrated_csv = str(integrated_csvs[0])
            substrate_ids = set()
            try:
                df = pd.read_csv(integrated_csv)
                if "locus_tag" in df.columns:
                    substrate_ids = set(df["locus_tag"].dropna().astype(str))
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
                    for rec in SeqIO.parse(str(pf), "fasta"):
                        if rec.id in substrate_ids and rec.id not in all_substrate_ids:
                            SeqIO.write(rec, out_f, "fasta")
                            all_substrate_ids[rec.id] = gdir
                            n_written += 1
                except Exception:
                    continue

    result["combined_fasta"] = combined_fasta
    result["n_proteins"] = n_written

    if n_written < 2:
        logger.info(f"Only {n_written} substrate(s) across all genomes -- skipping cross-genome ortholog grouping")
        return result

    if progress_callback:
        progress_callback(
            "Cross-genome orthologs",
            40,
            f"Running all-vs-all BLASTp on {n_written} substrates...",
        )

    # Step 2: Run ortholog grouping
    orthologs_csv = os.path.join(output_dir, "cross_genome_orthologs.csv")
    groups_csv = os.path.join(output_dir, "cross_genome_ortholog_groups.csv")

    rc, stdout, stderr = run_script(
        "run_ortholog_grouping.py",
        [
            "--substrates-fasta",
            combined_fasta,
            "--min-pident",
            str(min_pident),
            "--min-qcov",
            str(min_qcov),
            "--output",
            orthologs_csv,
            "--output-groups",
            groups_csv,
        ],
        timeout=3600,
    )

    if rc != 0:
        logger.error(f"Cross-genome ortholog grouping failed: {stderr[:500]}")
        return result

    result["orthologs_csv"] = orthologs_csv
    result["ortholog_groups_csv"] = groups_csv

    if progress_callback:
        progress_callback(
            "Cross-genome orthologs",
            70,
            "Merging ortholog groups into per-genome results...",
        )

    # Step 3: Read ortholog assignments and merge into each genome's integrated CSV
    try:
        df_og = pd.read_csv(orthologs_csv)
        result["n_groups"] = df_og["ortholog_group"].nunique()
    except Exception as e:
        logger.warning(f"Could not read ortholog results: {e}")
        return result

    # Prefix with "xg_" (cross-genome) to distinguish from per-genome groups
    df_og = df_og.rename(
        columns={
            "ortholog_group": "xg_ortholog_group",
            "og_n_members": "xg_og_n_members",
            "og_mean_pident": "xg_og_mean_pident",
        }
    )

    genomes_updated = 0
    for gdir in genome_outdirs:
        gpath = Path(gdir)
        integrated_csvs = list(gpath.glob("*_integrated.csv")) + list(gpath.glob("*integrated*.csv"))
        if not integrated_csvs:
            continue

        integrated_csv = str(integrated_csvs[0])
        try:
            df_int = pd.read_csv(integrated_csv)
            for col in ["xg_ortholog_group", "xg_og_n_members", "xg_og_mean_pident"]:
                if col in df_int.columns:
                    df_int = df_int.drop(columns=[col])

            df_merged = df_int.merge(df_og, on="locus_tag", how="left")
            df_merged.to_csv(integrated_csv, index=False)
            genomes_updated += 1
            logger.info(f"Merged cross-genome orthologs into {integrated_csv}")
        except Exception as e:
            logger.warning(f"Could not merge into {integrated_csv}: {e}")

    result["genomes_updated"] = genomes_updated

    if progress_callback:
        progress_callback(
            "Cross-genome orthologs",
            100,
            f"Done: {result['n_groups']} groups across {n_written} substrates from {genomes_updated} genomes",
        )

    return result


def pool_enrichment_stats(per_genome_tsvs: list[str], output_tsv: str, nulls_output: str = "") -> int:
    """Pool per-genome circular-shift enrichment outputs into one cross-genome view.

    Each genome's per-(type, tool) test carries an exact rotation null (dumped
    next to the stats TSV as ``*_enrichment_nulls.npz``). The single-genome nulls
    can't be enumerated as an exact product across genomes, so pooling Monte-Carlos
    it (the validated fleet approach): per (display type, tool), sum the observed
    counts across genomes and, for each of ``ENRICH_POOL_REPS`` replicates, add one
    randomly-drawn rotation from every contributing genome. fold = pooled-observed /
    pooled-null-mean; p = (#{null >= observed} + 1)/(reps + 1); BH across the pooled
    (type x tool) tests.

    A genome contributes to a (type, tool) cell only when its npz carries that
    cell's null array; genomes missing the npz are skipped for pooling. When
    ``nulls_output`` is given, the pooled null arrays are dumped there (same
    ``enrich_null_key`` scheme) so the figure script can render a combined
    cross-genome figure. Returns the number of pooled rows written.
    """
    import csv
    import re as _re
    import sys as _sys

    import numpy as np

    _scripts = os.path.join(os.path.dirname(__file__), "..", "scripts")
    if _scripts not in _sys.path:
        _sys.path.insert(0, _scripts)
    from enrichment_testing import OUT_FIELDS, bh_fdr_by_family

    from ssign_app.scripts.ssign_lib.constants import ENRICH_POOL_REPS, ENRICH_POOL_SEED, enrich_null_key

    # (ss_type, tool, mode) -> {observed, n_mask, nulls: [array, ...]}. Mode keeps a
    # T5a/c self and hitchhiker(window) result pooled separately, not merged.
    accum: dict[tuple[str, str, str], dict] = {}
    for tsv in per_genome_tsvs:
        if not tsv or not os.path.exists(tsv):
            continue
        npz_path = _re.sub(r"_enrichment_stats\.tsv$", "_enrichment_nulls.npz", tsv)
        if os.path.exists(npz_path):
            with np.load(npz_path) as npz:  # materialise + close the handle (don't leak fds across genomes)
                nulls = {k: npz[k] for k in npz.files}
        else:
            nulls = {}
        with open(tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                ss_type, tool = row.get("ss_type", ""), row.get("tool", "")
                mode = row.get("mode", "window") or "window"
                key_npz = enrich_null_key(ss_type, tool, mode)
                if not row.get("observed", "").strip() or key_npz not in nulls:
                    continue  # skipped row (e.g. DSE-T3SS) or no null to pool
                slot = accum.setdefault(
                    (ss_type, tool, mode),
                    {"observed": 0, "n_mask": 0, "nulls": []},
                )
                slot["observed"] += int(row["observed"])
                slot["n_mask"] += int(row.get("n_mask") or 0)
                slot["nulls"].append(nulls[key_npz])

    rng = np.random.default_rng(ENRICH_POOL_SEED)
    pooled = []
    pooled_nulls: dict[str, np.ndarray] = {}
    for (ss_type, tool, mode), s in accum.items():
        null_pool = np.zeros(ENRICH_POOL_REPS, dtype=np.int64)
        for arr in s["nulls"]:
            if arr.size:
                null_pool += arr[rng.integers(0, arr.size, size=ENRICH_POOL_REPS)]
        null_mean = float(null_pool.mean())
        obs = s["observed"]
        fold = round(obs / null_mean, 4) if null_mean > 0 else (float("inf") if obs > 0 else 0.0)
        p = (int(np.sum(null_pool >= obs)) + 1) / (ENRICH_POOL_REPS + 1)
        pooled.append(
            {
                "ss_type": ss_type,
                "tool": tool,
                "mode": mode,
                "observed": obs,
                "n_mask": s["n_mask"],
                "null_mean": round(null_mean, 4),
                "fold": fold,
                "p_perm": round(p, 6),
                "n_null": ENRICH_POOL_REPS,
            }
        )
        pooled_nulls[enrich_null_key(ss_type, tool, mode)] = null_pool

    bh_fdr_by_family(pooled)  # per-tool vs combined DLP-or-DSE families (shared rule)

    with open(output_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in pooled:
            writer.writerow({**r, "sample_id": "POOLED"})

    if nulls_output:
        np.savez_compressed(nulls_output, **pooled_nulls)

    return len(pooled)


def pool_and_plot_enrichment(per_genome_tsvs: list[str], out_dir: str) -> int:
    """Pool >=2 genomes' enrichment stats and render the combined figure.

    Shared by the CLI multi-genome path (MultiGenomeRunner) and the GUI so both
    emit the same cross-genome view. Writes pooled_enrichment_stats.tsv to out_dir
    and figures/pooled_enrichment_fold.png (the same combined fold-enrichment bar
    chart as a single genome, computed over all genomes' pooled per-type nulls).
    No-op (returns 0) when fewer than two genomes carry enrichment output. Each
    per-genome TSV must have its sibling *_enrichment_nulls.npz (pool_enrichment_stats
    finds it by name to recompute the pooled fold/p)."""
    tsvs = [p for p in per_genome_tsvs if p and os.path.exists(p)]
    if len(tsvs) < 2:
        return 0
    pooled_tsv = os.path.join(out_dir, "pooled_enrichment_stats.tsv")
    n = pool_enrichment_stats(tsvs, pooled_tsv)
    if n:
        fig_dir = os.path.join(out_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)
        out_png = os.path.join(fig_dir, "pooled_enrichment_fold.png")
        combined_png = os.path.join(fig_dir, "pooled_enrichment_fold_combined.png")
        rc, _out, err = run_script("run_enrichment_figure.py", ["--stats", pooled_tsv, "--out", out_png])
        if rc != 0:
            logger.warning("pooled enrichment figure failed: %s", err[:160])
        rc_c, _o, err_c = run_script(
            "run_enrichment_figure.py", ["--stats", pooled_tsv, "--out", combined_png, "--combined"]
        )
        if rc_c != 0:
            logger.warning("pooled combined enrichment figure failed: %s", err_c[:160])
    return n
