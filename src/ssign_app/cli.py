#!/usr/bin/env python3
"""CLI entry point for ssign.

Usage:
    ssign                           # launch Streamlit GUI
    ssign --version                 # print version
    ssign run input.gbff --outdir results [flags...]
                                    # run pipeline non-interactively (HPC,
                                    # scripting, batch use)
    ssign run --help                # full list of run-mode flags

Two execution modes:
  - GUI (default, no subcommand): launches Streamlit, suitable for desktop /
    Easy Mode use.
  - `run` subcommand: builds a PipelineConfig from CLI flags, drives
    PipelineRunner directly. Suitable for HPC, batch, and scripting.

The `run` subcommand exposes every PipelineConfig field as a flag, grouped
by phase. See `ssign run --help` for the full list.
"""

import argparse
import os
import re
import socket
import subprocess
import sys
import threading
from typing import TYPE_CHECKING

from ssign_app.scripts.ssign_lib.constants import DEFAULT_EXCLUDED_SYSTEMS  # noqa: F401  (argparse default below)

if TYPE_CHECKING:
    from ssign_app.core.runner import PipelineConfig

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


# ---------------------------------------------------------------------------
# `ssign run` subcommand
# ---------------------------------------------------------------------------


def _default_cpu_per_genome() -> int:
    """Default for --cpu-per-genome: cgroup-allocated count, never host total.

    Lazy import so `ssign --help` doesn't pay the import cost when the user
    isn't running the pipeline.
    """
    try:
        from ssign_app.scripts.ssign_lib.resources import effective_cpu_count

        return effective_cpu_count()
    except Exception:
        return os.cpu_count() or 4


def _add_run_parser(subparsers: argparse._SubParsersAction) -> None:
    """Build the `ssign run` subcommand parser.

    Every PipelineConfig field is exposed as a flag. Booleans use
    argparse.BooleanOptionalAction so each field also accepts its
    `--no-<flag>` inverse (e.g. `--no-skip-blastp`).
    """
    p = subparsers.add_parser(
        "run",
        help="Run the ssign pipeline non-interactively.",
        description=(
            "Run ssign on one or more input genomes (GenBank, GFF3, or FASTA). "
            "All PipelineConfig fields are exposed as flags. When N>1 genomes "
            "are passed, ssign pools predictions over neighborhoods and "
            "annotations over substrates so heavy startup costs (IPS JVM, "
            "EggNOG DB load, pLM-BLAST embeddings, PLM-Effector models) are "
            "paid once per batch rather than once per genome."
        ),
    )

    # ── Essentials ──────────────────────────────────────────────────────
    g = p.add_argument_group("essentials")
    g.add_argument(
        "input_path",
        nargs="+",
        help=(
            "One or more input genomes (GenBank .gbff/.gbk, GFF3 .gff, or "
            "FASTA). Pass multiple files to run them as a single batched job."
        ),
    )
    g.add_argument(
        "--outdir",
        default="./results",
        help=(
            "Output directory (default: ./results). For multi-genome runs, "
            "per-genome outputs land in <outdir>/<sample_id>/ and a "
            "combined_results.csv is written at the top level."
        ),
    )
    g.add_argument(
        "--sample-id",
        default="",
        help=(
            "Sample identifier used to prefix output files (single-genome only; "
            "for multi-genome runs the sample_id is derived per-genome from "
            "the input filename's stem)."
        ),
    )
    g.add_argument(
        "--combined-summary",
        dest="combined_summary",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Write a top-level combined_results.csv aggregating all genomes' "
            "substrates with a source_genome column (multi-genome only; "
            "default on). Flag is still named --combined-summary for back-compat."
        ),
    )
    g.add_argument(
        "--original-filename",
        default="",
        help="Original filename when input_path is a temp upload (informational).",
    )
    g.add_argument(
        "--resume",
        action="store_true",
        help="Skip steps that already have a successful entry in the progress "
        "manifest at <outdir>/.ssign/<sid>_progress.json.",
    )
    g.add_argument(
        "--tier",
        choices=("base", "extended", "full"),
        default=None,
        help=(
            "Install tier the run targets — sets each tool's default on/off "
            "state to what that tier ships. Leave unset to use what "
            "fetch_databases.sh recorded; defaults to 'extended'."
        ),
    )

    # ── Phase 2: SS detection ──────────────────────────────────────────
    g = p.add_argument_group("SS detection (MacSyFinder)")
    g.add_argument(
        "--wholeness-threshold", type=float, default=0.8, help="Minimum MacSyFinder system completeness (default: 0.8)."
    )
    g.add_argument(
        "--excluded-systems",
        nargs="+",
        default=list(DEFAULT_EXCLUDED_SYSTEMS),
        help=f"TXSScan models to exclude (default: {' '.join(DEFAULT_EXCLUDED_SYSTEMS)}). "
        "These are surface/uptake appendages, not protein secretion systems. "
        "T3SS is detected by default; DeepSecE is never trusted for T3SS calls. "
        "Append 'T3SS' to restore the old T3SS-excluded behaviour.",
    )
    g.add_argument(
        "--macsyfinder-db-type",
        choices=["ordered_replicon", "unordered"],
        default="ordered_replicon",
        help="MacSyFinder --db-type (default: ordered_replicon).",
    )
    g.add_argument(
        "--cpu-per-genome",
        type=int,
        default=_default_cpu_per_genome(),
        help="CPUs available to per-genome subtools (default: cgroup allocation, or all host CPUs).",
    )

    # ── Phase 3: Prediction thresholds ──────────────────────────────────
    g = p.add_argument_group("prediction thresholds")
    g.add_argument(
        "--conf-threshold", type=float, default=0.8, help="DeepLocPro extracellular probability minimum (default: 0.8)."
    )
    g.add_argument(
        "--proximity-window", type=int, default=3, help="+/-N genes per SS component for proximity (default: 3)."
    )
    g.add_argument(
        "--required-fraction-correct",
        type=float,
        default=0.8,
        help="Fraction of SS components correctly localized (default: 0.8).",
    )
    g.add_argument(
        "--dlp-confidence-threshold",
        type=float,
        default=0.8,
        help=(
            "Minimum DLP max-probability for an SS-machinery component to count in the "
            "localization-correctness gate (default: 0.8). Components below this are "
            "excluded from both numerator and denominator of fraction_correct."
        ),
    )
    g.add_argument(
        "--skip-localization-gate",
        action="store_true",
        help="Disable the literature-derived localization-correctness gate (debug / ad-hoc).",
    )
    g.add_argument(
        "--deepsece-min-prob", type=float, default=0.8, help="DeepSecE min probability to call secreted (default: 0.8)."
    )
    g.add_argument(
        "--signalp-min-prob",
        type=float,
        default=0.5,
        help="SignalP min probability for a signal peptide (default: 0.5).",
    )

    # ── Enrichment stats (opt-in) ───────────────────────────────────────
    g = p.add_argument_group("enrichment stats")
    g.add_argument(
        "--enrichment-stats",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Per-SS-type circular-shift enrichment test: emits fold (enrichment) + "
            "permutation p + BH q per system type for DeepLocPro, DeepSecE, and SignalP, "
            "plus enrichment figures. Forces whole-genome DeepLocPro + DeepSecE + SignalP "
            "(local; the rotation null needs every gene's positivity in gene order), which "
            "adds ~13 min/genome. Off by default."
        ),
    )

    # ── Phase 1: ORF prediction / annotation ────────────────────────────
    g = p.add_argument_group("ORF prediction + annotation")
    g.add_argument(
        "--use-input-annotations",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Trust input GenBank annotations (skip Bakta re-annotation).",
    )
    g.add_argument(
        "--run-bakta",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Run Bakta on FASTA contigs input (default: True per plan A.6). "
            "GenBank input is governed by --use-input-annotations instead."
        ),
    )
    g.add_argument("--bakta-db", default="", help="Path to Bakta database (required for Bakta runs).")
    g.add_argument(
        "--bakta-threads",
        type=int,
        default=0,
        help=(
            "Threads passed to Bakta (default: same as --cpu-per-genome, i.e. "
            "the cgroup-allocated count). Bakta enforces its own ceiling at the "
            "OS-visible CPU count and rejects values above it, so this knob is "
            "really only useful to set a lower-than-default thread count on "
            "shared machines."
        ),
    )

    # ── DTU prediction tools (DeepLocPro + SignalP) ─────────────────────
    g = p.add_argument_group("DTU prediction tools")
    g.add_argument(
        "--deeplocpro-mode",
        choices=["local", "remote"],
        default=None,
        help="DeepLocPro execution mode. Default: auto — local if 'deeplocpro' is "
        "on PATH or at --deeplocpro-path / $SSIGN_DEEPLOCPRO_PATH, else falls "
        "back to the DTU webserver with a warning. 'local' / 'remote' force "
        "the choice.",
    )
    g.add_argument(
        "--deeplocpro-path",
        default="",
        help="Path to local DeepLocPro install. Empty falls back to $SSIGN_DEEPLOCPRO_PATH, then PATH.",
    )
    g.add_argument(
        "--signalp-mode",
        choices=["local", "remote"],
        default=None,
        help="SignalP execution mode. Default: auto — local if 'signalp6' is "
        "on PATH or at --signalp-path / $SSIGN_SIGNALP_PATH, else falls back "
        "to the DTU webserver with a warning. 'local' / 'remote' force the choice.",
    )
    g.add_argument(
        "--signalp-path",
        default="",
        help="Path to local SignalP 6 install. Empty falls back to $SSIGN_SIGNALP_PATH, then PATH.",
    )
    g.add_argument(
        "--skip-deeplocpro",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip DeepLocPro step (overrides --tier default).",
    )
    g.add_argument(
        "--skip-signalp",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip SignalP step (overrides --tier default).",
    )
    g.add_argument(
        "--skip-deepsece",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip DeepSecE step (overrides --tier default).",
    )
    g.add_argument(
        "--dlp-whole-genome",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run DeepLocPro on every protein, not just the SS neighborhood.",
    )
    g.add_argument(
        "--dse-whole-genome",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run DeepSecE on every protein, not just the SS neighborhood.",
    )
    g.add_argument(
        "--sp-whole-genome",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run SignalP on every protein, not just the SS neighborhood.",
    )
    g.add_argument(
        "--plme-whole-genome",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run PLM-Effector on every protein, not just the SS neighborhood.",
    )
    g.add_argument(
        "--monitor-resources",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write outdir/runtime_data/{step_timings,resource_samples}.csv during a run. On by default.",
    )
    g.add_argument(
        "--monitor-interval-s",
        type=float,
        default=5.0,
        help="Sampling interval for resource_samples.csv (seconds, default 5).",
    )

    # ── Phase 5: Annotation tools ───────────────────────────────────────
    g = p.add_argument_group("BLASTp")
    g.add_argument(
        "--skip-blastp",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip BLASTp (overrides --tier default; on at tier=full only, defaults to Swiss-Prot).",
    )
    g.add_argument("--blastp-db", default="", help="Path to BLAST database (NR or Swiss-Prot).")
    g.add_argument("--blastp-exclude-taxid", default="", help="Comma-separated taxid(s) to exclude from BLASTp hits.")
    g.add_argument(
        "--blastp-min-pident", type=float, default=80.0, help="BLASTp percent identity floor (default: 80.0)."
    )
    g.add_argument("--blastp-min-qcov", type=float, default=80.0, help="BLASTp query coverage floor (default: 80.0).")
    g.add_argument("--blastp-evalue", type=float, default=1e-5, help="BLASTp e-value threshold (default: 1e-5).")

    g = p.add_argument_group("HH-suite")
    g.add_argument(
        "--skip-hhsuite",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip HH-suite step (overrides --tier default).",
    )
    g.add_argument(
        "--hhsuite-pfam-db", default="", help="Path to HH-suite Pfam database. Falls back to $SSIGN_HHSUITE_PFAM."
    )
    g.add_argument(
        "--hhsuite-pdb70-db", default="", help="Path to HH-suite PDB70 database. Falls back to $SSIGN_HHSUITE_PDB70."
    )
    g.add_argument(
        "--hhsuite-uniclust-db", default="", help="Path to UniClust DB. Falls back to $SSIGN_HHSUITE_UNICLUST."
    )
    g.add_argument(
        "--hhsuite-min-prob",
        type=float,
        help="HH-suite probability floor (default: ssign_lib.constants.HHSUITE_MIN_PROB).",
    )

    g = p.add_argument_group("InterProScan")
    g.add_argument(
        "--skip-interproscan",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip InterProScan (overrides --tier default).",
    )
    g.add_argument("--interproscan-db", default="", help="Path to InterProScan install dir.")
    g.add_argument(
        "--interproscan-min-evalue", type=float, default=1e-5, help="InterProScan e-value threshold (default: 1e-5)."
    )

    g = p.add_argument_group("pLM-BLAST")
    g.add_argument(
        "--skip-plmblast",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip pLM-BLAST step (overrides --tier default).",
    )
    g.add_argument("--plmblast-db", default="", help="Path to ECOD pLM-BLAST database (ECOD30 default).")
    g.add_argument(
        "--plmblast-cpc",
        type=int,
        default=90,
        help=(
            "pLM-BLAST cosine percentile cutoff (default: 90, the Kaminski 2023 "
            "paper setting). Drop to 70-80 for more permissive matching on "
            "short proteins (<200 aa, where cpc=90 often returns no hit) at the "
            "cost of longer search wallclock."
        ),
    )

    g = p.add_argument_group("EggNOG-mapper")
    g.add_argument(
        "--skip-eggnog",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip EggNOG-mapper step (overrides --tier default).",
    )
    g.add_argument("--eggnog-db", default="", help="Path to EggNOG database directory.")
    g.add_argument(
        "--eggnog-dbmem",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Pass --dbmem to emapper.py (default: enabled). Loads eggnog.db "
            "into RAM (~44 GB resident); required on NFS-backed cluster "
            "scratch where the SQLite mmap otherwise hangs. Disable only on "
            "RAM-constrained machines with the DB on local SSD."
        ),
    )

    g = p.add_argument_group("PLM-Effector")
    g.add_argument(
        "--skip-plm-effector",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "PLM-Effector is OFF by default (it over-predicts at genome scale; "
            "see docs). Opt in with --no-skip-plm-effector."
        ),
    )
    g.add_argument("--plm-effector-weights-dir", default="", help="Directory with PLM-Effector weights + ProtT5 cache.")
    g.add_argument(
        "--plm-effector-types",
        nargs="+",
        default=["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"],
        help="Secretion-system types to predict (default: T1SE T2SE T3SE T4SE T6SE).",
    )
    g.add_argument(
        "--plm-chunk-size",
        type=int,
        default=256,
        help=(
            "Proteins per PLM-Effector feature-extraction chunk (default: 256). "
            "Lower this if PLM-Effector OOM-kills on a large genome; raises only "
            "help when host RAM is abundant."
        ),
    )

    g = p.add_argument_group("misc annotation")
    g.add_argument(
        "--skip-protparam",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip ProtParam step (overrides --tier default).",
    )
    g.add_argument(
        "--t5ass-annotate-whole",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Run EggNOG / BLASTp / pLM-BLAST / HHsuite / ProtParam a "
            "SECOND time on T5aSS substrates with the full protein "
            "FASTA, emitting t5ass_whole_* columns alongside the "
            "default passenger-only annotations. Useful to compare "
            "functional (passenger) vs structural (whole-AT, "
            "β-barrel-dominated) annotation. IPS unchanged."
        ),
    )
    g.add_argument(
        "--filter-dse-type-mismatch",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Filter DSE-only substrates whose predicted SS type doesn't match the nearby MacSyFinder system.",
    )
    g.add_argument(
        "--ortholog-min-pident",
        type=float,
        default=40.0,
        help="Ortholog grouping percent-identity floor (default: 40.0).",
    )
    g.add_argument(
        "--ortholog-min-qcov", type=float, default=70.0, help="Ortholog grouping query-coverage floor (default: 70.0)."
    )

    # ── Figures ─────────────────────────────────────────────────────────
    g = p.add_argument_group("figures")
    g.add_argument("--dpi", type=int, default=300, help="Figure DPI (default: 300).")
    g.add_argument(
        "--fig-ss-comp",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="01 secreted proteins per genome, stacked by SS type.",
    )
    g.add_argument(
        "--fig-physicochemical",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="02 size & physicochemical properties by SS type (length + ProtParam when present).",
    )
    g.add_argument(
        "--fig-func-summary",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="03-06 functional categories by SS type (COG/KEGG/EggNOG/consensus).",
    )


def _config_from_args(
    args: argparse.Namespace,
    input_path: str,
    sample_id: str,
    outdir: str,
) -> "PipelineConfig":
    """Map argparse Namespace to a PipelineConfig for one genome.

    ``input_path``, ``sample_id``, and ``outdir`` are passed in explicitly so
    the same args Namespace can build N configs (one per genome) in a
    multi-genome run.
    """
    from ssign_app.core.runner import PipelineConfig

    cfg_kwargs = {
        "input_path": input_path,
        "original_filename": args.original_filename,
        "sample_id": sample_id,
        "outdir": outdir,
        "tier": args.tier,
        "wholeness_threshold": args.wholeness_threshold,
        "excluded_systems": list(args.excluded_systems),
        "macsyfinder_db_type": args.macsyfinder_db_type,
        "cpu_per_genome": args.cpu_per_genome,
        "conf_threshold": args.conf_threshold,
        "proximity_window": args.proximity_window,
        "required_fraction_correct": args.required_fraction_correct,
        "dlp_confidence_threshold": args.dlp_confidence_threshold,
        "skip_localization_gate": args.skip_localization_gate,
        "use_input_annotations": args.use_input_annotations,
        "run_bakta": args.run_bakta,
        "bakta_db": args.bakta_db,
        "bakta_threads": args.bakta_threads,
        "deeplocpro_mode": args.deeplocpro_mode,
        "deeplocpro_path": args.deeplocpro_path,
        "signalp_mode": args.signalp_mode,
        "signalp_path": args.signalp_path,
        "skip_deeplocpro": args.skip_deeplocpro,
        "skip_signalp": args.skip_signalp,
        "skip_deepsece": args.skip_deepsece,
        "dlp_whole_genome": args.dlp_whole_genome,
        "dse_whole_genome": args.dse_whole_genome,
        "sp_whole_genome": args.sp_whole_genome,
        "plme_whole_genome": args.plme_whole_genome,
        "monitor_resources": args.monitor_resources,
        "monitor_interval_s": args.monitor_interval_s,
        "enrichment_stats": args.enrichment_stats,
        "skip_blastp": args.skip_blastp,
        "blastp_db": args.blastp_db,
        "blastp_exclude_taxid": args.blastp_exclude_taxid,
        "blastp_min_pident": args.blastp_min_pident,
        "blastp_min_qcov": args.blastp_min_qcov,
        "blastp_evalue": args.blastp_evalue,
        "skip_hhsuite": args.skip_hhsuite,
        "hhsuite_pfam_db": args.hhsuite_pfam_db,
        "hhsuite_pdb70_db": args.hhsuite_pdb70_db,
        "hhsuite_uniclust_db": args.hhsuite_uniclust_db,
        "skip_interproscan": args.skip_interproscan,
        "interproscan_db": args.interproscan_db,
        "interproscan_min_evalue": args.interproscan_min_evalue,
        "skip_plmblast": args.skip_plmblast,
        "plmblast_db": args.plmblast_db,
        "plmblast_cpc": args.plmblast_cpc,
        "skip_eggnog": args.skip_eggnog,
        "eggnog_db": args.eggnog_db,
        "eggnog_dbmem": args.eggnog_dbmem,
        "skip_plm_effector": args.skip_plm_effector,
        "plm_effector_weights_dir": args.plm_effector_weights_dir,
        "plm_effector_types": list(args.plm_effector_types),
        "plm_chunk_size": args.plm_chunk_size,
        "skip_protparam": args.skip_protparam,
        "t5ass_annotate_whole": args.t5ass_annotate_whole,
        "filter_dse_type_mismatch": args.filter_dse_type_mismatch,
        "deepsece_min_prob": args.deepsece_min_prob,
        "signalp_min_prob": args.signalp_min_prob,
        "ortholog_min_pident": args.ortholog_min_pident,
        "ortholog_min_qcov": args.ortholog_min_qcov,
        "dpi": args.dpi,
        "fig_ss_comp": args.fig_ss_comp,
        "fig_physicochemical": args.fig_physicochemical,
        "fig_func_summary": args.fig_func_summary,
    }
    # hhsuite_min_prob is the only field with a non-trivial default
    # (constants.HHSUITE_MIN_PROB). argparse leaves it None when absent;
    # only override the dataclass default if the user supplied it explicitly.
    if args.hhsuite_min_prob is not None:
        cfg_kwargs["hhsuite_min_prob"] = args.hhsuite_min_prob

    return PipelineConfig(**cfg_kwargs)


def _run_pipeline(args: argparse.Namespace) -> int:
    """Execute the `ssign run` subcommand. Returns the process exit code."""
    from ssign_app.core.runner import PipelineRunner

    inputs: list[str] = list(args.input_path)
    for p in inputs:
        if not os.path.exists(p):
            print(f"Error: input file not found: {p}", file=sys.stderr)
            return 2

    def _terminal_progress(step: str, pct: int, msg: str) -> None:
        print(f"  [{pct:3d}%] {step} — {msg}", flush=True)

    if len(inputs) == 1:
        input_path = inputs[0]
        sample_id = args.sample_id or os.path.splitext(os.path.basename(input_path))[0]
        config = _config_from_args(args, input_path, sample_id, args.outdir)

        runner = PipelineRunner(config, progress_callback=_terminal_progress)
        print(f"ssign — running on {config.input_path}", flush=True)
        print(f"   outdir: {config.outdir}", flush=True)
        print(f"   sample_id: {config.sample_id}", flush=True)
        print(flush=True)

        try:
            results = runner.run(resume=args.resume)
        except KeyboardInterrupt:
            print("\nInterrupted.", file=sys.stderr)
            return 130

        return _report_single_genome(results)

    # Multi-genome path
    if args.sample_id:
        print(
            "Error: --sample-id is only valid for single-genome runs; "
            "per-genome sample_ids are derived from input filenames in "
            "multi-genome runs.",
            file=sys.stderr,
        )
        return 2

    from ssign_app.core.multi_runner import MultiGenomeRunner

    top_outdir = args.outdir
    configs = []
    seen_sids: set[str] = set()
    for input_path in inputs:
        sid = os.path.splitext(os.path.basename(input_path))[0]
        if sid in seen_sids:
            print(
                f"Error: duplicate sample_id {sid!r} derived from input "
                f"filenames; rename inputs so their basenames are distinct.",
                file=sys.stderr,
            )
            return 2
        seen_sids.add(sid)
        per_genome_outdir = os.path.join(top_outdir, sid)
        configs.append(_config_from_args(args, input_path, sid, per_genome_outdir))

    runner = MultiGenomeRunner(
        configs,
        progress_callback=_terminal_progress,
        write_combined_summary=args.combined_summary,
    )
    print(f"ssign — running on {len(inputs)} genome(s) (batched)", flush=True)
    print(f"   outdir: {top_outdir}", flush=True)
    print(f"   sample_ids: {', '.join(c.sample_id for c in configs)}", flush=True)
    print(flush=True)

    try:
        results_by_sid = runner.run(resume=args.resume)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130

    return _report_multi_genome(results_by_sid)


def _report_single_genome(results) -> int:
    n_success = sum(1 for r in results if r.success)
    n_total = len(results)
    print(flush=True)
    if n_success == n_total:
        print(f"Pipeline complete: {n_success}/{n_total} steps succeeded.", flush=True)
        return 0
    print(f"Pipeline finished with issues: {n_success}/{n_total} steps succeeded.", file=sys.stderr)
    for r in results:
        if not r.success:
            print(f"  - FAILED: {r.name} — {r.message[:200]}", file=sys.stderr)
    return 1


def _report_multi_genome(results_by_sid: dict) -> int:
    print(flush=True)
    any_failed = False
    for sid, results in results_by_sid.items():
        n_success = sum(1 for r in results if r.success)
        n_total = len(results)
        if n_success == n_total:
            print(f"  {sid}: {n_success}/{n_total} steps succeeded", flush=True)
        else:
            any_failed = True
            print(f"  {sid}: {n_success}/{n_total} steps succeeded (FAILED)", file=sys.stderr)
            for r in results:
                if not r.success:
                    print(f"      - {r.name} — {r.message[:200]}", file=sys.stderr)
    return 1 if any_failed else 0


# ---------------------------------------------------------------------------
# `ssign doctor` subcommand
# ---------------------------------------------------------------------------


def _add_doctor_parser(subparsers: argparse._SubParsersAction) -> None:
    """Build the `ssign doctor` subcommand parser.

    Implementation lives in ``ssign_app.scripts.doctor``; this stub just
    exposes the flags the runtime function consumes. Defaults are imported
    from there to avoid drift.
    """
    from ssign_app.scripts.doctor import DEFAULT_DATA_ROOT, DEFAULT_TIER

    p = subparsers.add_parser(
        "doctor",
        help="Verify the install: Python packages, external binaries, databases, model weights.",
        description=(
            "Check every dependency ssign needs and report what's missing with the exact "
            "fix command. Exit non-zero on any failure so you can chain `ssign doctor && "
            "ssign run ...` in scripts."
        ),
    )
    p.add_argument(
        "--tier",
        choices=("base", "extended", "full"),
        default=DEFAULT_TIER,
        help=f"Install tier to verify against (default: {DEFAULT_TIER}).",
    )
    p.add_argument(
        "--imports-only",
        action="store_true",
        help="Only check Python imports; skip binaries / DBs / weights (used by CI).",
    )
    p.add_argument(
        "--data-root",
        default=DEFAULT_DATA_ROOT,
        help=f"Root for databases + models (default: {DEFAULT_DATA_ROOT}). SSIGN_* env vars override per-DB paths.",
    )


# ---------------------------------------------------------------------------
# `ssign` (no subcommand) — Streamlit GUI launcher
# ---------------------------------------------------------------------------


def _launch_gui(args: argparse.Namespace) -> int:
    """Launch the Streamlit GUI. Preserves the historical `ssign` UX:
    free-port detection, banner, stderr filtering, and pass-through of
    [ssign] log lines on stdout."""
    print(BANNER, flush=True)

    def _port_free(p: int) -> bool:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("localhost", p))
                return True
            except OSError:
                return False

    port = args.port
    if not _port_free(port):
        for candidate in range(port + 1, port + 50):
            if _port_free(candidate):
                print(f"  Port {port} is in use, using {candidate} instead.", flush=True)
                port = candidate
                break
        else:
            print(f"Error: No free port found in range {port}-{port + 49}", file=sys.stderr)
            return 1

    app_dir = os.path.dirname(os.path.abspath(__file__))
    app_file = os.path.join(app_dir, "Home.py")
    config_dir = os.path.join(app_dir, ".streamlit")

    if not os.path.exists(app_file):
        print(f"Error: Could not find {app_file}", file=sys.stderr)
        return 1

    env = os.environ.copy()
    if os.path.isdir(config_dir):
        env["STREAMLIT_CONFIG_DIR"] = config_dir
    env["STREAMLIT_SERVER_MAX_UPLOAD_SIZE"] = "500"
    env["STREAMLIT_SERVER_MAX_MESSAGE_SIZE"] = "500"
    env["STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION"] = "false"
    env["STREAMLIT_SERVER_ENABLE_CORS"] = "false"
    env["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"

    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        app_file,
        "--server.port",
        str(port),
        "--server.headless",
        "true" if args.no_browser else "false",
        "--server.maxUploadSize",
        "500",
        "--server.maxMessageSize",
        "500",
        "--server.enableXsrfProtection",
        "false",
        "--server.enableCORS",
        "false",
    ]

    url = f"http://localhost:{port}"
    if args.no_browser:
        print(f"  Open in browser: {url}", flush=True)
    else:
        print("  Opening... If nothing automatically opens, try pasting", flush=True)
        print(f"  this into your browser: {url}", flush=True)
    print(flush=True)

    try:
        _banner_re = re.compile(r"You can now view|Network URL|External URL|Local URL")
        _ssign_line_re = re.compile(r"^\[ssign\]")

        proc = subprocess.Popen(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        def _filter_stderr() -> None:
            for line in proc.stderr:
                if _banner_re.search(line) or not line.strip():
                    continue
                sys.stderr.write(line)
                sys.stderr.flush()

        def _filter_stdout() -> None:
            for line in proc.stdout:
                if _ssign_line_re.match(line):
                    sys.stdout.write(f"  {line}")
                    sys.stdout.flush()

        t1 = threading.Thread(target=_filter_stderr, daemon=True)
        t2 = threading.Thread(target=_filter_stdout, daemon=True)
        t1.start()
        t2.start()
        proc.wait()
        t1.join(timeout=1)
        t2.join(timeout=1)
        return proc.returncode
    except KeyboardInterrupt:
        print("\nssign stopped.")
        return 130
    except FileNotFoundError:
        print("Error: Streamlit not found. Install with: pip install ssign", file=sys.stderr)
        return 1


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="ssign",
        description="Secretion-system Identification for Gram Negatives. "
        "Run with no subcommand to launch the GUI; use `ssign run` "
        "for non-interactive / HPC use.",
    )
    parser.add_argument("--version", action="store_true", help="Print version and exit.")
    parser.add_argument(
        "--no-browser", action="store_true", help="(GUI mode) Start the GUI server without opening a browser."
    )
    parser.add_argument(
        "--port", type=int, default=8501, help="(GUI mode) Port for the Streamlit server (default: 8501)."
    )

    subparsers = parser.add_subparsers(dest="subcommand")
    _add_run_parser(subparsers)
    _add_doctor_parser(subparsers)

    args = parser.parse_args()

    if args.version:
        from ssign_app import __version__

        print(f"ssign {__version__}")
        return 0

    if args.subcommand == "run":
        return _run_pipeline(args)

    if args.subcommand == "doctor":
        from ssign_app.scripts.doctor import run as doctor_run

        return doctor_run(
            tier=args.tier,
            imports_only=args.imports_only,
            data_root=args.data_root,
        )

    return _launch_gui(args)


if __name__ == "__main__":
    sys.exit(main())
