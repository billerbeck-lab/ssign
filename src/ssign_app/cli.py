#!/usr/bin/env python3
"""CLI entry point for ssign.

Usage:
    ssign                           # print a usage hint and exit
    ssign --version                 # print version
    ssign run input.gbff --outdir results [flags...]
                                    # run pipeline (HPC, scripting, batch use)
    ssign run --help                # full list of run-mode flags
    ssign doctor                    # verify the install
    ssign fetch-databases --tier T  # download reference databases

The `run` subcommand builds a PipelineConfig from CLI flags and drives
PipelineRunner directly. It exposes every PipelineConfig field as a flag,
grouped by phase. See `ssign run --help` for the full list.
"""

import argparse
import os
import subprocess
import sys
from typing import TYPE_CHECKING, Optional

from ssign_app.scripts.ssign_lib.constants import DEFAULT_EXCLUDED_SYSTEMS  # noqa: F401  (argparse default below)

if TYPE_CHECKING:
    from ssign_app.core.runner import PipelineConfig

GITHUB_URL = "https://github.com/billerbeck-lab/ssign"

BANNER = r"""
  _|_|_|       _|_|_|      _|        _|_|_|     _|_|_
_|_|         _|_|          _|      _|    _|     _|    _|
  _|           _|          _|      _|    _|     _|    _|
    _|_|         _|_|      _|      _|    _|     _|    _|
_|_|_|       _|_|_|        _|        _|_|_|     _|    _|
                                         _|
                                     _|_|
           .
       `    .                   .    ::         ::
     ,    .     *               .:   .:....::::=.
   .    *     ,    .     :.      ::.:-:---=--==:..
  `   .    *     .   *    -: ....:----:.      :+-.:
       .     .    .-.     .:.:-==-:.      ..   .+-.:....
         ,     .    -:....:-=-:::..   ..:...    -+ -
                      .:-=-:        ....::     :=:
               ..:..:-=-.      .::.::.:   .:::==: : ..:
                 : --:      ...:.:.     ::-=-::.-:
   ..            - +-    ... .      .-====:..
  =---           .:-=.   ...    ..-+*=-.....::
  -: .=         .--:-=:      .:-===-...-.    :-
   +  .-       .:. ::-=---:---=-:: .   .-.
   -:  :-     ..    ..::-::-*-   =:      .
    +   --   :.       .-    .+    -:
    -:   :-::.        .-
    ::     .
  .:-
 :
"""


def _print_banner(show_version: bool = True) -> None:
    """Print the ssign wordmark art at startup.

    Purely decorative, so it must never be fatal: if stdout is a redirected or
    non-UTF-8 stream, any failure is swallowed (skip the banner) rather than
    taking down a real run.
    """
    try:
        print(BANNER, flush=True)
        if show_version:
            from ssign_app import __version__

            print(f"  v{__version__}\n", flush=True)
    except Exception:
        pass


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
            "EggNOG DB load, pLM-BLAST embeddings) are paid once per batch "
            "rather than once per genome."
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
        "--scratch-dir",
        default="",
        help=(
            "Directory for scratch/temp files (tool working dirs). Default ('') "
            "auto-resolves: keep $TMPDIR if it has adequate free space, else fall "
            "back to a dir under --outdir. Set this when running in a container "
            "whose /tmp is a small tmpfs (avoids Bakta 'No space left on device')."
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
        help="DeepLocPro execution mode. Default: auto, local if 'deeplocpro' is "
        "on PATH or at --deeplocpro-path / $SSIGN_DEEPLOCPRO_PATH. If no local "
        "install is found, ssign STOPS with install instructions (it does not "
        "auto-upload to the DTU webserver); pass 'remote' to opt into the webserver.",
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
        help="SignalP execution mode. Default: auto, local if 'signalp6' is "
        "on PATH or at --signalp-path / $SSIGN_SIGNALP_PATH. If no local install "
        "is found, ssign STOPS with install instructions (it does not auto-upload "
        "to the DTU webserver); pass 'remote' to opt into the webserver.",
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
    g = p.add_argument_group("annotation (all tools)")
    g.add_argument(
        "--skip-annotation",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Skip every annotation tool at once (BLASTp, HH-suite, "
            "InterProScan, pLM-BLAST, EggNOG, ProtParam). Predictions and "
            "substrate calls still run; only the functional-annotation phase "
            "is dropped. A per-tool --no-skip-<tool> (e.g. --no-skip-eggnog) "
            "still overrides this to keep that one tool on."
        ),
    )

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
        default=None,
        help=(
            "Pass --dbmem to emapper.py (loads eggnog.db into RAM, ~44 GB "
            "resident). Default: auto — enabled only when the job's RAM share "
            "is >= 50 GB, else the on-disk SQLite is memory-mapped (the runner "
            "stages it to local scratch, so no NFS mmap stall). Force with "
            "--eggnog-dbmem / --no-eggnog-dbmem."
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

    # `--skip-annotation` is a master switch: force all six annotation-tool
    # skips on, except any the user explicitly re-enabled with `--no-skip-<tool>`
    # (that leaves the flag False, not None, so it is preserved here). A flag
    # left at None falls through to the tier default when skip-annotation is off.
    ann_skip = {
        name: getattr(args, name)
        for name in (
            "skip_blastp",
            "skip_hhsuite",
            "skip_interproscan",
            "skip_plmblast",
            "skip_eggnog",
            "skip_protparam",
        )
    }
    if args.skip_annotation:
        ann_skip = {name: (True if val is None else val) for name, val in ann_skip.items()}

    cfg_kwargs = {
        "input_path": input_path,
        "original_filename": args.original_filename,
        "sample_id": sample_id,
        "outdir": outdir,
        "tier": args.tier,
        "scratch_dir": args.scratch_dir,
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
        "monitor_resources": args.monitor_resources,
        "monitor_interval_s": args.monitor_interval_s,
        "enrichment_stats": args.enrichment_stats,
        "skip_blastp": ann_skip["skip_blastp"],
        "blastp_db": args.blastp_db,
        "blastp_exclude_taxid": args.blastp_exclude_taxid,
        "blastp_min_pident": args.blastp_min_pident,
        "blastp_min_qcov": args.blastp_min_qcov,
        "blastp_evalue": args.blastp_evalue,
        "skip_hhsuite": ann_skip["skip_hhsuite"],
        "hhsuite_pfam_db": args.hhsuite_pfam_db,
        "hhsuite_pdb70_db": args.hhsuite_pdb70_db,
        "hhsuite_uniclust_db": args.hhsuite_uniclust_db,
        "skip_interproscan": ann_skip["skip_interproscan"],
        "interproscan_db": args.interproscan_db,
        "interproscan_min_evalue": args.interproscan_min_evalue,
        "skip_plmblast": ann_skip["skip_plmblast"],
        "plmblast_db": args.plmblast_db,
        "plmblast_cpc": args.plmblast_cpc,
        "skip_eggnog": ann_skip["skip_eggnog"],
        "eggnog_db": args.eggnog_db,
        "eggnog_dbmem": args.eggnog_dbmem,
        "skip_protparam": ann_skip["skip_protparam"],
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

    _print_banner()

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

        rc = _report_single_genome(results)
        if rc == 0:
            _print_outputs_hint(config.outdir, config.sample_id)
        return rc

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

    rc = _report_multi_genome(results_by_sid)
    if rc == 0:
        _print_outputs_hint(top_outdir)
    return rc


def _print_outputs_hint(outdir: str, sample_id: Optional[str] = None) -> None:
    """Point the user at the key output files after a successful run."""
    outdir = os.path.abspath(outdir)
    print(f"\nResults written to {outdir}", flush=True)
    if sample_id:
        print(f"  {sample_id}_results.csv   secreted proteins + secretion systems (3 sections)", flush=True)
        print(f"  {sample_id}_summary.txt   readable summary + enrichment stats", flush=True)
        print(f"  figures/{sample_id}/      figures 01-06", flush=True)
    else:
        print("  combined_results.csv     merged results table across genomes", flush=True)
        print("  combined_summary.txt     aggregated summary across all genomes", flush=True)
        print("  <genome>/                per-genome results.csv, summary.txt, figures/", flush=True)


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
# `ssign fetch-databases`: download reference DBs via the bundled script
# ---------------------------------------------------------------------------


def _add_fetch_databases_parser(subparsers) -> None:
    p = subparsers.add_parser(
        "fetch-databases",
        help="Download the reference databases for a tier (wraps scripts/fetch_databases.sh).",
        description=(
            "Fetch ssign's reference databases (pinned versions, resumable download). This "
            "wraps the bundled fetch_databases.sh, so from the container you need no host "
            "tools: `apptainer run ssign.sif fetch-databases --tier extended --target /data/db` "
            "(bind /data/db). The DTU-licensed SignalP tool is separate; see ssign-setup-dtu."
        ),
    )
    p.add_argument(
        "--tier",
        choices=("base", "extended", "full"),
        required=True,
        help="Which tier's databases to download (base ~4 GB, extended ~100 GB, full ~500 GB).",
    )
    p.add_argument(
        "--target",
        default="",
        help="Destination directory (default: the script's own ~/.ssign/databases).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be downloaded without downloading anything.",
    )


def _find_fetch_script() -> "str | None":
    """Locate the bundled fetch_databases.sh: the container bake path first, then
    the repo root relative to this file (works for an editable install / a clone).
    A wheel-only install without the repo won't find it; the container is the
    supported path for that case."""
    here = os.path.dirname(os.path.abspath(__file__))  # src/ssign_app
    repo_root = os.path.dirname(os.path.dirname(here))  # repo root
    for cand in (
        "/opt/ssign/scripts/fetch_databases.sh",
        os.path.join(repo_root, "scripts", "fetch_databases.sh"),
    ):
        if os.path.isfile(cand):
            return cand
    return None


def _fetch_databases(args: argparse.Namespace) -> int:
    script = _find_fetch_script()
    if not script:
        print(
            "fetch_databases.sh not found. It ships with the container image "
            "(apptainer run ssign.sif fetch-databases ...) and the source repo. "
            "For a wheel-only install without the repo, clone it and run "
            "scripts/fetch_databases.sh directly.",
            file=sys.stderr,
        )
        return 1
    cmd = ["bash", script, "--tier", args.tier]
    if args.target:
        cmd += ["--target", args.target]
    if args.dry_run:
        cmd += ["--dry-run"]
    return subprocess.call(cmd)


# ---------------------------------------------------------------------------
# `ssign` (no subcommand) — usage hint
# ---------------------------------------------------------------------------


def _print_usage_hint() -> int:
    """No subcommand given: ssign is a command-line tool, so show how to run it
    and point at `ssign --help` for the full list."""
    _print_banner()
    print(
        "ssign is a command-line tool. Common commands:\n"
        "\n"
        "  ssign run <genome.gbff> --outdir <out>   run the pipeline\n"
        "  ssign doctor                             verify the install\n"
        "  ssign fetch-databases --tier <tier>      download reference databases\n"
        "\n"
        "  ssign --help        full usage\n"
        "  ssign run --help    every run-mode flag\n"
        "\n"
        f"Docs and source: {GITHUB_URL}",
        flush=True,
    )
    return 0


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="ssign",
        description="Secretion-system Identification for Gram Negatives. "
        "Use `ssign run` to run the pipeline; `ssign doctor` to check the "
        "install; `ssign fetch-databases` to download reference databases.",
        epilog=f"Docs and source: {GITHUB_URL}",
    )
    parser.add_argument("--version", action="store_true", help="Print version and exit.")

    subparsers = parser.add_subparsers(dest="subcommand")
    _add_run_parser(subparsers)
    _add_doctor_parser(subparsers)
    _add_fetch_databases_parser(subparsers)

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

    if args.subcommand == "fetch-databases":
        return _fetch_databases(args)

    return _print_usage_hint()


if __name__ == "__main__":
    sys.exit(main())
