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
            "Run ssign on one input genome (GenBank, GFF3, or FASTA). All "
            "PipelineConfig fields are exposed as flags. For batch / "
            "multi-genome runs, drive ssign from a shell loop or HPC array."
        ),
    )

    # ── Essentials ──────────────────────────────────────────────────────
    g = p.add_argument_group("essentials")
    g.add_argument(
        "input_path",
        help="Path to the input genome (GenBank .gbff/.gbk, GFF3 .gff, or FASTA).",
    )
    g.add_argument(
        "--outdir",
        default="./results",
        help="Output directory (default: ./results).",
    )
    g.add_argument(
        "--sample-id",
        default="",
        help="Sample identifier used to prefix output files. Defaults to the input filename's stem.",
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

    # ── Phase 2: SS detection ──────────────────────────────────────────
    g = p.add_argument_group("SS detection (MacSyFinder)")
    g.add_argument(
        "--wholeness-threshold", type=float, default=0.8, help="Minimum MacSyFinder system completeness (default: 0.8)."
    )
    g.add_argument(
        "--excluded-systems",
        nargs="+",
        default=["Flagellum", "Tad", "T3SS"],
        help="Secretion-system models to exclude (default: Flagellum Tad T3SS).",
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
        default=os.cpu_count() or 4,
        help="CPUs available to per-genome subtools (default: all).",
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
        "--deepsece-min-prob", type=float, default=0.8, help="DeepSecE min probability to call secreted (default: 0.8)."
    )
    g.add_argument(
        "--signalp-min-prob",
        type=float,
        default=0.5,
        help="SignalP min probability for a signal peptide (default: 0.5).",
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
        default=False,
        help="Run Bakta on FASTA input or to re-annotate GenBank.",
    )
    g.add_argument("--bakta-db", default="", help="Path to Bakta database (required when --run-bakta).")
    g.add_argument("--bakta-threads", type=int, default=4, help="Threads passed to Bakta (default: 4).")

    # ── DTU prediction tools (DeepLocPro + SignalP) ─────────────────────
    g = p.add_argument_group("DTU prediction tools")
    g.add_argument(
        "--deeplocpro-mode",
        choices=["local", "remote"],
        default="remote",
        help="DeepLocPro execution mode (default: remote = DTU web API).",
    )
    g.add_argument(
        "--deeplocpro-path",
        default="",
        help="Path to local DeepLocPro install (required when --deeplocpro-mode local).",
    )
    g.add_argument(
        "--signalp-mode",
        choices=["local", "remote"],
        default="remote",
        help="SignalP execution mode (default: remote = DTU web API).",
    )
    g.add_argument(
        "--signalp-path", default="", help="Path to local SignalP 6 install (required when --signalp-mode local)."
    )
    g.add_argument("--skip-signalp", action=argparse.BooleanOptionalAction, default=False, help="Skip SignalP step.")
    g.add_argument("--skip-deepsece", action=argparse.BooleanOptionalAction, default=False, help="Skip DeepSecE step.")
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

    # ── Phase 5: Annotation tools ───────────────────────────────────────
    g = p.add_argument_group("BLASTp")
    g.add_argument("--skip-blastp", action=argparse.BooleanOptionalAction, default=False)
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
        default=True,
        help="Skip HH-suite step (default: skip).",
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
    g.add_argument("--skip-interproscan", action=argparse.BooleanOptionalAction, default=False)
    g.add_argument("--interproscan-db", default="", help="Path to InterProScan install dir.")
    g.add_argument(
        "--interproscan-min-evalue", type=float, default=1e-5, help="InterProScan e-value threshold (default: 1e-5)."
    )

    g = p.add_argument_group("pLM-BLAST")
    g.add_argument(
        "--skip-plmblast",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip pLM-BLAST step (default: skip).",
    )
    g.add_argument("--plmblast-db", default="", help="Path to ECOD70 pLM-BLAST database.")

    g = p.add_argument_group("EggNOG-mapper")
    g.add_argument(
        "--skip-eggnog",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip EggNOG-mapper step (default: skip).",
    )
    g.add_argument("--eggnog-db", default="", help="Path to EggNOG database directory.")

    g = p.add_argument_group("PLM-Effector")
    g.add_argument(
        "--skip-plm-effector",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip PLM-Effector step (default: skip).",
    )
    g.add_argument("--plm-effector-weights-dir", default="", help="Directory with PLM-Effector weights + ProtT5 cache.")
    g.add_argument(
        "--plm-effector-types",
        nargs="+",
        default=["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"],
        help="Secretion-system types to predict (default: T1SE T2SE T3SE T4SE T6SE).",
    )

    g = p.add_argument_group("misc annotation")
    g.add_argument("--skip-protparam", action=argparse.BooleanOptionalAction, default=False)
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
        "--fig-category", action=argparse.BooleanOptionalAction, default=True, help="Render functional-category figure."
    )
    g.add_argument(
        "--fig-ss-comp",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render SS-component composition figure.",
    )
    g.add_argument(
        "--fig-tool-heatmap", action=argparse.BooleanOptionalAction, default=True, help="Render tool-coverage heatmap."
    )
    g.add_argument(
        "--fig-substrate-count",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render per-SS substrate count figure.",
    )
    g.add_argument(
        "--fig-func-summary",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render functional-summary figure.",
    )


def _config_from_args(args: argparse.Namespace) -> "PipelineConfig":
    """Map argparse Namespace to a PipelineConfig.

    Defaults are applied by argparse from PipelineConfig's documented defaults
    above; we just translate the argparse names back to dataclass field names
    (they match modulo `-` → `_`).
    """
    from ssign_app.core.runner import PipelineConfig

    # sample_id default: derive from input filename if not supplied
    sample_id = args.sample_id or os.path.splitext(os.path.basename(args.input_path))[0]

    cfg_kwargs = {
        "input_path": args.input_path,
        "original_filename": args.original_filename,
        "sample_id": sample_id,
        "outdir": args.outdir,
        "wholeness_threshold": args.wholeness_threshold,
        "excluded_systems": list(args.excluded_systems),
        "macsyfinder_db_type": args.macsyfinder_db_type,
        "cpu_per_genome": args.cpu_per_genome,
        "conf_threshold": args.conf_threshold,
        "proximity_window": args.proximity_window,
        "required_fraction_correct": args.required_fraction_correct,
        "use_input_annotations": args.use_input_annotations,
        "run_bakta": args.run_bakta,
        "bakta_db": args.bakta_db,
        "bakta_threads": args.bakta_threads,
        "deeplocpro_mode": args.deeplocpro_mode,
        "deeplocpro_path": args.deeplocpro_path,
        "signalp_mode": args.signalp_mode,
        "signalp_path": args.signalp_path,
        "skip_signalp": args.skip_signalp,
        "skip_deepsece": args.skip_deepsece,
        "dlp_whole_genome": args.dlp_whole_genome,
        "dse_whole_genome": args.dse_whole_genome,
        "sp_whole_genome": args.sp_whole_genome,
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
        "skip_eggnog": args.skip_eggnog,
        "eggnog_db": args.eggnog_db,
        "skip_plm_effector": args.skip_plm_effector,
        "plm_effector_weights_dir": args.plm_effector_weights_dir,
        "plm_effector_types": list(args.plm_effector_types),
        "skip_protparam": args.skip_protparam,
        "filter_dse_type_mismatch": args.filter_dse_type_mismatch,
        "deepsece_min_prob": args.deepsece_min_prob,
        "signalp_min_prob": args.signalp_min_prob,
        "ortholog_min_pident": args.ortholog_min_pident,
        "ortholog_min_qcov": args.ortholog_min_qcov,
        "dpi": args.dpi,
        "fig_category": args.fig_category,
        "fig_ss_comp": args.fig_ss_comp,
        "fig_tool_heatmap": args.fig_tool_heatmap,
        "fig_substrate_count": args.fig_substrate_count,
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

    if not os.path.exists(args.input_path):
        print(f"Error: input file not found: {args.input_path}", file=sys.stderr)
        return 2

    config = _config_from_args(args)

    def _terminal_progress(step: str, pct: int, msg: str) -> None:
        print(f"  [{pct:3d}%] {step} — {msg}", flush=True)

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

    args = parser.parse_args()

    if args.version:
        from ssign_app import __version__

        print(f"ssign {__version__}")
        return 0

    if args.subcommand == "run":
        return _run_pipeline(args)

    return _launch_gui(args)


if __name__ == "__main__":
    sys.exit(main())
