#!/usr/bin/env python3
"""Manual end-to-end pipeline runner for the T1SS fixture genome.

NOT a pytest test — this is a hand-run script that invokes the full
PipelineRunner on the bundled T1SS fixture and prints a human-readable
summary. The authoritative pytest integration test lives in
`tests/integration/test_pipeline_fixture.py`.

What runs:
- Input format detection + protein extraction (always runs)
- MacSyFinder secretion-system detection (requires MacSyFinder + TXSScan)
- DeepLocPro localisation prediction (currently runs via DTU's webface;
  will move to local install in Phase 2.3 once DTU redistribution is
  cleared)
- Proximity analysis + cross-validation + T5SS handling + reporting

Heavy annotation tools (BLASTp, HH-suite, InterProScan, Bakta, EggNOG,
PLM-Effector, pLM-BLAST) are **skipped** by default here — they all need
local databases and, in a few cases, a GPU. Turn any of them on by
providing the database path / weights dir plus `skip_*=False`.

Remote-mode tools that were removed in Phase 2:
- BLASTp NCBI remote fallback  (2.1b) — require --blastp-db
- HH-suite MPI Toolkit remote  (2.1c) — require local hh-suite install
- InterProScan EBI web service (2.1c) — require local InterProScan
- Foldseek scaffolding         (2.1a) — dropped entirely

Usage:
    python tests/test_pipeline_e2e.py [output_dir]

Output dir defaults to /tmp/ssign_e2e_test. Pipeline supports resume via
the persisted progress JSON.
"""

import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from ssign_app.core.runner import PipelineConfig, PipelineRunner

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> int:
    fixture = os.path.join(
        os.path.dirname(__file__),
        "fixtures",
        "Xanthobacter_tagetidis_ATCC_700314_contig_87.gbff",
    )
    if not os.path.exists(fixture):
        print(f"ERROR: T1SS fixture not found: {fixture}", file=sys.stderr)
        print(
            "Run this script from the repo root so the fixture path resolves.",
            file=sys.stderr,
        )
        return 1

    outdir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ssign_e2e_test"
    os.makedirs(outdir, exist_ok=True)

    config = PipelineConfig(
        input_path=fixture,
        sample_id="t1ss_fixture",
        outdir=outdir,
        # Core settings
        wholeness_threshold=0.8,
        excluded_systems=["Flagellum", "Tad", "T3SS"],
        conf_threshold=0.8,
        proximity_window=3,
        # Prediction tools
        deeplocpro_mode="remote",  # DTU webface; replaced by local in Phase 2.3
        skip_deepsece=True,  # ESM model is too large for a laptop run
        skip_signalp=True,  # Optional; skip by default
        # Annotation — all skipped; enable individually by providing a DB path
        skip_blastp=True,
        skip_hhsuite=True,
        skip_interproscan=True,
        skip_plmblast=True,
        skip_protparam=False,  # Fast and local, always run
        skip_structure=True,  # Foldseek was dropped in Phase 2
    )

    def progress_cb(step, pct, msg):
        print(f"[{pct:3d}%] {step}: {msg}")

    runner = PipelineRunner(config, progress_callback=progress_cb)
    results = runner.run(resume=True)

    print("\n" + "=" * 60)
    print("PIPELINE RESULTS")
    print("=" * 60)
    passed = failed = 0
    for r in results:
        status = "OK" if r.success else "FAIL"
        passed += int(r.success)
        failed += int(not r.success)
        print(f"  [{status:4s}] {r.name}: {r.message}")

    print(f"\n{passed} passed, {failed} failed out of {len(results)} steps")

    print("\n--- Output files ---")
    for entry in sorted(os.listdir(outdir)):
        path = os.path.join(outdir, entry)
        if os.path.isfile(path):
            print(f"  {entry}: {os.path.getsize(path):,} bytes")
        elif os.path.isdir(path):
            print(f"  {entry}/: {len(os.listdir(path))} files")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
