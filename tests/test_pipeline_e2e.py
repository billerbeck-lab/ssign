#!/usr/bin/env python3
"""End-to-end pipeline test using Xanthobacter genome.

Runs the full pipeline with:
- DeepSecE SKIPPED (requires 7.3GB ESM model)
- DeepLocPro in REMOTE mode (DTU web server)
- BLASTp in REMOTE mode (NCBI web)
- InterProScan in REMOTE mode (EBI web)
- ProtParam computed locally
- HH-suite SKIPPED (slow, optional)
- SignalP SKIPPED
"""
import logging
import os
import sys
import tempfile

# Add package to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ssign_app.core.runner import PipelineConfig, PipelineRunner

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def main():
    # Input: Xanthobacter genome (GenBank)
    test_data = os.path.join(os.path.dirname(__file__), 'data',
                             'Xanthobacter_tagetidis_TagT2C_genomic.gbff')
    if not os.path.exists(test_data):
        print(f"ERROR: Test genome not found: {test_data}")
        sys.exit(1)

    outdir = os.environ.get("SSIGN_TEST_OUTDIR", "/tmp/ssign_e2e_test")
    os.makedirs(outdir, exist_ok=True)

    config = PipelineConfig(
        input_path=test_data,
        sample_id="xantest",
        outdir=outdir,
        # Core settings
        wholeness_threshold=0.8,
        excluded_systems=["Flagellum", "Tad", "T3SS"],
        conf_threshold=0.8,
        proximity_window=3,
        # Skip heavy tools
        skip_deepsece=True,     # ESM model too large for test env
        skip_signalp=True,      # Optional
        skip_hhsuite=True,      # Optional, slow
        skip_foldseek=True,
        skip_plmblast=True,
        skip_structure=True,
        # Enable core annotation
        skip_blastp=False,
        blastp_mode="remote",
        skip_interproscan=False,
        interproscan_mode="remote",
        skip_protparam=False,
        # DeepLocPro
        deeplocpro_mode="remote",
    )

    def progress_cb(step, pct, msg):
        print(f"[{pct:3d}%] {step}: {msg}")

    runner = PipelineRunner(config, progress_callback=progress_cb)
    results = runner.run(resume=True)

    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE RESULTS")
    print("=" * 60)
    passed = 0
    failed = 0
    for r in results:
        status = "OK" if r.success else "FAIL"
        if r.success:
            passed += 1
        else:
            failed += 1
        print(f"  [{status:4s}] {r.name}: {r.message}")

    print(f"\n{passed} passed, {failed} failed out of {len(results)} steps")

    # Check key output files
    print("\n--- Output files ---")
    for f in sorted(os.listdir(outdir)):
        fpath = os.path.join(outdir, f)
        if os.path.isfile(fpath):
            size = os.path.getsize(fpath)
            print(f"  {f}: {size:,} bytes")
        elif os.path.isdir(fpath):
            n = len(os.listdir(fpath))
            print(f"  {f}/: {n} files")

    sys.exit(0 if failed == 0 else 1)

if __name__ == "__main__":
    main()
