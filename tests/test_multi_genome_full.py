#!/usr/bin/env python3
"""Multi-genome pipeline test with HHpred + ortholog grouping.

Runs 2 genomes, then cross-genome ortholog grouping.
"""

import logging
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from ssign_app.core.runner import (
    PipelineConfig,
    PipelineRunner,
    run_cross_genome_orthologs,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def run_genome(genome_path, sample_id, outdir, enable_hhpred=True):
    """Run full pipeline on one genome."""
    config = PipelineConfig(
        input_path=genome_path,
        sample_id=sample_id,
        outdir=outdir,
        # Core settings
        wholeness_threshold=0.8,
        excluded_systems=["Flagellum", "Tad", "T3SS"],
        conf_threshold=0.8,
        proximity_window=3,
        # Skip heavy tools
        skip_deepsece=True,
        skip_signalp=True,
        skip_plmblast=True,
        skip_structure=True,
        # Enable annotation tools
        skip_blastp=False,
        skip_hhsuite=not enable_hhpred,
        hhsuite_mode="remote",
        skip_interproscan=False,
        interproscan_mode="remote",
        skip_protparam=False,
        # DeepLocPro
        deeplocpro_mode="remote",
        # Ortholog thresholds
        ortholog_min_pident=40.0,
        ortholog_min_qcov=70.0,
    )

    def progress(step, pct, msg):
        logger.info(f"  [{pct:3d}%] {step}: {msg}")

    runner = PipelineRunner(config, progress_callback=progress)
    results = runner.run(resume=True)

    n_ok = sum(1 for r in results if r.success)
    n_total = len(results)
    logger.info(f"  => {sample_id}: {n_ok}/{n_total} steps succeeded")

    for r in results:
        icon = "OK" if r.success else "FAIL"
        logger.info(f"    [{icon}] {r.name}: {r.message}")

    return results


def main():
    test_dir = os.path.join(os.path.dirname(__file__), "data")
    genome1 = os.path.join(test_dir, "Xanthobacter_tagetidis_TagT2C_genomic.gbff")
    genome2 = os.path.join(
        test_dir, "Roseixanthobacter_finlandensis_VTT_E-85241_genomic.gbff"
    )

    for g in [genome1, genome2]:
        if not os.path.exists(g):
            print(f"ERROR: Genome not found: {g}")
            sys.exit(1)

    base_outdir = os.environ.get("SSIGN_TEST_OUTDIR", "/tmp/ssign_multi_full")

    genomes = [
        (genome1, "xanthobacter", os.path.join(base_outdir, "xanthobacter")),
        (genome2, "roseixanthobacter", os.path.join(base_outdir, "roseixanthobacter")),
    ]

    all_results = []
    genome_outdirs = []
    total_start = time.time()

    for genome_path, sample_id, outdir in genomes:
        os.makedirs(outdir, exist_ok=True)
        logger.info(f"\n{'=' * 60}")
        logger.info(f"GENOME: {sample_id} ({time.strftime('%H:%M:%S')})")
        logger.info(f"{'=' * 60}")

        start = time.time()
        results = run_genome(genome_path, sample_id, outdir, enable_hhpred=True)
        elapsed = time.time() - start

        all_results.extend(results)
        genome_outdirs.append(outdir)
        logger.info(f"  Time: {elapsed:.0f}s ({elapsed / 60:.1f}m)")

    # Cross-genome ortholog grouping
    logger.info(f"\n{'=' * 60}")
    logger.info(f"CROSS-GENOME ORTHOLOGS ({time.strftime('%H:%M:%S')})")
    logger.info(f"{'=' * 60}")

    start = time.time()
    xg_result = run_cross_genome_orthologs(
        genome_outdirs=genome_outdirs,
        output_dir=base_outdir,
        min_pident=40.0,
        min_qcov=70.0,
        progress_callback=lambda s, p, m: logger.info(f"  [{p:3d}%] {s}: {m}"),
    )
    elapsed = time.time() - start

    logger.info("\n=== CROSS-GENOME RESULTS ===")
    for k, v in xg_result.items():
        logger.info(f"  {k}: {v}")
    logger.info(f"  Time: {elapsed:.0f}s")

    # Final summary
    total = time.time() - total_start
    logger.info(f"\n{'=' * 60}")
    logger.info("FINAL SUMMARY")
    logger.info(f"{'=' * 60}")

    n_ok = sum(1 for r in all_results if r.success)
    n_total = len(all_results)
    logger.info(f"Steps: {n_ok}/{n_total} succeeded")
    logger.info(f"Cross-genome groups: {xg_result.get('n_groups', 0)}")
    logger.info(f"Total wall time: {total:.0f}s ({total / 60:.1f}m)")

    # List output files
    for gdir in genome_outdirs:
        name = os.path.basename(gdir)
        csvs = [f for f in os.listdir(gdir) if f.endswith(".csv")]
        logger.info(f"  {name}: {len(csvs)} CSVs")

    xg_files = [
        f for f in os.listdir(base_outdir) if f.endswith(".csv") or f.endswith(".faa")
    ]
    logger.info(f"  cross-genome: {xg_files}")

    # Check xg_ortholog_group column exists in integrated CSVs
    import pandas as pd

    for gdir in genome_outdirs:
        name = os.path.basename(gdir)
        int_csvs = [
            f for f in os.listdir(gdir) if "integrated" in f and f.endswith(".csv")
        ]
        for csv_f in int_csvs:
            df = pd.read_csv(os.path.join(gdir, csv_f))
            has_xg = "xg_ortholog_group" in df.columns
            logger.info(
                f"  {name}/{csv_f}: xg_ortholog_group={'YES' if has_xg else 'NO'}, "
                f"rows={len(df)}, cols={len(df.columns)}"
            )

    if n_ok < n_total:
        logger.warning(f"SOME STEPS FAILED: {n_total - n_ok} failures")
        sys.exit(1)
    else:
        logger.info("ALL STEPS PASSED")


if __name__ == "__main__":
    main()
