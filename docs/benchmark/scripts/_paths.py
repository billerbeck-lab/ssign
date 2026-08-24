"""Shared path resolution for the training-data overlap scripts.

Downloaded training data is large and stays out of git. Point SSIGN_OVERLAP_DIR
at a working directory holding `raw/` (downloaded training sets) and `work/`
(intermediate files); it defaults to the current working directory.
"""

import os

BASE = os.path.abspath(os.environ.get("SSIGN_OVERLAP_DIR", os.getcwd()))
RAW = os.path.join(BASE, "raw")
WORK = os.path.join(BASE, "work")


def repo_root():
    """Walk up from this file until the benchmarking list is found."""
    d = os.path.dirname(os.path.abspath(__file__))
    while d != os.path.dirname(d):
        if os.path.exists(os.path.join(d, "docs", "benchmark", "ssign_benchmarking_list.csv")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("could not locate the ssign repo root from " + __file__)


ROOT = repo_root()
BENCHMARK_CSV = os.path.join(ROOT, "docs", "benchmark", "ssign_benchmarking_list.csv")
OUT_CSV = os.path.join(ROOT, "docs", "benchmark", "training_data_overlap.csv")
os.makedirs(WORK, exist_ok=True)

# Verdict tiers. Written by match_training_set.py, read by aggregate_overlap.py;
# they live here so renaming one cannot silently reclassify the other's counts.
IN_TRAINING = "IN_TRAINING"
NEAR_IDENTICAL_HOMOLOG = "NEAR_IDENTICAL_HOMOLOG"
CLOSE_HOMOLOG = "CLOSE_HOMOLOG"
DISTANT_HOMOLOG = "DISTANT_HOMOLOG"
NO_MEANINGFUL_MATCH = "NO_MEANINGFUL_MATCH"
