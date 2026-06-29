#!/usr/bin/env python3
"""Per-system + per-broad-type binomial enrichment test.

For each MacSyFinder secretion system in a genome, compares the fraction
of DLP-positive (and DSE-positive) proteins in the system's +/-W
neighborhood against the genome-specific background rates ``p_DLP`` /
``p_DSE``, estimated from a random sample of non-neighborhood proteins
written by ``sample_null_proteins.py``. One-sided binomial test
(alternative is "enriched"); BH FDR across all (scope x tool) tests.

This replaces the earlier Fisher's-exact + circular-permutation path,
which had two bugs documented in NOTES.md: the permutation branch's
inner loop body was ``pass`` (dead code, never iterated), and the
Fisher contingency ``total`` counted ``(locus, ss_type)`` pairs rather
than unique substrates.

Positivity rules mirror proximity_analysis.py:162 minus its cross-genome
leakage guard (applying that guard asymmetrically to neighborhood vs
null would distort the background rate).
"""

import argparse
import csv
import logging
import os as _os
import re
import sys as _sys
from collections import defaultdict

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)

from extract_neighborhood import (  # noqa: E402
    get_neighborhood_proteins,
    load_gene_order,
)
from ssign_lib.tsv_io import load_tsv_by_key  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DSE_NEGATIVE = {"Non-secreted", "T3SS", ""}

OUT_FIELDS = [
    "sample_id",
    "scope_kind",
    "scope_id",
    "ss_type",
    "tool",
    "M",
    "k",
    "p_bg",
    "fold_enrich",
    "pvalue",
    "qvalue",
    "significant",
    "n_null",
]


def broad_type(ss_type: str) -> str:
    """Collapse subtype labels to TxSS for per-broad-type aggregation.

    T5aSS/T5bSS/T5cSS -> T5SS; T6SSi/T6SSii -> T6SS; pT4SSt -> T4SS;
    T1SS stays T1SS. Non-matching labels (Flagellum, Tad, ...) pass
    through unchanged.
    """
    m = re.match(r"p?(T\d+)[a-z]*SS", ss_type)
    return f"{m.group(1)}SS" if m else ss_type


def is_dlp_positive(row: dict, conf: float) -> bool:
    try:
        return float(row.get("dlp_extracellular_prob", row.get("extracellular_prob", 0))) >= conf
    except (ValueError, TypeError):
        return False


def is_dse_positive(row: dict, conf: float) -> bool:
    ss_type = row.get("dse_ss_type", "").strip()
    if ss_type in DSE_NEGATIVE:
        return False
    try:
        return float(row.get("dse_max_prob", 0)) >= conf
    except (ValueError, TypeError):
        return False


def load_predictions_keyed(path: str) -> dict:
    """Read a tool-output TSV indexed by locus_tag.

    Thin wrapper around ``load_tsv_by_key`` keyed strictly on
    ``locus_tag``; this script's inputs always come from the runner's
    prediction-tool wrappers and have the canonical column name.
    """
    # annotate the result so the flat ssign_lib import (Any under the mypy override) is narrowed to its real type
    keyed: dict[str, dict[str, str]] = load_tsv_by_key(path, key_columns=("locus_tag",))
    return keyed


def load_systems(ss_components_path: str):
    """Group components by sys_id; skip excluded rows (mirrors proximity_analysis).

    Returns (systems_by_sys_id: {sys_id: set(locus_tag)},
             ss_type_by_sys_id: {sys_id: ss_type}).
    """
    systems: dict[str, set] = defaultdict(set)
    ss_type_of_sys: dict[str, str] = {}
    with open(ss_components_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("excluded", "False").lower() == "true":
                continue
            sys_id = row.get("sys_id", "")
            if not sys_id:
                continue
            systems[sys_id].add(row["locus_tag"])
            ss_type_of_sys[sys_id] = row.get("ss_type", "")
    return dict(systems), ss_type_of_sys


def binom_pvalue(k: int, n: int, p: float) -> float:
    """One-sided binomial test ``P(X >= k | n, p)``.

    Returns 1.0 for degenerate inputs (n<=0, p<=0, p>=1) so the
    downstream BH step still ranks them.
    """
    if n <= 0 or p <= 0 or p >= 1:
        return 1.0
    from scipy.stats import binomtest

    return float(binomtest(k, n, p, alternative="greater").pvalue)


def bh_fdr(
    rows: list, pvalue_key: str = "pvalue", q_key: str = "qvalue", sig_key: str = "significant", alpha: float = 0.05
) -> None:
    """In-place BH FDR. Adds ``qvalue`` (monotone non-decreasing) and ``significant`` (q < alpha)."""
    if not rows:
        return
    indexed = sorted(enumerate(rows), key=lambda kv: (kv[1][pvalue_key], kv[1].get("scope_id", "")))
    n = len(indexed)
    q_raw = [r[pvalue_key] * n / (rank0 + 1) for rank0, (_, r) in enumerate(indexed)]
    # Enforce monotone non-decreasing q across ascending-p order: walk
    # backward keeping the running minimum.
    running_min = 1.0
    q_adj = q_raw[:]
    for i in range(n - 1, -1, -1):
        running_min = min(running_min, q_raw[i])
        q_adj[i] = min(running_min, 1.0)
    for rank0, (orig_idx, _) in enumerate(indexed):
        rows[orig_idx][q_key] = round(q_adj[rank0], 6)
        rows[orig_idx][sig_key] = q_adj[rank0] < alpha


def score_scope(scope_id, ss_type, scope_kind, neigh_loci, dlp, dse, p_dlp, p_dse, conf):
    """Build one row per tool (DLP, DSE) for a given scope (system or broad_type).

    PLM-Effector is deliberately NOT an enrichment predictor: at genome scale it
    calls a large fraction of the proteome positive (~25% on PAO1, recall-tuned +
    OR-of-5-ensembles), so its background swamps any neighborhood signal and it
    showed no reliable per-system enrichment. See the openspec change
    enrichment-background-and-plme-default-off.
    """
    if not neigh_loci:
        return []
    M = len(neigh_loci)
    tools = [
        ("DLP", sum(1 for L in neigh_loci if is_dlp_positive(dlp.get(L, {}), conf)), p_dlp),
        ("DSE", sum(1 for L in neigh_loci if is_dse_positive(dse.get(L, {}), conf)), p_dse),
    ]
    out = []
    for tool, k, p_bg in tools:
        fold = round((k / M) / p_bg, 4) if p_bg > 0 else ""
        out.append(
            {
                "scope_kind": scope_kind,
                "scope_id": scope_id,
                "ss_type": ss_type,
                "tool": tool,
                "M": M,
                "k": k,
                "p_bg": round(p_bg, 6),
                "fold_enrich": fold,
                "pvalue": round(binom_pvalue(k, M, p_bg), 6),
            }
        )
    return out


def write_rows(path: str, sample_id: str, rows: list, n_null: int) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t")
        writer.writeheader()
        for r in rows:
            r_out = dict(r)
            r_out["sample_id"] = sample_id
            r_out["n_null"] = n_null
            writer.writerow(r_out)


def main():
    parser = argparse.ArgumentParser(description="Per-system binomial enrichment test")
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--gene-order", required=True)
    parser.add_argument("--dlp", required=True, help="DeepLocPro output TSV")
    parser.add_argument("--dse", required=True, help="DeepSecE output TSV")
    parser.add_argument("--null-ids", required=True, help="One locus_tag per line: the null sample")
    parser.add_argument("--window", type=int, default=3)
    parser.add_argument("--conf-threshold", type=float, default=0.8)
    parser.add_argument("--sample", required=True, help="Sample / genome ID")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    null_ids = {line.strip() for line in open(args.null_ids) if line.strip()}
    dlp = load_predictions_keyed(args.dlp)
    dse = load_predictions_keyed(args.dse)

    n_null = len(null_ids)
    p_dlp = (
        sum(1 for nid in null_ids if is_dlp_positive(dlp.get(nid, {}), args.conf_threshold)) / n_null if n_null else 0.0
    )
    p_dse = (
        sum(1 for nid in null_ids if is_dse_positive(dse.get(nid, {}), args.conf_threshold)) / n_null if n_null else 0.0
    )
    logger.info("Null sample: %d proteins; p_DLP=%.4f, p_DSE=%.4f", n_null, p_dlp, p_dse)

    systems, ss_type_of_sys = load_systems(args.ss_components)
    if not systems:
        logger.warning("No SS components found; writing header-only output")
        write_rows(args.out, args.sample, [], n_null)
        return

    gene_order = load_gene_order(args.gene_order)
    all_components = set().union(*systems.values())

    rows = []
    # Per-system tests
    for sys_id, components in systems.items():
        ss_type = ss_type_of_sys.get(sys_id, "")
        neigh = get_neighborhood_proteins(gene_order, components, args.window) - all_components
        rows.extend(score_scope(sys_id, ss_type, "system", neigh, dlp, dse, p_dlp, p_dse, args.conf_threshold))

    # Per-broad-type aggregates (only when there are multiple systems of
    # that broad type; otherwise the aggregate equals the single system).
    type_to_sys_ids: dict[str, list] = defaultdict(list)
    for sys_id, ss_type in ss_type_of_sys.items():
        type_to_sys_ids[broad_type(ss_type)].append(sys_id)
    for bt, sys_ids in type_to_sys_ids.items():
        if len(sys_ids) <= 1:
            continue
        combined = set().union(*(systems[sid] for sid in sys_ids))
        neigh = get_neighborhood_proteins(gene_order, combined, args.window) - all_components
        rows.extend(score_scope(bt, bt, "broad_type", neigh, dlp, dse, p_dlp, p_dse, args.conf_threshold))

    bh_fdr(rows)
    write_rows(args.out, args.sample, rows, n_null)
    n_sig = sum(1 for r in rows if r.get("significant"))
    logger.info("Wrote %d (scope x tool) tests; %d significant at q < 0.05", len(rows), n_sig)


if __name__ == "__main__":
    main()
