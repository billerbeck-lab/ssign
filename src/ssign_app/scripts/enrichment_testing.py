#!/usr/bin/env python3
"""Per-SS-type circular-shift enrichment test.

For each secretion-system type in a genome, asks whether secreted-predicted
proteins cluster around that type's components more than a genome-structure-
preserving null. The null is the exact set of circular rotations of the
predictor's gene-ordered positivity vector: rotate the positives by 1, 2, ...,
n-1 genes and recount how many land in the type's windows. Offset 0 is the
observed value. The full all-rotations count is the circular cross-correlation
of the positivity vector and the window mask, computed in one FFT pass.

Three predictors (DLP, DSE, SignalP) are tested per type. Two type classes:
  - window types (T1SS, T2SS, T3SS, T4SS, pT4SSt, T5bSS, T6SSi, T6SSii): the
    window mask marks genes within +/-W of any component of that type (minus
    the component positions themselves), and the positivity vector is the
    standard secreted call (DLP extracellular, DSE secreted-type, SignalP
    Sec signal peptide).
  - autotransporter types (T5aSS, T5cSS): the component IS the substrate, so
    the mask marks the component positions themselves and the positivity vector
    is self-detection (DLP outer-membrane-or-extracellular; DSE secreted-type;
    SignalP Sec signal peptide). "Observed = how many of this type's components
    are self-secreted."

Per type x tool emits fold = observed / null-mean and an exact permutation
p-value; BH-corrects across the real (non-skipped) tests. DSE is not tested for
T3SS; SignalP is tested for every type (its low fold on Sec-bypassing T3SS/T6SS
effectors is an informative contrast, not a skip). A combined one-bar-per-type
track is "DLP or DSE" for window types and SignalP-alone for autotransporters.
PLM-Effector is deliberately not an enrichment predictor (it over-predicts at
genome scale). Replaces the earlier binomial test.

Approximation: multi-contig genomes are ordered contig-then-start and rotated as
one circular replicon (a rotation can wrap across a contig junction).
"""

import argparse
import csv
import logging
import os as _os
import re
import sys as _sys
from collections import defaultdict

import numpy as np

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)

from extract_neighborhood import load_gene_order  # noqa: E402
from ssign_lib.constants import (  # noqa: E402
    CONF_THRESHOLD,
    ENRICH_AUTOTRANSPORTER_TYPES,
    ENRICH_COMBINED_TOOL,
    ENRICH_DSE_NO_WINDOW,
    ENRICH_MAX_NULL,
    ENRICH_TOOLS,
    ENRICH_WINDOW_TYPES,
    PROXIMITY_WINDOW,
    display_type,
    enrich_null_key,
    is_sec_signal_peptide,
)
from ssign_lib.tsv_io import load_tsv_by_key  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Mirrors ENRICH_DSE_NO_WINDOW: DeepSecE can't reliably call T3SS (CLAUDE.md bug #4).
DSE_NEGATIVE = {"Non-secreted", "T3SS", ""}

OUT_FIELDS = [
    "sample_id",
    "ss_type",
    "tool",
    "mode",
    "observed",
    "n_mask",
    "null_mean",
    "fold",
    "p_perm",
    "qvalue",
    "significant",
    "n_rotations",
]


def broad_type(ss_type: str) -> str:
    """Collapse all subtype labels to TxSS (incl. T5a/T5b/T5c -> T5SS).

    Retained only for the validation/comparison scripts (e.g.
    fleet_67/03_permutation_refined.py); the shipped pipeline (per-run test and
    cross-genome pooling) uses ``display_type``, which keeps T5 subtypes distinct
    because autotransporters are tested differently.
    """
    m = re.match(r"p?(T\d+)[a-z]*SS", ss_type)
    return f"{m.group(1)}SS" if m else ss_type


def binom_pvalue(k: int, n: int, p: float) -> float:
    """One-sided binomial test ``P(X >= k | n, p)``.

    Retained as a utility for the validation/comparison scripts that contrast
    the binomial against the circular-shift null. The shipped enrichment test no
    longer uses it (the binomial under-states significance by ignoring secreted-
    gene clustering). Returns 1.0 for degenerate inputs.
    """
    if n <= 0 or p <= 0 or p >= 1:
        return 1.0
    from scipy.stats import binomtest

    return float(binomtest(k, n, p, alternative="greater").pvalue)


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


def is_signalp_positive(row: dict) -> bool:
    """Sec signal peptide present (shared rule: ``is_sec_signal_peptide``).

    No probability threshold: SignalP's own class call is the decision, matching
    the rest of the pipeline. Used in both the window (T2SS/T5bSS...) and
    autotransporter self-detection modes; SignalP positivity is the same vector
    either way (the mask differs, not the call).
    """
    return is_sec_signal_peptide(row.get("signalp_prediction", ""))


def is_dlp_self_positive(row: dict, conf: float) -> bool:
    """Autotransporter self-detection: outer-membrane OR extracellular >= conf.

    A T5aSS/T5cSS autotransporter threads its own passenger through an integral
    outer-membrane beta-barrel, so DeepLocPro should call it OM or extracellular;
    either counts as the component detecting itself.
    """
    try:
        om = float(row.get("outer_membrane_prob") or 0)
        extra = float(row.get("dlp_extracellular_prob", row.get("extracellular_prob", 0)) or 0)
    except (ValueError, TypeError):
        return False
    return om >= conf or extra >= conf


def load_predictions_keyed(path: str) -> dict:
    """Read a tool-output TSV indexed by locus_tag."""
    return load_tsv_by_key(path, key_columns=("locus_tag",))


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


def bh_fdr(
    rows: list, pvalue_key: str = "pvalue", q_key: str = "qvalue", sig_key: str = "significant", alpha: float = 0.05
) -> None:
    """In-place BH FDR. Adds ``qvalue`` (monotone non-decreasing) and ``significant`` (q < alpha)."""
    if not rows:
        return
    indexed = sorted(enumerate(rows), key=lambda kv: (kv[1][pvalue_key], kv[1].get("ss_type", "")))
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


def bh_fdr_by_family(rows: list) -> None:
    """BH-correct the per-tool (DLP/DSE) and combined (DLP-or-DSE) tests as two
    separate multiple-testing families. The combined track is derived from the
    per-tool vectors, so it gets its own BH rather than padding the per-tool
    denominator. Mutates ``rows`` in place; shared by the single-genome test and
    the cross-genome pooling so the family rule lives in one place."""
    bh_fdr([r for r in rows if r["tool"] in ENRICH_TOOLS], pvalue_key="p_perm")
    bh_fdr([r for r in rows if r["tool"] == ENRICH_COMBINED_TOOL], pvalue_key="p_perm")


# ── circular-shift core ──


def gene_order_flat(gene_order_path: str) -> list:
    """Flat list of locus_tags in circular gene order (contig-then-position)."""
    contigs = load_gene_order(gene_order_path)  # {contig: [(pos, locus), ...]} pre-sorted
    order = []
    for contig in sorted(contigs):
        order.extend(locus for _pos, locus in contigs[contig])
    return order


def positivity_vectors(order: list, dlp: dict, dse: dict, signalp: dict, conf: float):
    """Per-protein positivity arrays over the whole genome in gene order.

    Returns (idx, vecs) where idx maps locus_tag -> ordinal and vecs holds the
    gene-ordered positivity arrays keyed "dlp"/"dse"/"dlp_self"/"signalp" (the keys
    run_enrichment reads). A locus absent from a prediction TSV scores 0 (treated
    negative); whole-genome DLP/DSE/SignalP (forced when enrichment is on) means
    that is rare. SignalP uses one vector for both window and self modes (the
    Sec-signal call is the same; only the mask differs).
    """
    idx = {lt: i for i, lt in enumerate(order)}
    vecs = {
        "dlp": np.array([is_dlp_positive(dlp.get(lt, {}), conf) for lt in order], dtype=float),
        "dse": np.array([is_dse_positive(dse.get(lt, {}), conf) for lt in order], dtype=float),
        "dlp_self": np.array([is_dlp_self_positive(dlp.get(lt, {}), conf) for lt in order], dtype=float),
        "signalp": np.array([is_signalp_positive(signalp.get(lt, {})) for lt in order], dtype=float),
    }
    return idx, vecs


def window_mask(comp_idx, n: int, w: int, exclude=()) -> np.ndarray:
    """Length-n mask: 1 within +/-w of any component (circular), then 0 at every
    position in `exclude` (the components themselves, matching the production
    neighborhood which is scored minus its components)."""
    mask = np.zeros(n, dtype=float)
    for p in comp_idx:
        for d in range(-w, w + 1):
            mask[(p + d) % n] = 1.0
    for p in exclude:
        mask[p] = 0.0
    return mask


def self_mask(comp_idx, n: int) -> np.ndarray:
    """Length-n mask marking the component positions themselves (autotransporters)."""
    mask = np.zeros(n, dtype=float)
    for p in comp_idx:
        mask[p] = 1.0
    return mask


def rotation_counts(pos_vec: np.ndarray, win_mask: np.ndarray) -> np.ndarray:
    """All-rotations count array c[k] = # positives landing in the window after a
    circular shift of k. Circular cross-correlation via FFT; c[0] == observed."""
    n = len(pos_vec)
    c = np.fft.irfft(np.fft.rfft(win_mask) * np.conj(np.fft.rfft(pos_vec)), n)
    return np.rint(c).astype(int)


def single_test(pos_vec: np.ndarray, win_mask: np.ndarray, rng=None):
    """Observed + exact circular-shift null for one (positivity, mask) pair.

    Returns a dict with observed (c[0]), the null array (c[1:], optionally
    subsampled for storage when the genome is huge), null_mean, fold, exact
    permutation p = (#{null >= observed} + 1)/(len(null)+1), and n_rotations."""
    c = rotation_counts(pos_vec, win_mask)
    n = len(c)
    observed = int(c[0])
    null = c[1:]  # the n-1 non-identity rotations = the exact permutation null
    null_mean = float(null.mean()) if null.size else 0.0
    fold = observed / null_mean if null_mean > 0 else (float("inf") if observed > 0 else 0.0)
    p = (int(np.sum(null >= observed)) + 1) / (null.size + 1)
    stored = null
    if null.size > ENRICH_MAX_NULL:  # thin only the stored histogram; p/fold stay exact
        rng = rng or np.random.default_rng(0)
        stored = null[rng.integers(0, null.size, size=ENRICH_MAX_NULL)]
    return {
        "observed": observed,
        "n_mask": int(win_mask.sum()),
        "null": stored,
        "null_mean": round(null_mean, 4),
        "fold": round(fold, 4) if fold != float("inf") else fold,
        "p_perm": round(float(p), 6),
        "n_rotations": n,
    }


def components_by_display_type(systems: dict, ss_type_of_sys: dict):
    """{display_type: set(locus_tag)} pooling a genome's systems of each type, plus
    the full set of component loci (for window exclusion)."""
    by_type: dict[str, set] = defaultdict(set)
    all_components: set = set()
    for sys_id, loci in systems.items():
        dt = display_type(ss_type_of_sys.get(sys_id, ""))
        by_type[dt].update(loci)
        all_components.update(loci)
    return by_type, all_components


def run_enrichment(order, idx, vecs, by_type, all_comp_idx, window):
    """Build one row per (display type, tool); returns the rows list (pre-BH)."""
    n = len(order)
    rows = []
    for dt in sorted(by_type):
        is_auto = dt in ENRICH_AUTOTRANSPORTER_TYPES
        if not is_auto and dt not in ENRICH_WINDOW_TYPES:
            continue  # not an enrichment target (e.g. Flagellum, Tad slipping through)
        comp_idx = [idx[lt] for lt in by_type[dt] if lt in idx]
        if not comp_idx:
            continue
        mode = "self" if is_auto else "window"
        mask = self_mask(comp_idx, n) if is_auto else window_mask(comp_idx, n, window, exclude=all_comp_idx)
        dlp_track = vecs["dlp_self"] if is_auto else vecs["dlp"]
        # SignalP positivity is the same vector in both modes (the Sec-signal call
        # doesn't change; only the mask does), so it needs no per-mode track.
        tool_vec = {"DLP": dlp_track, "DSE": vecs["dse"], "SignalP": vecs["signalp"]}
        for tool in ENRICH_TOOLS:
            if tool == "DSE" and dt in ENRICH_DSE_NO_WINDOW:
                rows.append({"ss_type": dt, "tool": tool, "mode": mode, "skip": True})
                continue
            r = single_test(tool_vec[tool], mask)
            rows.append({"ss_type": dt, "tool": tool, "mode": mode, "skip": False, **r})
        # Combined one-bar-per-type track: "DLP OR DSE" for window types (DSE
        # dropped where unreliable, e.g. T3SS); SignalP-alone for autotransporters,
        # whose passengers are Sec-exported. SignalP's low genome background makes
        # its few-loci self-detection stronger and more specific than a 3-way union
        # (which raises the background and weakens the fold/p). openspec:
        # signalp-enrichment-track decision 4.
        if is_auto:
            combined_pos = vecs["signalp"]
        elif dt in ENRICH_DSE_NO_WINDOW:
            combined_pos = dlp_track
        else:
            combined_pos = np.maximum(dlp_track, vecs["dse"])
        r = single_test(combined_pos, mask)
        rows.append({"ss_type": dt, "tool": ENRICH_COMBINED_TOOL, "mode": mode, "skip": False, **r})
    return rows


def write_rows(path: str, sample_id: str, rows: list) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({**r, "sample_id": sample_id})


def write_nulls(path: str, rows: list) -> None:
    """Dump the per-(type,tool) null arrays for the figure step, keyed by
    ``enrich_null_key`` (shared with the pooling + figure readers)."""
    arrays = {
        enrich_null_key(r["ss_type"], r["tool"]): r["null"] for r in rows if not r["skip"] and r.get("null") is not None
    }
    np.savez_compressed(path, **arrays)


def main():
    parser = argparse.ArgumentParser(description="Per-SS-type circular-shift enrichment test")
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--gene-order", required=True)
    parser.add_argument("--dlp", required=True, help="DeepLocPro output TSV (whole-genome)")
    parser.add_argument("--dse", required=True, help="DeepSecE output TSV (whole-genome)")
    parser.add_argument("--signalp", required=True, help="SignalP output TSV (whole-genome)")
    parser.add_argument("--window", type=int, default=PROXIMITY_WINDOW)
    parser.add_argument("--conf-threshold", type=float, default=CONF_THRESHOLD)
    parser.add_argument("--sample", required=True, help="Sample / genome ID")
    parser.add_argument("--out", required=True, help="Output stats TSV path")
    parser.add_argument("--nulls-out", default="", help="Optional .npz dump of per-type null arrays (for the figure)")
    args = parser.parse_args()

    dlp = load_predictions_keyed(args.dlp)
    dse = load_predictions_keyed(args.dse)
    signalp = load_predictions_keyed(args.signalp)

    systems, ss_type_of_sys = load_systems(args.ss_components)
    if not systems:
        logger.warning("No SS components found; writing header-only output")
        write_rows(args.out, args.sample, [])
        if args.nulls_out:
            write_nulls(args.nulls_out, [])
        return

    order = gene_order_flat(args.gene_order)
    idx, vecs = positivity_vectors(order, dlp, dse, signalp, args.conf_threshold)

    by_type, all_components = components_by_display_type(systems, ss_type_of_sys)
    all_comp_idx = [idx[lt] for lt in all_components if lt in idx]

    rows = run_enrichment(order, idx, vecs, by_type, all_comp_idx, args.window)

    # BH only over the real tests; skipped (T3SS-DSE) rows must not pad the denominator.
    real = [r for r in rows if not r["skip"]]
    bh_fdr_by_family(real)
    for r in rows:
        if r["skip"]:
            r["qvalue"], r["significant"] = 1.0, False

    write_rows(args.out, args.sample, rows)
    if args.nulls_out:
        write_nulls(args.nulls_out, rows)
    n_sig = sum(1 for r in real if r.get("significant"))
    logger.info(
        "Wrote %d (type x tool) circular-shift tests across %d proteins; %d significant at q < 0.05",
        len(real),
        len(order),
        n_sig,
    )


if __name__ == "__main__":
    main()
