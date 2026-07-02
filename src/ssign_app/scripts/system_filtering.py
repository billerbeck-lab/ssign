#!/usr/bin/env python3
"""Apply system-level filtering and DSE cross-genome validation.

Merges proximity substrates with T5SS substrates, runs the
literature-derived localization-correctness gate over each
MacSyFinder-validated system, and produces filtered + unfiltered
substrate lists.

Preserves the dse_type_in_genome bug fix from Session 11.
"""

import argparse
import csv
import logging
import os as _os
import sys as _sys

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)

from ssign_lib.constants import (  # noqa: E402
    DEFAULT_EXCLUDED_SYSTEMS,  # noqa: F401  (default for --excluded-systems)
    SUBSTRATE_SOURCE_T5_SELF,
    is_sec_signal_peptide,
)
from ssign_lib.localization_gate import (  # noqa: E402
    aggregate_failed_ss_types,
    default_localization_table_path,
    evaluate_system,
    load_component_localizations,
)
from ssign_lib.tsv_io import load_tsv_by_key  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def _load_systems_with_components(ss_components_path: str) -> dict:
    """Group ss_components rows by sys_id.

    Returns ``{sys_id: {"ss_type": str, "components": [(locus_tag, gene_name), ...]}}``.
    Excluded rows (``excluded`` column truthy) are skipped because the
    downstream excluded-systems filter will drop them anyway, and they
    shouldn't influence the gate's per-system aggregation.
    """
    if not ss_components_path or not _os.path.exists(ss_components_path):
        return {}
    systems: dict[str, dict] = {}
    with open(ss_components_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if str(row.get("excluded", "False")).strip().lower() == "true":
                continue
            sys_id = (row.get("sys_id") or "").strip()
            locus = (row.get("locus_tag") or "").strip()
            gene = (row.get("gene_name") or "").strip()
            ss_type = (row.get("ss_type") or "").strip()
            if not sys_id or not locus:
                continue
            entry = systems.setdefault(sys_id, {"ss_type": ss_type, "components": []})
            entry["components"].append((locus, gene))
            if not entry["ss_type"]:
                entry["ss_type"] = ss_type
    return systems


def _t5_self_has_evidence(row: dict, conf_threshold: float) -> bool:
    """A T5SS-self component is a substrate iff it has localization evidence
    (DeepLocPro) OR a Sec signal peptide (SignalP).

    T5 passengers are Sec-exported, so a clear Sec signal calls the component even
    when DeepLocPro misses it; a component with neither is dropped (previously all
    detected T5 components were reported unconditionally). Preserves the
    per-component DeepLocPro rule (critical bug-fix #6): the T5bSS translocator
    (TpsB pore) passes on outer-membrane only, never extracellular; T5aSS/T5cSS
    pass on extracellular OR outer-membrane. openspec: signalp-t5ss-substrate-call.
    """
    ss_type = row.get("nearby_ss_types", "")

    def _p(col: str) -> float:
        try:
            return float(row.get(col, 0) or 0)
        except (ValueError, TypeError):
            return 0.0

    if ss_type == "T5bSS":  # translocator is OM-only
        dlp_pass = _p("outer_membrane_prob") >= conf_threshold
    else:  # T5aSS / T5cSS: extracellular or outer-membrane
        dlp_pass = max(_p("dlp_extracellular_prob"), _p("outer_membrane_prob")) >= conf_threshold
    return dlp_pass or is_sec_signal_peptide(row.get("signalp_prediction", ""))


def _run_localization_gate(
    ss_components_path: str,
    predictions_path: str,
    table_path: str,
    dlp_confidence_threshold: float,
    required_fraction_correct: float,
) -> tuple[set[str], list]:
    """Score every system; return (failed_ss_types, per-system verdicts).

    Empty inputs short-circuit to "no failures" so the gate is a no-op
    when ss_components.tsv or predictions.tsv is missing (e.g. genomes
    with no MacSyFinder hits at all).
    """
    systems = _load_systems_with_components(ss_components_path)
    if not systems:
        return set(), []
    rules = load_component_localizations(table_path) if table_path else load_component_localizations()
    # Reuses the shared TSV-to-dict loader (#100); missing_ok=True returns
    # {} for a missing/empty path which short-circuits the per-system scoring
    # to "no DLP info -> n_scored=0 -> fail-open" by design.
    predictions = load_tsv_by_key(predictions_path, key_columns=("locus_tag",))

    verdicts = []
    for sys_id, info in systems.items():
        verdict = evaluate_system(
            sys_id=sys_id,
            ss_type=info["ss_type"],
            components=info["components"],
            predictions=predictions,
            rules=rules,
            dlp_confidence_threshold=dlp_confidence_threshold,
            required_fraction_correct=required_fraction_correct,
        )
        verdicts.append(verdict)

    failed_types = aggregate_failed_ss_types(verdicts)
    n_fail = sum(1 for v in verdicts if not v.passed)
    logger.info(
        "Localization gate: %d/%d systems failed; %d ss_types where every system failed (%s)",
        n_fail,
        len(verdicts),
        len(failed_types),
        sorted(failed_types) if failed_types else "none",
    )
    return failed_types, verdicts


def _write_verdict_log(out_path: str, sample: str, verdicts: list) -> None:
    """Persist per-system gate verdicts for downstream reporting / debugging."""
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "sys_id", "ss_type", "n_scored", "n_correct", "fraction_correct", "passed"])
        for v in verdicts:
            writer.writerow(
                [
                    sample,
                    v.sys_id,
                    v.ss_type,
                    v.n_scored,
                    v.n_correct,
                    f"{v.fraction_correct:.3f}",
                    "True" if v.passed else "False",
                ]
            )


def main():
    parser = argparse.ArgumentParser(description="Filter substrates")
    parser.add_argument("--proximity-substrates", required=True)
    parser.add_argument("--t5ss-substrates", required=True)
    parser.add_argument("--valid-systems", required=True)
    parser.add_argument("--ss-components", default="", help="ss_components TSV for the localization gate")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--excluded-systems", default=",".join(DEFAULT_EXCLUDED_SYSTEMS))
    parser.add_argument("--required-fraction-correct", type=float, default=0.8)
    parser.add_argument(
        "--dlp-confidence-threshold",
        type=float,
        default=0.8,
        help="Minimum DLP max-probability for a system component to count in the gate's correctness calc",
    )
    parser.add_argument(
        "--component-localizations-tsv",
        default="",
        help="Override path to the literature-derived expected-localization table",
    )
    parser.add_argument(
        "--skip-localization-gate",
        action="store_true",
        help="Disable the literature-localization gate entirely (kept for ad-hoc debugging)",
    )
    parser.add_argument(
        "--filter-dse-type-mismatch",
        action="store_true",
        help="Remove DSE-only substrates where DSE predicted type doesn't match nearby MacSyFinder system type",
    )
    parser.add_argument("--out-filtered", required=True)
    parser.add_argument("--out-all", required=True)
    parser.add_argument("--out-gate-verdicts", default="", help="Optional path for per-system gate verdicts TSV")
    args = parser.parse_args()

    excluded = set(s.strip() for s in args.excluded_systems.split(",") if s.strip())

    # Run the localization-correctness gate (#58). Skipping the gate or
    # missing ss-components / predictions leaves failed_types empty, in
    # which case the rest of the function behaves exactly as before.
    failed_types: set[str] = set()
    verdicts: list = []
    if not args.skip_localization_gate:
        try:
            failed_types, verdicts = _run_localization_gate(
                ss_components_path=args.ss_components,
                predictions_path=args.predictions,
                table_path=args.component_localizations_tsv or str(default_localization_table_path()),
                dlp_confidence_threshold=args.dlp_confidence_threshold,
                required_fraction_correct=args.required_fraction_correct,
            )
        except (OSError, csv.Error, ValueError) as e:
            # Don't let a malformed table or missing column take down the
            # whole substrate-filtering step — log loudly and fail-open.
            logger.warning("Localization gate skipped due to error: %s", e)
        if args.out_gate_verdicts and verdicts:
            _write_verdict_log(args.out_gate_verdicts, args.sample, verdicts)

    # Drop systems where the localization-correctness gate decided "not real" —
    # at ss_type granularity (a type only fails if every system of that type
    # failed). These behave like additional excluded_systems for substrate filtering.
    excluded_or_failed = excluded | failed_types

    # Load proximity substrates
    substrates = []
    with open(args.proximity_substrates) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["substrate_source"] = "proximity"
            substrates.append(row)

    # Load T5SS substrates
    with open(args.t5ss_substrates) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["substrate_source"] = SUBSTRATE_SOURCE_T5_SELF
            substrates.append(row)

    # Write unfiltered (all substrates before exclusion). Union keys across
    # every row — proximity rows and T5SS-self rows carry different columns
    # (T5SS-self has t5_quality_flag, proximity has plm_effector_* etc.) and
    # taking substrates[0].keys() would silently drop one side's columns.
    if substrates:
        seen: dict[str, None] = {}
        for s in substrates:
            for k in s:
                seen.setdefault(k, None)
        fieldnames = list(seen)
    else:
        fieldnames = ["locus_tag", "sample_id", "tool", "nearby_ss_types", "substrate_source"]

    with open(args.out_all, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for s in substrates:
            writer.writerow(s)

    # Filter: remove substrates associated only with excluded or
    # gate-failed SS types. Also optionally filter DSE-type mismatches.
    filtered = []
    n_dse_mismatch_removed = 0
    n_gate_dropped = 0
    n_t5_no_evidence_dropped = 0
    for s in substrates:
        ss_types = set(s.get("nearby_ss_types", "").split(","))
        ss_types.discard("")

        # T5SS-self components must carry localization (DeepLocPro) or Sec-signal
        # (SignalP) evidence to be called substrates; a detected T5 component with
        # neither is dropped (previously every T5 component was kept
        # unconditionally). openspec: signalp-t5ss-substrate-call.
        if s.get("substrate_source") == SUBSTRATE_SOURCE_T5_SELF and not _t5_self_has_evidence(
            s, args.dlp_confidence_threshold
        ):
            n_t5_no_evidence_dropped += 1
            continue

        # Keep if any ss_type that's neither excluded nor gate-failed survives
        surviving = ss_types - excluded_or_failed
        if surviving or s.get("substrate_source") == SUBSTRATE_SOURCE_T5_SELF:
            # Trim the displayed nearby_ss_types to the surviving set
            if surviving:
                s["nearby_ss_types"] = ",".join(sorted(surviving))

            # DSE type-match filter: for DSE-only substrates, check if
            # DSE predicted type matches the nearby MacSyFinder system
            if args.filter_dse_type_mismatch:
                tool = s.get("tool", "")
                dse_match = s.get("dse_type_match", "True")
                # Only filter DSE-only substrates (not DLP or DLP+DSE)
                if tool == "DSE" and str(dse_match).lower() in ("false", "0", ""):
                    n_dse_mismatch_removed += 1
                    continue

            filtered.append(s)
        elif ss_types & failed_types:
            # Substrate has no surviving types and at least one of the failed
            # types is to blame; count it. Excluded-only drops are not counted
            # here (historical behaviour: only the new gate's drops are logged).
            n_gate_dropped += 1

    if n_dse_mismatch_removed:
        logger.info(
            f"DSE type-match filter removed {n_dse_mismatch_removed} DSE-only substrates "
            f"where predicted SS type didn't match nearby MacSyFinder system"
        )
    if n_gate_dropped:
        logger.info(
            "Localization gate removed %d substrate(s) tagged only with failed ss_types (%s)",
            n_gate_dropped,
            sorted(failed_types),
        )
    if n_t5_no_evidence_dropped:
        logger.info(
            "T5SS evidence gate removed %d T5 component(s) with neither DeepLocPro localization "
            "nor a Sec signal peptide",
            n_t5_no_evidence_dropped,
        )

    with open(args.out_filtered, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for s in filtered:
            writer.writerow(s)

    logger.info(
        f"{args.sample}: {len(substrates)} total substrates, {len(filtered)} after filtering "
        f"(excluded: {sorted(excluded)}, gate-failed types: {sorted(failed_types)})"
    )


if __name__ == "__main__":
    main()
