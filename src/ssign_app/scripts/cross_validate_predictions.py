#!/usr/bin/env python3
"""Cross-validate secreted-protein predictions across the core tool set.

Inputs (one row per protein each):
    --deeplocpro     DeepLocPro localisation scores (required)
    --deepsece       DeepSecE secretion-type probabilities (optional)
    --plm-effector   PLM-Effector combined-type predictions (optional,
                     one row per protein with `passes_threshold` reflecting
                     "secreted by at least one SS type" — the upstream runner
                     merges the five per-type PLM-Effector outputs before
                     calling cross_validate)
    --signalp        SignalP signal-peptide predictions (optional)
    --valid-systems  MacSyFinder-validated secretion systems for this genome
    --ss-components  Per-protein SS component table (locus_tag → ss_type).
                     Used to apply T5SS-specific localisation rule.

Rule (3.2.b):
    DeepLocPro, DeepSecE, and PLM-Effector are treated as equal
    secretion predictors. Any one flagging is enough to mark a protein
    as `is_secreted=True`; the count of agreeing tools is emitted as
    `n_prediction_tools_agreeing` (0-3). SignalP is evidence-only —
    its determination goes into `signalp_supports_secretion` but does
    not trip `is_secreted` on its own.

T5SS localisation rules are per-component, not per-subtype. TXSScan v2
only models T5aSS, T5bSS, and T5cSS (T5dSS/T5eSS exist in the literature
but lack TXSScan models, so MacSyFinder never emits them). Each modelled
subtype has exactly one HMM-modelled component:

  - T5aSS_PF03797 (autotransporter β-domain): passenger may be cleaved
    (Extracellular) or remain tethered (Outer membrane); max(ext, om) wins.
  - T5bSS_translocator (TpsB pore): Outer membrane only — DLP calling it
    Extracellular is biologically wrong (it's a transmembrane β-barrel).
    The TpsA passenger is NOT in the TXSScan model — proximity_analysis
    catches it as a neighborhood protein and the standard ext-only rule
    applies (correctly: TpsA should be Extracellular).
  - T5cSS_PF03895 (trimeric AT anchor): surface-displayed; either
    Extracellular or Outer membrane is valid.

The (ss_type, gene_name) → probability-columns map lives in
`ssign_lib.constants.T5SS_COMPONENT_RULES`. Unmapped components fall
back to the standard extracellular-only rule.

DSE T3SS reliability guard (preserved from earlier logic): DeepSecE
flags far more T3SS candidates than MacSyFinder validates
(~1800 vs 0 across a 74-genome benchmark). If the genome has no
MacSyFinder-validated T3SS, DSE T3SS calls are flagged
(`dse_T3SS_flagged=True`) and excluded from the DSE trigger count.
"""

import argparse
import csv
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = os.path.dirname(os.path.abspath(__file__))
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)
from ssign_lib.constants import T5SS_COMPONENT_RULES  # noqa: E402

# SignalP's "protein has no signal peptide" sentinels
_SP_NEGATIVE = {"OTHER", "", "No signal peptide"}
# DeepSecE's "protein is not a secretion substrate" sentinels
_DSE_NEGATIVE = {"Non-secreted", "", "OTHER"}


def _load_tsv_by_locus(path: str):
    """Return {locus_tag: row_dict} for a TSV keyed by `locus_tag`. Empty dict if path missing."""
    if not path or not os.path.exists(path):
        return {}
    rows = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tag = row.get("locus_tag") or row.get("protein_id") or row.get("seq_id")
            if tag:
                rows[tag] = row
    return rows


def _float_or_zero(value) -> float:
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0


def _genome_has_t3ss(valid_systems_path: str) -> bool:
    """True if any non-excluded valid system in the genome is a T3SS variant."""
    with open(valid_systems_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("excluded", "False").lower() == "true":
                continue
            if "T3SS" in row.get("ss_type", ""):
                return True
    return False


def _load_ss_component_info(path: str) -> dict:
    """Return {locus_tag: (ss_type, gene_name)} from validate_macsyfinder_systems.py output.

    Empty dict if the path is missing or empty (e.g. genomes with no
    validated systems). Excluded systems are skipped — those proteins
    aren't real components of the surviving SS calls. The gene_name is
    needed alongside ss_type for per-component rules (a T5bSS_translocator
    has different DLP expectations than a T5aSS_PF03797).
    """
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("excluded", "False").lower() == "true":
                continue
            tag = row.get("locus_tag", "").strip()
            ss_type = row.get("ss_type", "").strip()
            gene_name = row.get("gene_name", "").strip()
            if tag and ss_type:
                out[tag] = (ss_type, gene_name)
    return out


def _dlp_flag(
    dlp_row: dict, conf_threshold: float,
    ss_type: str = "", gene_name: str = "",
) -> tuple:
    """(is_secreted_by_dlp, extracellular_prob).

    Component-aware rule: when (ss_type, gene_name) maps to entry in
    T5SS_COMPONENT_RULES, the trigger is `max(probabilities for the
    listed columns) >= conf_threshold`. Otherwise (unmapped components,
    neighborhood proteins) it's the standard ext-only rule. The mapped
    columns encode per-T5-subtype biology — see constants.py.
    """
    ext_prob = _float_or_zero(dlp_row.get("extracellular_prob", 0))
    rule_columns = T5SS_COMPONENT_RULES.get((ss_type, gene_name))
    if rule_columns:
        max_prob = max(_float_or_zero(dlp_row.get(col, 0)) for col in rule_columns)
        return max_prob >= conf_threshold, ext_prob
    return ext_prob >= conf_threshold, ext_prob


def _dse_flag(dse_row: dict, has_t3ss: bool) -> tuple:
    """(is_secreted_by_dse, dse_ss_type, dse_max_prob, t3ss_flagged).

    T3SS flagging: if DeepSecE predicts T3SS but MacSyFinder did not
    validate a T3SS in this genome, mark the prediction as flagged
    (unreliable) and exclude it from the trigger count.
    """
    dse_type = dse_row.get("dse_ss_type", "Non-secreted") if dse_row else "Non-secreted"
    dse_max = _float_or_zero(dse_row.get("dse_max_prob", 0)) if dse_row else 0.0
    t3ss_flagged = dse_type == "T3SS" and not has_t3ss
    is_secreted = dse_type not in _DSE_NEGATIVE and not t3ss_flagged and dse_max > 0
    return is_secreted, dse_type, dse_max, t3ss_flagged


def _plm_effector_flag(plm_row: dict) -> bool:
    """True if PLM-Effector's `passes_threshold` is set for this protein.

    PLM-Effector emits one row per protein with `passes_threshold=1`
    if the ensemble called it a secreted effector for at least one
    secretion-system type. The runner (Phase 3.2.d) is responsible for
    merging the five per-type outputs into a single "secreted by at
    least one" summary row before handing to cross_validate.
    """
    if not plm_row:
        return False
    return str(plm_row.get("passes_threshold", "0")).strip() in ("1", "True", "true")


def _signalp_supports(sp_row: dict) -> tuple:
    """(supports_secretion, prediction, probability).

    SignalP is evidence-only under 3.2.b — `supports_secretion` does not
    contribute to `is_secreted` or `n_prediction_tools_agreeing`.
    """
    sp_pred = sp_row.get("signalp_prediction", "OTHER") if sp_row else "OTHER"
    sp_prob = _float_or_zero(sp_row.get("signalp_probability", 0)) if sp_row else 0.0
    return sp_pred not in _SP_NEGATIVE, sp_pred, sp_prob


FIELDNAMES = [
    "locus_tag",
    "sample_id",
    # DeepLocPro (all scores preserved)
    "predicted_localization",
    "dlp_extracellular_prob",
    "dlp_max_localization",
    "dlp_max_probability",
    "periplasmic_prob",
    "outer_membrane_prob",
    "cytoplasmic_prob",
    # DeepSecE
    "dse_ss_type",
    "dse_max_prob",
    "dse_T3SS_flagged",
    # PLM-Effector
    "plm_effector_secreted",
    "plm_effector_type",
    # SignalP (evidence-only in 3.2.b)
    "signalp_prediction",
    "signalp_probability",
    "signalp_cs_position",
    "signalp_supports_secretion",
    # Aggregate
    "is_secreted",
    "n_prediction_tools_agreeing",
    "secretion_evidence",
    "product",
]


def cross_validate(
    dlp_data: dict,
    dse_data: dict,
    plm_e_data: dict,
    sp_data: dict,
    sample_id: str,
    conf_threshold: float,
    has_t3ss: bool,
    ss_component_info: dict | None = None,
):
    """Yield one output dict per protein across the union of inputs.

    `ss_component_info` is an optional `{locus_tag: (ss_type, gene_name)}`
    map from validate_macsyfinder_systems.py — used to apply per-component
    DLP rules (e.g. T5bSS translocator must be OM, not extracellular).
    Non-component / neighborhood proteins pass through with the standard
    extracellular-only rule.

    Factored out as a pure function so it's directly unit-testable without
    touching the filesystem.
    """
    ss_component_info = ss_component_info or {}
    all_loci = sorted(
        set(dlp_data.keys())
        | set(dse_data.keys())
        | set(plm_e_data.keys())
        | set(sp_data.keys())
    )

    for locus in all_loci:
        dlp = dlp_data.get(locus, {})
        dse = dse_data.get(locus, {})
        plm_e = plm_e_data.get(locus, {})
        sp = sp_data.get(locus, {})

        ss_type, gene_name = ss_component_info.get(locus, ("", ""))
        dlp_secreted, ext_prob = _dlp_flag(dlp, conf_threshold, ss_type, gene_name)
        dse_secreted, dse_type, dse_max, t3ss_flagged = _dse_flag(dse, has_t3ss)
        plm_e_secreted = _plm_effector_flag(plm_e)
        sp_supports, sp_pred, sp_prob = _signalp_supports(sp)

        evidence = []
        if dlp_secreted:
            evidence.append("DeepLocPro")
        if dse_secreted:
            evidence.append("DeepSecE")
        if plm_e_secreted:
            evidence.append("PLM-Effector")

        n_agreeing = len(evidence)

        yield {
            "locus_tag": locus,
            "sample_id": sample_id,
            "predicted_localization": dlp.get("predicted_localization", ""),
            "dlp_extracellular_prob": ext_prob,
            "dlp_max_localization": dlp.get(
                "max_localization", dlp.get("predicted_localization", "")
            ),
            "dlp_max_probability": dlp.get("max_probability", ext_prob),
            "periplasmic_prob": dlp.get("periplasmic_prob", ""),
            "outer_membrane_prob": dlp.get("outer_membrane_prob", ""),
            "cytoplasmic_prob": dlp.get("cytoplasmic_prob", ""),
            "dse_ss_type": dse_type,
            "dse_max_prob": dse_max,
            "dse_T3SS_flagged": t3ss_flagged,
            "plm_effector_secreted": plm_e_secreted,
            "plm_effector_type": plm_e.get("effector_type", ""),
            "signalp_prediction": sp_pred,
            "signalp_probability": sp_prob,
            "signalp_cs_position": sp.get("signalp_cs_position", ""),
            "signalp_supports_secretion": sp_supports,
            "is_secreted": n_agreeing >= 1,
            "n_prediction_tools_agreeing": n_agreeing,
            "secretion_evidence": ",".join(evidence),
            "product": dlp.get("product", dse.get("product", "")),
        }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Cross-validate secretion predictions across DeepLocPro + "
            "DeepSecE + PLM-Effector (equal triggers) and SignalP (evidence-only)."
        )
    )
    parser.add_argument("--deeplocpro", required=True)
    parser.add_argument("--deepsece", default="")
    parser.add_argument(
        "--plm-effector",
        default="",
        help="PLM-Effector combined-type TSV (optional). Missing file == no PLM-E trigger.",
    )
    parser.add_argument("--signalp", default="")
    parser.add_argument("--valid-systems", required=True)
    parser.add_argument(
        "--ss-components",
        default="",
        help=(
            "Per-protein SS component table from validate_macsyfinder_systems.py. "
            "Used to apply the T5SS-specific 'Extracellular OR Outer membrane' "
            "DLP rule. Optional — without it, all proteins use the standard "
            "extracellular-only rule."
        ),
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--conf-threshold", type=float, default=0.8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    has_t3ss = _genome_has_t3ss(args.valid_systems)

    dlp_data = _load_tsv_by_locus(args.deeplocpro)
    dse_data = _load_tsv_by_locus(args.deepsece)
    plm_e_data = _load_tsv_by_locus(args.plm_effector)
    sp_data = _load_tsv_by_locus(args.signalp)
    ss_component_info = _load_ss_component_info(args.ss_components)

    if not dse_data:
        logger.info("DeepSecE not available — running without DSE trigger")
    if not plm_e_data:
        logger.info("PLM-Effector not available — running without PLM-E trigger")
    if not sp_data:
        logger.info("SignalP not available — no signal-peptide evidence")

    n_flagged_t3ss = 0
    n_rows = 0

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in cross_validate(
            dlp_data=dlp_data,
            dse_data=dse_data,
            plm_e_data=plm_e_data,
            sp_data=sp_data,
            sample_id=args.sample,
            conf_threshold=args.conf_threshold,
            has_t3ss=has_t3ss,
            ss_component_info=ss_component_info,
        ):
            if row["dse_T3SS_flagged"]:
                n_flagged_t3ss += 1
            writer.writerow(row)
            n_rows += 1

    if n_flagged_t3ss:
        logger.warning(
            f"Flagged {n_flagged_t3ss} DeepSecE T3SS predictions "
            f"(no MacSyFinder T3SS found in {args.sample})"
        )
    logger.info(f"Cross-validated {n_rows} proteins for {args.sample}")


if __name__ == "__main__":
    main()
