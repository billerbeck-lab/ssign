#!/usr/bin/env python3
"""Handle T5SS (Type V Secretion System) self-secreting autotransporters.

T5aSS proteins ARE their own substrates (self-secreting). Each MacSyFinder
T5aSS component is classified by a geometric Pfam-domain filter:

    passenger_length = (barrel_start - LINKER_LENGTH) - (sp_end + 1) + 1

where ``barrel_start`` is the envelope start of the PF03797 (autotransporter
β-barrel) hit on the protein, ``sp_end`` is the SignalP cleavage-site position,
and ``LINKER_LENGTH`` (30 aa) is the alpha-helical linker between passenger
and barrel. The HMM (PF03797 + the porin PF13505 for OMP discrimination) is
bundled with the package at ``src/ssign_app/data/t5aSS/``.

All five geometric classes are emitted as substrates. Lower-confidence ones
carry a ``t5_quality_flag`` so users can see them but the master CSV sorts
them to the bottom (see integrate_annotations.py).

Classification (``domain_group``) and quality flag:
- Classical AT          — PF03797, passenger >= 100 aa      → flag = ""
- Minimal passenger     — PF03797, passenger 1-99 aa        → flag = ""
- Barrel-only           — PF03797, passenger == 0 aa        → flag = "barrel_only"
- OMP/Porin             — PF13505 only, no PF03797          → flag = "omp_porin_no_at"
- Unclassified-AT       — MacSyFinder called T5aSS but neither HMM hit → flag = "unclassified"

When SignalP returned no cleavage site, an additional ``no_signalp`` flag
takes priority over a clean classification: T5SS substrates are Sec-
dependent so a missing signal peptide is a strong "not secreted" signal,
even if the geometry looks fine.

When SignalP was actually run and produced a prediction that is NOT a
valid Sec signal — SignalP 6 short labels ``TAT`` (Tat/SPI),
``TATLIPO`` (Tat/SPII), or ``PILIN`` (Sec/SPIII), see
``VALID_SEC_SIGNAL_TYPES`` — the row is flagged ``no_sec_signal``
instead. This is distinct from ``no_signalp``: that fires when no
cleavage site was reported at all (SignalP "OTHER" or empty), this
fires when SignalP positively predicted the wrong kind of signal.
Both flags mean "biologically incapable of T5SS" and sort the row
to the bottom equally. When SignalP was not run for the genome at
all (no row has a non-empty ``signalp_prediction``), neither flag
fires, preserving the pre-SignalP pathway.
"""

import argparse
import csv
import logging
import os as _os
import sys as _sys
from importlib.resources import files

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)

from ssign_lib.constants import (  # noqa: E402
    LINKER_LENGTH,
    MIN_PASSENGER_LENGTH,
    is_sec_signal_peptide,
)
from ssign_lib.parsing import parse_int_or_none  # noqa: E402

BUNDLED_HMMS = {"PF03797": "PF03797.hmm", "PF13505": "PF13505.hmm"}


def _hmm_path(filename: str) -> str:
    """Resolve a bundled HMM path via importlib.resources."""
    return str(files("ssign_app.data.t5aSS") / filename)


def _decode(val: bytes | str | None) -> str:
    """pyhmmer returns hit.name as bytes in <=0.10 and str in >=0.11."""
    if isinstance(val, bytes):
        return val.decode()
    return str(val) if val is not None else ""


def scan_bundled_pfams(proteins_fasta: str) -> dict[str, dict[str, tuple[int, int]]]:
    """Run pyhmmer for each bundled HMM. Returns {locus_tag: {pfam: (env_from, env_to)}}.

    Uses each HMM's gathering threshold (GA cutoff). Both PF03797 and PF13505
    have GA set in their headers (26 and 28.7 bits respectively).

    FRAGILE: pyhmmer.hmmsearch API. If pyhmmer changes the call signature or
    domain-attribute names, this breaks. Tested with pyhmmer 0.10.x. If this
    breaks: check src/ssign_app/shims/hmmsearch.py for the equivalent call
    pattern — it uses the same API and will fail in the same way.
    """
    import pyhmmer
    from pyhmmer.easel import Alphabet, SequenceFile
    from pyhmmer.plan7 import HMMFile

    alphabet = Alphabet.amino()
    with SequenceFile(proteins_fasta, digital=True, alphabet=alphabet) as sf:
        targets = sf.read_block()

    hits: dict[str, dict[str, tuple[int, int]]] = {}
    for pfam_id, hmm_filename in BUNDLED_HMMS.items():
        with HMMFile(_hmm_path(hmm_filename)) as hf:
            hmm = next(iter(hf))
        top_hits = next(pyhmmer.hmmsearch(hmm, targets, bit_cutoffs="gathering"))
        for hit in top_hits:
            name = _decode(hit.name)
            if not hit.domains:
                continue
            best = min(hit.domains, key=lambda d: d.i_evalue)
            hits.setdefault(name, {})[pfam_id] = (best.env_from, best.env_to)
    return hits


def classify_t5a(
    pfam_hits: dict[str, tuple[int, int]],
    sp_end: int | None,
) -> tuple[str, int, str]:
    """Apply geometric classifier. Returns (domain_group, passenger_length, quality_flag).

    Quality flags surface the per-row caveat: empty for a clean call, otherwise
    one of T5_QUALITY_FLAG_RANK's keys. Structural problems (no barrel, barrel-only) win
    over a missing signal peptide; `no_signalp` flags an otherwise-clean call
    where the Sec/Tat signal couldn't be detected. The passenger-length
    estimate falls back to sp_end=1 only to keep the group label informative.
    """
    has_barrel = "PF03797" in pfam_hits
    has_porin = "PF13505" in pfam_hits
    no_signalp = sp_end is None or sp_end <= 0

    if not has_barrel and has_porin:
        return "OMP/Porin (no AT barrel)", 0, "omp_porin_no_at"
    if not has_barrel:
        return "Unclassified-AT", 0, "unclassified"

    barrel_start, _barrel_end = pfam_hits["PF03797"]
    effective_sp_end = sp_end if not no_signalp else 1
    passenger_length = max(0, (barrel_start - LINKER_LENGTH) - (effective_sp_end + 1) + 1)

    if passenger_length >= MIN_PASSENGER_LENGTH:
        group = "Classical AT"
    elif passenger_length >= 1:
        group = "Minimal passenger"
    else:
        group = "Barrel-only"

    if group == "Barrel-only":
        flag = "barrel_only"
    elif no_signalp:
        flag = "no_signalp"
    else:
        flag = ""
    return group, passenger_length, flag


def _parse_sp_end(raw: str) -> int | None:
    """Parse SignalP CS position (formats vary: '22', '22-23', '22.0', '').

    Thin wrapper around ``parse_int_or_none(allow_range=True)``; kept
    locally so call sites read naturally as ``_parse_sp_end(...)``.
    """
    return parse_int_or_none(raw, allow_range=True)


def _signalp_was_run(predictions: dict[str, dict]) -> bool:
    """True iff at least one input row carries a non-empty ``signalp_prediction``.

    When SignalP is skipped for the genome the column is empty across
    every row; we use that as the signal to suppress Sec-gating flags.
    Auto-detect from data instead of plumbing a separate flag through
    the runner so wrapper-level invocations behave identically.
    """
    return any(p.get("signalp_prediction", "") for p in predictions.values())


def _apply_sec_signal_gate(existing_flag: str, signalp_prediction: str, signalp_was_run: bool) -> str:
    """Augment ``existing_flag`` with ``no_sec_signal`` if appropriate.

    Fires only when SignalP was run for the genome AND this protein has
    a non-empty, non-"OTHER" prediction that isn't in
    ``VALID_SEC_SIGNAL_TYPES`` (i.e. TAT / TATLIPO / PILIN in SignalP 6
    short labels). Empty / "OTHER" predictions fall through to the
    existing ``no_signalp`` logic. Structural flags (barrel_only,
    omp_porin_no_at, unclassified) win over the Sec-signal gate so the
    most informative failure mode is what the user sees.
    """
    if not signalp_was_run:
        return existing_flag
    if existing_flag not in ("", "no_signalp"):
        return existing_flag
    if not signalp_prediction or signalp_prediction == "OTHER":
        return existing_flag
    if is_sec_signal_peptide(signalp_prediction):
        return existing_flag
    return "no_sec_signal"


def main():
    parser = argparse.ArgumentParser(description="Handle T5SS autotransporters")
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--proteins", required=True, help="Proteins FASTA for Pfam scan")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-substrates", required=True)
    parser.add_argument("--out-domains", required=True)
    args = parser.parse_args()

    t5ss_components = []
    with open(args.ss_components) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("ss_type", "").startswith("T5"):
                t5ss_components.append(row)

    predictions = {}
    with open(args.predictions) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            predictions[row["locus_tag"]] = row

    pfam_hits = scan_bundled_pfams(args.proteins) if t5ss_components else {}
    signalp_was_run = _signalp_was_run(predictions)

    substrates = []
    classifications = []
    for comp in t5ss_components:
        locus = comp["locus_tag"]
        ss_type = comp.get("ss_type", "T5aSS")
        pred = predictions.get(locus, {})
        sp_end = _parse_sp_end(pred.get("signalp_cs_position", ""))

        # Geometric filter only applies to T5aSS (T5bSS/T5cSS use their own
        # MacSyFinder profiles and don't fragment the same way).
        if ss_type == "T5aSS":
            domain_group, passenger_length, quality_flag = classify_t5a(pfam_hits.get(locus, {}), sp_end)
        else:
            domain_group, passenger_length = f"{ss_type}-component", 0
            quality_flag = "no_signalp" if (sp_end is None or sp_end <= 0) else ""

        quality_flag = _apply_sec_signal_gate(quality_flag, pred.get("signalp_prediction", ""), signalp_was_run)

        classifications.append(
            {
                "locus_tag": locus,
                "sample_id": args.sample,
                "ss_type": ss_type,
                "domain_group": domain_group,
                "t5_quality_flag": quality_flag,
                "passenger_length": passenger_length,
                "sp_end": sp_end if sp_end is not None else "",
                "barrel_start": pfam_hits.get(locus, {}).get("PF03797", ("", ""))[0],
            }
        )

        try:
            dlp_prob = float(pred.get("dlp_extracellular_prob", pred.get("extracellular_prob", 0)))
        except (ValueError, TypeError):
            dlp_prob = 0.0

        try:
            om_prob = float(pred.get("outer_membrane_prob", 0) or 0)
        except (ValueError, TypeError):
            om_prob = 0.0

        substrates.append(
            {
                "locus_tag": locus,
                "sample_id": args.sample,
                "tool": "T5SS-self",
                "nearby_ss_types": ss_type,
                "t5_quality_flag": quality_flag,
                "dlp_extracellular_prob": dlp_prob,
                "outer_membrane_prob": om_prob,
                "predicted_localization": pred.get("predicted_localization", ""),
                "dlp_max_localization": pred.get("dlp_max_localization", ""),
                "dlp_max_probability": pred.get("dlp_max_probability", ""),
                "dse_ss_type": pred.get("dse_ss_type", ""),
                "dse_max_prob": pred.get("dse_max_prob", ""),
                "signalp_prediction": pred.get("signalp_prediction", ""),
                "signalp_probability": pred.get("signalp_probability", ""),
                "signalp_cs_position": pred.get("signalp_cs_position", ""),
                "product": pred.get("product", ""),
            }
        )

    sub_fields = [
        "locus_tag",
        "sample_id",
        "tool",
        "nearby_ss_types",
        "t5_quality_flag",
        "dlp_extracellular_prob",
        "outer_membrane_prob",
        "predicted_localization",
        "dlp_max_localization",
        "dlp_max_probability",
        "dse_ss_type",
        "dse_max_prob",
        "signalp_prediction",
        "signalp_probability",
        "signalp_cs_position",
        "product",
    ]
    with open(args.out_substrates, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=sub_fields, delimiter="\t")
        writer.writeheader()
        for s in substrates:
            writer.writerow(s)

    domain_fields = [
        "locus_tag",
        "sample_id",
        "ss_type",
        "domain_group",
        "t5_quality_flag",
        "passenger_length",
        "sp_end",
        "barrel_start",
    ]
    with open(args.out_domains, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=domain_fields, delimiter="\t")
        writer.writeheader()
        for c in classifications:
            writer.writerow(c)

    n_flagged = sum(1 for s in substrates if s["t5_quality_flag"])
    logger.info(
        "Found %d T5SS self-substrates in %s "
        "(%d flagged: barrel-only / OMP / unclassified / no-signalp / no-sec-signal)",
        len(substrates),
        args.sample,
        n_flagged,
    )


if __name__ == "__main__":
    main()
