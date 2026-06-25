#!/usr/bin/env python3
"""Score annotation correctness for every emitted effector and plot it.

Each emitted effector was annotated by up to four ssign tools (InterProScan, EggNOG,
pLM-BLAST/ECOD, Pfam). For each tool we hand-grade its call against the protein's known
function (corpus family/quote + UniProt + the verified backfill in
annotation_groundtruth_fill.tsv), as one of:

  c = correct   tool names the right family / enzyme class
  p = partial   right fold/superfamily or one aspect, but not the specific function
  w = wrong     names an unrelated function (usually a remote structural false hit)
  n = none      tool produced no call
  u = untestable  no gradeable ground truth (uncharacterised protein, or broken corpus entry)

overall = best call across the four tools (correct if any tool is correct, etc.).

The proximity-panel effectors are graded from the checked-in sheet (locus-keyed verdicts, V);
the self-secreting T5SS effectors are sourced fresh from the tier-2 rerun by coordinate join and
graded gene-keyed (T5SS_V), since Bakta renames their loci. Verdicts are authored by hand (no
LLM/keyword auto-scoring). This script writes the per-tool + overall verdicts back into
annotation_accuracy_sheet.tsv and renders summary/07_annotation_correctness.png + 08 (by SS type).

    .venv/bin/python annotation_correctness.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
SHEET = HERE / "annotation_accuracy_sheet.tsv"
FIGDIR = HERE.parent / "data" / "phase2" / "figures" / "summary"

sys.path.insert(0, str(HERE.parent / "scripts"))
from rerun_coords import RerunIndex  # noqa: E402  (needs scripts/ on sys.path above)

# Self-secreting T5SS autotransporters live on a separate self-detection path, not the proximity
# panel the sheet is built from, so they are sourced directly from the tier-2 rerun by COORDINATE
# join (Bakta renames every locus, so we bridge corpus->rerun on genome position; see rerun_coords).
# Every T5SS corpus effector the rerun EMITTED is graded; the rerun's Bakta emits no Pfam IDs of its
# own, so the Pfam channel reads EggNOG's Pfam-domain names (same short-name format the fleet used).
T5SS_CORPUS = HERE.parent / "data" / "dataset" / "t5ss_effectors.tsv"
T5SS_ANN = {
    "ssign_interpro": "interpro_descriptions",
    "ssign_eggnog": "eggnog_description",
    "ssign_plmblast_ecod": "ecod_top1_description",
    "ssign_pfam": "eggnog__pfam_ids",
}

# Proximity-panel verdicts, ssign_locus -> (interpro, eggnog, plmblast, pfam, overall, basis).
# (Self-secreting T5SS are graded separately in T5SS_V below, keyed by gene.)
V = {
    "CC_1007": ("p", "p", "c", "p", "correct", "S-layer/RTX protein; ECOD names S-layer, others capture RTX repeats"),
    "CO076_RS11485": ("c", "c", "c", "c", "correct", "Serralysin M10 metalloprotease"),
    "CO076_RS11480": ("c", "c", "c", "c", "correct", "Serralysin M10 metalloprotease"),
    "WCU84_RS14645": ("c", "c", "c", "c", "correct", "Serralysin M10 metalloprotease"),
    "plu0655": ("c", "w", "c", "c", "correct", "Serralysin; EggNOG wrongly calls Serine 3-dehydrogenase"),
    "AACX37_RS13005": ("c", "w", "c", "c", "correct", "Serralysin; EggNOG wrongly calls dehydrogenase"),
    "PA1249": ("c", "p", "c", "c", "correct", "Serralysin (AprA); EggNOG names only the Ca-binding repeat"),
    "LPG_RS09620": ("c", "c", "c", "c", "correct", "Cellulase GH5"),
    "STM2885": ("c", "c", "w", "c", "correct", "T3SS translocon SipB; ECOD remote false hit (bacteriocin AS-48)"),
    "STM2884": ("c", "c", "p", "c", "correct", "T3SS translocon SipC; ECOD hits adjacent SipD"),
    "YE_RS17825": ("u", "u", "u", "u", "untestable", "YspE uncharacterised in its primary ref - no ground truth"),
    "ECs_4562": ("c", "n", "c", "c", "correct", "Map WxxxE/IpgB-like GEF effector"),
    "ECs_4564": ("n", "n", "w", "n", "wrong", "EspH; only call (ECOD Colicin S4) is wrong"),
    "ECs_4571": ("n", "c", "w", "c", "correct", "EspZ/SepZ; EggNOG+Pfam name SepZ; ECOD wrong (hemolysin)"),
    "E2348C_RS21140": ("c", "n", "c", "c", "correct", "Map WxxxE/IpgB-like GEF effector"),
    "E2348C_RS21150": ("n", "n", "w", "n", "wrong", "EspH; only call (ECOD KilA-N/SMT3) is wrong"),
    "ROD_RS14635": ("n", "n", "n", "n", "no_call", "NleB (arginine-GlcNAc transferase) - no annotation produced"),
    "CP0181": ("c", "c", "c", "c", "correct", "VirA/EspG cysteine protease"),
    "BPS_RS27065": ("c", "c", "w", "c", "correct", "T3SS translocon BipB; ECOD remote false hit"),
    "BPS_RS27060": ("c", "n", "p", "n", "correct", "T3SS translocon BipC; ECOD hits adjacent SipD"),
    "RS_RS21300": ("n", "n", "n", "n", "no_call", "RipJ (YopJ acetyltransferase) - no annotation produced"),
    "CT_089": ("c", "c", "c", "c", "correct", "CopN T3SS gatekeeper (YopN/LcrE family)"),
    "CT_868": ("c", "c", "c", "c", "correct", "ChlaDub1 deubiquitinase (Peptidase C48)"),
    "PA1510": ("c", "c", "c", "c", "correct", "Tle4/TplE alpha/beta-hydrolase phospholipase"),
    "VC_RS06860": ("c", "c", "c", "c", "correct", "TseL lipase"),
    "STM14_RS02020": ("n", "n", "n", "n", "no_call", "Tae4 amidase - no annotation produced"),
    "STM14_RS02030": ("n", "c", "c", "c", "correct", "Tlde1 L,D-transpeptidase (DUF2778)"),
    "ATU_RS20360": ("c", "c", "c", "c", "correct", "Tae amidase effector (Tae4 family)"),
    "ATU_RS20375": ("c", "p", "c", "p", "correct", "Tde1 DNase (Ntox15/His-Me); EggNOG+Pfam get RHS/phage delivery"),
    "EC042_RS24190": ("n", "n", "n", "n", "no_call", "Tle1 phospholipase - no annotation produced"),
    "PSEEN1550": ("c", "w", "c", "c", "correct", "Serralysin (AprA); EggNOG wrongly calls dehydrogenase"),
    "D1108_RS06660": ("c", "c", "c", "c", "correct", "ApxIIIA RTX toxin"),
    "BB0324": (
        "c",
        "p",
        "c",
        "p",
        "correct",
        "CyaA adenylate-cyclase toxin; InterPro+ECOD name cyclase, others get RTX",
    ),
    "BP0760": ("c", "p", "c", "p", "correct", "CyaA adenylate-cyclase toxin"),
    "BVG90_RS01585": ("c", "c", "c", "c", "correct", "Hemophore HasA"),
    "ACJ7Z2_RS03035": ("c", "c", "p", "c", "correct", "Leukotoxin RTX; ECOD beta-roll mislabeled 'Lipase'"),
    "B824_RS12000": ("c", "c", "p", "c", "correct", "Leukotoxin RTX"),
    "KKKWG1_RS01520": ("c", "c", "p", "c", "correct", "RtxA cytolysin RTX"),
    "PA0688": ("c", "c", "c", "c", "correct", "LapA alkaline phosphatase / PstS family"),
    "LPG_RS04595": (
        "n",
        "n",
        "w",
        "n",
        "wrong",
        "PlaA phospholipase A; only call (ECOD SecE/PKG) is a spurious 50-aa hit",
    ),
    "IS_RS19530": (
        "u",
        "u",
        "u",
        "u",
        "untestable",
        "Corpus broken: Q3BPB1 deleted; 'Serralysin' conflicts with unanimous S8-subtilisin calls",
    ),
    "BTH_RS25940": (
        "c",
        "w",
        "w",
        "w",
        "correct",
        "TseA/VasX; InterPro names VasX; EggNOG/ECOD/Pfam all wrong (dehydrogenase/lactoferrin)",
    ),
    "BTH_RS32260": ("c", "c", "c", "c", "correct", "Tle1 phospholipase (DUF2235)"),
    "BTH_RS09640": ("n", "n", "n", "n", "no_call", "TseM Mn-scavenging effector - no annotation produced"),
    "BTH_RS09645": (
        "p",
        "n",
        "w",
        "n",
        "partial",
        "TseZ Zn-scavenging effector; InterPro gets DUF6277 family, ECOD wrong",
    ),
}

# T5SS verdicts, gene-keyed (the rerun's Bakta loci are volatile, so we resolve them at runtime via the
# coordinate join). Hand-graded from the tier-2 rerun annotation strings for every emitted T5SS corpus
# effector: (interpro, eggnog, plmblast, pfam, overall, basis). Pfam channel = EggNOG Pfam-domain names.
T5SS_V = {
    "espP": (
        "c",
        "c",
        "c",
        "c",
        "correct",
        "SPATE serine-protease AT; InterPro+ECOD+Pfam name Peptidase_S6 + Autotransporter",
    ),
    "pic": (
        "c",
        "c",
        "c",
        "c",
        "correct",
        "Pic SPATE serine-protease/mucinase AT; all four name SPATE/IgA-protease + S6",
    ),
    "flu": ("c", "c", "c", "c", "correct", "Antigen 43 self-associating AT adhesin; ECOD names Ag43, others AIDA/AT"),
    "sat": ("c", "c", "c", "c", "correct", "Sat SPATE toxin; InterPro+ECOD name Peptidase_S6 + autotransporter"),
    "iga": (
        "c",
        "c",
        "c",
        "c",
        "correct",
        "IgA1 protease (founding classical AT); InterPro+ECOD name IgA1-protease + S6",
    ),
    "icsA": (
        "c",
        "c",
        "c",
        "c",
        "correct",
        "IcsA/VirG AT; InterPro names VirG insertion domain, ECOD names IcsA autotransporter",
    ),
    "yadA": (
        "c",
        "c",
        "c",
        "w",
        "correct",
        "YadA trimeric AT adhesin; 3 tools name YadA TAA, Pfam wrong (DUF814/FbpA)",
    ),
    "nadA": (
        "c",
        "c",
        "c",
        "p",
        "correct",
        "NadA YadA-like TAA; InterPro/EggNOG/ECOD name YadA anchor, Pfam dilutes YadA_anchor with noise",
    ),
    "sadA": ("c", "c", "c", "c", "correct", "SadA TAA; all four name YadA-like head/stalk/anchor"),
    "cdrA": (
        "c",
        "c",
        "c",
        "c",
        "correct",
        "CdrA TpsA exoprotein; all four name TPS / filamentous-haemagglutinin (Haemagg_act)",
    ),
    "hxuA": (
        "c",
        "c",
        "c",
        "c",
        "correct",
        "HxuA TpsA heme-haemopexin binding; ECOD names heme/hemopexin protein, Haemagg_act",
    ),
    "hmw1A": ("c", "c", "c", "c", "correct", "HMW1A TpsA adhesin; all four name TPS / Haemagg_act"),
    "hmw2A": ("c", "c", "c", "c", "correct", "HMW2A TpsA adhesin (HMW1A paralog); all four name TPS / Haemagg_act"),
}

TOOLS = ["interpro", "eggnog", "plmblast", "pfam"]
TOOL_LABEL = {"interpro": "InterProScan", "eggnog": "EggNOG", "plmblast": "pLM-BLAST\n(ECOD)", "pfam": "Pfam"}
CODE_ORDER = ["c", "p", "w", "n"]
CODE_LABEL = {"c": "correct", "p": "partial", "w": "wrong", "n": "no call"}
COLOR = {"c": "#2E7D6F", "p": "#E0A33B", "w": "#C0504D", "n": "#B7BCC2", "u": "#E2E2E2"}
OVR_ORDER = ["correct", "partial", "wrong", "no_call", "untestable"]
OVR_COLOR = {
    "correct": "#2E7D6F",
    "partial": "#E0A33B",
    "wrong": "#C0504D",
    "no_call": "#B7BCC2",
    "untestable": "#E2E2E2",
}


def apply_to_sheet():
    rows = list(csv.DictReader(open(SHEET), delimiter="\t"))
    cols = list(rows[0].keys())
    new = [
        "interpro_correct",
        "eggnog_correct",
        "plmblast_correct",
        "pfam_correct",
        "overall_correct",
        "correctness_basis",
    ]
    for c in new:
        if c not in cols:
            cols.append(c)
    missing = []
    for r in rows:
        # T5SS verdicts are gene-keyed (Bakta loci vary per rerun); proximity verdicts are locus-keyed.
        v = T5SS_V.get(r["gene"]) if r.get("ss_type") == "T5SS" else V.get(r["ssign_locus"])
        if not v:
            missing.append(r["ssign_locus"])
            continue
        ip, eg, pb, pf, ov, basis = v
        r.update(
            interpro_correct=ip,
            eggnog_correct=eg,
            plmblast_correct=pb,
            pfam_correct=pf,
            overall_correct=ov,
            correctness_basis=basis,
        )
    with open(SHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    if missing:
        print(f"  WARNING: no verdict for loci: {missing}")
    return rows


def figure(rows):
    n = len(rows)
    testable = [r for r in rows if r["overall_correct"] != "untestable"]
    n_test = len(testable)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})

    # Panel A: overall per-effector correctness
    counts = {k: sum(1 for r in rows if r["overall_correct"] == k) for k in OVR_ORDER}
    bottom = 0
    for k in OVR_ORDER:
        if counts[k] == 0:
            continue
        axA.bar(
            0, counts[k], bottom=bottom, color=OVR_COLOR[k], width=0.6, label=f"{k.replace('_', ' ')} ({counts[k]})"
        )
        if counts[k] >= 2:
            axA.text(
                0,
                bottom + counts[k] / 2,
                str(counts[k]),
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
                fontsize=13,
            )
        bottom += counts[k]
    axA.set_xlim(-0.6, 0.6)
    axA.set_xticks([])
    axA.set_ylabel("secreted proteins")
    pct = 100 * counts["correct"] / n_test
    axA.set_title(
        f"Overall annotation correctness\n{counts['correct']}/{n_test} correct ({pct:.0f}%) of gradeable secreted proteins",
        fontsize=12,
        fontweight="bold",
    )
    axA.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.05))
    axA.spines[["top", "right"]].set_visible(False)

    # Panel B: per-tool correctness (over gradeable effectors)
    x = np.arange(len(TOOLS))
    bottoms = np.zeros(len(TOOLS))
    for code in CODE_ORDER:
        vals = np.array([sum(1 for r in testable if r[f"{t}_correct"] == code) for t in TOOLS])
        axB.bar(x, vals, bottom=bottoms, color=COLOR[code], width=0.62, label=CODE_LABEL[code])
        for xi, v, b in zip(x, vals, bottoms):
            if v >= 2:
                axB.text(xi, b + v / 2, str(v), ha="center", va="center", color="white", fontweight="bold", fontsize=10)
        bottoms += vals
    axB.set_xticks(x)
    axB.set_xticklabels([TOOL_LABEL[t] for t in TOOLS], fontsize=10)
    axB.set_ylabel("secreted proteins")
    axB.set_ylim(0, n_test * 1.04)
    # correct-rate annotation above each bar
    for xi, t in zip(x, TOOLS):
        c = sum(1 for r in testable if r[f"{t}_correct"] == "c")
        called = sum(1 for r in testable if r[f"{t}_correct"] in ("c", "p", "w"))
        axB.text(
            xi,
            n_test * 1.0,
            f"{c}/{called} called\ncorrect" if called else "",
            ha="center",
            va="bottom",
            fontsize=8,
            color="0.3",
        )
    axB.set_title(
        f"Per-tool annotation correctness  (n={n_test} gradeable secreted proteins)", fontsize=12, fontweight="bold"
    )
    axB.legend(frameon=False, fontsize=9, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.08))
    axB.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "ssign annotation-tool correctness on recovered secreted proteins",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / "07_annotation_correctness.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # text summary
    print(f"\noverall (n={n}): " + ", ".join(f"{k}={counts[k]}" for k in OVR_ORDER))
    print(f"correct rate (gradeable n={n_test}): {pct:.0f}%")
    print("per-tool correct / called:")
    for t in TOOLS:
        c = sum(1 for r in testable if r[f"{t}_correct"] == "c")
        p = sum(1 for r in testable if r[f"{t}_correct"] == "p")
        w = sum(1 for r in testable if r[f"{t}_correct"] == "w")
        no = sum(1 for r in testable if r[f"{t}_correct"] == "n")
        print(f"  {TOOL_LABEL[t].replace(chr(10), ' '):16s}: correct {c}, partial {p}, wrong {w}, no-call {no}")


SSTYPE_ORDER = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS"]


def _fold3(overall):
    """3-way fold for the by-type view: partial -> wrong, untestable -> no annotation."""
    if overall == "correct":
        return "correct"
    if overall in ("wrong", "partial"):
        return "wrong"
    return "none"  # no_call + untestable


def figure_by_sstype(rows, poster=True):
    """Annotation correctness per secretion-system type, styled like 06_recall_systems."""
    types = [t for t in SSTYPE_ORDER if any(r["ss_type"] == t for r in rows)]
    cats = ["correct", "wrong", "none"]
    # figure-6 palette (correct=green, wrong=amber, no-annotation=slate)
    col = {"correct": "#2E7D6F", "wrong": "#E0A33B", "none": "#6C8EAD"}
    lab = {"correct": "correctly annotated", "wrong": "wrong annotation", "none": "no annotation"}
    tab = {t: {c: 0 for c in cats} for t in types}
    for r in rows:
        if r["ss_type"] in tab:
            tab[r["ss_type"]][_fold3(r["overall_correct"])] += 1
    totals = {t: sum(tab[t].values()) for t in types}
    n = sum(totals.values())
    n_correct = sum(tab[t]["correct"] for t in types)
    ymax = max(totals.values())

    sz = (
        dict(fig=(14, 8.5), title=27, ylab=24, xtick=22, seg=21, ann=19, leg=17, barw=0.66)
        if poster
        else dict(fig=(9.8, 5.0), title=12, ylab=10, xtick=10, seg=9, ann=8.5, leg=8.5, barw=0.62)
    )
    fig, ax = plt.subplots(figsize=sz["fig"])
    x = np.arange(len(types))
    bottom = np.zeros(len(types))
    for c in cats:
        vals = np.array([tab[t][c] for t in types])
        ax.bar(x, vals, bottom=bottom, color=col[c], width=sz["barw"], label=lab[c])
        for xi, v, b in zip(x, vals, bottom):
            if v >= 1:
                ax.text(
                    xi,
                    b + v / 2,
                    str(v),
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=sz["seg"],
                    fontweight="bold",
                )
        bottom += vals
    for xi, t in zip(x, types):
        ax.text(
            xi,
            totals[t] + ymax * 0.02,
            f"{tab[t]['correct']}/{totals[t]} correct",
            ha="center",
            va="bottom",
            fontsize=sz["ann"],
        )
    ax.set_xticks(x)
    ax.set_xticklabels([f"{t}\n(n={totals[t]})" for t in types], fontsize=sz["xtick"])
    ax.tick_params(axis="y", labelsize=sz["xtick"])
    ax.set_ylabel("secreted proteins", fontsize=sz["ylab"])
    ax.set_ylim(0, ymax * 1.25)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=sz["leg"])
    ax.set_title(
        f"ssign correctly annotated {n_correct}/{n} secreted proteins",
        fontsize=sz["title"],
        fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    name = "08_annotation_correctness_by_type_poster.png" if poster else "08_annotation_correctness_by_type.png"
    fig.savefig(FIGDIR / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}  (" + ", ".join(f"{t} {tab[t]['correct']}/{totals[t]}" for t in types) + ")")


def ensure_t5ss_from_rerun():
    """Rebuild the sheet's T5SS rows from the tier-2 rerun (coordinate-joined).

    Self-secreting T5SS are absent from the proximity panel the sheet is built from, so we add every
    T5SS corpus effector the rerun EMITTED as secreted, matched corpus->rerun on genome coordinates
    (Bakta renames the loci). The rerun annotates only emitted proteins, so non-emitted T5SS effectors
    (TPS substrates ssign did not recover) carry no annotation and are left out, logged for the recall
    reconciliation. Idempotent: existing T5SS rows are dropped and rebuilt every run."""
    rows = list(csv.DictReader(open(SHEET), delimiter="\t"))
    cols = list(rows[0].keys())
    rows = [r for r in rows if r.get("ss_type") != "T5SS"]  # rebuild T5SS cleanly
    if not T5SS_CORPUS.exists():
        print("  (T5SS corpus unavailable; figure will omit T5SS)")
        return
    idx = RerunIndex()
    graded, not_recovered, no_unit = 0, [], []
    for e in csv.DictReader(open(T5SS_CORPUS), delimiter="\t"):
        if not e.get("contig") or not e.get("start"):
            no_unit.append(e["gene"])
            continue
        j = idx.join(e["contig"], int(e["start"]), int(e["stop"]))
        if j is None:
            no_unit.append(e["gene"])
            continue
        if not j["emitted"]:
            not_recovered.append(e["gene"])
            continue
        new = {c: "" for c in cols}
        new.update(
            uniprot=e.get("uniprot", ""),
            gene=e["gene"],
            organism=e.get("organism", ""),
            ss_type="T5SS",
            known_quote=(e.get("uniprot_note") or e.get("quote", ""))[:300],
            genome=j["unit"],
            ssign_locus=j["locus_tag"],
        )
        for out_col, raw_col in T5SS_ANN.items():
            new[out_col] = (j["raw"].get(raw_col, "") or "").strip()
        rows.append(new)
        graded += 1
    with open(SHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(
        f"  T5SS from rerun: {graded} emitted+graded; {len(not_recovered)} not recovered "
        f"({', '.join(not_recovered)}); {len(no_unit)} not in panel ({', '.join(no_unit)})"
    )


def main():
    ensure_t5ss_from_rerun()
    rows = apply_to_sheet()
    figure(rows)
    figure_by_sstype(rows, poster=False)
    figure_by_sstype(rows, poster=True)


if __name__ == "__main__":
    main()
