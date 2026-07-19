#!/usr/bin/env python3
"""Generate ssign pipeline summary report (HTML + text).

Summarizes: genome count, protein count, SS system counts, substrate counts,
annotation tool coverage, and key parameters.
"""

import argparse
import logging
import os
from collections import Counter
from datetime import datetime

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Annotation tools grouped by the column prefixes each one contributes. A
# secreted protein counts as "covered" by a tool if ANY of that tool's columns
# is filled. Tools whose columns are absent (e.g. BLASTp/HH-suite off-tier) are
# simply not listed. Order = display order.
_TOOL_COLUMN_PREFIXES = [
    ("SignalP", ("signalp_",)),
    ("InterProScan", ("interpro_",)),  # its Pfam col is interpro_pfam_ids
    ("pLM-BLAST (ECOD)", ("ecod_",)),
    ("EggNOG", ("eggnog_", "cog_", "kegg_", "pfam_ids")),  # bare pfam_ids is EggNOG's
    ("BLASTp", ("blastp_",)),
    ("HH-suite", ("pdb_top1_", "pfam_top1_")),  # parse_hhr emits {pfam,pdb}_top1_*
]


def _summarize_enrichment(path, n_genomes=1):
    """Distil the per-type enrichment test into a short human summary.

    The full per-(type, tool, mode) table is kept verbatim in the standalone
    ``*_enrichment_stats.tsv``; this reduces it to the headline COMBINED-predictor
    track -- one fold + BH-q line per secretion-system type/mode, most significant
    first. Returns report lines (empty if the file can't be summarised).
    """
    import sys

    _d = os.path.dirname(os.path.abspath(__file__))
    if _d not in sys.path:
        sys.path.insert(0, _d)
    from ssign_lib.constants import ENRICH_COMBINED_TOOL, ENRICH_MODE_SELF, T5_HITCH_TAG, T5_SELF_TAG

    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return []
    if df.empty or "tool" not in df.columns or "ss_type" not in df.columns:
        return []
    comb = df[df["tool"] == ENRICH_COMBINED_TOOL]
    if comb.empty:
        return []

    has_mode = "mode" in comb.columns
    types_with_self = set(comb.loc[comb["mode"] == ENRICH_MODE_SELF, "ss_type"].astype(str)) if has_mode else set()

    rows = []
    for _, r in comb.iterrows():
        try:
            fold, q = float(r["fold"]), float(r["qvalue"])
        except (KeyError, ValueError, TypeError):
            continue
        st = str(r["ss_type"])
        mode = str(r["mode"]) if has_mode else "window"
        if mode == ENRICH_MODE_SELF:
            lbl = f"{st} {T5_SELF_TAG}"
        elif st in types_with_self:
            lbl = f"{st} {T5_HITCH_TAG}"
        else:
            lbl = st
        rows.append((lbl, fold, q, q < 0.05))
    if not rows:
        return []

    rows.sort(key=lambda t: (t[2], -t[1]))  # BH q ascending, then fold descending
    width = max(len(lbl) for lbl, *_ in rows)
    n_sig = sum(1 for *_, sig in rows if sig)

    out = ["Secretion-system enrichment (per-type permutation test, COMBINED predictor):"]
    for lbl, fold, q, sig in rows:
        q_str = "q<0.001" if q < 0.001 else f"q={q:.3f}"
        out.append(f"  {lbl:<{width}}  {fold:>5.1f}x   {q_str}{'  *' if sig else ''}")
    if n_sig:
        out.append("  * enriched at Benjamini-Hochberg q < 0.05")
    elif n_genomes <= 1:
        out.append("  (nothing reached q < 0.05; single-genome power is low -- pool genomes)")
    else:
        out.append("  (nothing reached q < 0.05)")
    return out


def generate_text_report(master_csvs, enrichment_file, output_path, tier=""):
    """Generate plain text summary report."""
    # Load and combine all master CSVs
    dfs = []
    for f in master_csvs:
        try:
            dfs.append(pd.read_csv(f))
        except Exception as e:
            logger.warning(f"Could not read {f}: {e}")

    if not dfs:
        with open(output_path, "w") as f:
            f.write("ssign Report\n============\nNo results to report.\n")
        return

    df = pd.concat(dfs, ignore_index=True)
    n = len(df)

    # ── Run metadata line ──
    try:
        from ssign_app import __version__
    except Exception:
        __version__ = "unknown"
    genomes = sorted(df["sample_id"].dropna().astype(str).unique()) if "sample_id" in df.columns else []
    meta_bits = [f"ssign v{__version__}"]
    if tier:
        meta_bits.append(f"{tier} tier")
    meta_bits.append(datetime.now().strftime("%Y-%m-%d %H:%M"))

    lines = [
        "=" * 60,
        "ssign — Secretion-System Identification for Gram Negatives",
        "=" * 60,
        "",
        "  |  ".join(meta_bits),
    ]
    if genomes:
        shown = ", ".join(genomes[:6]) + (f" (+{len(genomes) - 6} more)" if len(genomes) > 6 else "")
        lines.append(f"Genomes ({len(genomes)}): {shown}")
    lines += ["", f"Secreted proteins found: {n}", ""]

    # SS type distribution
    if "nearby_ss_types" in df.columns:
        ss_counts: Counter[str] = Counter()
        for val in df["nearby_ss_types"].dropna():
            for ss in str(val).split(","):
                ss = ss.strip()
                if ss:
                    ss_counts[ss] += 1

        if ss_counts:
            lines.append("Secreted proteins per secretion system type:")
            for ss, count in ss_counts.most_common():
                lines.append(f"  {ss}: {count}")
            lines.append("")

    # Annotation coverage, grouped by tool (only tools whose columns are present)
    cov = []
    for tool, prefixes in _TOOL_COLUMN_PREFIXES:
        cols = [c for c in df.columns if c.startswith(prefixes)]
        if cols:
            covered = int(df[cols].notna().any(axis=1).sum())
            cov.append((tool, covered))
    if cov:
        lines.append("Annotation coverage:")
        width = max(len(t) for t, _ in cov)
        for tool, covered in cov:
            pct = 100 * covered / max(n, 1)
            lines.append(f"  {tool:<{width}}  {covered}/{n} ({pct:.0f}%)")
        lines.append("")

    # Enrichment summary -- distilled headline; the full per-tool table stays in
    # the standalone *_enrichment_stats.tsv.
    if enrichment_file and os.path.exists(enrichment_file):
        enrich_lines = _summarize_enrichment(enrichment_file, n_genomes=len(genomes))
        if enrich_lines:
            lines += enrich_lines + [""]

    lines.append("=" * 60)

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    logger.info(f"Wrote text report to {output_path}")


def generate_html_report(master_csvs, enrichment_file, output_path):
    """Generate HTML summary report."""
    # Simple HTML wrapping the text report
    text_path = output_path.replace(".html", ".txt")
    if os.path.exists(text_path):
        with open(text_path) as f:
            text_content = f.read()
    else:
        text_content = "No results."

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>ssign Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, sans-serif; max-width: 900px; margin: 40px auto; padding: 0 20px; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        pre {{ background: #f8f9fa; padding: 20px; border-radius: 5px; overflow-x: auto; }}
        .footer {{ color: #7f8c8d; font-size: 0.85em; margin-top: 40px; }}
    </style>
</head>
<body>
    <h1>ssign Report</h1>
    <p>Secretion-System Identification for Gram Negatives</p>
    <pre>{text_content}</pre>
    <div class="footer">
        <p>Generated by ssign pipeline</p>
    </div>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(html)

    logger.info(f"Wrote HTML report to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate ssign report")
    parser.add_argument("--master-csvs", nargs="+", required=True)
    parser.add_argument("--enrichment", default="")
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--out-txt", required=True)
    parser.add_argument("--tier", default="", help="Install tier (shown in the report metadata line).")
    args = parser.parse_args()

    generate_text_report(args.master_csvs, args.enrichment, args.out_txt, tier=args.tier)
    generate_html_report(args.master_csvs, args.enrichment, args.out_html)


if __name__ == "__main__":
    main()
