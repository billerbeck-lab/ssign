#!/usr/bin/env python3
"""ssign — Streamlit GUI Home page.

Launch with: ssign (after pip install) or streamlit run Home.py
"""

import os
import shutil
import tempfile
from pathlib import Path

import streamlit as st

from ssign_app.core.runner import PipelineConfig, PipelineRunner, StepResult

# ─────────────────────────────────────────────────────────────────────
# Page config
# ─────────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="ssign",
    page_icon="\U0001f9ec",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get help": "https://github.com/billerbeck-lab/ssign",
        "Report a Bug": "https://github.com/billerbeck-lab/ssign/issues",
        "About": "ssign — Secretion-system Identification for Gram Negatives",
    },
)

# Dark-mode compatible styling — uses Streamlit CSS variables
st.markdown(
    """
<style>
    .stDeployButton { display: none; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        padding: 8px 20px;
        border-radius: 6px 6px 0 0;
    }
    div[data-testid="stMetric"] {
        background-color: var(--secondary-background-color, #F0F5F8);
        padding: 12px 16px;
        border-radius: 8px;
        border-left: 4px solid #1B6B93;
    }
    .tool-section {
        padding: 12px 0;
        border-bottom: 1px solid var(--secondary-background-color, #E8EEF2);
    }
</style>
""",
    unsafe_allow_html=True,
)

# Custom connection error — use setInterval to periodically check for and
# replace Streamlit's generic error message. MutationObserver can fail if the
# iframe or observer gets destroyed during connection loss.
import streamlit.components.v1 as components

components.html(
    """
<script>
(function() {
    var doc;
    try { doc = window.parent.document; } catch(e) { doc = document; }
    setInterval(function() {
        try {
            var walker = doc.createTreeWalker(doc.body, NodeFilter.SHOW_TEXT);
            while (walker.nextNode()) {
                var n = walker.currentNode;
                if (n.nodeValue.indexOf('Is Streamlit still running') > -1) {
                    n.nodeValue = 'ssign server stopped. To restart, run: ssign';
                }
                if (n.nodeValue.indexOf('streamlit run') > -1) {
                    n.nodeValue = n.nodeValue.replace(/streamlit run [^\\s]+/g, 'ssign');
                }
            }
        } catch(e) {}
    }, 500);
})();
</script>
""",
    height=0,
)

# ─────────────────────────────────────────────────────────────────────
# Session state defaults
# ─────────────────────────────────────────────────────────────────────

if "running" not in st.session_state:
    st.session_state.running = False
if "results" not in st.session_state:
    st.session_state.results = []
if "uploaded_files_data" not in st.session_state:
    st.session_state.uploaded_files_data = []


def _has_previous_progress(outdir_val: str) -> bool:
    """Check if any progress files exist in the output directory."""
    if not outdir_val or not os.path.isdir(outdir_val):
        return False
    import glob

    ssign_dir = os.path.join(outdir_val, ".ssign")
    return bool(
        glob.glob(os.path.join(ssign_dir, "*_progress.json"))
        # Legacy locations
        or glob.glob(os.path.join(outdir_val, "*_progress.json"))
        or os.path.exists(os.path.join(outdir_val, "ssign_progress.json"))
    )


def _needs_run_mode_gate() -> bool:
    """Return True if a previous run exists but the user hasn't chosen a run mode yet."""
    outdir_val = st.session_state.get(
        "outdir_input",
        os.path.join(os.path.expanduser("~"), "ssign_results"),
    )
    if _has_previous_progress(outdir_val):
        return not bool(st.session_state.get("run_mode_choice"))
    return False


_GATE_MSG = (
    "A previous run was detected in the output directory. "
    "Please go to the **Upload & Configure** tab and choose how to proceed "
    "(Resume / Start fresh / Selective rerun) before continuing."
)


def _merge_genome_outputs(outdir: str, sample_names: list[str]):
    """Merge per-genome output files into single combined files.

    Produces:
      ssign_results.csv     — Normal (chunked: secreted proteins, their SS, other SS)
      ssign_results_raw.csv — Raw (all tool data, no filtering)
      ssign_summary.txt     — Combined summary
      figures/              — Combined figures
    """
    import pandas as pd

    # ── 1. Merge normal results CSVs (chunked format) ──
    # The per-genome CSVs have comment-header chunks separated by blank lines
    all_chunks = {"proteins": [], "ss_with": [], "ss_other": []}
    for sid in sample_names:
        csv_path = os.path.join(outdir, f"{sid}_results.csv")
        if not os.path.exists(csv_path):
            continue

        with open(csv_path) as f:
            content = f.read()

        # Parse chunks by comment headers
        current_chunk = None
        current_lines = []
        for line in content.split("\n"):
            if line.startswith("# Secreted Proteins"):
                if current_chunk and current_lines:
                    all_chunks[current_chunk].append("\n".join(current_lines))
                current_chunk = "proteins"
                current_lines = []
            elif line.startswith("# Secretion Systems (with"):
                if current_chunk and current_lines:
                    all_chunks[current_chunk].append("\n".join(current_lines))
                current_chunk = "ss_with"
                current_lines = []
            elif line.startswith("# Secretion Systems (other"):
                if current_chunk and current_lines:
                    all_chunks[current_chunk].append("\n".join(current_lines))
                current_chunk = "ss_other"
                current_lines = []
            elif line.strip():
                current_lines.append(line)
        if current_chunk and current_lines:
            all_chunks[current_chunk].append("\n".join(current_lines))

        os.remove(csv_path)

    # Concatenate each chunk type across genomes
    import io as _io

    merged = {}
    for key, blocks in all_chunks.items():
        dfs = []
        for block in blocks:
            if block.strip():
                try:
                    dfs.append(pd.read_csv(_io.StringIO(block)))
                except Exception:
                    pass
        if dfs:
            merged[key] = pd.concat(dfs, ignore_index=True)

    combined_csv = os.path.join(outdir, "ssign_results.csv")
    with open(combined_csv, "w", newline="") as f:
        written = False
        if "proteins" in merged and not merged["proteins"].empty:
            f.write("# Secreted Proteins\n")
            merged["proteins"].to_csv(f, index=False)
            written = True
        if "ss_with" in merged and not merged["ss_with"].empty:
            if written:
                f.write("\n")
            f.write("# Secretion Systems (with secreted proteins)\n")
            merged["ss_with"].to_csv(f, index=False)
            written = True
        if "ss_other" in merged and not merged["ss_other"].empty:
            if written:
                f.write("\n")
            f.write("# Secretion Systems (other)\n")
            merged["ss_other"].to_csv(f, index=False)

    # ── 2. Merge raw results CSVs ──
    raw_dfs = []
    for sid in sample_names:
        raw_path = os.path.join(outdir, f"{sid}_results_raw.csv")
        if os.path.exists(raw_path):
            try:
                raw_dfs.append(pd.read_csv(raw_path))
            except Exception:
                pass
            os.remove(raw_path)
    if raw_dfs:
        pd.concat(raw_dfs, ignore_index=True).to_csv(
            os.path.join(outdir, "ssign_results_raw.csv"), index=False
        )

    # ── 2. Merge summary texts ──
    summary_parts = []
    for sid in sample_names:
        txt_path = os.path.join(outdir, f"{sid}_summary.txt")
        if not os.path.exists(txt_path):
            continue
        with open(txt_path) as f:
            content = f.read()
        if len(sample_names) > 1:
            summary_parts.append(f"\n{'=' * 60}\n  {sid}\n{'=' * 60}\n\n{content}")
        else:
            summary_parts.append(content)
        os.remove(txt_path)

    if summary_parts:
        with open(os.path.join(outdir, "ssign_summary.txt"), "w") as f:
            # Cross-genome summary header
            if (
                len(sample_names) > 1
                and "proteins" in merged
                and not merged["proteins"].empty
            ):
                df_all = merged["proteins"]
                n_genomes = (
                    df_all["sample_id"].nunique()
                    if "sample_id" in df_all.columns
                    else len(sample_names)
                )
                n_total = len(df_all)
                avg_per_genome = n_total / n_genomes if n_genomes > 0 else 0

                f.write("=" * 60 + "\n")
                f.write("  CROSS-GENOME SUMMARY\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Genomes analyzed: {n_genomes}\n")
                f.write(f"Total secreted proteins: {n_total}\n")
                f.write(f"Average per genome: {avg_per_genome:.1f}\n\n")

                if "nearby_ss_types" in df_all.columns:
                    f.write("Secreted proteins by SS type (total across genomes):\n")
                    for ss, count in df_all["nearby_ss_types"].value_counts().items():
                        f.write(f"  {ss}: {count}\n")
                    f.write("\n")

                if "broad_annotation" in df_all.columns:
                    non_empty = df_all["broad_annotation"].dropna()
                    non_empty = non_empty[non_empty != ""]
                    if not non_empty.empty:
                        f.write("Most common annotations:\n")
                        for ann, count in non_empty.value_counts().head(10).items():
                            f.write(f"  {ann}: {count}\n")
                        f.write("\n")

                f.write("\n")

            f.write("\n".join(summary_parts))

    # ── 3. Regenerate combined figures from merged data ──
    # Remove per-genome figure subdirectories
    fig_base = os.path.join(outdir, "figures")
    if os.path.isdir(fig_base):
        shutil.rmtree(fig_base)

    # Use the raw CSV (clean tabular data, no comment headers) for figure generation
    raw_csv = os.path.join(outdir, "ssign_results_raw.csv")
    if os.path.exists(raw_csv):
        try:
            from ssign_app.core.runner import run_script

            run_script(
                "generate_figures.py",
                [
                    "--master-csvs",
                    raw_csv,
                    "--outdir",
                    os.path.join(outdir, "figures"),
                    "--dpi",
                    "300",
                ],
            )
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────
# Header
# ─────────────────────────────────────────────────────────────────────

st.markdown(
    '<h1 style="margin-bottom: 0;">'
    '<span style="color: #1B6B93;">S</span>ecretion-<span style="color: #1B6B93;">s</span>ystem '
    '<span style="color: #1B6B93;">I</span>dentification for '
    '<span style="color: #1B6B93;">G</span>ram <span style="color: #1B6B93;">N</span>egatives</h1>',
    unsafe_allow_html=True,
)
st.divider()


# ─────────────────────────────────────────────────────────────────────
# Sidebar
# ─────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown(
        '<h2><span style="color: #1B6B93;">ss</span>i<span style="color: #1B6B93;">g</span>'
        '<span style="color: #1B6B93;">n</span></h2>',
        unsafe_allow_html=True,
    )

    st.markdown(
        "ssign links a variety of bioinformatic tools together primarily to identify "
        "secretion systems & secreted protein pairs in gram-negative bacteria."
    )
    st.markdown(
        "It can be run on single genomes or groupings of organisms and optional "
        "functional annotation tools can be run along with premade figure generators "
        "to more easily understand the data."
    )
    st.markdown(
        "All tools being used are accredited below and should be cited along with "
        "this tool if used."
    )

    st.divider()

    st.subheader("Pipeline Overview")

    st.markdown("**1. Secretion System Identification**")
    st.caption(
        "Detects secretion systems and their components in the genome using "
        "MacSyFinder v2 with TXSScan models."
    )

    st.markdown("**2. Secreted Protein Identification**")
    st.caption(
        "Identifies candidate secreted proteins using guilt by association — "
        "proteins genomically proximal to secretion system components are evaluated. "
        "DeepLocPro predicts localization of secretion system components (secondary "
        "verification) and nearby proteins for possible secretion. DeepSecE provides "
        "a secondary prediction of secreted proteins, good for enrichment and verification."
    )

    st.markdown("**3. Improve Secreted Protein Annotations**")
    st.caption(
        "Many secreted proteins in lesser-studied organisms are classified as "
        "hypothetical or domain of unknown function (DUF). A suite of optional "
        "annotation tools (BLASTp, HHpred, InterProScan, ProtParam) "
        "can help resolve some of these."
    )

    st.markdown("**4. Generate Data & Figures**")
    st.caption(
        "Produces summary figures and result tables for the identified "
        "secreted proteins and their annotations."
    )

    st.divider()

    st.info(
        "Annotation tools run locally against databases fetched from Zenodo "
        "(see `scripts/fetch_databases.sh`). Pick an install tier — **base** "
        "(~17 GB) for prediction only, **extended** (~130 GB) adds EggNOG + "
        "HH-suite + InterProScan + pLM-BLAST, **full** (~630 GB) adds BLAST NR. "
        "See the README for details."
    )

    st.divider()

    # ── Citations ──
    with st.expander("Tool Citations", expanded=False):
        st.markdown(
            "**MacSyFinder v2:** Neron, B., et al. (2023). MacSyFinder v2: Improved "
            "modelling and search engine to identify molecular systems in genomes. "
            "*Peer Community Journal*, 3, e28. doi:10.24072/pcjournal.250\n\n"
            "**TXSScan:** Abby, S.S., et al. (2016). Identification of protein secretion "
            "systems in bacterial genomes. *Scientific Reports*, 6, 23080. "
            "doi:10.1038/srep23080\n\n"
            "**DeepLocPro:** Moreno, J., et al. (2024). Predicting the subcellular "
            "location of prokaryotic proteins with DeepLocPro. *Bioinformatics*, "
            "40(12). doi:10.1093/bioinformatics/btae677\n\n"
            "**DeepSecE:** Zhang, Y., et al. (2023). DeepSecE: A Deep-Learning-Based "
            "Framework for Multiclass Prediction of Secreted Proteins in Gram-Negative "
            "Bacteria. *Research*, 6, 0258. doi:10.34133/research.0258\n\n"
            "**SignalP 6.0:** Teufel, F., et al. (2022). SignalP 6.0 predicts all five "
            "types of signal peptides using protein language models. *Nature "
            "Biotechnology*, 40(7), 1023-1025. doi:10.1038/s41587-021-01156-3\n\n"
            "**BLAST+:** Camacho, C., et al. (2009). BLAST+: architecture and "
            "applications. *BMC Bioinformatics*, 10, 421. doi:10.1186/1471-2105-10-421\n\n"
            "**HH-suite3:** Steinegger, M., et al. (2019). HH-suite3 for fast remote "
            "homology detection and deep protein annotation. *BMC Bioinformatics*, "
            "20(1), 473. doi:10.1186/s12859-019-3019-7\n\n"
            "**InterProScan 5:** Jones, P., et al. (2014). InterProScan 5: genome-scale "
            "protein function classification. *Bioinformatics*, 30(9), 1236-1240. "
            "doi:10.1093/bioinformatics/btu031\n\n"
            "**Pyrodigal:** Larralde, M. (2022). Pyrodigal: Python bindings and "
            "interface to Prodigal, an efficient method for gene prediction in "
            "prokaryotes. *Journal of Open Source Software*, 7(72), 4296. "
            "doi:10.21105/joss.04296\n\n"
            "**Bakta:** Schwengers, O., et al. (2021). Bakta: rapid and standardized "
            "annotation of bacterial genomes via alignment-free sequence identification. "
            "*Microbial Genomics*, 7(11). doi:10.1099/mgen.0.000685\n\n"
            "**Biopython:** Cock, P.J.A., et al. (2009). Biopython: freely available "
            "Python tools for computational molecular biology and bioinformatics. "
            "*Bioinformatics*, 25(11), 1422-1423. doi:10.1093/bioinformatics/btp163"
        )

    st.divider()
    from ssign_app import __version__

    st.caption(f"ssign v{__version__} | GPLv3 | Billerbeck Lab")
    st.markdown(
        "[GitHub](https://github.com/billerbeck-lab/ssign) · "
        "[Report a bug](https://github.com/billerbeck-lab/ssign/issues)"
    )


# ─────────────────────────────────────────────────────────────────────
# Tabs
# ─────────────────────────────────────────────────────────────────────

tab_upload, tab_pipeline, tab_run = st.tabs(
    ["Upload & Configure", "Pipeline Overview", "Run & Results"]
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TAB 1: Upload & Configure
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

with tab_upload:
    st.subheader("Upload Genome Files")

    uploaded_files = st.file_uploader(
        "Choose genome file(s)",
        type=["gbff", "gbk", "gb", "fasta", "fna", "fa", "faa"],
        help=(
            "GenBank (.gbff/.gbk/.gb), FASTA contigs (.fasta/.fna/.fa), or "
            "protein FASTA (.faa). Upload multiple files for batch processing. "
            "GenBank files are recommended as they include gene names and "
            "functional annotations. GFF3 files require a companion FASTA "
            "and are only supported via the CLI (`ssign run`)."
        ),
        accept_multiple_files=True,
    )

    if uploaded_files:
        # Clear resolved taxonomy if files changed
        current_names = sorted(f.name for f in uploaded_files)
        prev_names = st.session_state.get("_last_uploaded_names", [])
        if current_names != prev_names:
            st.session_state.pop("per_genome_taxonomy", None)
            st.session_state.pop("organism_resolved", None)
            st.session_state.pop("detected_organism", None)
            st.session_state["_last_uploaded_names"] = current_names

        st.session_state.uploaded_files_data = uploaded_files

        # Derive sample names internally (not displayed)
        sample_names = []
        for f in uploaded_files:
            sn = f.name
            for suffix in [".gbff", ".gbk", ".gb", ".fasta", ".fna", ".fa", ".faa"]:
                sn = sn.replace(suffix, "")
            sn = sn.replace("_genomic", "")
            sample_names.append(sn)
        if len(uploaded_files) == 1:
            st.session_state.sample_name = sample_names[0]
        st.session_state.sample_names = sample_names

        st.success(
            f"**{len(uploaded_files)} genome(s) loaded:** "
            + ", ".join(f.name for f in uploaded_files)
        )

        # Extract organism names from GenBank files (local parsing only, no API)
        # Taxonomy resolution (NCBI API) is deferred to when the user configures
        # BLASTp exclusion, to avoid blocking the GUI on upload.
        if "per_genome_organisms" not in st.session_state:
            per_genome_orgs = {}
            for f in uploaded_files:
                ext_lower = Path(f.name).suffix.lower()
                if ext_lower not in (".gbff", ".gbk", ".gb"):
                    continue
                try:
                    import io

                    from Bio import SeqIO

                    f.seek(0)
                    content = f.read().decode("utf-8", errors="replace")
                    f.seek(0)
                    record = next(SeqIO.parse(io.StringIO(content), "genbank"))
                    org_name = record.annotations.get("organism", "").strip()
                    if not org_name:
                        org_name = record.annotations.get("source", "").strip()
                    if not org_name or len(org_name.split()) < 2:
                        for feat in record.features:
                            if feat.type == "source":
                                src_org = feat.qualifiers.get("organism", [""])[
                                    0
                                ].strip()
                                if src_org and len(src_org.split()) >= 2:
                                    org_name = src_org
                                break
                    if not org_name or len(org_name.split()) < 2:
                        stem = Path(f.name).stem
                        for sfx in ("_genomic", "_protein", "_cds", "_rna"):
                            stem = stem.replace(sfx, "")
                        parts = stem.replace("_", " ").split()
                        if (
                            len(parts) >= 2
                            and parts[0][0].isupper()
                            and parts[1][0].islower()
                            and parts[1].isalpha()
                        ):
                            org_name = f"{parts[0]} {parts[1]}"
                    if org_name:
                        per_genome_orgs[f.name] = org_name
                except Exception:
                    pass
            st.session_state.per_genome_organisms = per_genome_orgs

    st.divider()

    col1, col2 = st.columns(2)
    with col1:
        st.info(
            "**Recommended:** Use GenBank (.gbff) files, which include gene names and "
            "functional annotations from tools like Bakta.\n\n"
            "Raw FASTA contigs (.fasta/.fna/.fa) will use Prodigal for gene prediction "
            "(no additional setup) or Bakta (richer annotation, requires database download). "
            "Protein FASTA (.faa) files are read directly.\n\n"
            "**Note:** GFF3/GTF files are not supported in the GUI because they require "
            "a companion FASTA file. Use GenBank format instead, or the CLI "
            "(`ssign run genome.gff --fasta genome.fna ...`) for GFF3 input."
        )
        bakta_available = shutil.which("bakta") is not None
        use_bakta = st.checkbox(
            "Use Bakta for ORF prediction (raw FASTA input only)",
            value=False,
            key="use_bakta",
            disabled=not bakta_available,
            help="Bakta provides richer genome annotation than Prodigal, including "
            "gene names and functional descriptions. Only needed for raw FASTA input.",
        )
        # Phase 3.3.c: GenBank input is re-annotated through Bakta by
        # default for uniform annotation across a cohort. Original
        # annotations are kept as `gbff_annotation` for consensus voting.
        st.checkbox(
            "Use input GenBank annotations as-is (skip Bakta re-annotation)",
            value=False,
            key="use_input_annotations",
            disabled=not bakta_available,
            help="By default, GenBank input is re-annotated by Bakta so all genomes "
            "in a cohort share the same annotation pipeline. Original product strings "
            "are preserved as a `gbff_annotation` column. Tick this box to skip "
            "re-annotation when you trust the input annotations (e.g. hand-curated GenBank).",
        )
        if bakta_available:
            if use_bakta or not st.session_state.get("use_input_annotations", False):
                st.text_input(
                    "Bakta database path",
                    key="bakta_db_path",
                    placeholder="/path/to/bakta_db",
                )
        else:
            st.info(
                "**Bakta is not installed.** Required for GenBank re-annotation "
                "(default) and optional for raw FASTA input.\n\n"
                "Install (~2 GB): `pip install ssign[bakta]` then download the database.\n\n"
                "See the [installation guide](https://github.com/billerbeck-lab/ssign/blob/main/docs/optional_tools.md) "
                "for detailed instructions."
            )

    with col2:
        outdir_default = os.path.join(os.path.expanduser("~"), "ssign_results")
        outdir = st.text_input(
            "Output directory",
            value=outdir_default,
            key="outdir_input",
            help="All results, figures, and intermediate files will be saved here. "
            "The directory will be created if it doesn't exist.",
        )
        # Output directory validation
        if outdir:
            expanded = os.path.expanduser(outdir)
            if os.path.isdir(expanded):
                st.caption(f"Directory exists: `{expanded}`")
            else:
                parent = os.path.dirname(expanded)
                if os.path.isdir(parent):
                    st.caption(f"Will be created: `{expanded}`")
                else:
                    st.warning(f"Parent directory not found: `{parent}`")

    # ── Previous run detection ──
    outdir_val = st.session_state.get("outdir_input", "")
    if outdir_val and _has_previous_progress(outdir_val):
        try:
            import glob as _glob
            import json

            # Collect per-genome progress files from .ssign/ subdirectory
            ssign_dir = os.path.join(outdir_val, ".ssign")
            prog_files = sorted(_glob.glob(os.path.join(ssign_dir, "*_progress.json")))
            # Legacy: progress files in outdir root
            if not prog_files:
                prog_files = sorted(
                    _glob.glob(os.path.join(outdir_val, "*_progress.json"))
                )
            if not prog_files:
                legacy = os.path.join(outdir_val, "ssign_progress.json")
                if os.path.exists(legacy):
                    prog_files = [legacy]

            total_done = 0
            total_steps = 0
            genome_summaries = []
            latest_time = "unknown"
            saved_config = {}
            for pf_path in prog_files:
                with open(pf_path) as pf:
                    prev = json.load(pf)
                prev_steps = prev.get("steps", [])
                n_done = sum(1 for s in prev_steps if s.get("success"))
                total_done += n_done
                total_steps += len(prev_steps)
                sid = prev.get("sample_id", "unknown")
                latest_time = prev.get("timestamp", latest_time)
                status = (
                    "complete"
                    if n_done == len(prev_steps)
                    else f"{n_done}/{len(prev_steps)} steps"
                )
                genome_summaries.append(f"- **{sid}**: {status}")
                # Keep the first genome's config as representative
                if not saved_config and prev.get("config"):
                    saved_config = prev["config"]

            st.divider()
            st.markdown("#### Previous Run Detected")
            col_info, col_action = st.columns([3, 2])
            with col_info:
                n_genomes = len(prog_files)
                summary = "\n".join(genome_summaries[:10])
                if len(genome_summaries) > 10:
                    summary += f"\n- ...and {len(genome_summaries) - 10} more"
                st.info(
                    f"**Genomes:** {n_genomes}  \n"
                    f"**Progress:** {total_done}/{total_steps} steps completed  \n"
                    f"**Last run:** {latest_time}\n\n"
                    f"{summary}"
                )
            with col_action:
                run_mode = st.radio(
                    "How to proceed?",
                    [
                        "Resume (skip completed steps)",
                        "Start fresh (rerun everything)",
                        "Selective rerun (choose steps)",
                    ],
                    key="run_mode_choice",
                    index=0,
                )

            # Restore saved tool settings when resuming
            if saved_config and run_mode and "Resume" in run_mode:
                _cfg_to_session = {
                    "skip_hhsuite": ("run_hh", lambda v: not v),
                    "skip_blastp": ("run_blastp", lambda v: not v),
                    "skip_interproscan": ("run_iprs", lambda v: not v),
                    "skip_signalp": ("run_signalp", lambda v: not v),
                    "skip_deepsece": ("run_deepsece", lambda v: not v),
                    "skip_protparam": ("run_protparam", lambda v: not v),
                    "conf_threshold": ("conf_threshold", None),
                    "wholeness_threshold": ("wholeness_threshold", None),
                    "proximity_window": ("proximity_window", None),
                    "blastp_min_pident": ("blastp_min_pident", None),
                    "blastp_evalue": ("blastp_evalue", None),
                }
                restored = []
                for cfg_key, (session_key, transform) in _cfg_to_session.items():
                    if cfg_key in saved_config:
                        val = saved_config[cfg_key]
                        if transform:
                            val = transform(val)
                        st.session_state[session_key] = val
                        restored.append(session_key)
                if restored:
                    st.caption(
                        f"Restored settings from previous run: {', '.join(restored)}"
                    )
            if run_mode and "Selective" in run_mode:
                st.markdown(
                    "**Select steps to rerun** (unchecked = keep previous result):"
                )
                cols = st.columns(3)
                for i, step in enumerate(prev_steps):
                    with cols[i % 3]:
                        st.checkbox(
                            f"{step['name']}",
                            value=not step.get("success", False),
                            key=f"rerun_{step['name']}",
                        )
        except Exception:
            pass


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TAB 2: Pipeline Overview
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

with tab_pipeline:
    if _needs_run_mode_gate():
        st.warning(_GATE_MSG)
    else:
        # ──────────────────────────────────────────────────────────────
        # Stage 1: Secretion System Identification
        # ──────────────────────────────────────────────────────────────

        st.subheader("Stage 1: Secretion System Identification")
        st.markdown(
            "Detects secretion systems in the genome using **MacSyFinder v2** with "
            "TXSScan models. Proteins are extracted from the input genome, then "
            "MacSyFinder searches for conserved secretion system components using "
            "HMM profiles. Systems that meet the completeness threshold are retained.",
            help="MacSyFinder identifies macromolecular systems by searching for sets "
            "of co-occurring protein components using HMM profiles and genetic "
            "organization rules defined in TXSScan models.",
        )

        col1, col2 = st.columns(2)
        with col1:
            st.slider(
                "System completeness threshold",
                0.0,
                1.0,
                0.8,
                0.05,
                help="MacSyFinder wholeness score minimum. 0.8 = at least 80% of "
                "expected components must be present for a system to be called.",
                key="wholeness",
            )
        with col2:
            all_system_types = [
                "T1SS",
                "T2SS",
                "T3SS",
                "T4SS",
                "T5aSS",
                "T5bSS",
                "T5cSS",
                "T6SSi",
                "T6SSii",
                "T6SSiii",
                "T9SS",
                "Flagellum",
                "Tad",
                "pT4SSt",
            ]
            st.multiselect(
                "Exclude these system types",
                all_system_types,
                default=["Flagellum", "Tad", "T3SS"],
                key="excluded",
                help="Flagellum and Tad are excluded because they are not true secretion "
                "systems. T3SS is excluded by default because DeepSecE T3SS "
                "predictions are unreliable (mostly flagellar misclassification).",
            )

        st.divider()

        # ──────────────────────────────────────────────────────────────
        # Stage 2: Secreted Protein Identification
        # ──────────────────────────────────────────────────────────────

        st.subheader("Stage 2: Secreted Protein Identification")
        st.markdown(
            "Identifies candidate secreted proteins using **guilt by association** — "
            "proteins genomically proximal to secretion system components are evaluated. "
            "Localization predictions verify that detected systems are correctly assembled "
            "and flag nearby proteins as potential secretion candidates.",
            help="Proteins within a configurable window (default: +/- 3 genes) of each "
            "secretion system component are extracted and evaluated using "
            "localization prediction tools.",
        )

        # Proximity & thresholds
        with st.expander("Detection Parameters", expanded=False):
            pc1, pc2, pc3 = st.columns(3)
            with pc1:
                st.slider(
                    "Proximity window (genes)",
                    1,
                    15,
                    3,
                    help="How many genes upstream/downstream of each secretion system "
                    "component to search for candidate secreted proteins.",
                    key="window",
                )
            with pc2:
                st.slider(
                    "DeepLocPro extracellular threshold",
                    0.0,
                    1.0,
                    0.8,
                    0.05,
                    help="Minimum probability for a protein to be called extracellular "
                    "by DeepLocPro.",
                    key="conf",
                )
            with pc3:
                st.slider(
                    "Required fraction correctly localized",
                    0.0,
                    1.0,
                    0.8,
                    0.05,
                    help="Fraction of secretion system components that must have correct "
                    "predicted localization for the system to be considered valid.",
                    key="frac",
                )

        # ── Whole-genome prediction option ──
        with st.expander(
            "Advanced: Run predictions on entire proteome", expanded=False
        ):
            st.warning(
                "By default, predictions run only on proteins near detected secretion "
                "systems (typically 50-200 proteins per genome). Enabling whole-genome "
                "mode runs predictions on **all proteins** in the genome (typically "
                "3,000-6,000 proteins), which will **significantly increase runtime** "
                "(10-50x longer for cloud APIs). Only enable if you need proteome-wide "
                "predictions for downstream analysis."
            )
            wg1, wg2, wg3 = st.columns(3)
            with wg1:
                st.checkbox(
                    "DeepLocPro (whole genome)",
                    value=False,
                    key="dlp_whole_genome",
                    help="Run DeepLocPro on all proteins, not just those near "
                    "secretion systems. Significantly increases runtime.",
                )
            with wg2:
                st.checkbox(
                    "DeepSecE (whole genome)",
                    value=False,
                    key="dse_whole_genome",
                    help="Run DeepSecE on all proteins, not just those near "
                    "secretion systems. Significantly increases runtime.",
                )
            with wg3:
                st.checkbox(
                    "SignalP (whole genome)",
                    value=False,
                    key="sp_whole_genome",
                    help="Run SignalP on all proteins, not just those near "
                    "secretion systems. Significantly increases runtime.",
                )

        st.markdown("---")

        # ── DeepLocPro ──
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(
                "**DeepLocPro** — Subcellular Localization",
                help="Predicts protein subcellular localization using protein language "
                "models. Used to verify secretion system components localize as "
                "expected, and to identify nearby extracellular proteins as "
                "secretion candidates.",
            )
            st.caption(
                "Predicts localization of secretion system components (secondary "
                "verification) and nearby proteins for possible secretion."
            )
            dlp_mode = st.radio(
                "Running mode",
                ["BioLib cloud (no install needed)", "Local install (DTU license)"],
                key="dlp_mode",
                help="BioLib cloud: free, no license required, ~5-10 min per genome. "
                "Local: faster, requires free DTU academic license + ~5 GB model.",
            )
            if "Local" in dlp_mode:
                dlp_found = shutil.which("deeplocpro") is not None
                dlp_path = st.text_input("DeepLocPro install path", key="dlp_path")
                if not dlp_found and not dlp_path:
                    st.info(
                        "**DeepLocPro is not installed locally.** Requires DTU academic "
                        "license (~5 GB, GPU recommended). See the "
                        "[installation guide](https://github.com/billerbeck-lab/ssign/blob/main/docs/optional_tools.md). "
                        "Or use **cloud mode** (no install needed)."
                    )
            else:
                st.caption("Cloud mode (~5-10 min per genome, no setup needed).")
                dlp_path = ""

        with col2:
            st.markdown(
                "**SignalP 6.0** — Signal Peptides",
                help="Predicts N-terminal signal peptides (Sec/SPI, Sec/SPII, Tat/SPI, "
                "Tat/SPII, PILIN) recognised by signal peptidases. Only useful for "
                "T2SS and T5SS substrates which use the Sec pathway to cross the "
                "inner membrane. T1SS, T3SS, T4SS, T6SS substrates do NOT have "
                "N-terminal signal peptides — SignalP will return 'OTHER' for them.",
            )
            st.caption(
                "Detects Sec/Tat signal peptides — most useful for T2SS and T5SS. "
                "T1SS/T3SS/T4SS/T6SS substrates use different secretion mechanisms "
                "and will be predicted as 'OTHER' (not a tool failure)."
            )
            run_signalp = st.checkbox(
                "Enable SignalP 6.0",
                value=False,
                key="run_signalp",
                help="Adds signal peptide predictions as an additional layer of evidence "
                "for secretion.",
            )
            if run_signalp:
                sp_mode = st.radio(
                    "Running mode",
                    ["BioLib cloud (no install needed)", "Local install (DTU license)"],
                    key="sp_mode",
                    help="BioLib cloud: free, ~2-5 min per genome. "
                    "Local: faster, requires free DTU academic license.",
                )
                if "Local" in sp_mode:
                    sp_found = shutil.which("signalp6") is not None
                    sp_path = st.text_input("SignalP install path", key="sp_path")
                    if not sp_found and not sp_path:
                        st.info(
                            "**SignalP 6.0 is not installed locally.** Requires DTU academic "
                            "license (~1 GB). See the "
                            "[installation guide](https://github.com/billerbeck-lab/ssign/blob/main/docs/optional_tools.md). "
                            "Or use **cloud mode** (no install needed)."
                        )
                else:
                    st.caption("Cloud mode (~2-5 min per genome, no setup needed).")
                    sp_path = ""
                with st.expander("SignalP threshold", expanded=False):
                    st.slider(
                        "Min. probability",
                        0.0,
                        1.0,
                        0.5,
                        0.05,
                        key="sp_min_prob",
                        help="Minimum SignalP probability to call a signal peptide.",
                    )
            else:
                sp_mode = "BioLib cloud"
                sp_path = ""

        st.markdown("---")

        # ── DeepSecE ──
        st.markdown(
            "**DeepSecE** — Secretion Type Prediction",
            help="A deep-learning model that predicts which secretion system type "
            "each protein is secreted by. Cross-checks DeepLocPro localization "
            "predictions and can identify additional secreted proteins that "
            "DeepLocPro may miss. Does not work for T5SS autotransporters.",
        )
        st.caption(
            "Secondary prediction of secreted proteins — cross-checks DeepLocPro "
            "predictions and provides additional candidates, good for enrichment "
            "and verification."
        )

        deepsece_available = False
        deepsece_checkpoint_ok = False
        try:
            import DeepSecE  # noqa: F401 -- import IS the availability check; ImportError caught below

            deepsece_available = True
            # Check if the model checkpoint exists (~2.5 GB)
            _dse_ckpt = os.path.join(
                os.path.expanduser("~"), ".ssign", "models", "deepsece_checkpoint.pt"
            )
            deepsece_checkpoint_ok = (
                os.path.isfile(_dse_ckpt) and os.path.getsize(_dse_ckpt) > 100_000_000
            )
        except ImportError:
            pass

        if deepsece_available and deepsece_checkpoint_ok:
            st.success("DeepSecE is installed and ready.")
            run_deepsece = st.checkbox(
                "Enable DeepSecE",
                value=True,
                key="run_deepsece",
            )
            if run_deepsece:
                with st.expander("DeepSecE threshold", expanded=False):
                    st.slider(
                        "Min. probability",
                        0.0,
                        1.0,
                        0.8,
                        0.05,
                        key="dse_min_prob",
                        help="Minimum DeepSecE probability to call a protein as secreted.",
                    )
        elif deepsece_available and not deepsece_checkpoint_ok:
            run_deepsece = st.checkbox(
                "Enable DeepSecE",
                value=False,
                key="run_deepsece",
                disabled=True,
            )
            _dse_path = os.path.join("~", ".ssign", "models", "deepsece_checkpoint.pt")
            st.warning(
                "**DeepSecE package is installed but the model checkpoint is missing** "
                f"(~2.5 GB, expected at `{_dse_path}`).\n\n"
                "Download it manually (primary Zenodo mirror; URL will be populated at v1.0.0 release):\n"
                "```bash\n"
                "mkdir -p ~/.ssign/models\n"
                "wget -c --tries=5 --timeout=60 \\\n"
                "  https://zenodo.org/records/PLACEHOLDER/files/deepsece_checkpoint.pt \\\n"
                "  -O ~/.ssign/models/deepsece_checkpoint.pt\n"
                "```\n"
                "Fallback — origin server (SJTU, China; known-unreliable):\n"
                "```bash\n"
                "wget -c --tries=5 --timeout=60 \\\n"
                "  https://tool2-mml.sjtu.edu.cn/DeepSecE/checkpoint.pt \\\n"
                "  -O ~/.ssign/models/deepsece_checkpoint.pt\n"
                "```\n"
                "Or set `SSIGN_DEEPSECE_CHECKPOINT_URL` to an institutional mirror, then "
                "restart ssign — it will auto-download on the next run."
            )
        else:
            run_deepsece = st.checkbox(
                "Enable DeepSecE",
                value=False,
                key="run_deepsece",
                disabled=True,
            )
            st.info(
                "**DeepSecE is not installed** but is recommended. It provides an "
                "independent cross-check against DeepLocPro.\n\n"
                "Install (~7.3 GB): `pip install ssign[deepsece]` (run inside your ssign venv)\n\n"
                "See the [installation guide](https://github.com/billerbeck-lab/ssign/blob/main/docs/optional_tools.md) "
                "for detailed instructions."
            )

        st.divider()

        # ──────────────────────────────────────────────────────────────
        # Stage 3: Improve Secreted Protein Annotations
        # ──────────────────────────────────────────────────────────────

        st.subheader("Stage 3: Improve Secreted Protein Annotations")
        st.markdown(
            "Many secreted proteins in lesser-studied organisms are classified as "
            "hypothetical or domain of unknown function (DUF). These optional annotation "
            "tools can help resolve protein function. All run via **cloud APIs** by "
            "default — no local install needed. Uncheck to skip.",
            help="Each tool queries a different database or algorithm to add functional "
            "annotations to your secreted protein candidates. Running more tools "
            "takes longer but provides richer annotations.",
        )
        st.info(
            "If your genome is already well-annotated (e.g. from Bakta, Prokka, or "
            "NCBI RefSeq), you can skip these tools to save significant time. The "
            "existing annotations from your input file are already included in the "
            "results. You can always re-run with annotation tools enabled later using "
            "the **Resume** feature."
        )

        st.markdown("---")

        # ── BLASTp ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_blastp = st.checkbox(
                "BLASTp",
                value=True,
                key="run_blastp",
                help="Searches the NCBI nr (non-redundant) protein database for "
                "homologous proteins. Returns functional annotations, gene names, "
                "and descriptions from characterized homologs.",
            )
        with col_info:
            if run_blastp:
                st.caption("NCBI web API | ~30-60 min per genome")

                with st.expander("Taxonomy exclusion", expanded=False):
                    st.markdown(
                        "If your organism is already in NCBI nr, BLASTp returns your own "
                        "proteins as top hits, duplicating existing annotations instead of "
                        "providing new functional information.\n\n"
                        "- **Exclude species** — removes only your exact species, keeping "
                        "closely related hits from the same genus.\n"
                        "- **Exclude genus** — removes the entire genus, finding homologs "
                        "from more distant organisms. Useful when many genomes from your "
                        "genus are in nr.\n"
                        "- **Custom taxonomy ID(s)** — manually specify NCBI taxonomy IDs "
                        "to exclude. Enter one or more comma-separated taxids (look up at "
                        "[NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)). Useful "
                        "for excluding specific clades or strains.\n"
                        "- **No exclusion** — keeps all hits. Use if your organism is novel "
                        "or not yet in NCBI databases."
                    )

                # Show detected organisms (from local GenBank parsing, no API call)
                per_genome_orgs = st.session_state.get("per_genome_organisms", {})
                detected_organisms = set(per_genome_orgs.values())

                if detected_organisms:
                    if len(detected_organisms) == 1:
                        st.info(
                            f"Detected organism: **{next(iter(detected_organisms))}**"
                        )
                    else:
                        org_list = ", ".join(
                            f"*{o}*" for o in sorted(detected_organisms)
                        )
                        st.info(
                            f"Detected organisms: {org_list}\n\n"
                            "Taxonomy exclusion will be applied **per genome** — each "
                            "genome's own species/genus is excluded from its BLASTp search. "
                            "Taxonomy IDs are resolved automatically when the pipeline runs."
                        )

                # Exclusion level options (taxonomy IDs resolved at run time)
                exclusion_options = []
                if detected_organisms:
                    exclusion_options.append("Exclude input species (per genome)")
                    exclusion_options.append("Exclude input genus (per genome)")
                exclusion_options.append("Custom taxonomy ID(s)")
                exclusion_options.append("No exclusion (include all hits)")

                exclusion_choice = st.radio(
                    "BLASTp taxonomy exclusion",
                    exclusion_options,
                    index=0,
                    key="blastp_exclusion_mode",
                )

                if "species" in exclusion_choice.lower():
                    st.session_state.blastp_exclusion_level = "species"
                elif "genus" in exclusion_choice.lower():
                    st.session_state.blastp_exclusion_level = "genus"
                elif "Custom" in exclusion_choice:
                    st.session_state.blastp_exclusion_level = "custom"
                    st.text_input(
                        "Enter NCBI taxonomy ID(s)",
                        key="blastp_taxid",
                        placeholder="e.g. 339 or 339,340,338",
                        help="Comma-separate multiple taxonomy IDs to exclude "
                        "several organisms.",
                    )
                else:
                    st.session_state.blastp_exclusion_level = "none"
                    st.session_state.blastp_taxid = ""

                if not detected_organisms:
                    st.caption(
                        "No organisms auto-detected (FASTA input has no organism metadata). "
                        "Select **Custom** to enter NCBI taxonomy IDs, or **No exclusion** "
                        "if your organisms are not yet in NCBI databases."
                    )

                # BLASTp result filters
                st.markdown("**Result filters:**")
                bc1, bc2, bc3 = st.columns(3)
                with bc1:
                    st.slider(
                        "Min. % identity",
                        0,
                        100,
                        80,
                        5,
                        key="blastp_pident",
                        help="Minimum percent identity to keep a BLASTp hit.",
                    )
                with bc2:
                    st.slider(
                        "Min. query coverage (%)",
                        0,
                        100,
                        80,
                        5,
                        key="blastp_qcov",
                        help="Minimum query coverage to keep a BLASTp hit.",
                    )
                with bc3:
                    st.number_input(
                        "E-value threshold",
                        value=1e-5,
                        format="%.0e",
                        key="blastp_evalue",
                        help="Maximum e-value for BLASTp hits.",
                    )

        # ── HHpred ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_hh = st.checkbox(
                "HHpred (Pfam + PDB)",
                value=True,
                key="run_hh",
                help="Detects remote homology using profile-profile comparison via the "
                "MPI Bioinformatics Toolkit. Searches Pfam and PDB databases to "
                "find structural and functional domains even in highly diverged proteins.",
            )
        with col_info:
            if run_hh:
                st.caption("MPI Toolkit API | ~45-90 min per genome")
                st.slider(
                    "Min. probability (%)",
                    0,
                    100,
                    40,
                    5,
                    key="hhpred_min_prob",
                    help="Minimum HHpred probability to keep a hit. "
                    "Default 40% balances sensitivity and specificity.",
                )

        # ── InterProScan ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_iprs = st.checkbox(
                "InterProScan",
                value=True,
                key="run_iprs",
                help="Scans proteins against InterPro's consortium of protein signature "
                "databases (Pfam, SMART, CDD, PANTHER, etc.) to identify domains, "
                "families, and Gene Ontology (GO) terms.",
            )
        with col_info:
            if run_iprs:
                st.caption("EBI REST API | ~20-40 min per genome")
                st.number_input(
                    "E-value threshold",
                    value=1e-5,
                    format="%.0e",
                    key="iprs_evalue",
                    help="Maximum e-value for InterProScan domain hits.",
                )

        # ── ProtParam ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_pp = st.checkbox(
                "ProtParam",
                value=True,
                key="run_pp",
                help="Computes physicochemical properties including molecular weight, "
                "isoelectric point (pI), GRAVY hydropathicity, instability index, "
                "and amino acid composition. Runs locally via BioPython — instant.",
            )
        with col_info:
            if run_pp:
                st.caption("Local (BioPython) | No download needed | Instant")

        st.divider()

        st.markdown("**Coming in v1.0.0**")
        st.caption(
            "The following tool will be integrated as a first-class "
            "annotation source in v1.0.0. Disabled for now."
        )

        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            st.checkbox("pLM-BLAST (ECOD70)", value=False, key="run_plm", disabled=True)
        with col_info:
            st.caption(
                "Protein language model-based remote homology detection against the "
                "ECOD70 structural database (~20 GB). Included in the **extended** "
                "and **full** install tiers."
            )

        st.divider()

        # ──────────────────────────────────────────────────────────────
        # Stage 4: Generate Data & Figures
        # ──────────────────────────────────────────────────────────────

        st.subheader("Stage 4: Generate Data & Figures")
        st.markdown(
            "Produces result tables and summary figures for the identified "
            "secreted proteins and their annotations.",
            help="Figures are saved as SVG and PNG to the output directory. "
            "All results are also exported as CSV tables.",
        )

        st.markdown("**Figures:**")
        fc1, fc2 = st.columns(2)
        with fc1:
            st.checkbox("Category distribution", value=True, key="fig_category")
            st.checkbox("SS composition", value=True, key="fig_ss_comp")
            st.checkbox("Tool coverage heatmap", value=True, key="fig_tool_heatmap")
        with fc2:
            st.checkbox(
                "Secreted protein count per genome",
                value=True,
                key="fig_substrate_count",
            )
            st.checkbox(
                "Functional annotation summary", value=True, key="fig_func_summary"
            )

        st.markdown("---")

        # ── Ortholog Group Assignment ──
        st.markdown("**Ortholog Group Assignment**")

        blastp_available = shutil.which("blastp") is not None

        if blastp_available:
            st.success("BLAST+ is installed and ready for ortholog grouping.")
            st.caption(
                "Secreted proteins are grouped into ortholog groups using all-vs-all "
                "BLASTp with single-linkage clustering."
            )
            col1, col2 = st.columns(2)
            with col1:
                st.slider(
                    "Min. % identity for ortholog grouping",
                    0,
                    100,
                    40,
                    5,
                    key="og_min_pident",
                    help="Minimum percent identity between two proteins to "
                    "consider them orthologs. Default 40% is permissive; "
                    "increase for stricter groups.",
                )
            with col2:
                st.slider(
                    "Min. query coverage (%) for ortholog grouping",
                    0,
                    100,
                    70,
                    5,
                    key="og_min_qcov",
                    help="Minimum query coverage to consider a BLASTp hit "
                    "as an ortholog relationship. Default 70%.",
                )
        else:
            st.info(
                "**BLAST+ is not installed** (needed for ortholog grouping).\n\n"
                "Ortholog grouping clusters secreted proteins across genomes into "
                "groups of related proteins, useful for comparative analysis.\n\n"
                "To install (~200 MB):\n"
                "```\nsudo apt install ncbi-blast+\n```\n"
                "Or: `conda install -c bioconda blast`\n\n"
                "ssign works without BLAST+ — ortholog grouping will be skipped."
            )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TAB 3: Run & Results
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

with tab_run:
    if _needs_run_mode_gate():
        st.warning(_GATE_MSG)
    else:
        st.subheader("Run Pipeline")

        uploaded_files = st.session_state.get("uploaded_files_data", [])

        # Validation
        issues = []
        if not uploaded_files:
            issues.append("No genome file uploaded. Go to **Upload & Configure** tab.")

        can_run = len(issues) == 0

        if can_run:
            st.success(
                f"Ready to run on **{len(uploaded_files)} genome(s)**: "
                + ", ".join(uf.name for uf in uploaded_files)
            )
        else:
            for issue in issues:
                st.warning(issue)

        # Estimated time — based on measured per-protein API timings:
        #   DeepLocPro:    ~6s per neighborhood protein (DTU batch, ~100 proteins/11 min)
        #   SignalP:       ~3s per neighborhood protein (DTU batch)
        #   BLASTp:        ~30s per secreted protein (NCBI, 15s between batches)
        #   HHpred:        ~260s per secreted protein (MPI, 20s submit + ~200s processing)
        #   InterProScan:  ~15s per secreted protein (EBI, 1s submit + polling)
        #   ProtParam:     instant (local)
        #   MacSyFinder:   ~2-5 min per genome (local)
        #
        # Typical genome: ~100 neighborhood proteins, ~5-10 secreted proteins
        st.divider()
        n_genomes = len(uploaded_files) or 1
        # Assume typical counts (refined once pipeline runs)
        est_neighborhood = 100  # typical neighborhood proteins per genome
        est_secreted = 8  # typical secreted proteins per genome

        st.markdown("**Estimated run time:**")
        st.caption(
            f"Based on typical genome (~{est_neighborhood} neighborhood proteins, "
            f"~{est_secreted} secreted proteins). Actual time depends on genome "
            f"content and API server load."
        )

        # Per-genome time breakdown
        base_min = 5  # MacSyFinder + local steps
        dlp_min = max(5, est_neighborhood * 6 / 60)  # ~6s per protein
        sp_min = (
            max(3, est_neighborhood * 3 / 60)
            if st.session_state.get("run_signalp")
            else 0
        )
        dse_min = 3 if st.session_state.get("run_deepsece") else 0
        # Prediction tools run in parallel — take the max
        prediction_min = max(dlp_min, sp_min, dse_min)

        blastp_min = est_secreted * 30 / 60 if st.session_state.get("run_blastp") else 0
        hh_min = est_secreted * 260 / 60 if st.session_state.get("run_hh") else 0
        iprs_min = est_secreted * 15 / 60 if st.session_state.get("run_iprs") else 0
        # Annotation tools run in parallel — take the max
        annotation_min = max(blastp_min, hh_min, iprs_min)

        per_genome_min = base_min + prediction_min + annotation_min

        time_parts = [f"Local steps (MacSyFinder, filtering): ~{base_min} min"]
        time_parts.append(
            f"Predictions (DeepLocPro{' + SignalP' if sp_min else ''}{' + DeepSecE' if dse_min else ''}): ~{prediction_min:.0f} min"
        )
        if annotation_min > 0:
            ann_tools = []
            if blastp_min:
                ann_tools.append(f"BLASTp ~{blastp_min:.0f} min")
            if hh_min:
                ann_tools.append(f"HHpred ~{hh_min:.0f} min")
            if iprs_min:
                ann_tools.append(f"InterProScan ~{iprs_min:.0f} min")
            time_parts.append(
                f"Annotations (parallel: {', '.join(ann_tools)}): ~{annotation_min:.0f} min"
            )

        for t in time_parts:
            st.markdown(f"- {t}")

        # Total estimate accounting for parallelism
        if n_genomes == 1:
            total_min = per_genome_min
        else:
            # With semaphore-based scheduling:
            # - Local steps: all parallel (fast)
            # - DTU predictions: 5 concurrent → n_genomes/5 batches
            # - HHpred: sequential (1 at a time) → bottleneck
            # - Other annotations: 5 concurrent
            dtu_batches = max(1, n_genomes / 5)
            hh_sequential = n_genomes  # 1 at a time
            local_total = base_min  # all parallel, ~same time
            pred_total = prediction_min * dtu_batches
            if hh_min > 0:
                ann_total = hh_min * hh_sequential  # HHpred is the bottleneck
            else:
                ann_batches = max(1, n_genomes / 5)
                ann_total = annotation_min * ann_batches
            total_min = local_total + pred_total + ann_total

        total_max = total_min * 1.5  # API variability
        genome_note = f" for {n_genomes} genome(s)" if n_genomes > 1 else ""
        if total_min >= 100:
            est_low = f"{total_min / 60:.1f}"
            est_high = f"{total_max / 60:.1f}"
            st.info(
                f"**Estimated total: ~{est_low}-{est_high} hours{genome_note}** (rough estimate, depends on API load)"
            )
        else:
            st.info(
                f"**Estimated total: ~{total_min:.0f}-{total_max:.0f} minutes{genome_note}** (rough estimate, depends on API load)"
            )

        if n_genomes > 1 and hh_min > 0:
            st.caption(
                f"HHpred is the bottleneck for multi-genome runs (~{hh_min:.0f} min per genome, "
                f"runs one genome at a time due to MPI rate limits)."
            )

        st.divider()

        st.warning(
            "**Before you run:**\n"
            "- The time estimate above is an **initial estimate** based on typical genomes. "
            "It will update with actual numbers once the pipeline discovers how many "
            "secreted proteins each genome has.\n"
            "- Closing the **browser tab** is safe — the pipeline continues in the "
            "background and results are saved to the output directory.\n"
            "- Closing the **terminal** will stop the pipeline. If interrupted, re-run "
            "with the same output directory and **Resume** to continue where you left off.\n"
            "- An internet disconnection will cause active API calls to fail, but the "
            "retry logic will attempt to recover. If retries are exhausted, the failed "
            "step can be resumed later."
        )

        # Resume — use run mode from Upload tab if set, otherwise detect here
        run_mode = st.session_state.get("run_mode_choice", "")
        outdir_val = st.session_state.get("outdir_input", "")
        progress_exists = _has_previous_progress(outdir_val) if outdir_val else False

        if progress_exists and not run_mode:
            st.checkbox(
                "Resume from previous run (skip completed steps)",
                value=True,
                key="resume_run",
                help="Steps that completed successfully will be skipped.",
            )
            st.caption(
                "Previous progress detected. Configure run mode in the Upload & Configure tab for more options."
            )
        elif "Resume" in run_mode:
            st.info("Resuming from previous run — completed steps will be skipped.")
            st.session_state.resume_run = True
        elif "fresh" in run_mode.lower() if run_mode else False:
            st.info("Starting fresh — all steps will rerun.")
            st.session_state.resume_run = False
        elif "Selective" in run_mode:
            st.info(
                "Selective rerun — only checked steps from the Upload & Configure tab will rerun."
            )
            st.session_state.resume_run = True
        else:
            st.session_state.resume_run = False

        # ── Helper: resolve per-genome BLASTp taxid ──

        def _resolve_blastp_taxid(filename: str) -> str:
            """Resolve the BLASTp exclusion taxid for a specific genome file.

            Taxonomy resolution (NCBI API) happens here at run time, not during
            upload, to avoid blocking the GUI. Results are cached in session state.
            """
            level = st.session_state.get("blastp_exclusion_level", "species")
            if level == "none":
                return ""
            if level == "custom":
                return st.session_state.get("blastp_taxid", "")

            # Get organism name (extracted locally during upload)
            per_genome_orgs = st.session_state.get("per_genome_organisms", {})
            org_name = per_genome_orgs.get(filename, "")
            if not org_name:
                return ""

            # Resolve taxonomy via NCBI (cached after first call)
            cache_key = "_tax_cache"
            if cache_key not in st.session_state:
                st.session_state[cache_key] = {}
            cache = st.session_state[cache_key]

            if org_name not in cache:
                try:
                    from ssign_app.scripts.resolve_taxonomy import resolve_organism

                    cache[org_name] = resolve_organism(org_name)
                except Exception:
                    cache[org_name] = {"species": None, "genus": None}

            tax_info = cache[org_name]
            tax_entry = tax_info.get(level)
            if tax_entry and isinstance(tax_entry, dict):
                return str(tax_entry.get("taxid", ""))
            return ""

        # Run button
        if st.button(
            "Run ssign",
            type="primary",
            disabled=not can_run or st.session_state.running,
            use_container_width=True,
        ):
            st.session_state.running = True
            st.session_state.results = []

            # Write uploaded files to temp directory (preserve original filenames)
            tmpdir = tempfile.mkdtemp(prefix="ssign_input_")
            input_paths = []
            original_filenames = []
            for uf in uploaded_files:
                dest = os.path.join(tmpdir, uf.name)
                uf.seek(0)
                with open(dest, "wb") as out:
                    out.write(uf.read())
                input_paths.append(dest)
                original_filenames.append(uf.name)

            # ── Build configs and progress UI for each genome ──
            n_genomes_to_run = len(input_paths)
            genome_configs = []
            genome_progress = []  # (progress_bar, status_text) per genome
            # Divide local CPU budget across concurrent genomes — without this,
            # macsyfinder -w {cpu_count} × N genomes oversubscribes cores.
            _cpu_per_genome = max(1, (os.cpu_count() or 4) // max(1, n_genomes_to_run))

            for file_idx, input_path in enumerate(input_paths):
                if n_genomes_to_run > 1:
                    sample_id = (
                        st.session_state.get("sample_names", ["sample"])[file_idx]
                        if file_idx < len(st.session_state.get("sample_names", []))
                        else f"sample_{file_idx + 1}"
                    )
                else:
                    sample_id = st.session_state.get("sample_name", "sample")

                orig_fname = (
                    original_filenames[file_idx]
                    if file_idx < len(original_filenames)
                    else ""
                )
                genome_taxid = _resolve_blastp_taxid(orig_fname)

                config = PipelineConfig(
                    input_path=input_path,
                    original_filename=orig_fname,
                    sample_id=sample_id,
                    outdir=st.session_state.get("outdir_input", "./results"),
                    run_bakta=st.session_state.get("use_bakta", False),
                    bakta_db=st.session_state.get("bakta_db_path", ""),
                    use_input_annotations=st.session_state.get(
                        "use_input_annotations", False
                    ),
                    wholeness_threshold=st.session_state.get("wholeness", 0.8),
                    excluded_systems=st.session_state.get(
                        "excluded", ["Flagellum", "Tad", "T3SS"]
                    ),
                    conf_threshold=st.session_state.get("conf", 0.8),
                    proximity_window=st.session_state.get("window", 3),
                    required_fraction_correct=st.session_state.get("frac", 0.8),
                    deeplocpro_mode="local"
                    if "Local" in st.session_state.get("dlp_mode", "")
                    else "remote",
                    deeplocpro_path=st.session_state.get("dlp_path", ""),
                    skip_deepsece=not st.session_state.get("run_deepsece", False),
                    skip_signalp=not st.session_state.get("run_signalp", False),
                    signalp_mode="local"
                    if "Local" in st.session_state.get("sp_mode", "")
                    else "remote",
                    signalp_path=st.session_state.get("sp_path", ""),
                    dlp_whole_genome=st.session_state.get("dlp_whole_genome", False),
                    dse_whole_genome=st.session_state.get("dse_whole_genome", False),
                    sp_whole_genome=st.session_state.get("sp_whole_genome", False),
                    skip_blastp=not st.session_state.get("run_blastp", True),
                    blastp_db=st.session_state.get("blastp_db", ""),
                    blastp_exclude_taxid=genome_taxid,
                    blastp_min_pident=float(st.session_state.get("blastp_pident", 80)),
                    blastp_min_qcov=float(st.session_state.get("blastp_qcov", 80)),
                    blastp_evalue=float(st.session_state.get("blastp_evalue", 1e-5)),
                    skip_hhsuite=not st.session_state.get("run_hh", False),
                    hhsuite_pfam_db="",
                    hhsuite_pdb70_db="",
                    hhpred_min_probability=float(
                        st.session_state.get("hhpred_min_prob", 40)
                    ),
                    skip_interproscan=not st.session_state.get("run_iprs", True),
                    interproscan_db="",
                    skip_plmblast=not st.session_state.get("run_plm", False),
                    plmblast_db=st.session_state.get("plm_db", ""),
                    skip_protparam=not st.session_state.get("run_pp", True),
                    fig_category=st.session_state.get("fig_category", True),
                    fig_ss_comp=st.session_state.get("fig_ss_comp", True),
                    fig_tool_heatmap=st.session_state.get("fig_tool_heatmap", True),
                    fig_substrate_count=st.session_state.get(
                        "fig_substrate_count", True
                    ),
                    fig_func_summary=st.session_state.get("fig_func_summary", True),
                    interproscan_evalue=float(
                        st.session_state.get("iprs_evalue", 1e-5)
                    ),
                    deepsece_min_prob=float(st.session_state.get("dse_min_prob", 0.8)),
                    signalp_min_prob=float(st.session_state.get("sp_min_prob", 0.5)),
                    ortholog_min_pident=float(
                        st.session_state.get("og_min_pident", 40)
                    ),
                    ortholog_min_qcov=float(st.session_state.get("og_min_qcov", 70)),
                    cpu_per_genome=_cpu_per_genome,
                )
                genome_configs.append(config)

                # Create per-genome progress UI
                if n_genomes_to_run > 1:
                    st.markdown(
                        f"**Genome {file_idx + 1}/{n_genomes_to_run}: {sample_id}**"
                    )
                progress_bar = st.progress(0)
                status_text = st.empty()
                genome_progress.append((progress_bar, status_text))

            # ── Dynamic time estimate placeholder ──
            estimate_box = st.empty()

            # ── Run genomes (parallel if multiple, sequential if single) ──
            all_results = []
            resume = st.session_state.get("resume_run", False)

            if n_genomes_to_run > 1:
                import re as _re
                import threading
                import time as _time
                from concurrent.futures import ThreadPoolExecutor, as_completed

                from streamlit.runtime.scriptrunner import (
                    add_script_run_ctx,
                    get_script_run_ctx,
                )

                _ctx = get_script_run_ctx()
                _start_time = _time.monotonic()

                # Shared counters for dynamic estimate
                _counts_lock = threading.Lock()
                _genome_counts = {}  # {idx: {'neighborhood': N, 'secreted': N, 'done': bool}}

                # Per-protein timing constants (measured)
                _SEC_PER_NEIGH_DLP = 6  # DeepLocPro: ~6s per neighborhood protein
                _SEC_PER_SEC_HH = 260  # HHpred: ~260s per secreted protein
                _SEC_PER_SEC_BLAST = 30  # BLASTp: ~30s per secreted protein
                _SEC_PER_SEC_IPRS = 15  # InterProScan: ~15s per secreted protein

                def _update_estimate():
                    """Recalculate and display time estimate based on discovered protein counts."""
                    with _counts_lock:
                        counted = [
                            c
                            for c in _genome_counts.values()
                            if c.get("secreted") is not None
                        ]
                        done = sum(1 for c in _genome_counts.values() if c.get("done"))
                        remaining = n_genomes_to_run - done

                    if not counted or remaining == 0:
                        return

                    # Average from discovered genomes
                    avg_sec = sum(c["secreted"] for c in counted) / len(counted)
                    total_sec_remaining = avg_sec * remaining

                    # HHpred is sequential (sem=1), others are parallel
                    hh_enabled = st.session_state.get("run_hh", False)
                    if hh_enabled:
                        hh_remaining_min = total_sec_remaining * _SEC_PER_SEC_HH / 60
                    else:
                        hh_remaining_min = 0

                    blastp_enabled = st.session_state.get("run_blastp", False)
                    iprs_enabled = st.session_state.get("run_iprs", False)
                    other_ann_min = max(
                        total_sec_remaining * _SEC_PER_SEC_BLAST / 60 / 5
                        if blastp_enabled
                        else 0,
                        total_sec_remaining * _SEC_PER_SEC_IPRS / 60 / 5
                        if iprs_enabled
                        else 0,
                    )

                    elapsed = (_time.monotonic() - _start_time) / 60
                    est_remaining = max(0, hh_remaining_min + other_ann_min)

                    if est_remaining >= 60:
                        est_str = f"~{est_remaining / 60:.1f} hours"
                    else:
                        est_str = f"~{est_remaining:.0f} minutes"

                    estimate_box.info(
                        f"**{done}/{n_genomes_to_run} genomes complete** | "
                        f"Elapsed: {elapsed:.0f} min | "
                        f"Est. remaining: {est_str} | "
                        f"Avg {avg_sec:.0f} secreted proteins/genome"
                    )

                # DTU BioLib is the only remaining rate-limited remote service
                # (used by DeepLocPro + SignalP in remote mode). The ncbi/mpi/
                # ebi semaphores were removed with Phase 2.1 (BLAST/HH-suite/
                # InterProScan are local-only now).
                _api_sem = {
                    "dtu": threading.Semaphore(5),
                }

                def _run_one_genome(idx):
                    add_script_run_ctx(threading.current_thread(), _ctx)
                    cfg = genome_configs[idx]
                    bar, status = genome_progress[idx]

                    def _update(step, pct, msg):
                        bar.progress(min(pct, 100) / 100)
                        status.markdown(f"**{step}** \u2014 {msg}")

                        # Parse protein counts from step messages
                        try:
                            m = _re.search(r"(\d+)\s+neighborhood proteins", msg)
                            if m:
                                with _counts_lock:
                                    _genome_counts.setdefault(idx, {})[
                                        "neighborhood"
                                    ] = int(m.group(1))

                            m = _re.search(r"(\d+)\s+secreted proteins", msg)
                            if m:
                                with _counts_lock:
                                    _genome_counts.setdefault(idx, {})["secreted"] = (
                                        int(m.group(1))
                                    )
                                _update_estimate()
                        except Exception:
                            pass

                    runner = PipelineRunner(
                        cfg, progress_callback=_update, api_semaphores=_api_sem
                    )
                    result = runner.run(resume=resume)

                    with _counts_lock:
                        _genome_counts.setdefault(idx, {})["done"] = True
                    _update_estimate()

                    return result

                with ThreadPoolExecutor(max_workers=n_genomes_to_run) as executor:
                    futures = {
                        executor.submit(_run_one_genome, i): i
                        for i in range(n_genomes_to_run)
                    }
                    for future in as_completed(futures):
                        idx = futures[future]
                        try:
                            results = future.result()
                            all_results.extend(results)
                        except Exception as e:
                            bar, status = genome_progress[idx]
                            status.markdown(f":red[**Error**] \u2014 {e}")
                            all_results.append(StepResult("pipeline", False, str(e)))
            else:
                # Single genome — run directly
                import re as _re
                import time as _time

                _start_time = _time.monotonic()
                bar, status = genome_progress[0]

                def _update_single(step, pct, msg):
                    bar.progress(min(pct, 100) / 100)
                    status.markdown(f"**{step}** \u2014 {msg}")
                    # Update estimate when protein counts are discovered
                    try:
                        m = _re.search(r"(\d+)\s+secreted proteins", msg)
                        if m:
                            n_sec = int(m.group(1))
                            elapsed = (_time.monotonic() - _start_time) / 60
                            hh_min = (
                                n_sec * 260 / 60
                                if st.session_state.get("run_hh")
                                else 0
                            )
                            blast_min = (
                                n_sec * 30 / 60
                                if st.session_state.get("run_blastp")
                                else 0
                            )
                            iprs_min = (
                                n_sec * 15 / 60
                                if st.session_state.get("run_iprs")
                                else 0
                            )
                            ann_min = max(hh_min, blast_min, iprs_min)
                            est_str = (
                                f"~{ann_min:.0f} min"
                                if ann_min < 60
                                else f"~{ann_min / 60:.1f} hr"
                            )
                            estimate_box.info(
                                f"**{n_sec} secreted proteins found** | "
                                f"Elapsed: {elapsed:.0f} min | "
                                f"Annotation est: {est_str}"
                            )
                    except Exception:
                        pass

                runner = PipelineRunner(
                    genome_configs[0], progress_callback=_update_single
                )
                with st.spinner("Running ssign pipeline..."):
                    results = runner.run(resume=resume)
                all_results.extend(results)

            st.session_state.results = all_results
            st.session_state.running = False
            outdir_final = st.session_state.get("outdir_input", "./results")
            st.session_state.output_dir = outdir_final

            # ── Merge per-genome outputs into combined files ──
            sample_names_run = [c.sample_id for c in genome_configs]
            _merge_genome_outputs(outdir_final, sample_names_run)

            # Show results summary
            n_success = sum(1 for r in all_results if r.success)
            n_total = len(all_results)
            n_genomes_run = len(input_paths)

            if n_success == n_total:
                st.success(
                    f"Pipeline completed successfully for {n_genomes_run} genome(s)! ({n_success}/{n_total} steps)"
                )
            else:
                st.warning(
                    f"Pipeline finished with issues for {n_genomes_run} genome(s) ({n_success}/{n_total} steps succeeded)"
                )

        # Show progress if running
        if st.session_state.running:
            st.info(f"Running: {st.session_state.get('progress_step', '')} ...")

        # ── Results section ──

        st.divider()
        st.subheader("Results")

        if st.session_state.results:
            results = st.session_state.results
            output_dir = st.session_state.get("output_dir", "")

            # Summary metrics
            n_ok = sum(1 for r in results if r.success)
            n_fail = sum(1 for r in results if not r.success)
            mc1, mc2, mc3 = st.columns(3)
            with mc1:
                st.metric("Steps completed", f"{n_ok}/{len(results)}")
            with mc2:
                st.metric("Failed", n_fail)
            with mc3:
                n_secreted = "\u2014"
                if output_dir:
                    results_csv = os.path.join(output_dir, "ssign_results.csv")
                    if os.path.exists(results_csv):
                        try:
                            with open(results_csv) as _f:
                                lines = _f.readlines()
                            # Count data rows in the Secreted Proteins chunk
                            in_proteins = False
                            count = 0
                            for line in lines:
                                if "# Secreted Proteins" in line:
                                    in_proteins = True
                                    continue
                                if line.startswith("#") or (
                                    not line.strip() and in_proteins
                                ):
                                    if in_proteins:
                                        break
                                    continue
                                if (
                                    in_proteins
                                    and line.strip()
                                    and not line.startswith("locus_tag")
                                ):
                                    count += 1
                            n_secreted = str(count) if count > 0 else "\u2014"
                        except Exception:
                            pass
                st.metric("Secreted proteins found", n_secreted)

            st.divider()

            # Show step results
            for r in results:
                if r.success:
                    st.markdown(f":green[**\u25cf**] **{r.name}** \u2014 {r.message}")
                else:
                    st.markdown(f":red[**\u25cf**] **{r.name}** \u2014 {r.message}")

            # Download buttons for output files
            if output_dir and os.path.isdir(output_dir):
                st.divider()
                st.info(f"All output files are saved to: `{output_dir}`")

                # ZIP download
                import io
                import zipfile

                zip_buf = io.BytesIO()
                with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                    for fname in sorted(os.listdir(output_dir)):
                        fpath = os.path.join(output_dir, fname)
                        if os.path.isfile(fpath) and not fname.startswith("."):
                            zf.write(fpath, fname)
                    fig_dir = os.path.join(output_dir, "figures")
                    if os.path.isdir(fig_dir):
                        for root, dirs, files in os.walk(fig_dir):
                            for fname in sorted(files):
                                fpath = os.path.join(root, fname)
                                arcname = os.path.relpath(fpath, output_dir)
                                zf.write(fpath, arcname)

                st.download_button(
                    "\U0001f4e6 Download all results (ZIP)",
                    zip_buf.getvalue(),
                    file_name="ssign_results.zip",
                    mime="application/zip",
                    use_container_width=True,
                    key="dl_zip",
                )

                st.divider()
                st.markdown("**Individual files:**")
                for fname in sorted(os.listdir(output_dir)):
                    fpath = os.path.join(output_dir, fname)
                    if os.path.isfile(fpath) and not fname.startswith("."):
                        with open(fpath, "rb") as f:
                            st.download_button(
                                f"\U0001f4e5 {fname}",
                                f.read(),
                                file_name=fname,
                                key=f"dl_{fname}",
                            )
        else:
            st.info(
                "No results yet. Configure your pipeline in the **Pipeline Overview** tab and click **Run ssign** above."
            )
