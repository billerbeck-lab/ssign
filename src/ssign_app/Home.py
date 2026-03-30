#!/usr/bin/env python3
"""ssign — Streamlit GUI Home page.

Launch with: ssign (after pip install) or streamlit run Home.py
"""

import os
import tempfile
from pathlib import Path
import shutil

import streamlit as st
import streamlit.components.v1 as components

from ssign_app.core.runner import PipelineConfig, PipelineRunner, StepResult

# ─────────────────────────────────────────────────────────────────────
# Page config — hide deploy button
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

# Hide deploy button + custom styling
st.markdown("""
<style>
    .stDeployButton { display: none; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        padding: 8px 20px;
        border-radius: 6px 6px 0 0;
    }
    div[data-testid="stMetric"] {
        background-color: #F0F5F8;
        padding: 12px 16px;
        border-radius: 8px;
        border-left: 4px solid #1B6B93;
    }
    .tool-section {
        padding: 12px 0;
        border-bottom: 1px solid #E8EEF2;
    }
</style>
""", unsafe_allow_html=True)

# Custom connection error popup via components.html (st.markdown sanitizes <script>)
components.html('''
<script>
const observer = new MutationObserver(function(mutations) {
    document.querySelectorAll('[data-testid="stConnectionStatus"] span, .stException span').forEach(function(el) {
        if (el.textContent.includes('Is Streamlit still running')) {
            el.textContent = 'ssign server disconnected. To restart, run: ssign';
        }
    });
});
observer.observe(document.body, {childList: true, subtree: true});
</script>
''', height=0)

# ─────────────────────────────────────────────────────────────────────
# Session state defaults
# ─────────────────────────────────────────────────────────────────────

if "running" not in st.session_state:
    st.session_state.running = False
if "results" not in st.session_state:
    st.session_state.results = []
if "uploaded_files_data" not in st.session_state:
    st.session_state.uploaded_files_data = []


def _needs_run_mode_gate() -> bool:
    """Return True if a previous run exists but the user hasn't chosen a run mode yet."""
    outdir_val = st.session_state.get(
        "outdir_input",
        os.path.join(os.path.expanduser("~"), "ssign_results"),
    )
    if not outdir_val:
        return False
    progress_path = os.path.join(outdir_val, "ssign_progress.json")
    if os.path.exists(progress_path):
        return not bool(st.session_state.get("run_mode_choice"))
    return False


_GATE_MSG = (
    "A previous run was detected in the output directory. "
    "Please go to the **1. Upload** tab and choose how to proceed "
    "(Resume / Start fresh / Selective rerun) before continuing."
)

# ─────────────────────────────────────────────────────────────────────
# Header — make acronym obvious
# ─────────────────────────────────────────────────────────────────────

st.markdown(
    '<h1 style="margin-bottom: 0;">'
    '<span style="color: #1B6B93;">S</span>ecretion-<span style="color: #1B6B93;">s</span>ystem '
    '<span style="color: #1B6B93;">I</span>dentification for '
    '<span style="color: #1B6B93;">G</span>ram <span style="color: #1B6B93;">N</span>egatives</h1>',
    unsafe_allow_html=True,
)
st.caption("ssign v0.1.0 | Upload a genome, configure tools, identify secretion system substrates.")
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
        "Identifies secretion system substrates in Gram-negative "
        "bacterial genomes using secretion system detection, "
        "localization prediction, and functional annotation."
    )

    st.divider()

    st.subheader("Pipeline Steps")
    steps = [
        "1. Extract proteins from genome",
        "2. Detect secretion systems (MacSyFinder v2)",
        "3. Predict localization (DeepLocPro)",
        "4. Predict secretion type (DeepSecE)*",
        "5. Identify substrates (proximity analysis)",
        "6. Annotate substrates (optional tools)",
        "7. Generate report & figures",
    ]
    for s in steps:
        st.markdown(s)
    st.caption("*Optional — requires `pip install ssign[full]`")

    st.divider()

    st.markdown(
        "> **How it works:** ssign identifies putative secretion system "
        "substrates by *guilt by association* — proteins are flagged as "
        "candidates if they (1) are predicted to be extracellularly "
        "localized by DeepLocPro, AND (2) are genomically proximal to "
        "secretion system components detected by MacSyFinder. This does "
        "not constitute experimental evidence of secretion."
    )

    st.divider()

    st.info(
        "**ssign lite** runs annotation tools via cloud APIs.\n\n"
        "For local/HPC execution with full databases, install "
        "`ssign-full` and use the command line interface."
    )

    st.divider()
    from ssign_app import __version__
    st.caption(f"ssign v{__version__} | GPLv3 | Billerbeck Lab")


# ─────────────────────────────────────────────────────────────────────
# Tabs (all always rendered — prevents widget re-creation issues)
# ─────────────────────────────────────────────────────────────────────

tab_upload, tab_settings, tab_annotation, tab_run = st.tabs([
    "1. Upload", "2. Settings", "3. Annotation Tools", "4. Run & Results"
])


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TAB 1: Upload
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

with tab_upload:
    st.subheader("Upload Genome File")

    uploaded_files = st.file_uploader(
        "Choose genome file(s)",
        type=["gbff", "gbk", "gb", "gff", "gff3", "fasta", "fna", "fa"],
        help="GenBank (.gbff/.gbk), GFF3 (.gff3), or FASTA contigs (.fasta/.fna). Upload multiple files for batch processing.",
        accept_multiple_files=True,
    )

    if uploaded_files:
        # Clear resolved taxonomy if files changed
        current_names = sorted(f.name for f in uploaded_files)
        prev_names = st.session_state.get("_last_uploaded_names", [])
        if current_names != prev_names:
            st.session_state.pop("organism_resolved", None)
            st.session_state.pop("detected_organism", None)
            st.session_state["_last_uploaded_names"] = current_names

        st.session_state.uploaded_files_data = uploaded_files

        # Sample names
        if len(uploaded_files) == 1:
            sample_name = uploaded_files[0].name
            for suffix in ['.gbff', '.gbk', '.gb', '.gff3', '.gff', '.fasta', '.fna', '.fa']:
                sample_name = sample_name.replace(suffix, '')
            sample_name = sample_name.replace('_genomic', '')
            st.text_input("Sample name (auto-detected, editable)", value=sample_name,
                           key="sample_name")
        else:
            st.markdown("**Sample names** (auto-derived from filenames):")
            sample_names = []
            for f in uploaded_files:
                sn = f.name
                for suffix in ['.gbff', '.gbk', '.gb', '.gff3', '.gff', '.fasta', '.fna', '.fa']:
                    sn = sn.replace(suffix, '')
                sn = sn.replace('_genomic', '')
                sample_names.append(sn)
                st.caption(f"- {sn}")
            st.session_state.sample_names = sample_names

        # Auto-detect organism from first GenBank file
        first_gbk = None
        for f in uploaded_files:
            ext_lower = Path(f.name).suffix.lower()
            if ext_lower in ('.gbff', '.gbk', '.gb'):
                first_gbk = f
                break

        # Re-resolve if not cached, or if previous resolve failed (species=None)
        _prev = st.session_state.get("organism_resolved")
        _needs_resolve = _prev is None or (_prev and _prev.get("species") is None)
        if first_gbk and _needs_resolve:
            try:
                from Bio import SeqIO
                import io
                first_gbk.seek(0)
                content = first_gbk.read().decode("utf-8", errors="replace")
                first_gbk.seek(0)
                record = next(SeqIO.parse(io.StringIO(content), "genbank"))
                org_name = record.annotations.get("organism", "").strip()
                if not org_name:
                    org_name = record.annotations.get("source", "").strip()
                # Check source feature /organism qualifier for fuller name
                if not org_name or len(org_name.split()) < 2:
                    for feat in record.features:
                        if feat.type == "source":
                            src_org = feat.qualifiers.get("organism", [""])[0].strip()
                            if src_org and len(src_org.split()) >= 2:
                                org_name = src_org
                            break
                # Fallback: infer from filename
                if not org_name or len(org_name.split()) < 2:
                    stem = Path(first_gbk.name).stem
                    for sfx in ('_genomic', '_protein', '_cds', '_rna'):
                        stem = stem.replace(sfx, '')
                    parts = stem.replace('_', ' ').split()
                    if (len(parts) >= 2
                            and parts[0][0].isupper()
                            and parts[1][0].islower()
                            and parts[1].isalpha()):
                        inferred = f"{parts[0]} {parts[1]}"
                        if not org_name or len(org_name.split()) < 2:
                            org_name = inferred
                if org_name:
                    st.session_state.detected_organism = org_name
                    try:
                        from ssign_app.scripts.resolve_taxonomy import resolve_organism
                        tax_info = resolve_organism(org_name)
                        st.session_state.organism_resolved = tax_info
                    except Exception:
                        st.session_state.organism_resolved = {"species": None, "genus": None}
            except Exception:
                pass

    col1, col2 = st.columns(2)
    with col1:
        st.info(
            "Recommended to use GenBank (.gbff) files, which generally include "
            "gene names and functional annotations from tools like Bakta.\n\n"
            "Raw FASTA contigs will use Prodigal (no additional setup) or "
            "Bakta (richer annotation, requires database download ~2GB)."
        )
        bakta_available = shutil.which("bakta") is not None
        use_bakta = st.checkbox(
            "Use Bakta for ORF prediction (raw FASTA input only)",
            value=False, key="use_bakta",
            disabled=not bakta_available,
        )
        if bakta_available:
            if use_bakta:
                st.text_input("Bakta database path", key="bakta_db_path",
                              placeholder="/path/to/bakta_db")
        else:
            st.caption(
                "Bakta not detected. Install with:\n"
                "`pip install bakta && bakta_db download --output /path/db --type light`"
            )
    with col2:
        outdir = st.text_input(
            "Output directory",
            value=os.path.join(os.path.expanduser("~"), "ssign_results"),
            key="outdir_input",
        )

    # ── Previous run detection ──
    outdir_val = st.session_state.get("outdir_input", "")
    progress_path = os.path.join(outdir_val, "ssign_progress.json") if outdir_val else ""
    if progress_path and os.path.exists(progress_path):
        try:
            import json
            with open(progress_path) as pf:
                prev = json.load(pf)
            prev_steps = prev.get("steps", [])
            n_done = sum(1 for s in prev_steps if s.get("success"))
            n_total = len(prev_steps)
            prev_time = prev.get("timestamp", "unknown")
            prev_sample = prev.get("sample_id", "unknown")

            st.divider()
            st.markdown("#### Previous Run Detected")
            col_info, col_action = st.columns([3, 2])
            with col_info:
                st.info(
                    f"**Sample:** {prev_sample}  \n"
                    f"**Progress:** {n_done}/{n_total} steps completed  \n"
                    f"**Last run:** {prev_time}"
                )
            with col_action:
                run_mode = st.radio(
                    "How to proceed?",
                    ["Resume (skip completed steps)",
                     "Start fresh (rerun everything)",
                     "Selective rerun (choose steps)"],
                    key="run_mode_choice",
                    index=0,
                )
            if run_mode and "Selective" in run_mode:
                st.markdown("**Select steps to rerun** (unchecked = keep previous result):")
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
# TAB 2: Settings
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

with tab_settings:
    if _needs_run_mode_gate():
        st.warning(_GATE_MSG)
    else:
        st.subheader("Detection Parameters")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Secretion System Detection (MacSyFinder v2)**")
            wholeness = st.slider(
                "System completeness threshold",
                0.0, 1.0, 0.8, 0.05,
                help="MacSyFinder wholeness score minimum. 0.8 = at least 80% of "
                     "expected components must be present.",
                key="wholeness",
            )
            window = st.slider(
                "Proximity window (genes)",
                1, 15, 3,
                help="How many genes upstream/downstream of each SS component to "
                     "search for putative substrates.",
                key="window",
            )

        with col2:
            st.markdown("**Prediction Thresholds**")
            conf_threshold = st.slider(
                "DeepLocPro extracellular threshold",
                0.0, 1.0, 0.8, 0.05,
                help="Minimum probability for a protein to be called extracellular.",
                key="conf",
            )
            fraction_correct = st.slider(
                "Required fraction correctly localized",
                0.0, 1.0, 0.8, 0.05,
                help="Fraction of SS components that must have correct predicted "
                     "localization for the system to be considered valid.",
                key="frac",
            )

        st.divider()

        st.subheader("System Filtering")
        st.markdown(
            "Select which secretion system types to **exclude** from analysis. "
            "T3SS is excluded by default because DeepSecE T3SS predictions are "
            "unreliable (mostly flagellar misclassification)."
        )

        # All MacSyFinder system types
        all_system_types = [
            "T1SS", "T2SS", "T3SS", "T4SS", "T5aSS", "T5bSS", "T5cSS",
            "T6SSi", "T6SSii", "T6SSiii", "T9SS",
            "Flagellum", "Tad", "pT4SSt",
        ]
        excluded = st.multiselect(
            "Exclude these system types",
            all_system_types,
            default=["Flagellum", "Tad", "T3SS"],
            key="excluded",
        )

        st.divider()

        # ── DeepLocPro & SignalP ──

        st.subheader("Localization & Signal Peptide Prediction")

        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**DeepLocPro** (subcellular localization)")
            dlp_mode = st.radio(
                "DeepLocPro mode",
                ["BioLib cloud (no install needed)", "Local install (DTU license)"],
                key="dlp_mode",
                help="BioLib cloud: free, no license required. "
                     "Local: faster, requires free DTU academic license.",
            )
            if "Local" in dlp_mode:
                dlp_found = shutil.which("deeplocpro") is not None
                dlp_path = st.text_input("DeepLocPro install path", key="dlp_path")
                if not dlp_found and not dlp_path:
                    st.warning(
                        "deeplocpro not found in PATH. Local mode requires a DTU academic "
                        "license and ~5GB model download + GPU recommended."
                    )
            else:
                st.caption(
                    "Cloud recommended (~5-10 min per genome, no setup needed)."
                )
                dlp_path = ""

        with col2:
            st.markdown("**SignalP 6.0** (signal peptides)")
            run_signalp = st.checkbox("Enable SignalP 6.0", value=False, key="run_signalp",
                                       help="Predicts Sec/SPI, Sec/SPII, Tat/SPI, Tat/SPII signal peptides")
            if run_signalp:
                sp_mode = st.radio(
                    "SignalP mode",
                    ["BioLib cloud (no install needed)", "Local install (DTU license)"],
                    key="sp_mode",
                )
                if "Local" in sp_mode:
                    sp_found = shutil.which("signalp6") is not None
                    sp_path = st.text_input("SignalP install path", key="sp_path")
                    if not sp_found and not sp_path:
                        st.warning(
                            "signalp6 not found in PATH. Local mode requires a DTU academic license."
                        )
                else:
                    st.caption("Cloud recommended (~2-5 min per genome, no setup needed).")
                    sp_path = ""
                with st.expander("SignalP threshold", expanded=False):
                    st.slider("Min. probability", 0.0, 1.0, 0.5, 0.05,
                              key="sp_min_prob",
                              help="Minimum SignalP probability to call a signal peptide")
            else:
                sp_mode = "BioLib cloud"
                sp_path = ""

        st.divider()

        # ── DeepSecE ──

        st.subheader("Secretion Type Prediction (DeepSecE)")

        deepsece_available = False
        try:
            import DeepSecE
            deepsece_available = True
        except ImportError:
            pass

        if deepsece_available:
            run_deepsece = st.checkbox(
                "Enable DeepSecE", value=True, key="run_deepsece",
                help="Predicts if proteins are secreted, and by which SS type. "
                     "Does not work for autotransporters. Helps to cross-validate "
                     "DeepLocPro and MacSyFinder results, may also yield more hits.",
            )
            st.success("DeepSecE is installed and ready.")
            if run_deepsece:
                with st.expander("DeepSecE threshold", expanded=False):
                    st.slider("Min. probability", 0.0, 1.0, 0.8, 0.05,
                              key="dse_min_prob",
                              help="Minimum DeepSecE probability to call a protein as secreted")
        else:
            run_deepsece = st.checkbox(
                "Enable DeepSecE", value=False, key="run_deepsece",
                disabled=True,
                help="Predicts if proteins are secreted, and by which SS type. "
                     "Does not work for autotransporters. Helps to cross-validate "
                     "DeepLocPro and MacSyFinder results, may also yield more hits.",
            )
            st.info(
                "**DeepSecE is not installed** but is recommended. It predicts which "
                "secretion system type each protein is secreted by, providing an "
                "independent cross-check against MacSyFinder. This improves confidence "
                "in substrate identification and helps catch false positives.\n\n"
                "To install (adds ~2 GB for PyTorch):\n"
                "```\npip install ssign[full]\n```\n\n"
                "ssign works without DeepSecE \u2014 it will use DeepLocPro alone for "
                "substrate identification, which is still effective but loses the "
                "cross-validation benefit."
            )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TAB 3: Annotation Tools
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

with tab_annotation:
    if _needs_run_mode_gate():
        st.warning(_GATE_MSG)
    else:
        st.subheader("Optional Annotation Tools")
        st.markdown(
            "These tools annotate your substrate candidates with functional information. "
            "All run via **cloud APIs** by default (no install needed, slower). "
            "Uncheck to skip. All are optional."
        )

        st.divider()

        # ── BLASTp ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_blastp = st.checkbox("BLASTp", value=True, key="run_blastp",
                                      help="Search NCBI nr for homologous proteins")
        with col_info:
            if run_blastp:
                st.caption("NCBI web API | Rate limit: 3 req/sec | ~30-60 min per genome")

                # Taxonomy exclusion — explanation
                with st.expander("Why exclude self-hits?", expanded=False):
                    st.markdown(
                        "If your organism is already in the NCBI nr database, BLASTp will "
                        "return your own proteins (or close relatives) as the top hits. "
                        "This duplicates the annotation you already have from Prokka/Bakta "
                        "instead of providing **new functional information** from more "
                        "distant homologs.\n\n"
                        "- **Exclude species** removes only your exact organism, so you still "
                        "see hits from closely related species in the same genus.\n"
                        "- **Exclude genus** removes the entire genus, forcing BLASTp to find "
                        "homologs from more distantly related organisms. This is useful when "
                        "many genomes from your genus are in nr.\n"
                        "- **No exclusion** keeps all hits — use this if your organism is novel "
                        "or not yet in NCBI databases."
                    )

                # Taxonomy exclusion options
                tax_info = st.session_state.get("organism_resolved", {})
                species_info = tax_info.get("species") if tax_info else None
                genus_info = tax_info.get("genus") if tax_info else None
                detected_org = st.session_state.get("detected_organism", "")

                if detected_org:
                    st.info(f"Detected organism: **{detected_org}**")

                # Build exclusion radio options
                exclusion_options = []
                option_taxids = {}

                if species_info:
                    sp_label = (
                        f"Exclude input species: *{species_info['name']}* "
                        f"(taxid {species_info['taxid']})"
                    )
                    exclusion_options.append(sp_label)
                    option_taxids[sp_label] = species_info["taxid"]

                if genus_info:
                    gen_label = (
                        f"Exclude input genus: *{genus_info['name']}* "
                        f"(taxid {genus_info['taxid']})"
                    )
                    exclusion_options.append(gen_label)
                    option_taxids[gen_label] = genus_info["taxid"]

                exclusion_options.append("Custom taxonomy ID(s)")
                option_taxids["Custom taxonomy ID(s)"] = "__custom__"
                exclusion_options.append("No exclusion (include all hits)")
                option_taxids["No exclusion (include all hits)"] = ""

                # Default to species exclusion if available
                exclusion_choice = st.radio(
                    "BLASTp taxonomy exclusion",
                    exclusion_options,
                    index=0,
                    key="blastp_exclusion_mode",
                )

                resolved_taxid = option_taxids.get(exclusion_choice, "")

                if resolved_taxid == "__custom__":
                    st.text_input(
                        "Enter NCBI taxonomy ID(s)",
                        key="blastp_taxid",
                        placeholder="e.g. 339 or 339,340,338",
                        help="Comma-separate multiple taxonomy IDs to exclude "
                             "several organisms. Look up taxids at "
                             "https://www.ncbi.nlm.nih.gov/taxonomy",
                    )
                elif resolved_taxid:
                    st.session_state.blastp_taxid = resolved_taxid
                else:
                    st.session_state.blastp_taxid = ""

                if not detected_org and not species_info:
                    st.caption(
                        "No organism auto-detected (FASTA input has no organism metadata). "
                        "Select **Custom** to enter your organism's NCBI taxonomy ID, "
                        "or **No exclusion** if your organism is not yet in NCBI databases."
                    )

                # BLASTp threshold sliders (visible, not hidden in expander)
                st.markdown("**BLASTp result filters:**")
                bc1, bc2, bc3 = st.columns(3)
                with bc1:
                    st.slider("Min. % identity", 0, 100, 80, 5,
                              key="blastp_pident",
                              help="Minimum percent identity to keep a BLASTp hit")
                with bc2:
                    st.slider("Min. query coverage (%)", 0, 100, 80, 5,
                              key="blastp_qcov",
                              help="Minimum query coverage to keep a BLASTp hit")
                with bc3:
                    st.number_input("E-value threshold", value=1e-5,
                                    format="%.0e", key="blastp_evalue",
                                    help="Maximum e-value for BLASTp hits")

        # ── HHpred ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_hh = st.checkbox("HHpred (Pfam + PDB)", value=True, key="run_hh",
                                  help="Remote homology detection via MPI Toolkit")
        with col_info:
            if run_hh:
                st.caption("MPI Toolkit API | Rate limit: 200 jobs/hr | ~45-90 min per genome")
                st.slider("Min. probability (%)", 0, 100, 40, 5,
                          key="hhpred_min_prob",
                          help="Minimum HHpred probability to keep a hit. "
                               "Default 40% balances sensitivity and specificity.")

        # ── InterProScan ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_iprs = st.checkbox("InterProScan", value=True, key="run_iprs",
                                    help="Domain and GO annotation via EBI")
        with col_info:
            if run_iprs:
                st.caption("EBI REST API | Rate limit: 25k seq/day | ~20-40 min per genome")
                st.number_input("E-value threshold", value=1e-5,
                                format="%.0e", key="iprs_evalue",
                                help="Maximum e-value for InterProScan domain hits")

        # ── ProtParam ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_pp = st.checkbox("ProtParam", value=True, key="run_pp",
                                  help="Physicochemical properties (MW, pI, GRAVY, etc.)")
        with col_info:
            if run_pp:
                st.caption("Local (BioPython) | No database needed | Instant")

        st.divider()

        st.subheader("Advanced / Large Database Tools")
        st.markdown(
            "These require significant local storage or are experimental. "
            "Disabled by default in ssign lite."
        )

        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            run_foldseek = st.checkbox("Foldseek", value=False, key="run_fs",
                                        help="Structural homology search")
        with col_info:
            if run_foldseek:
                fs_mode = st.radio("Foldseek mode", ["Web API", "Local"], key="fs_mode")
                if fs_mode == "Local":
                    st.text_input("Foldseek DB path (~10GB)", key="fs_db")
                else:
                    st.caption("Foldseek web API | Requires AlphaFold DB structures")
                fc1, fc2 = st.columns(2)
                with fc1:
                    st.number_input("E-value threshold", value=1e-3,
                                    format="%.0e", key="fs_evalue",
                                    help="Maximum e-value for Foldseek hits")
                with fc2:
                    st.slider("Min. TM-score", 0.0, 1.0, 0.5, 0.05,
                              key="fs_tmscore",
                              help="Minimum TM-score to keep a Foldseek structural hit")

        # ── pLM-BLAST ──
        col_check, col_info = st.columns([1.5, 3.5])
        with col_check:
            st.checkbox("pLM-BLAST (ECOD70)", value=False, key="run_plm",
                         disabled=True,
                         help="Structure-based remote homology detection")
        with col_info:
            st.caption(
                "Local only — MPI web API currently unavailable. "
                "Requires local pLM-BLAST installation + ECOD70 database (~10GB). "
                "Not yet integrated into ssign pipeline."
            )

        st.divider()

        # ── Figures ──
        st.subheader("Figures")
        st.markdown("Select which figures to generate after the pipeline finishes.")

        st.checkbox("Category distribution", value=True, key="fig_category")
        st.checkbox("SS composition", value=True, key="fig_ss_comp")
        st.checkbox("Tool coverage heatmap", value=True, key="fig_tool_heatmap")
        st.checkbox("Substrate count per genome", value=True, key="fig_substrate_count")
        st.checkbox("Functional annotation summary", value=True, key="fig_func_summary")

        st.divider()

        # ── Ortholog Group Assignment ──
        st.subheader("Ortholog Group Assignment")

        blastp_available = shutil.which("blastp") is not None

        if blastp_available:
            st.success("BLAST+ is installed and ready for ortholog grouping.")
            st.markdown(
                "Substrates are grouped into ortholog groups using all-vs-all "
                "BLASTp with single-linkage clustering."
            )
            col1, col2 = st.columns(2)
            with col1:
                st.slider("Min. % identity for ortholog grouping",
                          0, 100, 40, 5,
                          key="og_min_pident",
                          help="Minimum percent identity between two proteins to "
                               "consider them orthologs. Default 40% is permissive; "
                               "increase for stricter groups.")
            with col2:
                st.slider("Min. query coverage (%) for ortholog grouping",
                          0, 100, 70, 5,
                          key="og_min_qcov",
                          help="Minimum query coverage to consider a BLASTp hit "
                               "as an ortholog relationship. Default 70%.")
        else:
            st.info(
                "**BLAST+ is not installed** (needed for ortholog grouping).\n\n"
                "Ortholog grouping clusters substrate proteins across genomes into "
                "groups of related proteins, useful for comparative analysis.\n\n"
                "To install:\n"
                "```\nsudo apt install ncbi-blast+\n```\n"
                "Or: `conda install -c bioconda blast`\n\n"
                "ssign works without BLAST+ — ortholog grouping will be skipped."
            )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TAB 4: Run & Results
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
            issues.append("No genome file uploaded. Go back to Step 1.")

        can_run = len(issues) == 0

        if can_run:
            st.success(
                f"Ready to run on **{len(uploaded_files)} file(s)**: "
                + ", ".join(uf.name for uf in uploaded_files)
            )
        else:
            for issue in issues:
                st.warning(issue)

        # Estimated time (rough, scales per genome)
        st.divider()
        n_genomes = len(uploaded_files) or 1
        st.markdown("**Estimated run time** (per genome):")

        time_parts = ["MacSyFinder v2: ~2-5 min", "DeepLocPro (BioLib): ~5-10 min"]
        if st.session_state.get("run_deepsece"):
            time_parts.append("DeepSecE: ~2-5 min")
        if st.session_state.get("run_blastp"):
            time_parts.append("BLASTp (NCBI API): ~30-60 min")
        if st.session_state.get("run_hh"):
            time_parts.append("HHpred (MPI API): ~45-90 min")
        if st.session_state.get("run_iprs"):
            time_parts.append("InterProScan (EBI API): ~20-40 min")

        for t in time_parts:
            st.markdown(f"- {t}")

        total_min = 15  # base per genome
        if st.session_state.get("run_blastp"):
            total_min += 45
        if st.session_state.get("run_hh"):
            total_min += 60
        if st.session_state.get("run_iprs"):
            total_min += 30
        total_min *= n_genomes
        genome_note = f" for {n_genomes} genome(s)" if n_genomes > 1 else ""
        st.info(f"**Estimated total: ~{total_min}-{total_min*2} minutes{genome_note}** (rough estimate, depends on API load)")

        st.divider()

        st.warning(
            "**Before you run:**\n"
            "- Closing the **browser tab** is safe — the pipeline continues in the background and results are saved.\n"
            "- Closing the **terminal** will stop the pipeline. Partial progress is saved to the output directory.\n"
            "- If interrupted, re-run with the same output directory and **Resume** enabled to continue.\n"
            "- A full run with all tools enabled takes **1-3 hours** depending on genome size and API load."
        )

        # Resume — use run mode from Upload tab if set, otherwise detect here
        run_mode = st.session_state.get("run_mode_choice", "")
        outdir_val = st.session_state.get("outdir_input", "")
        progress_exists = os.path.exists(os.path.join(outdir_val, "ssign_progress.json")) if outdir_val else False

        if progress_exists and not run_mode:
            resume_enabled = st.checkbox(
                "Resume from previous run (skip completed steps)",
                value=True,
                key="resume_run",
                help="Steps that completed successfully will be skipped.",
            )
            st.caption("Previous progress detected. Configure run mode in the Upload tab for more options.")
        elif "Resume" in run_mode:
            st.info("Resuming from previous run — completed steps will be skipped.")
            st.session_state.resume_run = True
        elif "fresh" in run_mode.lower() if run_mode else False:
            st.info("Starting fresh — all steps will rerun.")
            st.session_state.resume_run = False
        elif "Selective" in run_mode:
            st.info("Selective rerun — only checked steps from the Upload tab will rerun.")
            st.session_state.resume_run = True  # resume base, but clear selected steps
        else:
            st.session_state.resume_run = False

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
                with open(dest, 'wb') as out:
                    out.write(uf.read())
                input_paths.append(dest)
                original_filenames.append(uf.name)

            # Process each genome (sequentially for now)
            all_results = []
            for file_idx, input_path in enumerate(input_paths):
                if len(input_paths) > 1:
                    sample_id = st.session_state.get("sample_names", ["sample"])[file_idx] if file_idx < len(st.session_state.get("sample_names", [])) else f"sample_{file_idx+1}"
                    st.subheader(f"Running genome {file_idx+1}/{len(input_paths)}: {sample_id}")
                else:
                    sample_id = st.session_state.get("sample_name", "sample")

                # Build config
                orig_fname = original_filenames[file_idx] if file_idx < len(original_filenames) else ""
                config = PipelineConfig(
                    input_path=input_path,
                    original_filename=orig_fname,
                    sample_id=sample_id,
                    outdir=st.session_state.get("outdir_input", "./results"),
                    run_bakta=st.session_state.get("use_bakta", False),
                    bakta_db=st.session_state.get("bakta_db_path", ""),
                    wholeness_threshold=st.session_state.get("wholeness", 0.8),
                    excluded_systems=st.session_state.get("excluded", ["Flagellum", "Tad", "T3SS"]),
                    conf_threshold=st.session_state.get("conf", 0.8),
                    proximity_window=st.session_state.get("window", 3),
                    required_fraction_correct=st.session_state.get("frac", 0.8),
                    deeplocpro_mode="local" if "Local" in st.session_state.get("dlp_mode", "") else "remote",
                    deeplocpro_path=st.session_state.get("dlp_path", ""),
                    skip_deepsece=not st.session_state.get("run_deepsece", False),
                    skip_signalp=not st.session_state.get("run_signalp", False),
                    signalp_mode="local" if "Local" in st.session_state.get("sp_mode", "") else "remote",
                    signalp_path=st.session_state.get("sp_path", ""),
                    skip_blastp=not st.session_state.get("run_blastp", True),
                    blastp_mode="remote",  # Always remote in lite
                    blastp_db="",
                    blastp_exclude_taxid=st.session_state.get("blastp_taxid", ""),
                    blastp_min_pident=float(st.session_state.get("blastp_pident", 80)),
                    blastp_min_qcov=float(st.session_state.get("blastp_qcov", 80)),
                    blastp_evalue=float(st.session_state.get("blastp_evalue", 1e-5)),
                    skip_hhsuite=not st.session_state.get("run_hh", False),
                    hhsuite_mode="remote",
                    hhsuite_pfam_db="",
                    hhsuite_pdb70_db="",
                    hhpred_min_probability=float(st.session_state.get("hhpred_min_prob", 40)),
                    skip_interproscan=not st.session_state.get("run_iprs", True),
                    interproscan_mode="remote",
                    interproscan_db="",
                    skip_foldseek=not st.session_state.get("run_fs", False),
                    foldseek_db=st.session_state.get("fs_db", ""),
                    skip_plmblast=not st.session_state.get("run_plm", False),
                    plmblast_db=st.session_state.get("plm_db", ""),
                    skip_protparam=not st.session_state.get("run_pp", True),
                    # Figure toggles
                    fig_category=st.session_state.get("fig_category", True),
                    fig_ss_comp=st.session_state.get("fig_ss_comp", True),
                    fig_tool_heatmap=st.session_state.get("fig_tool_heatmap", True),
                    fig_substrate_count=st.session_state.get("fig_substrate_count", True),
                    fig_func_summary=st.session_state.get("fig_func_summary", True),
                    # Tool thresholds
                    interproscan_evalue=float(st.session_state.get("iprs_evalue", 1e-5)),
                    foldseek_evalue=float(st.session_state.get("fs_evalue", 1e-3)),
                    foldseek_min_tmscore=float(st.session_state.get("fs_tmscore", 0.5)),
                    deepsece_min_prob=float(st.session_state.get("dse_min_prob", 0.8)),
                    signalp_min_prob=float(st.session_state.get("sp_min_prob", 0.5)),
                    ortholog_min_pident=float(st.session_state.get("og_min_pident", 40)),
                    ortholog_min_qcov=float(st.session_state.get("og_min_qcov", 70)),
                )

                # Progress display
                progress_bar = st.progress(0)
                status_text = st.empty()

                def update_progress(step, pct, msg):
                    progress_bar.progress(pct / 100)
                    status_text.markdown(f"**{step}** \u2014 {msg}")

                # Run pipeline
                runner = PipelineRunner(config, progress_callback=update_progress)

                with st.spinner(f"Running ssign pipeline on {sample_id}..."):
                    results = runner.run(resume=st.session_state.get("resume_run", False))

                all_results.extend(results)

            st.session_state.results = all_results
            st.session_state.running = False
            st.session_state.output_dir = st.session_state.get("outdir_input", "./results")

            # Show results summary
            n_success = sum(1 for r in all_results if r.success)
            n_total = len(all_results)
            n_genomes = len(input_paths)

            if n_success == n_total:
                st.success(f"Pipeline completed successfully for {n_genomes} genome(s)! ({n_success}/{n_total} steps)")
            else:
                st.warning(f"Pipeline finished with issues for {n_genomes} genome(s) ({n_success}/{n_total} steps succeeded)")

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
                # Count substrates from results CSV (second section, after blank row)
                n_substrates = "—"
                if output_dir:
                    import glob as _glob
                    csvs = _glob.glob(os.path.join(output_dir, "*_results.csv"))
                    if csvs:
                        try:
                            with open(csvs[0]) as _f:
                                lines = _f.readlines()
                            # Find blank separator row; substrates are after it
                            for idx, line in enumerate(lines):
                                if line.strip() == '':
                                    # Substrates = lines after blank row, minus header
                                    sub_lines = [l for l in lines[idx+2:] if l.strip()]
                                    n_substrates = str(len(sub_lines))
                                    break
                            else:
                                # No blank row = only substrates (no systems section)
                                n_substrates = str(max(0, len(lines) - 1))
                        except Exception:
                            pass
                st.metric("Substrates found", n_substrates)

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
                import zipfile, io
                zip_buf = io.BytesIO()
                with zipfile.ZipFile(zip_buf, 'w', zipfile.ZIP_DEFLATED) as zf:
                    for fname in sorted(os.listdir(output_dir)):
                        fpath = os.path.join(output_dir, fname)
                        if os.path.isfile(fpath) and not fname.startswith('.'):
                            zf.write(fpath, fname)
                    # Include figures subdirectory
                    fig_dir = os.path.join(output_dir, "figures")
                    if os.path.isdir(fig_dir):
                        for fname in sorted(os.listdir(fig_dir)):
                            fpath = os.path.join(fig_dir, fname)
                            if os.path.isfile(fpath):
                                zf.write(fpath, f"figures/{fname}")

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
                    if os.path.isfile(fpath) and not fname.startswith('.'):
                        with open(fpath, 'rb') as f:
                            st.download_button(
                                f"\U0001f4e5 {fname}",
                                f.read(),
                                file_name=fname,
                                key=f"dl_{fname}",
                            )
            else:
                st.info("No results yet. Configure your run and click **Run ssign** above.")
