"""Central constants and default thresholds for the ssign pipeline.

All configurable thresholds are defined here so that (a) the Python scripts
have sensible defaults when run standalone and (b) the Nextflow params
override them at invocation time via CLI arguments.
"""

# --- DeepLocPro ---
CONF_THRESHOLD = 0.8  # Minimum extracellular probability

# --- MacSyFinder ---
WHOLENESS_THRESHOLD = 0.8  # Minimum system completeness
REQUIRED_FRACTION_CORRECT = 0.8  # Fraction of SS components correctly localized

# --- Proximity ---
PROXIMITY_WINDOW = 3  # +/- N genes per SS component

# --- Structure quality ---
PLDDT_THRESHOLD = 70  # Minimum mean pLDDT for structure acceptance

# --- BLASTp ---
BLASTP_MIN_PIDENT = 80  # Minimum percent identity
BLASTP_MIN_QCOV = 80  # Minimum query coverage

# --- Tool subprocess timeout ---
# Generic upper bound for any external bioinformatics tool ssign shells
# out to (BLAST, IPS, EggNOG, Bakta, pLM-BLAST, DLP, SignalP). 4h is
# deliberately generous — pathological inputs (giant proteome, huge DB,
# slow disk) can stretch any of these tools well past the typical
# 10-60 min runtime. A timeout firing means the tool is wedged or the
# input is unrealistically large; give up rather than block the rest of
# the pipeline forever. Bump if you see legitimate timeouts.
# (HH-suite has its own per-stage timeouts a few sections below — those
# are tighter because they fire per-protein, not per-batch.)
TOOL_TIMEOUT_S = 14400  # 4h

# --- HH-suite ---
# HHR Prob (0-100) cutoff for keeping the top-1 hit per DB. Söding-lab
# guidance: ≥95 near-certain homolog, ≥50 worth considering. ssign default
# of 80 is permissive-but-meaningful. Filtering on Prob, not E-value.
HHSUITE_MIN_PROB = 80.0
# Three iterations balance sensitivity for remote homologs against
# runtime; ssign exists to surface distant homologs of secretion substrates.
HHBLITS_ITERATIONS = 3
# Per-protein subprocess timeouts. Empirical wall-time for a 1500-aa
# autotransporter against UniRef30 (~204 GB) with -n 3 iterations on a
# 4-thread CPU run is 30-60 min for hhblits; hhsearch against Pfam/PDB70
# is faster (~5-15 min) but capped at half the hhblits budget so total
# wall-time per protein stays bounded.
HHBLITS_TIMEOUT_S = 3600
HHSEARCH_TIMEOUT_S = 1800

# --- T5aSS domain classification ---
MIN_PASSENGER_LENGTH = 100  # aa — below this = "minimal passenger"
LINKER_LENGTH = 30  # aa — alpha-helix linker between passenger and barrel

# --- System filtering ---
DEFAULT_EXCLUDED_SYSTEMS = ["Flagellum", "Tad", "T3SS"]

# --- T5SS per-component DLP rules ---
# Per-component biology rationale lives in cross_validate_predictions.py's
# module docstring (the consumer); this dict is the data.
#
# FRAGILE: gene_names are TXSScan v2.1 model IDs. If TXSScan renames a
# component (Pfam succession, model rebuild), this dict goes stale
# silently — affected subtypes revert to the strict ext-only rule and
# under-call passengers. If this breaks: pin TXSScan or update keys to
# match the new model names.
T5SS_COMPONENT_RULES = {
    ("T5aSS", "T5aSS_PF03797"): ("extracellular_prob", "outer_membrane_prob"),
    ("T5bSS", "T5bSS_translocator"): ("outer_membrane_prob",),
    ("T5cSS", "T5cSS_PF03895"): ("extracellular_prob", "outer_membrane_prob"),
}

# --- DeepSecE to MacSyFinder SS type mapping ---
# DeepSecE predicts broad types; MacSyFinder uses specific names.
DSE_TO_MACSYFINDER = {
    "T1SS": ["T1SS"],
    "T2SS": ["T2SS"],
    "T3SS": ["T3SS"],
    "T4SS": ["pT4SSt", "T4SS"],
    "T6SS": ["T6SSi", "T6SS"],
}
