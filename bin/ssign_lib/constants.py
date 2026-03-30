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
BLASTP_MIN_QCOV = 80    # Minimum query coverage

# --- Foldseek ---
FOLDSEEK_EVALUE = 0.001
FOLDSEEK_TMSCORE = 0.8

# --- T5aSS domain classification ---
MIN_PASSENGER_LENGTH = 100  # aa — below this = "minimal passenger"
LINKER_LENGTH = 30          # aa — alpha-helix linker between passenger and barrel

# --- System filtering ---
DEFAULT_EXCLUDED_SYSTEMS = ["Flagellum", "Tad", "T3SS"]

# --- DeepSecE to MacSyFinder SS type mapping ---
# DeepSecE predicts broad types; MacSyFinder uses specific names.
DSE_TO_MACSYFINDER = {
    "T1SS": ["T1SS"],
    "T2SS": ["T2SS"],
    "T3SS": ["T3SS"],
    "T4SS": ["pT4SSt", "T4SS"],
    "T6SS": ["T6SSi", "T6SS"],
}
