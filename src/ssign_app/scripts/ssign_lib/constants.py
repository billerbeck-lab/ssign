"""Central constants and default thresholds for the ssign pipeline.

All configurable thresholds are defined here so that (a) the Python scripts
have sensible defaults when run standalone and (b) the Nextflow params
override them at invocation time via CLI arguments.
"""

import os
import re  # used by display_type()

# --- DeepLocPro ---
CONF_THRESHOLD = 0.8  # Minimum extracellular probability

# Max protein length DeepLocPro will be asked to embed. Its ESM-style transformer
# has ~O(L^2) attention memory, so a single mega-protein OOMs the GPU and, since
# DeepLocPro is a core step, crashes the whole genome. The motivating case:
# plu2670 / kolossin synthetase = 16,367 aa, the largest known prokaryotic protein
# (UniProt Q7N3P5). Sequences longer than this are withheld from the model and
# marked "not predicted" — they are cytoplasmic NRPS/PKS megasynthases and giant
# assembly lines, not secretion substrates. Override with SSIGN_DEEPLOCPRO_MAX_AA.
DEEPLOCPRO_MAX_AA = int(os.environ.get("SSIGN_DEEPLOCPRO_MAX_AA", "5000"))

# --- MacSyFinder ---
WHOLENESS_THRESHOLD = 0.8  # Minimum system completeness
REQUIRED_FRACTION_CORRECT = 0.8  # Fraction of SS components correctly localized

# --- Proximity ---
PROXIMITY_WINDOW = 3  # +/- N genes per SS component

# --- Enrichment (circular-shift permutation) ---
# The per-run enrichment test (opt-in via --enrichment-stats) asks, per SS type,
# whether secreted-predicted proteins cluster around that type's components more
# than a genome-structure-preserving null. The null is the exact set of circular
# rotations of the predictor's gene-ordered positivity vector (offset 0 = observed),
# so no Monte-Carlo sampling is needed for a single genome. Window size and the
# positivity threshold reuse PROXIMITY_WINDOW / CONF_THRESHOLD (single source of truth).
ENRICH_TOOLS = ("DLP", "DSE", "SignalP")  # the enrichment predictors (PLM-Effector excluded)
# Pseudo-predictor: "secreted by DLP OR DSE" (only the valid predictors per type;
# DSE dropped where unreliable). One combined fold/significance bar per SS type;
# emitted alongside DLP/DSE in the stats TSV, BH-corrected as its own family.
ENRICH_COMBINED_TOOL = "COMBINED"
# Above this genome length, subsample rotation offsets for the stored null (the
# fold + p stay exact; only the histogram dump is thinned). Bacterial genomes
# (~5k genes) never hit this, so in practice the full exact null is always used.
ENRICH_MAX_NULL = 200_000
# SS types (display form, after subtype collapse) tested with a +/-W window around
# each component. Autotransporters are the complement: the component IS the substrate
# (self-secreting), so they're tested by self-detection. A display type in neither set
# (e.g. Flagellum, Tad) is not an enrichment target and is skipped.
ENRICH_WINDOW_TYPES = frozenset({"T1SS", "T2SS", "T3SS", "T4SS", "T5bSS", "T6SS"})
ENRICH_AUTOTRANSPORTER_TYPES = frozenset({"T5aSS", "T5cSS"})
# All T5SS subtypes. Their substrates are Sec-dependent (signal-peptide-bearing),
# so the combined enrichment track pairs DLP with SignalP for these (not DSE, which
# is unreliable for T5); every other type pairs DLP with DSE. openspec:
# signalp-enrichment-track.
ENRICH_T5SS_TYPES = frozenset({"T5aSS", "T5bSS", "T5cSS"})
# DSE cannot call T3SS (DeepSecE T3SS unreliable, CLAUDE.md bug #4; mirrors DSE_NEGATIVE
# in enrichment_testing.py), so T3SS gets no DSE test.
ENRICH_DSE_NO_WINDOW = frozenset({"T3SS"})

# Canonical display-type order (after subtype collapse). The single source of truth
# for "which SS types exist" used by the enrichment test and the figures; must equal
# the union of the window + autotransporter sets above.
SS_TYPE_DISPLAY_ORDER = ("T1SS", "T2SS", "T3SS", "T4SS", "T5aSS", "T5bSS", "T5cSS", "T6SS")
assert set(SS_TYPE_DISPLAY_ORDER) == (ENRICH_WINDOW_TYPES | ENRICH_AUTOTRANSPORTER_TYPES)


def display_type(ss_type: str) -> str:
    """Per-type label: keep T5 subtypes distinct (they behave differently), but
    collapse T6SSi/T6SSii/T6SSiii -> T6SS and pT4SSt -> T4SS. T1SS stays T1SS."""
    if ss_type.startswith("T5"):
        return ss_type
    m = re.match(r"p?(T\d+)[a-z]*SS", ss_type)
    return f"{m.group(1)}SS" if m else ss_type


# Cross-genome pooling (multi-genome runs) can't enumerate the exact product of
# every genome's rotations, so it Monte-Carlos: draw one random rotation per
# genome per replicate and sum. Single-genome runs stay exact (no sampling).
ENRICH_POOL_REPS = 10_000
ENRICH_POOL_SEED = 42


def enrich_null_key(ss_type: str, tool: str) -> str:
    """Key for a (SS type, predictor) rotation-null array in *_enrichment_nulls.npz.

    Shared by the writer (enrichment_testing) and the readers (cross-genome pooling +
    the figure) so the format can't drift across files. Must stay a valid identifier
    because np.savez takes it as a kwarg name."""
    return f"{ss_type}__{tool}"


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

# --- pLM-BLAST ---
# Per-batch timeout for the embedding + search subprocesses. pLM-BLAST
# scales roughly linearly with query count at fixed cpc — empirically on
# CX3 RTX6000 ~17 s per query at cpc=90, so a 184-substrate pooled run
# (4-genome batched) takes ~50 min and a single PAO1 with ~92 substrates
# took 85 min on 2026-06-04. A T1SS-rich or larger-N genome could push
# over the 4 h generic TOOL_TIMEOUT_S, so this step gets its own 24 h cap.
# If you see a real timeout: increase substrate filtering (--conf-threshold)
# or split the run into smaller batches; don't keep bumping the cap.
PLMBLAST_TIMEOUT_S = 86400  # 24h

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

# --- DTU API HTTP timeouts (run_deeplocpro.py + run_signalp.py) ---
# Three tiers matching what the DTU webserver actually does at each call.
# Submit is a multipart upload that the server queues — generous budget
# accounts for the upload itself plus the server's `?ajax=1` redirect dance.
# Status polls are a single JSON GET against the queue, so much shorter.
# Result download is a tar fetch, in between.
DTU_API_SUBMIT_TIMEOUT_S = 60
DTU_API_STATUS_TIMEOUT_S = 15
DTU_API_DOWNLOAD_TIMEOUT_S = 30

# --- map_gbff_to_bakta_cds.py ---
# Minimum coordinate-overlap fraction used when matching a Bakta CDS back
# to its GenBank counterpart. 0.8 catches near-identical predictions while
# rejecting spurious matches caused by gene-prediction tool disagreement.
MAP_GBFF_BAKTA_MIN_OVERLAP = 0.8

# --- T5aSS domain classification ---
MIN_PASSENGER_LENGTH = 100  # aa — below this = "minimal passenger"
LINKER_LENGTH = 30  # aa — alpha-helix linker between passenger and barrel

# Minimum passenger length for which we route annotation tools to the
# passenger-only sequence (EggNOG / BLASTp / pLM-BLAST / HHsearch /
# ProtParam). Below this, the full protein is used so the cpc-90 pre-screen
# in pLM-BLAST and the homology searches have enough residues to work with.
# IPS always sees the full protein regardless.
MIN_PASSENGER_FOR_ANNOTATION = 25  # aa

# Per-row quality flags emitted by t5ss_handler.py, ordered worst → best.
# Lower numbers in T5_QUALITY_FLAG_RANK sort toward the top of the master CSV;
# missing or empty flag = clean call, ranked 0. `no_signalp` and
# `no_sec_signal` share rank 1 — both mean "no valid Sec-pathway export
# signal" and either is sufficient to push the row to the bottom.
T5_QUALITY_FLAG_RANK = {
    "": 0,
    "no_signalp": 1,
    "no_sec_signal": 1,
    "unclassified": 2,
    "barrel_only": 3,
    "omp_porin_no_at": 4,
}

# T5SS components are exported via the Sec pathway across the inner
# membrane before reaching the OM secretion machinery. Sec/SPI (classic
# signal peptide) and Sec/SPII (lipoprotein) are biologically valid for
# T5SS substrates and machinery; Tat-type, PILIN, or OTHER predictions
# are not. Used by t5ss_handler.py to flag rows whose SignalP prediction
# is a non-Sec signal as `no_sec_signal`.
#
# SignalP 6.0 emits SHORT labels in its "Prediction" column (see
# installation_instructions.md upstream): "SP" = Sec/SPI, "LIPO" =
# Sec/SPII, "TAT" = Tat/SPI, "TATLIPO" = Tat/SPII, "PILIN" = Sec/SPIII,
# "OTHER" = no signal peptide. We match the short forms here because
# that's what run_signalp.py passes through into signalp_prediction.
VALID_SEC_SIGNAL_TYPES = ("SP", "LIPO")


def is_sec_signal_peptide(prediction: str) -> bool:
    """True iff a SignalP call is a Sec signal peptide (Sec/SPI or Sec/SPII).

    Single source of truth for "this signal peptide supports Sec secretion", shared
    by the T5SS Sec-signal gate (t5ss_handler) and the enrichment SignalP track.
    Tat/TATLIPO/PILIN/OTHER are not Sec signals. SignalP 6 emits the short labels in
    VALID_SEC_SIGNAL_TYPES; strip first so stray whitespace never flips the call.
    """
    return str(prediction).strip() in VALID_SEC_SIGNAL_TYPES


# --- System filtering ---
# Default-excluded TXSScan models are surface/uptake appendages, NOT protein
# secretion systems, so the proteins near them must not be substrate-called:
#   Flagellum, Tad             - flagellum and tight-adherence pilus
#   T4aP, T4bP, Archaeal-T4P   - type IV pili (a/b and archaeal)
#   MSH                        - MSHA pilus
#   ComM                       - DNA-uptake competence machinery
# Genuine secretion systems stay: T1-T6SS, T9SS, pT4SSi, pT4SSt. Labels are the
# TXSScan model_fqn names (the strings that reach `ss_type`). T3SS is detected by
# default; DeepSecE is excluded from T3SS calls (see cross_validate_predictions),
# so T3SS is NOT in this list.
DEFAULT_EXCLUDED_SYSTEMS = ["Flagellum", "Tad", "T4aP", "T4bP", "MSH", "ComM", "Archaeal-T4P"]

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

# --- PLM-Effector type to MacSyFinder SS type mapping ---
# PLM-Effector reports per-effector-type predictions (T1SE = T1-secreted-
# effector, etc.); the gate at proximity time needs to map those back to
# MacSyFinder SS-system names to verify the system actually exists in the
# genome (same cross-genome leakage class of bug DSE_TO_MACSYFINDER fixes).
PLME_TO_MACSYFINDER = {
    "T1SE": ["T1SS"],
    "T2SE": ["T2SS"],
    "T3SE": ["T3SS"],
    "T4SE": ["pT4SSt", "T4SS"],
    "T6SE": ["T6SSi", "T6SS"],
}


# --- Install-tier defaults ---
# Which optional tools default on at each install tier. The user's tier
# comes from --tier on the CLI, or ~/.ssign/tier (written by
# fetch_databases.sh at the end of a successful fetch), or falls back to
# "extended". Per-tool overrides (--skip-X / --no-skip-X) still win when
# the user is explicit.
#
# Mapping reflects "what does this tier ship that's also actually usable?":
#   - base: Bakta light DB + DeepSecE/PLM-Effector weights → those tools
#     default on. Annotation tools that need extended DBs (EggNOG ~50 GB,
#     IPS ~24 GB, ECOD30 ~11 GB) default off.
#   - extended: adds those three DBs → EggNOG/IPS/pLM-BLAST default on.
#     Pfam + PDB70 also get fetched, but HH-suite stays default-off (see
#     _EXTENDED_ADDS comment) until the wrapper degrades gracefully
#     without UniRef30.
#   - full: adds BLAST NR (390 GB) + HH-suite UniRef30 → BLASTp and
#     HH-suite default on at the full tier.
#
# `run_bakta` is governed by --use-input-annotations + Bakta-DB-present
# invariant separately; not in this table.
# Each tool listed exactly once at the tier it first becomes available.
# The full per-tier on/off map is built below.
_BASE_ENABLED = frozenset({"deeplocpro", "signalp", "deepsece", "protparam"})
# PLM-Effector is installable at every tier but OFF by default: at genome scale it
# over-predicts heavily (called ~25% of the PAO1 proteome; the model is validated only
# on balanced sets, see openspec change enrichment-background-and-plme-default-off and
# the PLM-E over-prediction memo). Opt in with --no-skip-plm-effector. Listed here (not
# in any tier's enabled set) so the tier resolver emits a definite skip=True default
# rather than leaving skip_plm_effector unresolved (which would default it back ON).
_OPT_IN = frozenset({"plm_effector"})
# HH-suite is NOT in extended because its hhblits MSA step needs UniRef30
# (~25 GB), which only ships with --tier full. The wrapper currently
# aborts when UniRef30 is missing rather than degrading to hhsearch-only;
# until that lands, default-on at extended would mislead users into
# thinking HH-suite is available when it can't actually run. Pfam +
# PDB70 still get fetched at extended for when the user wants to enable
# HH-suite manually (--no-skip-hhsuite) on a node with UniRef30 present.
_EXTENDED_ADDS = frozenset({"interproscan", "eggnog", "plmblast"})
_FULL_ADDS = frozenset({"blastp", "hhsuite"})

_TIER_ENABLED = {
    "base": _BASE_ENABLED,
    "extended": _BASE_ENABLED | _EXTENDED_ADDS,
    "full": _BASE_ENABLED | _EXTENDED_ADDS | _FULL_ADDS,
}
_ALL_TOOLS = _TIER_ENABLED["full"] | _OPT_IN

TIER_TOOL_DEFAULTS = {tier: {tool: (tool in enabled) for tool in _ALL_TOOLS} for tier, enabled in _TIER_ENABLED.items()}
DEFAULT_TIER = "extended"
