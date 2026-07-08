#!/usr/bin/env bash
# scripts/fetch_databases.sh — download reference databases for ssign.
#
# Usage:
#   bash scripts/fetch_databases.sh --tier {base,extended,full} [--target DIR] [--dry-run]
#
# Tier sizes (post-extraction):
#   base       ~22 GB   NCBI taxdump + Bakta light + PLM-Effector weights
#   extended   ~100 GB  + EggNOG + InterProScan + ECOD30 (pLM-BLAST)
#   full       ~500 GB  + Bakta full + HH-suite (Pfam + PDB70 + UniRef30) + BLASTp-vs-Swiss-Prot
#              (HH-suite extracted dominates; NR is opt-in, not fetched by default — see fetch_blast_nr)
# (sizes above are rough — see task #221 for the audit.)
#
# Default target: ~/.ssign/databases (override with --target /path).
# Resumes interrupted downloads (wget -c). Skips items already extracted.
#
# Required:
#   wget, tar, unzip
#
# Per-tier tool dependencies (the script will tell you which to install):
#   all tiers         hf                          pip install 'huggingface_hub[cli]'
#   base | full       bakta_db, amrfinder         pip install ssign[bakta] + mamba ncbi-amrfinderplus
#   full              update_blastdb.pl           OS package: ncbi-blast+

set -euo pipefail

# ---------------------------------------------------------------------------
# Pinned versions / sources
# ---------------------------------------------------------------------------
# These are the upstream URLs as of 2026-04-30. v1.0.0 freezes against these
# pins; newer versions land in v1.1.0+. Mirror-fallback chains documented in
# the project plan (Addendum D.5 for HH-suite).

NCBI_TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# HH-suite — Tübingen MPI is the canonical fresher source per Söding lab
# issue #382. GWDG is the older mirror; used only when Tübingen is down.
HHSUITE_TUEBINGEN_BASE="http://ftp.tuebingen.mpg.de/pub/ebio/protevo/toolkit/databases/hhsuite_dbs"
HHSUITE_GWDG_BASE="https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs"

HHSUITE_PFAM_TUEBINGEN="${HHSUITE_TUEBINGEN_BASE}/PfamA_v38_2.tar.gz"
HHSUITE_PFAM_GWDG="${HHSUITE_GWDG_BASE}/pfamA_35.0.tar.gz"

HHSUITE_PDB70_TUEBINGEN="${HHSUITE_TUEBINGEN_BASE}/pdb70_from_mmcif_2026-02-20.tar.gz"
HHSUITE_PDB70_GWDG="${HHSUITE_GWDG_BASE}/pdb70_from_mmcif_latest.tar.gz"

# UniRef30 only exists at GWDG; no Tübingen mirror as of plan D.5.
HHSUITE_UNIREF30_GWDG="https://wwwuser.gwdg.de/~compbiol/uniclust/2023_02/UniRef30_2023_02_hhsuite.tar.gz"

# pLM-BLAST — ECOD30 prebuilt embedding database. Switched from ECOD70
# (21 GB / 2024-04-17) on 2026-06-04: ECOD30 is half the size (10 GB),
# represents every ECOD F-group, and the paper's published benchmarks
# already use this cluster level. The other levels (ECOD50/70/90) are
# also hosted at the same FTP path if a user wants more redundancy.
ECOD_URL="http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs/ecod30db_20240417.tar.gz"

# EggNOG — current host. The legacy hostname `eggnogdb.embl.de` was retired;
# eggnog-mapper 2.1.13 (latest on bioconda as of 2026-05) still hardcodes the
# dead host and produces 0-byte files with exit-code 0. 2.1.14 fixed the URL
# on GitHub but never reached PyPI. We wget directly from the live host
# instead of relying on download_eggnog_data.py.
EGGNOG_BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

# PLM-Effector — trained model weights from upstream (slow Chinese academic
# mirror at ~2 MB/s; ~11 min wall for the 1.5 GB sourcecode.zip). The PLMs
# themselves come from HuggingFace (~17 GB across four repos).
PLM_EFFECTOR_SOURCECODE_URL="https://www.mgc.ac.cn/PLM-Effector/download/sourcecode.zip"

# InterProScan — pin version explicitly. Bump together with the
# `interproscan-*-bin.tar.gz` checksum file from the EBI release page.
IPS_VERSION="5.77-108.0"
IPS_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${IPS_VERSION}/interproscan-${IPS_VERSION}-64-bit.tar.gz"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

TIER=""
TARGET="${HOME}/.ssign/databases"
DRY_RUN=0

usage() {
    cat <<'EOF'
fetch_databases.sh — download reference databases for ssign.

Usage:
  bash scripts/fetch_databases.sh --tier {base,extended,full} [--target DIR] [--dry-run]

Tier sizes (post-extraction):
  base       ~22 GB   NCBI taxdump + Bakta light + PLM-Effector weights
  extended   ~100 GB  + EggNOG + InterProScan + ECOD30 (pLM-BLAST)
  full       ~500 GB  + Bakta full + HH-suite (Pfam + PDB70 + UniRef30) + BLASTp-vs-Swiss-Prot

Default target: ~/.ssign/databases (override with --target /path).
Resumes interrupted downloads (wget -c). Skips items already extracted.
EOF
    exit "${1:-1}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tier)        TIER="$2"; shift 2 ;;
        --target)      TARGET="$2"; shift 2 ;;
        --dry-run)     DRY_RUN=1; shift ;;
        -h|--help)     usage 0 ;;
        *)             echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

case "$TIER" in
    base|extended|full) ;;
    "") echo "Error: --tier is required" >&2; usage 1 ;;
    *)  echo "Error: --tier must be one of base, extended, full (got: $TIER)" >&2; exit 1 ;;
esac

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_log() {
    printf '[fetch_databases] %s\n' "$*"
}

_run() {
    # Echo the command, then either run it or skip if --dry-run.
    _log "+ $*"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "$@"
    fi
}

# Install hints for the external tools a fetch needs. Defined once so the
# up-front pre-flight (_preflight_tier) and the per-fetch _require_command
# calls quote the same fix-it text and can't drift.
_HINT_BAKTA_DB="pip install ssign[bakta]"
_HINT_AMRFINDER="mamba install -c bioconda ncbi-amrfinderplus  (or 'mamba create -n bakta-deps -c bioconda ncbi-amrfinderplus -y' then 'export PATH=~/.conda/envs/bakta-deps/bin:\$PATH')"
_HINT_HF="pip install 'huggingface_hub[cli]'"
_HINT_UNZIP="OS package: unzip"
_HINT_UPDATE_BLASTDB="sudo apt install ncbi-blast+ (or brew/conda)"

_require_command() {
    # Require a command on PATH; print a fix-it message if missing.
    # In dry-run mode this is informational rather than fatal, so the user
    # can preview what a full run would do without first installing every
    # tier's tooling.
    local cmd="$1"; shift
    local install_hint="$*"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        if [[ "$DRY_RUN" -eq 1 ]]; then
            _log "(would require '$cmd' — install with: $install_hint)"
            return 0
        fi
        echo "Error: '$cmd' not found on PATH." >&2
        echo "    Install with: $install_hint" >&2
        exit 1
    fi
}

_wget_with_fallback() {
    # Try the primary URL, fall back to the mirror on failure.
    # $1 = output path, $2 = primary URL, $3 = mirror URL (optional).
    # Note: wget -c is itself idempotent — it returns immediately if the
    # remote file is already fully fetched, and resumes from byte offset
    # if partial. No early "skip if present" check is needed and adding
    # one would treat a killed-mid-download partial file as complete.
    local out="$1" primary="$2" mirror="${3:-}"

    # Some callers (HH-suite) download to a path whose parent dir isn't
    # created by the per-fetcher mkdir. Cover that here so callers don't have to.
    _run mkdir -p "$(dirname "$out")"

    _log "Fetching $primary -> $out"
    if [[ "$DRY_RUN" -eq 1 ]]; then
        _log "(dry-run; not downloading)"
        return 0
    fi

    if wget -c -O "$out" "$primary"; then
        return 0
    fi

    if [[ -n "$mirror" ]]; then
        _log "Primary failed; trying mirror $mirror"
        # wget -c on a partial file with a different URL is risky if the file
        # contents differ; remove the partial first so resume can't corrupt.
        rm -f "$out"
        if wget -c -O "$out" "$mirror"; then
            return 0
        fi
    fi

    # Distinguish "version rotated off server" (404) from a transient failure.
    # IPS, NCBI, and EMBL all silently remove old version dirs when newer ones
    # ship — surface that as a specific pointer rather than a generic wget error.
    if wget --server-response --spider "$primary" 2>&1 | grep -q "404 Not Found"; then
        echo "Error: $primary returned 404." >&2
        echo "  The pinned version may have been rotated off the server." >&2
        echo "  Check the index page and update the matching *_VERSION constant" >&2
        echo "  at the top of scripts/fetch_databases.sh." >&2
    else
        echo "Error: download failed (not a 404): $primary" >&2
        echo "  URL is reachable; the failure was likely transient — retry the fetch." >&2
    fi
    return 1
}

# Sentinel filename written into a destination directory after a successful
# extraction. Re-runs check for this rather than the directory itself, so a
# half-extracted dir (e.g. interrupted by Ctrl-C, OOM, disk-full) is correctly
# detected as incomplete and re-extracted.
_FETCH_DONE=".ssign_fetch_complete"

_extract_and_finalize() {
    # $1 = tarball, $2 = destination dir.
    # Extracts $tarball into $dest, writes the completion sentinel, and
    # removes the tarball to reclaim disk space (full tier saves ~75 GB).
    local tarball="$1" dest="$2"
    _run mkdir -p "$dest"
    _run tar -xzf "$tarball" -C "$dest"
    _run touch "$dest/$_FETCH_DONE"
    _run rm -f "$tarball"
}

# ---------------------------------------------------------------------------
# Per-database fetchers
# ---------------------------------------------------------------------------

_already_done() {
    # $1 = directory expected to contain $_FETCH_DONE after a previous successful run.
    # Returns 0 (skip) if the marker is present, 1 (proceed) otherwise.
    local dir="$1"
    if [[ -f "$dir/$_FETCH_DONE" ]]; then
        _log "Skipping ($dir already complete)"
        return 0
    fi
    return 1
}

_external_db_exists() {
    # Honor a user-set env var pointing at an externally-managed database.
    # On shared HPC the user often has DBs at /rds/.../databases/ and only
    # wants to fetch the subset they're missing; without this check the
    # script tries to re-bootstrap everything under $TARGET and pulls in
    # `_require_command` checks (e.g. bakta_db CLI) the user shouldn't
    # need to satisfy when their DB is already populated elsewhere.
    #
    # Returns 0 (skip the fetch) iff the env var is set AND its dir
    # contains the verification glob. Otherwise returns 1; a mismatched
    # env var only warns — it doesn't abort — so the user can still get
    # the normal fetch behaviour with a noisy hint.
    #
    # $1 = env var name (e.g. "BAKTA_DB")
    # $2 = verification glob relative to the dir (e.g. "version.json")
    local var_name="$1"
    local check_glob="$2"
    local path="${!var_name:-}"
    [[ -z "$path" ]] && return 1
    if [[ ! -d "$path" ]]; then
        _log "Note: \$$var_name=$path is set but not a directory; ignoring and fetching to \$TARGET."
        return 1
    fi
    if ! compgen -G "$path/$check_glob" >/dev/null; then
        _log "Note: \$$var_name=$path is set but missing $check_glob; ignoring and fetching to \$TARGET."
        return 1
    fi
    _log "Skipping (\$$var_name=$path already populated)"
    return 0
}

_preflight_tier() {
    # Fail fast on missing external tools BEFORE any (possibly hours-long)
    # download starts, instead of only when the fetch reaches that DB's step.
    # Motivating case: update_blastdb.pl (BLAST NR) is checked at the LAST
    # step, so without this a full-tier fetch would pull Bakta + all HH-suite
    # (~135 GB, hours) and only THEN die on a missing NR tool. This mirrors the
    # per-fetch _require_command checks exactly (same tools, same guards), so
    # it never fails on a tool the actual fetch wouldn't have needed. Honors
    # DRY_RUN via _require_command (informational, not fatal).
    local tier="$1"
    # Bakta (every tier fetches it) — skipped when an external Bakta DB is set.
    if ! _external_db_exists BAKTA_DB "db*/version.json"; then
        _require_command bakta_db "$_HINT_BAKTA_DB"
        _require_command amrfinder "$_HINT_AMRFINDER"
    fi
    # PLM-Effector weights (every tier).
    _require_command hf "$_HINT_HF"
    _require_command unzip "$_HINT_UNZIP"
    # BLAST NR (full tier only).
    if [[ "$tier" == "full" ]]; then
        _require_command update_blastdb.pl "$_HINT_UPDATE_BLASTDB"
    fi
}

fetch_taxdump() {
    _log "==> NCBI taxdump (~1.5 GB)"
    local dir="$TARGET/taxdump"
    local tarball="$dir/taxdump.tar.gz"

    _external_db_exists SSIGN_TAXDUMP_DIR "nodes.dmp" && return 0
    _already_done "$dir" && return 0

    _run mkdir -p "$dir"
    _wget_with_fallback "$tarball" "$NCBI_TAXDUMP_URL"
    # Selective extract — the full taxdump tarball is ~400 MB compressed but
    # ssign only needs nodes.dmp + names.dmp via taxopy.
    _run tar -xzf "$tarball" -C "$dir" nodes.dmp names.dmp
    _run touch "$dir/$_FETCH_DONE"
    _run rm -f "$tarball"
    _log "OK — set SSIGN_TAXDUMP_DIR=$dir"
}

fetch_bakta() {
    # $1 = "light" or "full"
    local variant="$1"
    local size_hint
    if [[ "$variant" == "light" ]]; then size_hint="~2 GB"; else size_hint="~84 GB"; fi
    _log "==> Bakta DB ($variant, $size_hint)"

    # Check externally-managed DB before requiring bakta_db on PATH —
    # users who already have a Bakta DB at $BAKTA_DB shouldn't need the
    # bakta_db CLI installed just to run this fetcher for OTHER databases.
    _external_db_exists BAKTA_DB "db*/version.json" && return 0

    _require_command bakta_db "$_HINT_BAKTA_DB"
    # bakta_db's startup runs the same dependency check as `bakta` itself,
    # which fails fast if AMRFinderPlus is not on PATH — even though the
    # tool isn't actually needed for downloading the DB. `pip install
    # ssign[bakta]` only pulls the Python wrapper; AMRFinderPlus is a
    # separate NCBI binary. Catch it here with a useful pointer instead
    # of letting bakta_db crash mid-run.
    _require_command amrfinder "$_HINT_AMRFINDER"

    local dir="$TARGET/bakta"
    # bakta_db creates either db/ or db-light/ inside --output depending on
    # version + --type. Glob both rather than guessing the exact convention.
    if compgen -G "$dir/db*/version.json" >/dev/null; then
        _log "Skipping (Bakta DB already at $dir)"
        return 0
    fi

    _run mkdir -p "$dir"
    _run bakta_db download --output "$dir" --type "$variant"
    # bakta_db creates `db/` (full) or `db-light/` (light) inside --output.
    # Glob both rather than guessing; fall back to a sensible default if the
    # glob comes up empty so the log line doesn't print `BAKTA_DB=`.
    local bakta_subdir
    bakta_subdir=$(compgen -G "$dir/db*" | head -n 1)
    : "${bakta_subdir:=$dir/db-${variant}}"
    _log "OK — Bakta $variant DB ready (set BAKTA_DB=$bakta_subdir)"
}

fetch_eggnog() {
    # Three files match what download_eggnog_data.py would fetch for the
    # default (non-HMMER, non-novel-families, non-MMseqs) install path —
    # the same defaults ssign's eggnog wrapper relies on at runtime.
    _log "==> EggNOG database (~25 GB extracted; emapperdb v5.0.2)"
    local dir="$TARGET/eggnog"

    _external_db_exists EGGNOG_DATA_DIR "eggnog.db" && return 0

    if [[ -f "$dir/eggnog.db" && -f "$dir/eggnog.taxa.db" && -f "$dir/eggnog_proteins.dmnd" ]]; then
        _log "Skipping (eggnog.{db,taxa.db,_proteins.dmnd} already at $dir)"
        return 0
    fi

    _run mkdir -p "$dir"

    if [[ ! -f "$dir/eggnog.db" ]]; then
        _wget_with_fallback "$dir/eggnog.db.gz" "$EGGNOG_BASE_URL/eggnog.db.gz"
        _run gunzip -f "$dir/eggnog.db.gz"
    fi

    if [[ ! -f "$dir/eggnog.taxa.db" ]]; then
        _wget_with_fallback "$dir/eggnog.taxa.tar.gz" "$EGGNOG_BASE_URL/eggnog.taxa.tar.gz"
        _run tar -zxf "$dir/eggnog.taxa.tar.gz" -C "$dir"
        _run rm -f "$dir/eggnog.taxa.tar.gz"
    fi

    if [[ ! -f "$dir/eggnog_proteins.dmnd" ]]; then
        _wget_with_fallback "$dir/eggnog_proteins.dmnd.gz" "$EGGNOG_BASE_URL/eggnog_proteins.dmnd.gz"
        _run gunzip -f "$dir/eggnog_proteins.dmnd.gz"
    fi

    _log "OK — EggNOG DB ready (set EGGNOG_DATA_DIR=$dir)"
}

fetch_interproscan() {
    _log "==> InterProScan ${IPS_VERSION} (~24 GB)"
    local dir="$TARGET/interproscan"
    local tarball="$dir/interproscan-${IPS_VERSION}-64-bit.tar.gz"
    local extracted="$dir/interproscan-${IPS_VERSION}"

    _external_db_exists SSIGN_INTERPROSCAN_PATH "interproscan.sh" && return 0
    _already_done "$extracted" && return 0

    _run mkdir -p "$dir"
    _wget_with_fallback "$tarball" "$IPS_URL"
    # IPS tarball contains its own wrapper dir; extract to parent.
    _run tar -xzf "$tarball" -C "$dir"
    _run touch "$extracted/$_FETCH_DONE"
    _run rm -f "$tarball"
    _log "OK — InterProScan at $extracted (set SSIGN_INTERPROSCAN_PATH=$extracted)"
}

fetch_hhsuite_pfam() {
    _log "==> HH-suite Pfam (~3 GB; Tübingen v38, GWDG v35 fallback)"
    local dir="$TARGET/hhsuite/pfam"
    local tarball="$TARGET/hhsuite/PfamA.tar.gz"

    _already_done "$dir" && return 0

    _wget_with_fallback "$tarball" "$HHSUITE_PFAM_TUEBINGEN" "$HHSUITE_PFAM_GWDG"
    _extract_and_finalize "$tarball" "$dir"
    _log "OK — set SSIGN_HHSUITE_PFAM=$dir"
}

fetch_hhsuite_pdb70() {
    _log "==> HH-suite PDB70 (~23 GB; Tübingen 2026-02 build, GWDG fallback)"
    local dir="$TARGET/hhsuite/pdb70"
    local tarball="$TARGET/hhsuite/pdb70.tar.gz"

    _already_done "$dir" && return 0

    _wget_with_fallback "$tarball" "$HHSUITE_PDB70_TUEBINGEN" "$HHSUITE_PDB70_GWDG"
    _extract_and_finalize "$tarball" "$dir"
    _log "OK — set SSIGN_HHSUITE_PDB70=$dir"
}

fetch_hhsuite_uniref30() {
    _log "==> HH-suite UniRef30 (~25 GB; GWDG only)"
    local dir="$TARGET/hhsuite/uniref30"
    local tarball="$TARGET/hhsuite/UniRef30_2023_02_hhsuite.tar.gz"

    _already_done "$dir" && return 0

    _wget_with_fallback "$tarball" "$HHSUITE_UNIREF30_GWDG"
    _extract_and_finalize "$tarball" "$dir"
    _log "OK — set SSIGN_HHSUITE_UNICLUST=$dir"
}

fetch_ecod() {
    _log "==> pLM-BLAST ECOD30 (~11 GB extracted)"
    local dir="$TARGET/plm_blast"
    local tarball="$dir/ecod30db_20240417.tar.gz"
    local extracted="$dir/ECOD30"

    # SSIGN_ECOD_DB is the parent dir holding ECOD30/ + the ECOD30.csv
    # sidecar. The sibling .csv is required by pLM-BLAST's DataObject loader.
    _external_db_exists SSIGN_ECOD_DB "ECOD30.csv" && return 0
    _already_done "$extracted" && return 0

    _run mkdir -p "$dir"
    _wget_with_fallback "$tarball" "$ECOD_URL"
    # ECOD30 tarball contains its own ECOD30/ wrapper dir; extract to parent.
    _run tar -xzf "$tarball" -C "$dir"
    _run touch "$extracted/$_FETCH_DONE"
    _run rm -f "$tarball"
    _log "OK — set SSIGN_ECOD_DB=$extracted"
}

fetch_plm_effector_weights() {
    # PLM-Effector ships its own trained_models bundle plus four pretrained
    # protein language models from HuggingFace. ssign reads them via
    # SSIGN_PLM_EFFECTOR_WEIGHTS. GPU strongly recommended at runtime
    # (~100x speedup over CPU).
    _log "==> PLM-Effector weights (~19 GB total; mgc.ac.cn + HuggingFace)"
    _require_command hf "$_HINT_HF"
    _require_command unzip "$_HINT_UNZIP"

    local dir="$TARGET/plm_effector_weights"
    local hf_dir="$dir/transformers_pretrained"
    _run mkdir -p "$dir" "$hf_dir"

    # Part 1: trained_models (~1.7 GB, slow mirror).
    if [[ ! -d "$dir/trained_models" ]]; then
        local zip="$dir/sourcecode.zip"
        _wget_with_fallback "$zip" "$PLM_EFFECTOR_SOURCECODE_URL"
        _run unzip -q -d "$dir" "$zip"
        _run mv "$dir/sourcecode/trained_models" "$dir/trained_models"
        _run rm -rf "$dir/sourcecode" "$zip"
    else
        _log "Skipping trained_models (already present)"
    fi

    # Part 2: four pretrained PLMs from HuggingFace (~17 GB).
    # prot_t5_xl_uniref50 carries a TF checkpoint we don't need; include-list
    # trims the download to the PyTorch + tokenizer files.
    if [[ ! -d "$hf_dir/prot_t5_xl_uniref50" ]]; then
        _run hf download Rostlab/prot_t5_xl_uniref50 \
            --include "config.json" "spiece.model" "tokenizer_config.json" \
                "special_tokens_map.json" "pytorch_model.bin" \
            --local-dir "$hf_dir/prot_t5_xl_uniref50"
    else
        _log "Skipping prot_t5_xl_uniref50 (already present)"
    fi
    if [[ ! -d "$hf_dir/esm1b_t33_650M_UR50S" ]]; then
        _run hf download facebook/esm1b_t33_650M_UR50S \
            --local-dir "$hf_dir/esm1b_t33_650M_UR50S"
    else
        _log "Skipping esm1b_t33_650M_UR50S (already present)"
    fi
    if [[ ! -d "$hf_dir/esm2_t33_650M_UR50D" ]]; then
        _run hf download facebook/esm2_t33_650M_UR50D \
            --local-dir "$hf_dir/esm2_t33_650M_UR50D"
    else
        _log "Skipping esm2_t33_650M_UR50D (already present)"
    fi
    if [[ ! -d "$hf_dir/prot_bert" ]]; then
        _run hf download Rostlab/prot_bert \
            --local-dir "$hf_dir/prot_bert"
    else
        _log "Skipping prot_bert (already present)"
    fi

    _log "OK — set SSIGN_PLM_EFFECTOR_WEIGHTS=$dir"
}

fetch_blast_swissprot() {
    # Swiss-Prot is the full-tier default BLASTp DB: ~300 MB curated, each
    # entry carrying a reviewed function name — far more useful for "what is
    # this substrate" than NR's redundant hits, and it searches in minutes
    # instead of the day+ full NR takes on a real substrate set. NR remains
    # available for the exhaustive case via `fetch_blast_nr` + `--blastp-db`.
    _log "==> BLAST Swiss-Prot (~300 MB; via update_blastdb.pl)"
    _require_command update_blastdb.pl "$_HINT_UPDATE_BLASTDB"

    local dir="$TARGET/blast_swissprot"
    _run mkdir -p "$dir"
    _run bash -c "cd \"$dir\" && update_blastdb.pl --decompress swissprot"
    _log "OK — BLAST Swiss-Prot at $dir"
}

fetch_blast_nr() {
    # NOT fetched by any tier by default (see fetch_blast_swissprot). Kept for
    # the operator who wants exhaustive NR homology: run this function's body
    # manually, then pass `--blastp-db <dir>/nr`. ~800 GB decompressed.
    _log "==> BLAST NR (~800 GB; via update_blastdb.pl)"
    _require_command update_blastdb.pl "$_HINT_UPDATE_BLASTDB"

    local dir="$TARGET/blast_nr"
    _run mkdir -p "$dir"

    # update_blastdb.pl is idempotent — it checks each volume's MD5 and only
    # downloads new/changed ones. Safe to re-run. Wrapped in a subshell so
    # the cd stays local to this call.
    _run bash -c "cd \"$dir\" && update_blastdb.pl --decompress nr"
    _log "OK — BLAST NR at $dir (set BLASTDB=$dir)"
}

# ---------------------------------------------------------------------------
# Tier dispatch
# ---------------------------------------------------------------------------

run_base() {
    fetch_taxdump
    fetch_bakta light
    fetch_plm_effector_weights
}

run_extended() {
    run_base
    fetch_eggnog
    fetch_interproscan
    fetch_ecod
}

run_full() {
    # Full = extended set + Bakta full + HH-suite (Pfam/PDB70/UniRef30) +
    # BLASTp-vs-Swiss-Prot. NR is deliberately NOT fetched here (it's ~800 GB
    # and too slow to blast a real substrate set — see fetch_blast_swissprot);
    # the operator opts into NR manually. full uses Bakta full instead of
    # light, so we don't call run_extended directly.
    fetch_taxdump
    fetch_bakta full
    fetch_plm_effector_weights
    fetch_eggnog
    fetch_hhsuite_pfam
    fetch_hhsuite_pdb70
    fetch_hhsuite_uniref30
    fetch_interproscan
    fetch_ecod
    fetch_blast_swissprot
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_require_command wget "OS package: wget"
_require_command tar  "OS package: tar"

_log "Tier:   $TIER"
_log "Target: $TARGET"
[[ "$DRY_RUN" -eq 1 ]] && _log "(DRY-RUN; no downloads will occur)"

# Check every external tool this tier needs up front, so a missing one fails
# now rather than after hours of downloading (see _preflight_tier).
_preflight_tier "$TIER"

_run mkdir -p "$TARGET"

case "$TIER" in
    base)     run_base ;;
    extended) run_extended ;;
    full)     run_full ;;
esac

# Record the data root so `ssign doctor` can find these databases without
# the user having to export an env var per DB. Doctor reads this file first
# (overriding the ~/.ssign/databases default) and resolves every sub-path
# beneath it via the same layout fetch_databases.sh just wrote.
if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "${HOME}/.ssign"
    printf '%s\n' "$(cd "$TARGET" && pwd)" > "${HOME}/.ssign/db_root"
    printf '%s\n' "$TIER" > "${HOME}/.ssign/tier"
    _log "Recorded data root at ~/.ssign/db_root → $TARGET"
    _log "Recorded tier      at ~/.ssign/tier     → $TIER"
fi

_log "Done. ssign doctor will read ~/.ssign/db_root automatically; SSIGN_* env vars override per-DB if you need them."
