#!/usr/bin/env bash
# scripts/fetch_databases.sh — download reference databases for ssign.
#
# Usage:
#   bash scripts/fetch_databases.sh --tier {base,extended,full} [--target DIR] [--dry-run]
#
# Tier sizes (post-extraction):
#   base       ~3 GB    NCBI taxdump + Bakta light
#   extended   ~150 GB  + EggNOG + HH-suite (Pfam + PDB70) + InterProScan + ECOD70
#   full       ~630 GB  + BLAST NR + Bakta full + HH-suite UniRef30
#
# Default target: ~/.ssign/databases (override with --target /path).
# Resumes interrupted downloads (wget -c). Skips items already extracted.
#
# Required:
#   wget, tar
#
# Per-tier tool dependencies (the script will tell you which to install):
#   base | full       bakta_db                    pip install ssign[bakta]
#   extended | full   download_eggnog_data.py     conda install -c bioconda eggnog-mapper
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

# pLM-BLAST — ECOD70 prebuilt embedding database.
ECOD70_URL="http://ftp.tuebingen.mpg.de/ebio/protevo/toolkit/databases/plmblast_dbs/ecod70db_20240417.tar.gz"

# InterProScan — pin version explicitly. Bump together with the
# `interproscan-*-bin.tar.gz` checksum file from the EBI release page.
IPS_VERSION="5.71-103.0"
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
  base       ~3 GB    NCBI taxdump + Bakta light
  extended   ~150 GB  + EggNOG + HH-suite (Pfam + PDB70) + InterProScan + ECOD70
  full       ~630 GB  + BLAST NR + Bakta full + HH-suite UniRef30

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
        wget -c -O "$out" "$mirror"
    else
        echo "Error: download failed and no mirror available: $primary" >&2
        return 1
    fi
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

fetch_taxdump() {
    _log "==> NCBI taxdump (~1.5 GB)"
    local dir="$TARGET/taxdump"
    local tarball="$dir/taxdump.tar.gz"

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

    _require_command bakta_db "pip install ssign[bakta]"

    local dir="$TARGET/bakta"
    # bakta_db creates either db/ or db-light/ inside --output depending on
    # version + --type. Glob both rather than guessing the exact convention.
    if compgen -G "$dir/db*/version.json" >/dev/null; then
        _log "Skipping (Bakta DB already at $dir)"
        return 0
    fi

    _run mkdir -p "$dir"
    _run bakta_db download --output "$dir" --type "$variant"
    _log "OK — Bakta $variant DB ready"
}

fetch_eggnog() {
    _log "==> EggNOG database (~50 GB; eggnog-mapper 2.1.13 + EggNOG v6.0)"
    _require_command download_eggnog_data.py "conda install -c bioconda eggnog-mapper"

    local dir="$TARGET/eggnog"
    if [[ -f "$dir/eggnog.db" ]]; then
        _log "Skipping (eggnog.db already at $dir)"
        return 0
    fi

    _run mkdir -p "$dir"
    _run download_eggnog_data.py -y --data_dir "$dir"
    _log "OK — EggNOG DB ready (set EGGNOG_DATA_DIR=$dir)"
}

fetch_interproscan() {
    _log "==> InterProScan ${IPS_VERSION} (~24 GB)"
    local dir="$TARGET/interproscan"
    local tarball="$dir/interproscan-${IPS_VERSION}-64-bit.tar.gz"
    local extracted="$dir/interproscan-${IPS_VERSION}"

    _already_done "$extracted" && return 0

    _run mkdir -p "$dir"
    _wget_with_fallback "$tarball" "$IPS_URL"
    # IPS tarball contains its own wrapper dir; extract to parent.
    _run tar -xzf "$tarball" -C "$dir"
    _run touch "$extracted/$_FETCH_DONE"
    _run rm -f "$tarball"
    _log "OK — InterProScan at $extracted"
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

fetch_ecod70() {
    _log "==> pLM-BLAST ECOD70 (~24 GB extracted)"
    local dir="$TARGET/plm_blast"
    local tarball="$dir/ecod70db_20240417.tar.gz"
    local extracted="$dir/ECOD70"

    _already_done "$extracted" && return 0

    _run mkdir -p "$dir"
    _wget_with_fallback "$tarball" "$ECOD70_URL"
    # ECOD70 tarball contains its own ECOD70/ wrapper dir; extract to parent.
    _run tar -xzf "$tarball" -C "$dir"
    _run touch "$extracted/$_FETCH_DONE"
    _run rm -f "$tarball"
    _log "OK — set SSIGN_ECOD70_DB=$extracted"
}

fetch_blast_nr() {
    _log "==> BLAST NR (~390 GB; via update_blastdb.pl)"
    _require_command update_blastdb.pl "sudo apt install ncbi-blast+ (or brew/conda)"

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
}

run_extended() {
    run_base
    fetch_eggnog
    fetch_hhsuite_pfam
    fetch_hhsuite_pdb70
    fetch_interproscan
    fetch_ecod70
}

run_full() {
    # Full = extended set, plus the big-three: NR, Bakta full, UniRef30.
    # We start with base + the extended additions, then the full additions.
    # Note: full uses Bakta full instead of Bakta light, so we don't call
    # run_extended directly.
    fetch_taxdump
    fetch_bakta full
    fetch_eggnog
    fetch_hhsuite_pfam
    fetch_hhsuite_pdb70
    fetch_hhsuite_uniref30
    fetch_interproscan
    fetch_ecod70
    fetch_blast_nr
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_require_command wget "OS package: wget"
_require_command tar  "OS package: tar"

_log "Tier:   $TIER"
_log "Target: $TARGET"
[[ "$DRY_RUN" -eq 1 ]] && _log "(DRY-RUN; no downloads will occur)"

_run mkdir -p "$TARGET"

case "$TIER" in
    base)     run_base ;;
    extended) run_extended ;;
    full)     run_full ;;
esac

_log "Done. Set the SSIGN_* env vars listed above to point ssign at this directory."
