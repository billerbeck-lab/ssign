#!/usr/bin/env bash
# scripts/fetch_weights.sh — download model weights for ssign.
#
# Usage:
#   bash scripts/fetch_weights.sh [--target DIR] [--dry-run]
#
# Downloads:
#   DeepSecE checkpoint (~2.5 GB)        Zenodo, SJTU fallback
#   PLM-Effector trained_models (~1.7 GB) MGC China (sourcecode.zip)
#   ProtT5 weights (~2.5 GB)             HuggingFace (Rostlab/prot_t5_xl_uniref50)
#   ESM-1b weights (~7 GB)               HuggingFace (facebook/esm1b_t33_650M_UR50S)
#   ESM-2 weights (~3 GB)                HuggingFace (facebook/esm2_t33_650M_UR50D)
#   ProtBert weights (~1.6 GB)           HuggingFace (Rostlab/prot_bert)
#
# Total download: ~18 GB.
# Default target: ~/.ssign/models (override with --target /path).
# Resumes interrupted downloads (wget -c, hf download is itself resumable).
# Skips items already complete.
#
# Required:
#   wget, unzip, hf (or huggingface-cli)
#
# Tool dependency notes:
#   `hf` is the modern command from huggingface_hub >= 0.21; older installs
#   expose `huggingface-cli` instead. The script picks whichever is present.
#   Both come from `pip install huggingface_hub` (transitively pulled in
#   by `pip install ssign[extended]`).

set -euo pipefail

# ---------------------------------------------------------------------------
# Pinned sources
# ---------------------------------------------------------------------------
# These are the upstream URLs as of 2026-04-30. v1.0.0 freezes against these
# pins. Zenodo URLs are placeholders until the deposits land in Phase 8.

# DeepSecE — Zenodo (primary, post-Phase 8) + SJTU origin (current fallback;
# server is known-unreliable, see the longevity-commitment memory).
# DeepSecE checkpoint URLs are duplicated from src/ssign_app/scripts/run_deepsece.py;
# keep both sources in sync when changing.
DEEPSECE_URL_ZENODO="https://zenodo.org/records/PLACEHOLDER/files/deepsece_checkpoint.pt"
DEEPSECE_URL_SJTU="https://tool2-mml.sjtu.edu.cn/DeepSecE/checkpoint.pt"
# Reject DeepSecE files smaller than this as truncated (matches
# MIN_CHECKPOINT_BYTES in run_deepsece.py).
DEEPSECE_MIN_BYTES=$((100 * 1024 * 1024))

# PLM-Effector trained models — bundled inside sourcecode.zip on the MGC
# (Chinese Academy of Medical Sciences) server. Mirroring to Zenodo is part
# of the longevity stack but pending until Phase 8.
PLM_EFFECTOR_SOURCECODE_URL="https://www.mgc.ac.cn/PLM-Effector/download/sourcecode.zip"

# HuggingFace model IDs. ProtT5 is loaded only with the files PLM-Effector
# actually needs (saves ~10 GB vs the full repo) — this is the same include
# list as docs/optional_tools.md.
HF_PROTT5_ID="Rostlab/prot_t5_xl_uniref50"
HF_PROTT5_FILES=(config.json spiece.model tokenizer_config.json special_tokens_map.json pytorch_model.bin)
HF_ESM1B_ID="facebook/esm1b_t33_650M_UR50S"
HF_ESM2_ID="facebook/esm2_t33_650M_UR50D"
HF_PROTBERT_ID="Rostlab/prot_bert"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

TARGET="${HOME}/.ssign/models"
DRY_RUN=0

usage() {
    cat <<'EOF'
fetch_weights.sh — download model weights for ssign.

Usage:
  bash scripts/fetch_weights.sh [--target DIR] [--dry-run]

Downloads (~18 GB total):
  DeepSecE checkpoint        ~2.5 GB    Zenodo + SJTU fallback
  PLM-Effector trained models ~1.7 GB   MGC China
  ProtT5 / ESM-1b / ESM-2 / ProtBert (~14 GB total) HuggingFace

Default target: ~/.ssign/models (override with --target /path).
Resumes interrupted downloads. Skips items already complete.
EOF
    exit "${1:-1}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --target)      TARGET="$2"; shift 2 ;;
        --dry-run)     DRY_RUN=1; shift ;;
        -h|--help)     usage 0 ;;
        *)             echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Helpers (kept in sync with fetch_databases.sh — see deferral E.8 in plan
# about consolidating these into a shared helper)
# ---------------------------------------------------------------------------

_log() {
    printf '[fetch_weights] %s\n' "$*"
}

_run() {
    _log "+ $*"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "$@"
    fi
}

_require_command() {
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
    # See fetch_databases.sh for the full comment. Short version: rely on
    # wget -c for resume/skip semantics rather than checking file presence
    # first (which would treat a partial file as complete).
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
        rm -f "$out"
        wget -c -O "$out" "$mirror"
    else
        echo "Error: download failed and no mirror available: $primary" >&2
        return 1
    fi
}

_FETCH_DONE=".ssign_fetch_complete"

_already_done() {
    local dir="$1"
    if [[ -f "$dir/$_FETCH_DONE" ]]; then
        _log "Skipping ($dir already complete)"
        return 0
    fi
    return 1
}

# Pick whichever HuggingFace CLI binary is on PATH.
_hf_cli() {
    if command -v hf >/dev/null 2>&1; then
        echo "hf"
    elif command -v huggingface-cli >/dev/null 2>&1; then
        echo "huggingface-cli"
    else
        # Caller will hit _require_command before reaching here; this branch
        # only runs in --dry-run on a system without HF tools installed.
        echo "hf"
    fi
}

# ---------------------------------------------------------------------------
# Per-weight fetchers
# ---------------------------------------------------------------------------

fetch_deepsece() {
    _log "==> DeepSecE checkpoint (~2.5 GB)"
    # Land the checkpoint at the location run_deepsece.py looks for it by
    # default (DEFAULT_CHECKPOINT in run_deepsece.py:51). With the default
    # --target=~/.ssign/models, this matches; with a custom --target, the
    # user must symlink or pass --checkpoint to the DeepSecE step.
    local out="$TARGET/deepsece_checkpoint.pt"

    # No sentinel here — the file is a single .pt with a known minimum size,
    # so wget -c idempotency + the size check below is enough. Avoids
    # inventing a per-file marker that would mismatch run_deepsece.py's own
    # _validate_checkpoint logic.
    _run mkdir -p "$TARGET"
    _wget_with_fallback "$out" "$DEEPSECE_URL_ZENODO" "$DEEPSECE_URL_SJTU"

    if [[ "$DRY_RUN" -eq 0 && -f "$out" ]]; then
        local size
        size=$(stat -c%s "$out" 2>/dev/null || stat -f%z "$out")
        if (( size < DEEPSECE_MIN_BYTES )); then
            echo "Error: DeepSecE checkpoint too small (${size} bytes < ${DEEPSECE_MIN_BYTES})." >&2
            echo "    Likely truncated; remove $out and re-run." >&2
            exit 1
        fi
    fi

    _log "OK — DeepSecE checkpoint at $out"
}

fetch_plm_effector_trained_models() {
    _log "==> PLM-Effector trained_models (~1.7 GB)"
    local dir="$TARGET/plm_effector"
    local extracted="$dir/trained_models"
    local zipfile="$dir/sourcecode.zip"
    local stage="$dir/_sourcecode_stage"

    _already_done "$extracted" && return 0

    _require_command unzip "OS package: unzip"

    _run mkdir -p "$dir"
    _wget_with_fallback "$zipfile" "$PLM_EFFECTOR_SOURCECODE_URL"

    # Validate the zip before trusting it — `wget -c` can leave a partial
    # file looking complete if the prior wget process was killed mid-download.
    if [[ "$DRY_RUN" -eq 0 && -f "$zipfile" ]]; then
        if ! unzip -tq "$zipfile" >/dev/null; then
            echo "Error: $zipfile failed integrity check; removing and retrying." >&2
            rm -f "$zipfile"
            _wget_with_fallback "$zipfile" "$PLM_EFFECTOR_SOURCECODE_URL"
            unzip -tq "$zipfile" >/dev/null
        fi
    fi

    # If a previous run died after extracting partway, $extracted exists but
    # has no sentinel — clean it before re-extracting to avoid mixed state.
    _run rm -rf "$extracted"
    # The zip contains sourcecode/trained_models/ plus all the (unneeded)
    # training/inference Python source. Extract to a stage dir, move just
    # trained_models/ into place, drop the rest.
    _run rm -rf "$stage"
    _run mkdir -p "$stage"
    _run unzip -q "$zipfile" -d "$stage"
    _run mv "$stage/sourcecode/trained_models" "$extracted"
    _run rm -rf "$stage"
    _run rm -f "$zipfile"
    _run touch "$extracted/$_FETCH_DONE"
    _log "OK — PLM-Effector trained_models at $extracted"
}

fetch_hf_model() {
    # $1 = HF model ID, $2 = local subdir name (relative to PLM-Effector
    # transformers_pretrained/), $3+ = optional --include patterns.
    local model_id="$1" subdir="$2"; shift 2
    local include_args=("$@")

    local local_dir="$TARGET/plm_effector/transformers_pretrained/$subdir"
    _log "==> HuggingFace: $model_id"

    _already_done "$local_dir" && return 0

    local hf_bin
    hf_bin="$(_hf_cli)"
    _require_command "$hf_bin" "pip install ssign[extended]"

    _run mkdir -p "$local_dir"

    # --local-dir-use-symlinks False writes real files into local_dir without
    # also caching them under ~/.cache/huggingface — saves ~13 GB across the
    # four HF models.
    local cmd=("$hf_bin" download "$model_id" --local-dir "$local_dir" --local-dir-use-symlinks False)
    for f in "${include_args[@]}"; do
        cmd+=(--include "$f")
    done
    _run "${cmd[@]}"

    _run touch "$local_dir/$_FETCH_DONE"
    _log "OK — $model_id at $local_dir"
}

fetch_prott5() {
    fetch_hf_model "$HF_PROTT5_ID" "prot_t5_xl_uniref50" "${HF_PROTT5_FILES[@]}"
}

fetch_esm1b() {
    fetch_hf_model "$HF_ESM1B_ID" "esm1b_t33_650M_UR50S"
}

fetch_esm2() {
    fetch_hf_model "$HF_ESM2_ID" "esm2_t33_650M_UR50D"
}

fetch_protbert() {
    fetch_hf_model "$HF_PROTBERT_ID" "prot_bert"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_require_command wget "OS package: wget"

_log "Target: $TARGET"
[[ "$DRY_RUN" -eq 1 ]] && _log "(DRY-RUN; no downloads will occur)"

_run mkdir -p "$TARGET"

fetch_deepsece
fetch_plm_effector_trained_models
fetch_prott5
fetch_esm1b
fetch_esm2
fetch_protbert

_log "Done."
_log "Set:"
_log "  SSIGN_PLM_EFFECTOR_WEIGHTS=$TARGET/plm_effector"
_log "  (DeepSecE checkpoint at $TARGET/deepsece_checkpoint.pt is auto-detected.)"
