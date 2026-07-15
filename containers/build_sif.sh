#!/usr/bin/env bash
# Build the ssign Singularity/Apptainer image and smoke-test it.
#
# Run from anywhere; must be on Linux with apptainer installed
# (Arch: `sudo pacman -S apptainer`; most HPCs already provide it).
#
# Usage:
#   containers/build_sif.sh [/path/to/dir/containing/deeplocpro]
#
# The DeepLocPro directory (arg, or $SSIGN_DEEPLOCPRO_PATH) enables the full
# offline golden run; without it the script does build + doctor smoke only
# (DeepLocPro is DTU-licensed and not bundled).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
SIF="${SSIGN_SIF:-$HOME/ssign.sif}"
DLP_PATH="${1:-${SSIGN_DEEPLOCPRO_PATH:-}}"
FIXTURE="tests/fixtures/Xanthobacter_T5aSS_minimal.gbff"
EXPECTED_SUBSTRATE="BIMENO_04457"

# The build sandbox + torch install need many GB. If TMPDIR is a small RAM-backed
# tmpfs (common for /tmp), the build dies with "disk quota exceeded". Force the
# Apptainer build tmp/cache onto real disk under $HOME.
BUILDTMP="${SSIGN_BUILD_TMPDIR:-$HOME/.ssign_sif_build}"
mkdir -p "$BUILDTMP/tmp" "$BUILDTMP/cache"
export APPTAINER_TMPDIR="$BUILDTMP/tmp" APPTAINER_CACHEDIR="$BUILDTMP/cache" TMPDIR="$BUILDTMP/tmp"

# Persistent build-weight cache (#36): the def bakes ~10 GB of model weights
# (ESM2, ESM1b, DeepSecE checkpoint, ProtT5). Bind this host dir to /build_cache
# so the FIRST build populates it and every rebuild after reuses it instead of
# re-downloading. Survives across builds (lives under $HOME, not the build tmp).
SSIGN_BUILD_CACHE="${SSIGN_BUILD_CACHE:-$HOME/.ssign_build_weights}"
mkdir -p "$SSIGN_BUILD_CACHE"

echo "== 1/3  ensure uv.lock is current with pyproject.toml =="
if ! uv lock --check 2>/dev/null; then
    echo "   uv.lock stale -> running uv lock"
    uv lock
fi

echo "== 2/3  apptainer build (fakeroot) -> $SIF =="
echo "   build-weight cache: $SSIGN_BUILD_CACHE ($(du -sh "$SSIGN_BUILD_CACHE" 2>/dev/null | cut -f1))"
apptainer build --fakeroot --bind "$SSIGN_BUILD_CACHE:/build_cache" "$SIF" containers/ssign.def

echo "== 3/3  smoke test (offline, --containall) =="
apptainer run --containall "$SIF" --version
apptainer run --writable-tmpfs --containall "$SIF" doctor --tier base | sed -n '1,40p'

if [ -n "$DLP_PATH" ] && [ -d "$DLP_PATH" ]; then
    echo "-- offline golden run (DeepLocPro bind-mounted from $DLP_PATH) --"
    OUT="$(mktemp -d)"; chmod 777 "$OUT"
    apptainer run --writable-tmpfs --containall \
        -B "$REPO/$FIXTURE:/work/in.gbff:ro" \
        -B "$DLP_PATH:/opt/deeplocpro:ro" \
        -B "$OUT:/work/out" \
        "$SIF" run /work/in.gbff --outdir /work/out \
            --use-input-annotations \
            --deeplocpro-mode local --deeplocpro-path /opt/deeplocpro \
            --skip-deepsece --skip-signalp --skip-blastp \
            --skip-hhsuite --skip-interproscan --skip-plmblast --skip-protparam
    if grep -q "$EXPECTED_SUBSTRATE" "$OUT"/*results*.csv 2>/dev/null; then
        echo "PASS: golden substrate $EXPECTED_SUBSTRATE found"
    else
        echo "FAIL: $EXPECTED_SUBSTRATE not found in $OUT"; exit 1
    fi
else
    echo "(no DeepLocPro path -> golden genome run skipped; pass the dir containing"
    echo " the deeplocpro executable to enable the full offline check)"
fi
echo "DONE.  image: $SIF"
