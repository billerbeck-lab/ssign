#!/usr/bin/env bash
# Build the ssign release image and run the offline golden smoke test.
#
# Run on a CUDA-capable Docker host (the laptop can't build the CUDA image).
#
# Usage:
#   containers/build_and_test.sh [/path/to/dir/containing/deeplocpro]
#
# The DeepLocPro directory (arg, or $SSIGN_DEEPLOCPRO_PATH) enables the full
# offline golden run; without it the script does a container-starts smoke only
# (DeepLocPro is DTU-licensed and not bundled in the image — see the Dockerfile).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
IMAGE="${SSIGN_IMAGE:-ssign:1.0.0}"
DLP_PATH="${1:-${SSIGN_DEEPLOCPRO_PATH:-}}"
FIXTURE="tests/fixtures/Xanthobacter_T5aSS_minimal.gbff"
EXPECTED_SUBSTRATE="BIMENO_04457"

echo "== 1/3  ensure uv.lock is current with pyproject.toml =="
# --frozen export in the Dockerfile silently omits deps missing from a stale
# lock; refresh so the image gets the full [extended] set.
if ! uv lock --check 2>/dev/null; then
    echo "   uv.lock stale -> running uv lock"
    uv lock
fi

echo "== 2/3  docker build $IMAGE =="
docker build -f containers/Dockerfile -t "$IMAGE" .

echo "== 3/3  smoke test =="
echo "-- container starts + ssign runs --"
docker run --rm --network none "$IMAGE" --version

if [ -n "$DLP_PATH" ] && [ -d "$DLP_PATH" ]; then
    echo "-- offline golden run (DeepLocPro bind-mounted from $DLP_PATH) --"
    # The container runs as uid 1000 (USER ssign); make the output dir writable
    # by it regardless of the host user's uid (else the run fails on permission,
    # not logic). Overriding --user instead would break the baked-in
    # /home/ssign/.macsyfinder TXSScan profiles, so widen the temp dir instead.
    OUT="$(mktemp -d)"; chmod 777 "$OUT"
    docker run --rm --network none \
        -v "$REPO/$FIXTURE:/work/in.gbff:ro" \
        -v "$DLP_PATH:/opt/deeplocpro:ro" \
        -v "$OUT:/work/out" \
        "$IMAGE" run /work/in.gbff --outdir /work/out \
            --use-input-annotations \
            --deeplocpro-mode local --deeplocpro-path /opt/deeplocpro \
            --skip-deepsece --skip-signalp --skip-plm-effector --skip-blastp \
            --skip-hhsuite --skip-interproscan --skip-plmblast --skip-protparam
    if grep -q "$EXPECTED_SUBSTRATE" "$OUT"/*results*.csv 2>/dev/null; then
        echo "PASS: golden substrate $EXPECTED_SUBSTRATE found in output"
    else
        echo "FAIL: $EXPECTED_SUBSTRATE not found in $OUT"
        exit 1
    fi
else
    echo "(no DeepLocPro path -> golden run skipped; pass the dir containing the"
    echo " deeplocpro executable to enable the full offline check)"
fi
echo "DONE."
