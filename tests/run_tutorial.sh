#!/usr/bin/env bash
# tests/run_tutorial.sh — end-to-end smoke test on E. coli K-12.
#
# Runs ssign on a real K-12 genome via the DTU webserver fallback (so the
# smoke test does not require local DTU installs in CI), checks expected
# outputs exist + are non-empty, and spot-checks that the T2SS gsp operon
# shows up in the secretion-systems section. Use this as a smoke test
# before tagging a release, or after touching anything in the run-pipeline
# path.
#
# Usage:
#     bash tests/run_tutorial.sh            # ~30 min, internet required
#     bash tests/run_tutorial.sh --force    # re-run even if output exists
#
# Requires: ssign on PATH (pip install ssign), wget, gunzip, python3.
# Internet required for: K-12 GenBank download + DTU webserver calls.

set -euo pipefail

WORK_DIR="${HOME}/ssign-tutorial-test"
OUTDIR="${WORK_DIR}/ecoli_results"
GENBANK="${WORK_DIR}/ecoli_k12.gbff"
GENBANK_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

# ── Sanity checks ────────────────────────────────────────────────────────
command -v ssign >/dev/null || { echo "FAIL: ssign not on PATH. Install with 'pip install ssign'." >&2; exit 2; }
command -v wget  >/dev/null || { echo "FAIL: wget not installed." >&2; exit 2; }
command -v python3 >/dev/null || { echo "FAIL: python3 not installed." >&2; exit 2; }

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Download K-12 (~13 MB) ──────────────────────────────────────────────
if [[ ! -f "$GENBANK" ]]; then
    echo ">> Downloading E. coli K-12 GenBank from NCBI..."
    wget -q "$GENBANK_URL" -O ecoli_k12.gbff.gz
    gunzip -f ecoli_k12.gbff.gz
    [[ -f "$GENBANK" ]] || { echo "FAIL: download did not produce $GENBANK" >&2; exit 1; }
fi
echo ">> K-12 GenBank: $(du -h "$GENBANK" | cut -f1)"

# ── Clean output ─────────────────────────────────────────────────────────
if [[ -d "$OUTDIR" ]]; then
    if (( FORCE )); then
        echo ">> --force: clearing $OUTDIR"
        rm -rf "$OUTDIR"
    else
        echo "FAIL: $OUTDIR exists. Pass --force to re-run, or delete manually." >&2
        exit 1
    fi
fi

# ── Run ssign (base-tier command, webserver DTU fallback) ────────────────
echo ">> Running ssign — DLP + SignalP via DTU webserver (opt-in fallback used here so CI does not need a DTU install), ~30 min total..."
START=$(date +%s)
ssign run "$GENBANK" --outdir "$OUTDIR" \
    --use-input-annotations \
    --signalp-mode remote \
    --deeplocpro-mode remote \
    --skip-blastp \
    --skip-eggnog \
    --skip-hhsuite \
    --skip-interproscan \
    --skip-plmblast
ELAPSED=$(( $(date +%s) - START ))
echo ">> ssign finished in ${ELAPSED}s"

# ── Verify outputs ───────────────────────────────────────────────────────
echo ">> Checking output files..."
EXPECTED=(
    "$OUTDIR/ecoli_k12_results.csv"
    "$OUTDIR/ecoli_k12_results_raw.csv"
    "$OUTDIR/ecoli_k12_summary.txt"
    "$OUTDIR/figures/ecoli_k12"
)
for f in "${EXPECTED[@]}"; do
    if [[ ! -e "$f" ]]; then
        echo "FAIL: missing $f" >&2
        exit 1
    fi
    if [[ -f "$f" && ! -s "$f" ]]; then
        echo "FAIL: $f is empty" >&2
        exit 1
    fi
done
echo "   all expected outputs present + non-empty"

# ── Spot-check: T2SS gsp operon must appear in the K-12 results ─────────
echo ">> Spot-checking K-12 biology..."
RESULTS="$OUTDIR/ecoli_k12_results.csv"
if ! grep -q "^# Secretion Systems (with secreted proteins)" "$RESULTS"; then
    echo "FAIL: no '# Secretion Systems (with secreted proteins)' chunk header in $RESULTS" >&2
    exit 1
fi
T2SS_COUNT=$(grep -c -E ",T2SS,|T2SS_gsp" "$RESULTS" || true)
if (( T2SS_COUNT == 0 )); then
    echo "FAIL: no T2SS rows in K-12 results — gsp operon not detected?" >&2
    exit 1
fi
echo "   T2SS rows found: $T2SS_COUNT"

echo
echo "=========================================="
echo "TUTORIAL VERIFICATION PASSED in ${ELAPSED}s"
echo "=========================================="
echo "Outputs: $OUTDIR"
