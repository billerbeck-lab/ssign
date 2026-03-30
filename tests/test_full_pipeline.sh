#!/bin/bash
# Full pipeline end-to-end test with neighborhood optimization
set -e

cd /mnt/d/ssign_package
source ~/ssign_test/bin/activate

GENOME="/mnt/d/ssign_package/tests/data/Xanthobacter_tagetidis_TagT2C_genomic.gbff"
OUTDIR="/tmp/ssign_full_test"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "=========================================="
echo "FULL PIPELINE TEST — Xanthobacter TagT2C"
echo "=========================================="
echo ""

# Step 1: Detect format
echo "--- Step 1: Detect input format ---"
python src/ssign_app/scripts/detect_input_format.py "$GENOME"
echo ""

# Step 2: Extract proteins
echo "--- Step 2: Extract proteins ---"
python src/ssign_app/scripts/extract_proteins.py \
  --input "$GENOME" \
  --sample test \
  --out-proteins "$OUTDIR/proteins.faa" \
  --out-gene-info "$OUTDIR/gene_info.tsv" \
  --out-metadata "$OUTDIR/metadata.json"
N_PROTS=$(grep -c '^>' "$OUTDIR/proteins.faa")
echo "Extracted: $N_PROTS proteins"
echo ""

# Step 2b: Gene order
echo "--- Step 2b: Extract gene order ---"
python src/ssign_app/scripts/extract_gene_order.py \
  --gene-info "$OUTDIR/gene_info.tsv" \
  --output "$OUTDIR/gene_order.tsv"
echo "Gene order lines: $(wc -l < "$OUTDIR/gene_order.tsv")"
echo ""

# Step 3: MacSyFinder
echo "--- Step 3: MacSyFinder ---"
rm -rf "$OUTDIR/msf_out"
macsyfinder \
  --sequence-db "$OUTDIR/proteins.faa" \
  --db-type ordered_replicon \
  --models TXSScan all \
  --out-dir "$OUTDIR/msf_out" \
  --mute 2>&1
echo "MacSyFinder done"
echo ""

# Step 4: Validate systems
echo "--- Step 4: Validate systems ---"
python src/ssign_app/scripts/validate_macsyfinder_systems.py \
  --msf-dir "$OUTDIR/msf_out" \
  --gene-info "$OUTDIR/gene_info.tsv" \
  --sample test \
  --wholeness-threshold 0.8 \
  --excluded-systems "Flagellum,Tad,T3SS" \
  --out-components "$OUTDIR/ss_components.tsv" \
  --out-systems "$OUTDIR/valid_systems.tsv"
N_COMP=$(tail -n +2 "$OUTDIR/ss_components.tsv" | wc -l)
echo "SS components: $N_COMP"
echo ""

# Step 5: Extract neighborhood (NEW — optimization!)
echo "--- Step 5: Extract neighborhood ---"
START_N=$(date +%s)
python src/ssign_app/scripts/extract_neighborhood.py \
  --gene-order "$OUTDIR/gene_order.tsv" \
  --ss-components "$OUTDIR/ss_components.tsv" \
  --proteins "$OUTDIR/proteins.faa" \
  --window 3 \
  --output "$OUTDIR/neighborhood.faa"
END_N=$(date +%s)
N_NBHD=$(grep -c '^>' "$OUTDIR/neighborhood.faa")
echo "Neighborhood: $N_NBHD / $N_PROTS proteins ($(echo "scale=1; 100*$N_NBHD/$N_PROTS" | bc)%)"
echo "Time: $((END_N - START_N))s"
echo ""

# Step 6: DeepLocPro (on neighborhood only!)
echo "--- Step 6: DeepLocPro (on $N_NBHD neighborhood proteins) ---"
START_DLP=$(date +%s)
python src/ssign_app/scripts/run_deeplocpro.py \
  --input "$OUTDIR/neighborhood.faa" \
  --sample test \
  --output "$OUTDIR/deeplocpro.tsv" \
  --mode remote
END_DLP=$(date +%s)
DLP_TIME=$((END_DLP - START_DLP))
N_DLP=$(tail -n +2 "$OUTDIR/deeplocpro.tsv" | wc -l)
echo "DLP: $N_DLP predictions in ${DLP_TIME}s"
echo ""

# Step 7: Cross-validate
echo "--- Step 7: Cross-validate ---"
python src/ssign_app/scripts/cross_validate_predictions.py \
  --deeplocpro "$OUTDIR/deeplocpro.tsv" \
  --valid-systems "$OUTDIR/valid_systems.tsv" \
  --sample test \
  --conf-threshold 0.8 \
  --output "$OUTDIR/predictions.tsv"
echo ""

# Step 8: Proximity analysis
echo "--- Step 8: Proximity analysis ---"
python src/ssign_app/scripts/proximity_analysis.py \
  --gene-order "$OUTDIR/gene_order.tsv" \
  --ss-components "$OUTDIR/ss_components.tsv" \
  --predictions "$OUTDIR/predictions.tsv" \
  --sample test \
  --window 3 \
  --conf-threshold 0.8 \
  --output "$OUTDIR/substrates.tsv"
N_SUB=$(tail -n +2 "$OUTDIR/substrates.tsv" | wc -l)
echo "Substrates found: $N_SUB"
echo ""

# Step 9: T5SS handler
echo "--- Step 9: T5SS handler ---"
python src/ssign_app/scripts/t5ss_handler.py \
  --ss-components "$OUTDIR/ss_components.tsv" \
  --predictions "$OUTDIR/predictions.tsv" \
  --sample test \
  --out-substrates "$OUTDIR/t5ss_substrates.tsv" \
  --out-domains "$OUTDIR/t5ss_domains.tsv"
echo ""

# Step 10: System filtering
echo "--- Step 10: System filtering ---"
python src/ssign_app/scripts/system_filtering.py \
  --proximity-substrates "$OUTDIR/substrates.tsv" \
  --t5ss-substrates "$OUTDIR/t5ss_substrates.tsv" \
  --valid-systems "$OUTDIR/valid_systems.tsv" \
  --predictions "$OUTDIR/predictions.tsv" \
  --sample test \
  --excluded-systems "Flagellum,Tad,T3SS" \
  --filter-dse-type-mismatch \
  --out-filtered "$OUTDIR/substrates_filtered.tsv" \
  --out-all "$OUTDIR/substrates_all.tsv"
N_FILT=$(tail -n +2 "$OUTDIR/substrates_filtered.tsv" | wc -l)
echo "Filtered substrates: $N_FILT"
echo ""

# Step 11: BLASTp (on filtered substrates only)
echo "--- Step 11: BLASTp (on $N_FILT substrates) ---"
START_BP=$(date +%s)
python src/ssign_app/scripts/run_blastp.py \
  --mode remote \
  --substrates "$OUTDIR/substrates_filtered.tsv" \
  --proteins "$OUTDIR/proteins.faa" \
  --sample test \
  --output "$OUTDIR/blastp.csv" \
  --min-pident 80 \
  --min-qcov 80 \
  --evalue 1e-5 2>&1
END_BP=$(date +%s)
echo "BLASTp: $((END_BP - START_BP))s"
echo ""

# Step 12: InterProScan (on filtered substrates only)
echo "--- Step 12: InterProScan (on $N_FILT substrates) ---"
START_IPRS=$(date +%s)
python src/ssign_app/scripts/run_interproscan.py \
  --mode remote \
  --substrates "$OUTDIR/substrates_filtered.tsv" \
  --proteins "$OUTDIR/proteins.faa" \
  --sample test \
  --output "$OUTDIR/interproscan.csv" 2>&1
END_IPRS=$(date +%s)
echo "InterProScan: $((END_IPRS - START_IPRS))s"
echo ""

# Step 13: Integrate
echo "--- Step 13: Integrate annotations ---"
python src/ssign_app/scripts/integrate_annotations.py \
  --substrates-filtered "$OUTDIR/substrates_filtered.tsv" \
  --substrates-all "$OUTDIR/substrates_all.tsv" \
  --sample test \
  --output "$OUTDIR/integrated.csv" \
  --annotations "$OUTDIR/blastp.csv" "$OUTDIR/interproscan.csv" 2>&1
echo ""

# Step 14: Figures (with toggles)
echo "--- Step 14: Generate figures ---"
mkdir -p "$OUTDIR/figures"
python src/ssign_app/scripts/generate_figures.py \
  --master-csvs "$OUTDIR/integrated.csv" \
  --outdir "$OUTDIR/figures" \
  --dpi 150 2>&1
echo ""

echo "=========================================="
echo "PIPELINE COMPLETE"
echo "=========================================="
echo ""
echo "Output files:"
ls -la "$OUTDIR/"*.tsv "$OUTDIR/"*.csv 2>/dev/null
echo ""
echo "Figures:"
ls -la "$OUTDIR/figures/" 2>/dev/null
echo ""
echo "KEY STATS:"
echo "  Total proteins: $N_PROTS"
echo "  Neighborhood:   $N_NBHD ($(echo "scale=1; 100*$N_NBHD/$N_PROTS" | bc)%)"
echo "  DLP time:       ${DLP_TIME}s (on $N_NBHD instead of $N_PROTS proteins)"
echo "  Substrates:     $N_FILT"
