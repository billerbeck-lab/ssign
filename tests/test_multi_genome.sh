#!/bin/bash
# Multi-genome pipeline test with HHpred and ortholog grouping
set -e

cd /mnt/d/ssign_package
source ~/ssign_test/bin/activate

GENOME1="/mnt/d/ssign_package/tests/data/Xanthobacter_tagetidis_TagT2C_genomic.gbff"
GENOME2="/mnt/d/ssign_package/tests/data/Roseixanthobacter_finlandensis_VTT_E-85241_genomic.gbff"
BASE="/tmp/ssign_multi_test"
rm -rf "$BASE"

echo "============================================="
echo "MULTI-GENOME TEST — 2 genomes + HHpred + OG"
echo "============================================="

# Run each genome through core pipeline (up to filtering)
for IDX in 1 2; do
    if [ "$IDX" = "1" ]; then
        GENOME="$GENOME1"
        SAMPLE="xanthobacter"
    else
        GENOME="$GENOME2"
        SAMPLE="roseixanthobacter"
    fi
    OUTDIR="$BASE/$SAMPLE"
    mkdir -p "$OUTDIR"

    echo ""
    echo "====== Genome $IDX: $SAMPLE ======"
    echo ""

    # Step 1-2: Extract
    echo "--- Detect + Extract ---"
    python src/ssign_app/scripts/detect_input_format.py "$GENOME"
    python src/ssign_app/scripts/extract_proteins.py \
      --input "$GENOME" --sample "$SAMPLE" \
      --out-proteins "$OUTDIR/proteins.faa" \
      --out-gene-info "$OUTDIR/gene_info.tsv" \
      --out-metadata "$OUTDIR/metadata.json"
    N=$(grep -c '^>' "$OUTDIR/proteins.faa")
    echo "Proteins: $N"

    # Gene order
    python src/ssign_app/scripts/extract_gene_order.py \
      --gene-info "$OUTDIR/gene_info.tsv" --output "$OUTDIR/gene_order.tsv"

    # Step 3: MacSyFinder
    echo "--- MacSyFinder ---"
    rm -rf "$OUTDIR/msf_out"
    macsyfinder --sequence-db "$OUTDIR/proteins.faa" --db-type ordered_replicon \
      --models TXSScan all --out-dir "$OUTDIR/msf_out" --mute

    # Step 4: Validate
    echo "--- Validate ---"
    python src/ssign_app/scripts/validate_macsyfinder_systems.py \
      --msf-dir "$OUTDIR/msf_out" --gene-info "$OUTDIR/gene_info.tsv" \
      --sample "$SAMPLE" --wholeness-threshold 0.8 \
      --excluded-systems "Flagellum,Tad,T3SS" \
      --out-components "$OUTDIR/ss_components.tsv" \
      --out-systems "$OUTDIR/valid_systems.tsv"
    echo "Components: $(tail -n +2 "$OUTDIR/ss_components.tsv" | wc -l)"

    # Step 5: Neighborhood
    echo "--- Neighborhood ---"
    python src/ssign_app/scripts/extract_neighborhood.py \
      --gene-order "$OUTDIR/gene_order.tsv" \
      --ss-components "$OUTDIR/ss_components.tsv" \
      --proteins "$OUTDIR/proteins.faa" \
      --window 3 --output "$OUTDIR/neighborhood.faa"
    N_NBHD=$(grep -c '^>' "$OUTDIR/neighborhood.faa")
    echo "Neighborhood: $N_NBHD / $N proteins"

    # Step 6: DLP
    echo "--- DeepLocPro (on $N_NBHD proteins) ---"
    START=$(date +%s)
    python src/ssign_app/scripts/run_deeplocpro.py \
      --input "$OUTDIR/neighborhood.faa" --sample "$SAMPLE" \
      --output "$OUTDIR/deeplocpro.tsv" --mode remote
    END=$(date +%s)
    echo "DLP: $((END-START))s"

    # Step 7: Cross-validate
    echo "--- Cross-validate ---"
    python src/ssign_app/scripts/cross_validate_predictions.py \
      --deeplocpro "$OUTDIR/deeplocpro.tsv" \
      --valid-systems "$OUTDIR/valid_systems.tsv" \
      --sample "$SAMPLE" --conf-threshold 0.8 \
      --output "$OUTDIR/predictions.tsv"

    # Step 8-10: Proximity + T5SS + Filtering
    echo "--- Proximity + T5SS + Filtering ---"
    python src/ssign_app/scripts/proximity_analysis.py \
      --gene-order "$OUTDIR/gene_order.tsv" \
      --ss-components "$OUTDIR/ss_components.tsv" \
      --predictions "$OUTDIR/predictions.tsv" \
      --sample "$SAMPLE" --window 3 --conf-threshold 0.8 \
      --output "$OUTDIR/substrates.tsv"

    python src/ssign_app/scripts/t5ss_handler.py \
      --ss-components "$OUTDIR/ss_components.tsv" \
      --predictions "$OUTDIR/predictions.tsv" \
      --sample "$SAMPLE" \
      --out-substrates "$OUTDIR/t5ss_substrates.tsv" \
      --out-domains "$OUTDIR/t5ss_domains.tsv"

    python src/ssign_app/scripts/system_filtering.py \
      --proximity-substrates "$OUTDIR/substrates.tsv" \
      --t5ss-substrates "$OUTDIR/t5ss_substrates.tsv" \
      --valid-systems "$OUTDIR/valid_systems.tsv" \
      --predictions "$OUTDIR/predictions.tsv" \
      --sample "$SAMPLE" --excluded-systems "Flagellum,Tad,T3SS" \
      --filter-dse-type-mismatch \
      --out-filtered "$OUTDIR/substrates_filtered.tsv" \
      --out-all "$OUTDIR/substrates_all.tsv"
    N_FILT=$(tail -n +2 "$OUTDIR/substrates_filtered.tsv" | wc -l)
    echo "Filtered substrates: $N_FILT"

    # Step 11: BLASTp
    echo "--- BLASTp ---"
    START=$(date +%s)
    python src/ssign_app/scripts/run_blastp.py \
      --mode remote --substrates "$OUTDIR/substrates_filtered.tsv" \
      --proteins "$OUTDIR/proteins.faa" --sample "$SAMPLE" \
      --output "$OUTDIR/blastp.csv" --min-pident 80 --min-qcov 80 --evalue 1e-5
    END=$(date +%s)
    echo "BLASTp: $((END-START))s"

    # Step 12: HHpred (Pfam + PDB)
    echo "--- HHpred ---"
    START=$(date +%s)
    python src/ssign_app/scripts/run_hhsuite.py \
      --mode remote --substrates "$OUTDIR/substrates_filtered.tsv" \
      --proteins "$OUTDIR/proteins.faa" --sample "$SAMPLE" \
      --output "$OUTDIR/hhsuite.csv"
    END=$(date +%s)
    echo "HHpred: $((END-START))s"

    # Step 13: InterProScan
    echo "--- InterProScan ---"
    START=$(date +%s)
    python src/ssign_app/scripts/run_interproscan.py \
      --mode remote --substrates "$OUTDIR/substrates_filtered.tsv" \
      --proteins "$OUTDIR/proteins.faa" --sample "$SAMPLE" \
      --output "$OUTDIR/interproscan.csv"
    END=$(date +%s)
    echo "InterProScan: $((END-START))s"

    # Step 14: Integrate
    echo "--- Integrate ---"
    ANNOT_FILES=""
    for AF in "$OUTDIR/blastp.csv" "$OUTDIR/hhsuite.csv" "$OUTDIR/interproscan.csv"; do
        if [ -f "$AF" ]; then
            ANNOT_FILES="$ANNOT_FILES $AF"
        fi
    done
    python src/ssign_app/scripts/integrate_annotations.py \
      --substrates-filtered "$OUTDIR/substrates_filtered.tsv" \
      --substrates-all "$OUTDIR/substrates_all.tsv" \
      --sample "$SAMPLE" --output "$OUTDIR/integrated.csv" \
      --annotations $ANNOT_FILES

    echo "Genome $IDX done!"
done

# Step 15: Ortholog grouping (across both genomes)
echo ""
echo "====== ORTHOLOG GROUPING (cross-genome) ======"
echo ""

# Combine all substrate sequences
COMBINED_FASTA="$BASE/all_substrates.faa"
> "$COMBINED_FASTA"
for SAMPLE in xanthobacter roseixanthobacter; do
    OUTDIR="$BASE/$SAMPLE"
    # Extract substrate IDs
    SUBS=$(tail -n +2 "$OUTDIR/substrates_filtered.tsv" | cut -f1)
    for SID in $SUBS; do
        python3 -c "
from Bio import SeqIO
for rec in SeqIO.parse('$OUTDIR/proteins.faa', 'fasta'):
    if rec.id == '$SID':
        SeqIO.write(rec, open('$COMBINED_FASTA', 'a'), 'fasta')
        break
"
    done
done
N_COMBINED=$(grep -c '^>' "$COMBINED_FASTA")
echo "Combined substrates: $N_COMBINED"

# Run ortholog grouping
echo "--- Ortholog grouping ---"
python src/ssign_app/scripts/run_ortholog_grouping.py \
  --substrates-fasta "$COMBINED_FASTA" \
  --min-pident 40 --min-qcov 70 \
  --output "$BASE/ortholog_assignments.csv" \
  --output-groups "$BASE/ortholog_groups.csv"

echo ""
echo "--- Ortholog assignments ---"
cat "$BASE/ortholog_assignments.csv"
echo ""
echo "--- Ortholog groups ---"
cat "$BASE/ortholog_groups.csv"

# Figures
echo ""
echo "====== FIGURES ======"
mkdir -p "$BASE/figures"
# Combine integrated CSVs
cat "$BASE/xanthobacter/integrated.csv" > "$BASE/combined_integrated.csv"
tail -n +2 "$BASE/roseixanthobacter/integrated.csv" >> "$BASE/combined_integrated.csv"

python src/ssign_app/scripts/generate_figures.py \
  --master-csvs "$BASE/combined_integrated.csv" \
  --outdir "$BASE/figures" --dpi 150

echo ""
echo "============================================="
echo "MULTI-GENOME TEST COMPLETE"
echo "============================================="
echo ""
echo "Output files:"
for SAMPLE in xanthobacter roseixanthobacter; do
    echo "  $SAMPLE:"
    ls "$BASE/$SAMPLE/"*.csv "$BASE/$SAMPLE/"*.tsv 2>/dev/null | while read f; do
        echo "    $(basename $f)"
    done
done
echo ""
echo "Cross-genome:"
ls "$BASE/"*.csv 2>/dev/null | while read f; do echo "  $(basename $f)"; done
echo ""
echo "Figures:"
ls "$BASE/figures/" 2>/dev/null
