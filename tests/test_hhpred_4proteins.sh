#!/bin/bash
# Test HHpred on the 4 Xanthobacter substrate proteins
set -e
cd /mnt/d/ssign_package
source ~/ssign_test/bin/activate

# Recreate the needed files
OUTDIR="/tmp/ssign_hhpred_test"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

GENOME="/mnt/d/ssign_package/tests/data/Xanthobacter_tagetidis_TagT2C_genomic.gbff"

echo "=== Setting up test data ==="
python src/ssign_app/scripts/extract_proteins.py \
  --input "$GENOME" --sample test \
  --out-proteins "$OUTDIR/proteins.faa" \
  --out-gene-info "$OUTDIR/gene_info.tsv" \
  --out-metadata "$OUTDIR/metadata.json"

python src/ssign_app/scripts/extract_gene_order.py \
  --gene-info "$OUTDIR/gene_info.tsv" --output "$OUTDIR/gene_order.tsv"

rm -rf "$OUTDIR/msf_out"
macsyfinder --sequence-db "$OUTDIR/proteins.faa" --db-type ordered_replicon \
  --models TXSScan all --out-dir "$OUTDIR/msf_out" --mute

python src/ssign_app/scripts/validate_macsyfinder_systems.py \
  --msf-dir "$OUTDIR/msf_out" --gene-info "$OUTDIR/gene_info.tsv" \
  --sample test --wholeness-threshold 0.8 \
  --excluded-systems "Flagellum,Tad,T3SS" \
  --out-components "$OUTDIR/ss_components.tsv" \
  --out-systems "$OUTDIR/valid_systems.tsv"

python src/ssign_app/scripts/extract_neighborhood.py \
  --gene-order "$OUTDIR/gene_order.tsv" \
  --ss-components "$OUTDIR/ss_components.tsv" \
  --proteins "$OUTDIR/proteins.faa" \
  --window 3 --output "$OUTDIR/neighborhood.faa"

python src/ssign_app/scripts/run_deeplocpro.py \
  --input "$OUTDIR/neighborhood.faa" --sample test \
  --output "$OUTDIR/deeplocpro.tsv" --mode remote

python src/ssign_app/scripts/cross_validate_predictions.py \
  --deeplocpro "$OUTDIR/deeplocpro.tsv" \
  --valid-systems "$OUTDIR/valid_systems.tsv" \
  --sample test --conf-threshold 0.8 --output "$OUTDIR/predictions.tsv"

python src/ssign_app/scripts/proximity_analysis.py \
  --gene-order "$OUTDIR/gene_order.tsv" \
  --ss-components "$OUTDIR/ss_components.tsv" \
  --predictions "$OUTDIR/predictions.tsv" \
  --sample test --window 3 --conf-threshold 0.8 \
  --output "$OUTDIR/substrates.tsv"

python src/ssign_app/scripts/t5ss_handler.py \
  --ss-components "$OUTDIR/ss_components.tsv" \
  --predictions "$OUTDIR/predictions.tsv" --sample test \
  --out-substrates "$OUTDIR/t5ss_substrates.tsv" \
  --out-domains "$OUTDIR/t5ss_domains.tsv"

python src/ssign_app/scripts/system_filtering.py \
  --proximity-substrates "$OUTDIR/substrates.tsv" \
  --t5ss-substrates "$OUTDIR/t5ss_substrates.tsv" \
  --valid-systems "$OUTDIR/valid_systems.tsv" \
  --predictions "$OUTDIR/predictions.tsv" \
  --sample test --excluded-systems "Flagellum,Tad,T3SS" \
  --filter-dse-type-mismatch \
  --out-filtered "$OUTDIR/substrates_filtered.tsv" \
  --out-all "$OUTDIR/substrates_all.tsv"

N_SUB=$(tail -n +2 "$OUTDIR/substrates_filtered.tsv" | wc -l)
echo "Substrates: $N_SUB"

echo ""
echo "=== Testing HHpred (Pfam + PDB) on $N_SUB substrates ==="
START=$(date +%s)
python -u src/ssign_app/scripts/run_hhsuite.py \
  --mode remote \
  --substrates "$OUTDIR/substrates_filtered.tsv" \
  --proteins "$OUTDIR/proteins.faa" \
  --sample test \
  --output "$OUTDIR/hhsuite.csv"
END=$(date +%s)
echo ""
echo "=== HHpred DONE in $((END-START))s ==="
echo ""
echo "=== Results ==="
cat "$OUTDIR/hhsuite.csv"
echo ""
N_HITS=$(tail -n +2 "$OUTDIR/hhsuite.csv" | wc -l)
echo "HHpred hits: $N_HITS / $N_SUB substrates"
