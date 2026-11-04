#!/bin/bash
#SBATCH -p Cnode
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -J mag_abun
#SBATCH -o %j.log
#SBATCH -e %j.err

export LC_ALL=C


cd ~
cd ~/LWZ


MAGs_DIR="MAGs"
QC_DIR="QC"
HOST_DIR="host"
OUT_DIR="Abundance_Result"


GENE_PRED_DIR="${OUT_DIR}/gene_predictions"
mkdir -p "$GENE_PRED_DIR"

source ~/.bashrc

echo ">>> Phase 1: Activating prodigal_env for gene prediction..."

source activate prodigal_env

if ! command -v prodigal &> /dev/null; then
    echo "Error: prodigal not found. Environment activation failed."
    exit 1
fi

MERGED_REF="${OUT_DIR}/all_MAG_genes_prefixed.fna"
> "$MERGED_REF" 

echo "Predicting genes from MAGs..."

shopt -s nullglob
for mag_file in "$MAGs_DIR"/*.{fa,fasta}; do
    [ -f "$mag_file" ] || continue

    basename=$(basename "$mag_file")
    mag_name=${basename%.*} 
    
    prodigal -i "$mag_file" -d "$GENE_PRED_DIR/${mag_name}.genes.fna" -p meta -q
    
    sed "s/^>/>${mag_name}_/" "$GENE_PRED_DIR/${mag_name}.genes.fna" >> "$MERGED_REF"
done
shopt -u nullglob

echo "Gene prediction finished. Merged reference: $MERGED_REF"


echo ">>> Switching environment..."

source deactivate

source activate coverm

if ! command -v coverm &> /dev/null; then
    echo "Error: coverm not found after switching environment."
    exit 1
fi

echo "Environment switched to: $CONDA_DEFAULT_ENV"


echo ">>> Phase 2: Collecting reads and running CoverM..."

declare -a ALL_READS


for r1 in "$HOST_DIR"/*_clean_1.fq.gz; do
    if [ -f "$r1" ]; then
        r2="${r1/_clean_1.fq.gz/_clean_2.fq.gz}"
        if [ -f "$r2" ]; then
            ALL_READS+=("$r1" "$r2")
        fi
    fi
done


for r1 in "$QC_DIR"/*-S*.R1.raw_val_1.fq.gz; do
    if [ -f "$r1" ]; then
        # Replacement rule: R1.raw_val_1 -> R2.raw_val_2
        r2="${r1/R1.raw_val_1.fq.gz/R2.raw_val_2.fq.gz}"
        
        if [ -f "$r2" ]; then
            ALL_READS+=("$r1" "$r2")
        else
            echo "Warning: Missing R2 for Soil sample $(basename "$r1")"
        fi
    fi
done

num_files=${#ALL_READS[@]}
echo "Total paired files collected (R1+R2): $num_files"

if [ "$num_files" -eq 0 ]; then
    echo "Error: No valid reads found in host/ or QC/ directories."
    exit 1
fi

echo "Running CoverM mapping..."

coverm contig \
    --reference "$MERGED_REF" \
    --threads 24 \
    --methods tpm count mean \
    --min-read-percent-identity 95 \
    --coupled "${ALL_READS[@]}" \
    --output-file "$OUT_DIR/Final_Gene_Abundance.tsv"

echo "All analysis finished at $(date)"