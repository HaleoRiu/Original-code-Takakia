#!/bin/bash
#SBATCH -p Cnode
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J eggnog_fix
#SBATCH -o %j.log
#SBATCH -e %j.err

export LC_ALL=C


cd ~
cd LWZ 

source ~/.bashrc
conda activate eggnog-envcwj

GENE_FASTA="Abundance_Result/all_MAG_genes_prefixed.fna"
ABUNDANCE_TABLE="Abundance_Result/Final_Gene_Abundance.tsv"

OUT_DIR="EggNOG_Result"
mkdir -p "$OUT_DIR"

echo "Step 1: Running EggNOG-mapper..."

emapper.py \
    -i "$GENE_FASTA" \
    --itype CDS \
    --output "mag_genes" \
    --output_dir "$OUT_DIR" \
    --cpu 32 \
    -m diamond \
    --evalue 0.00001 \
    --translate \
    --override \
    --no_file_comments


ANNOTATION_FILE="${OUT_DIR}/mag_genes.emapper.annotations"

if [ ! -s "$ANNOTATION_FILE" ]; then
    echo "Error: EggNOG annotation failed or output is empty."
    echo "Please check the log above for diamond errors."
    exit 1
fi

echo "EggNOG annotation successful. Proceeding to calculation..."


echo "Step 2: Calculating KO abundance (Weighted logic)..."

OUTPUT_KO_TABLE="${OUT_DIR}/Final_KO_Abundance_Weighted.tsv"

awk '
BEGIN {
    FS = "\t"; OFS = "\t"
}

FNR == NR {
    if ($0 ~ /^#/) next;
    
    gene = $1
    pathways = $12
    
    if (pathways != "-" && pathways != "") {
        gsub("ko:", "", pathways)
        split(pathways, pathway_array, ",")
        
        num_pathways = 0
        delete temp_path_list 
        
        for (p in pathway_array) {
            # 筛选 K 号
            if (pathway_array[p] ~ /^K[0-9]+$/) {
                temp_path_list[pathway_array[p]] = 1
                num_pathways++
            }
        }
        
        if (num_pathways > 0) {
            gene_pathway_count[gene] = num_pathways
            for (kp in temp_path_list) {
                map_gene_ko_str[gene] = map_gene_ko_str[gene] ? map_gene_ko_str[gene] "," kp : kp
            }
        }
    }
    next
}

{
    if (FNR == 1) {
        for (i = 2; i <= NF; i++) {
            samples[i] = $i
        }
        num_samples = NF
        next
    }
    
    gene = $1
    if (gene in gene_pathway_count) {
        N = gene_pathway_count[gene]
        split(map_gene_ko_str[gene], ko_list, ",")
        
        for (i = 2; i <= num_samples; i++) {
            # 核心：丰度平分给 N 个功能
            weighted_abundance = $i / N
            for (k in ko_list) {
                ko_id = ko_list[k]
                sum_pathway[ko_id, i] += weighted_abundance
                ko_exists[ko_id] = 1
            }
        }
    }
}

END {
    printf "KO_ID" > "'"$OUTPUT_KO_TABLE"'"
    for (i = 2; i <= num_samples; i++) {
        printf "\t%s", samples[i] > "'"$OUTPUT_KO_TABLE"'"
    }
    print "" > "'"$OUTPUT_KO_TABLE"'"
    
    for (ko in ko_exists) {
        printf "%s", ko > "'"$OUTPUT_KO_TABLE"'"
        for (i = 2; i <= num_samples; i++) {
            printf "\t%.4f", sum_pathway[ko, i] > "'"$OUTPUT_KO_TABLE"'"
        }
        print "" > "'"$OUTPUT_KO_TABLE"'"
    }
}
' "$ANNOTATION_FILE" "$ABUNDANCE_TABLE"

echo "Calculation finished. Result: $OUTPUT_KO_TABLE"
echo "Job finished at $(date)"