#!/bin/bash
#SBATCH -p Cnode
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -J bowtie_host
#SBATCH -o %j.log
#SBATCH -e %j.err


cd ~
cd LWZ

source ~/.bashrc
source activate Bowtie2

QC_DIR="QC"
HOST_REF="host.fna"
OUT_DIR="host"
HOST_INDEX_PREFIX="host"


echo "Current working directory: $(pwd)"


if [ ! -f "${HOST_INDEX_PREFIX}.1.bt2" ] && [ ! -f "${HOST_INDEX_PREFIX}.1.bt2l" ]; then
    echo "Index not found. Building Bowtie2 index..."
    bowtie2-build --threads 16 "$HOST_REF" "$HOST_INDEX_PREFIX"
else
    echo "Host index found."
fi


echo "Scanning for Plant samples in $QC_DIR..."

for r1_file in "$QC_DIR"/*-P*.R1.raw_val_1.fq.gz; do
    
    if [ ! -f "$r1_file" ]; then
        echo "No Plant (P) samples found matching pattern."
        exit 1
    fi

    filename=$(basename "$r1_file")
    sample_name=${filename%%.R1.raw_val_1.fq.gz}

    r2_file="${r1_file/R1.raw_val_1.fq.gz/R2.raw_val_2.fq.gz}"

    echo "-------------------------------------------"
    echo "Task: Host Removal for $sample_name"
    echo "Input R1: $r1_file" 

    if [ -f "$r2_file" ]; then
        
        out_r1="${OUT_DIR}/${sample_name}_clean_1.fq.gz"
        out_r2="${OUT_DIR}/${sample_name}_clean_2.fq.gz"

        (
            bowtie2 -x "$HOST_INDEX_PREFIX" \
                -1 "$r1_file" \
                -2 "$r2_file" \
                --very-sensitive-local \
                -p 16 \
                2> "${OUT_DIR}/${sample_name}_bowtie2.log" | \
            samtools view -b -f 4 -@ 16 - | \
            samtools sort -n -@ 16 - | \
            samtools fastq -@ 16 \
                -1 "$out_r1" \
                -2 "$out_r2" \
                -0 /dev/null -s /dev/null -n -
        )

        if [ $? -eq 0 ]; then
            echo "Success: $sample_name"
        else
            echo "Error: Failed to process $sample_name"
        fi

    else
        echo "Error: R2 file not found: $r2_file"
    fi
done

echo "All tasks finished at $(date)"