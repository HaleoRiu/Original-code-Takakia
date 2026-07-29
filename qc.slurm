#!/bin/bash
#SBATCH -p Cnode
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J trim_qc
#SBATCH -o %j.log
#SBATCH -e %j.err

cd ~
cd ~/LWZ

source ~/.bashrc
source activate trim_env

input_dir="rawdata" 
output_dir="QC"


for r1_file in "$input_dir"/*.R1.raw.fastq.gz; do

    [ -f "$r1_file" ] || { echo "No R1 files found in $input_dir"; exit 1; }
    

    echo "Processing R1: $r1_file"
    r2_file=$(echo "$r1_file" | sed 's/R1\.raw\.fastq\.gz/R2.raw.fastq.gz/')
    echo "Inferred R2: $r2_file"
    

    if [ -f "$r2_file" ]; then
        # 运行 trim_galore
        trim_galore --quality 30 --length 50 --paired --cores 4 --output_dir "$output_dir" "$r1_file" "$r2_file"
    else
        echo "Warning: Corresponding R2 file not found for $r1_file"
    fi
done

echo "Job finished at $(date)"