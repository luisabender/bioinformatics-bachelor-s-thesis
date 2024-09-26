#!/bin/bash

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"
output_fastqc="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_results_trimmed"

mkdir "$output_fastqc"

for dir in "$base_dir"/*/"fastq"/; do
    echo "Processing directory: $dir"

    
    for fastq_file in "$dir"*.fastq-cut; do
        if [ -f "$fastq_file" ]; then
            echo "Running FastQC on $fastq_file"
            fastqc "$fastq_file" --outdir="$output_fastqc"
        fi
    done
done

# move trimmed fastq files in the right directory
for dir in "$base_dir"/"samples"/*/"fastq"/; do
    
    for fastq_file in "$dir"*.fastq-cut; do
        mv $fastq_file "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/sequences"
    done
done