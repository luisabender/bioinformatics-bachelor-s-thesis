#!/bin/bash
#### fastqc after trimming ####

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

output_fastqc_trimmed="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_trimmed"

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

#mkdir "$output_fastqc"

for fastq_file in "$base_dir"/"trimmed_files"/"failed_multiqc"/*; do
        if [ -f "$fastq_file" ]; then
            echo "Running FastQC on $fastq_file"
            fastqc "$fastq_file" --outdir="$output_fastqc_trimmed"
        fi
done
