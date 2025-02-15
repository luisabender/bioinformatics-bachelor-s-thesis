#!/bin/bash
#### fastqc after trimming ####

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

output_fastqc_trimmed="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_trimmed/added_dataset"
outputdir_trimmed="$base_dir/trimmed_files/added_dataset"

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

#mkdir "$output_fastqc"

for sample in "$outputdir_trimmed"/*; do
    if [ -f "$sample" ]; then
        echo "Running FastQC on $sample"
        fastqc "$sample" --outdir="$output_fastqc_trimmed"
    fi
done

#### multiqc ####
echo "Running MultiQC"
cd $output_fastqc_trimmed
multiqc .
cd $base_dir