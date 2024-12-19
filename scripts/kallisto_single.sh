#!/bin/bash

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"
#id=0

sample_dir=$base_dir/trimmed_files/sample_single.paths

for read in $(cat $sample_dir); do
    #id=$id+1
    #reads=$(cat $base_dir/trimmed_files/sample_single.paths | head -n $id | tail -n1)


    outputfile=$base_dir/kallisto_single/$(basename $read .fastq)
    mkdir -p $( dirname $outputfile)
    

    error_file="${outputfile/.error.log}"
    log_file="${outputfile/.output.log}"

    # submit job to slurm
    sbatch -M cm4 -p cm4_tiny -q cm4_tiny -N 1 -c 112 -n 1 -J kallisto -e "$error_file" -o "$log_file" --mem 55G \
        --wrap="kallisto quant -i \"$base_dir/genome/Arabidopsis_thaliana_transcript/athaliana.idx\" -o \"$outputfile\" --single -l 116.41 -s 25.45 \"$read\""

done




