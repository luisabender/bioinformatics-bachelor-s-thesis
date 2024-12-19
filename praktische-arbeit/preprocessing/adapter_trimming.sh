#!/bin/bash
# use trimmomatic for basic adapter removal
base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"


for dir in "$base_dir"/*/"fastq"/; do
    echo "Processing directory: $dir"

    for fastq_file in "$dir"*.fastq; do
        if [ -f "$fastq_file" ]; then
            
            #running trimmomatic
            echo "Running trimmomatic on $fastq_file"
            cd $base_dir/Trimmomatic-0.39
            java -jar trimmomatic-0.39.jar SE -phred33 $fastq_file $fastq_file-trim ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
            
            #running cutadapt to remove polyA tails
            if [ -f "$fastq_file-trim" ]; then
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate cutadapt
            echo "Running cutadapt on $fastq_file-trim"
            cutadapt --poly-a -o $fastq_file-cut $fastq_file-trim
            conda deactivate
            
            fi
        fi
    done
done
