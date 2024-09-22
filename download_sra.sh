#!/bin/bash

# SRA Accessions List as input file: SRA_Acc_List.txt

# check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <accession_list_file>"
    exit 1
fi

# read in the accession list
ACCESSION_LIST="$1"

# Check if file exists
if [ ! -f "$ACCESSION_LIST" ]; then
    echo "File not found: $ACCESSION_LIST"
    exit 1
fi

while IFS= read -r ACC; do

    if [ -n "$ACC" ]; then 
        echo "Downloading $ACC"
    
        # download the SRA file
        prefetch "$ACC"
    
        if [ $? -eq 0 ]; then
            echo "$ACC downloaded successfully."
            echo "Extracting $ACC"
        
            # Extract fastq file from .sra file
            cd "$ACC"
            fasterq-dump --outdir fastq "$ACC" 
            
            if [ $? -eq 0 ]; then
                echo "$ACC extracted successfully."
                cd ".."
            else
                echo "Error extracting $ACC."
            fi
        else
            echo "Error downloading $ACC."
        fi
    fi
done < "$ACCESSION_LIST"

###### fastqc #######

#perform fastqc on every extracted .fastq file

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"
output_fastqc="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_results"

mkdir "$output_fastqc"

for dir in "$base_dir"/*/"fastq"/; do
    echo "Processing directory: $dir"

    # find all fastq files
    for fastq_file in "$dir"*.fastq; do
        if [ -f "$fastq_file" ]; then
            echo "Running FastQC on $fastq_file"
            fastqc "$fastq_file" --outdir="$output_fastqc"
        fi
    done
done

### adapter trimming with fastp ####


for dir in "$base_dir"/*/"fastq"/; do
    echo "Processing directory: $dir"

    for fastq_file in "$dir"*.fastq; do
        if [ -f "$fastq_file" ]; then
            echo "Running fastp on $fastq_file"
            fastp -i "$fastq_file" -o "$fastq_file"_trimmed
        fi
    done
done


