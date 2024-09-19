#!/bin/bash

# List of SRA accession numbers
ACCESSIONS=("SRR13823992")


for ACC in "${ACCESSIONS[@]}"; do
    echo "Downloading $ACC"
    
    # download the SRA file
    prefetch "$ACC"
    
    if [ $? -eq 0 ]; then
        echo "$ACC downloaded successfully."
        echo "Extracting $ACC"
    
        # Extract fastq file from .sra file
        cd "$ACC"
        fasterq-dump "$ACC" --outdir fastq
        
        if [ $? -eq 0 ]; then
            echo "$ACC extracted successfully."
        else
            echo "Error extracting $ACC."
        fi
    else
        echo "Error downloading $ACC."
    fi
done

