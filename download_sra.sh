#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <accession_list_file>"
    exit 1
fi

# Read the input file
ACCESSION_FILE="$1"

# Check if the file exists
if [ ! -f "$ACCESSION_FILE" ]; then
    echo "File not found: $ACCESSION_FILE"
    exit 1
fi

while IFS= read -r ACC; do

    # if line is not empty
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
done < "$ACCESSION_FILE"


