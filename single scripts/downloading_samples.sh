#!/bin/bash

set -e
set -u
set -o pipefail


base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

# SRA Accessions List as input file: SRA_Acc_List.txt

# check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <accession_list_file>"
    exit 1
fi

# read in the accession list
ACCESSION_LIST="$1"
echo $ACCESSION_LIST

# Check if file exists
if [ ! -f "$ACCESSION_LIST" ]; then
    echo "File not found: $ACCESSION_LIST"
    exit 1
fi

if [ ! -d "$base_dir/raw_files/" ]; then
    mkdir $base_dir/raw_files/
    cd $base_dir/raw_files/
    echo "created samples file and changed to directory"
else 
    cd $base_dir/raw_files/
    echo "changed to directory samples"
fi


while IFS= read -r ACC || [ -n "$ACC" ]; do

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
                #cd ".."
            fi
        else
            echo "Error downloading $ACC."
        fi
    fi
done < "$ACCESSION_LIST"

cd $base_dir