#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=55G
#SBATCH --time=24:00:00
#SBATCH --job-name=downloading_sra

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_sra/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_sra/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny


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

# create file for raw samples if not already provided
if [ ! -d "$base_dir/raw_files/" ]; then
    mkdir $base_dir/raw_files/
    cd $base_dir/raw_files/
    echo "changed to directory raw_files"
else 
    cd $base_dir/raw_files/
    echo "changed to directory raw_files"
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