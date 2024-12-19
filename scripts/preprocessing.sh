#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=55G
#SBATCH --time=48:00:00
#SBATCH --job-name=preprocessing

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_preprocessing/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_preprocessing/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=serial
#SBATCH --partition=serial_long

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


#### adapter trimming ####
echo "Beginning adapter trimming"

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
outputdir_trimmed="$base_dir/trimmed_files"

# create file for trimmed samples if not already provided
if [ ! -d "$outputdir_trimmed" ]; then
    mkdir $outputdir_trimmed
fi


for sample_dir in $base_dir/"raw_files"/*; do
    # Check if it's a directory
    if [ -d "$sample_dir/fastq" ]; then
        # Get the sample name (folder name)
        sample_name=$(basename "$sample_dir")
        
        # Check for paired-end or single-end files
        paired_end_1=$(find "$sample_dir/fastq" -name "*_1.fastq*")
        paired_end_2=$(find "$sample_dir/fastq" -name "*_2.fastq*")
        single_end=$(find "$sample_dir/fastq" -name "*.fastq*" ! -name "*_1.fastq*" ! -name "*_2.fastq*")
        
        # Paired-end trimming
        if [ -n "$paired_end_1" ] && [ -n "$paired_end_2" ]; then
            echo "Processing paired-end sample: $sample_name"
            fastp \
                -i "$paired_end_1" \
                -I "$paired_end_2" \
                -o "$outputdir_trimmed/${sample_name}_R1_trimmed.fastq" \
                -O "$outputdir_trimmed/${sample_name}_R2_trimmed.fastq" 
        
        # Single-end trimming
        elif [ -n "$single_end" ]; then
            echo "Processing single-end sample: $sample_name"
            fastp \
                -i "$single_end" \
                -o "$outputdir_trimmed/${sample_name}_trimmed.fastq" 
        
        # No suitable files found
        else
            echo "No fastq files found for sample: $sample_name"
        fi
    fi
done

#### fastqc after trimming ####
output_fastqc_trimmed="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_trimmed"

if [ ! -d "$output_fastqc_trimmed" ]; then
    mkdir "$output_fastqc_trimmed"
fi


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
