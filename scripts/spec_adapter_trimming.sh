#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=55G
#SBATCH --time=01:00:00
#SBATCH --job-name=fastqc_adapt
#SBATCH --mail-user=luisa.bender@tum.de
#SBATCH --mail-type=ALL
 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_adapt/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_adapt/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=serial
#SBATCH --partition=serial_long

#### adapter trimming ####

set -e
set -u
set -o pipefail

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

echo "Beginning adapter trimming"

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
outputdir_trimmed="$base_dir/trimmed_files"

# create file for trimmed samples if not already provided
if [ ! -d "$outputdir_trimmed" ]; then
    mkdir $outputdir_trimmed
fi


for sample_dir in $base_dir/"raw_files"/"failed_multiqc"/*; do
    # Check if its a directory
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
    
            fastp --cut_front --qualified_quality_phred 30 --cut_mean_quality 30 --trim_poly_g --n_base_limit 3 \
                -i "$paired_end_1" \
                -I "$paired_end_2" \
                -o "$outputdir_trimmed/${sample_name}_R1_trimmed.fastq" \
                -O "$outputdir_trimmed/${sample_name}_R2_trimmed.fastq" \
                --html report.html --json report.json

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
: 'output_fastqc_trimmed="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_trimmed"
outputdir_trimmed="$base_dir/trimmed_files/failed_multiqc"
if [ ! -d "$output_fastqc_trimmed" ]; then
    mkdir "$output_fastqc_trimmed"
fi

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
for sample in "$outputdir_trimmed"/*; do
    if [ -f "$sample" ]; then
        echo "Running FastQC on $sample"
        fastqc "$sample" --outdir="$output_fastqc_trimmed"
    fi
done
'