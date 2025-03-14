#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=55G
#SBATCH --time=24:00:00
#SBATCH --job-name=adapter_trimming

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_adapt/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_adapt/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#### adapter trimming ####

set -e
set -u
set -o pipefail

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

echo "Beginning adapter trimming"

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
outputdir_trimmed="$base_dir/trimmed_files/added_dataset"

# create file for trimmed samples if not already provided
if [ ! -d "$outputdir_trimmed" ]; then
    mkdir $outputdir_trimmed
fi


for sample in $base_dir/"raw_files/Stegmann_RNAseq/RNAseq_WT"/*; do
    
    : '
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
    fi'

    #specific adapter trimming for sample with bad quality scores
    sample_name=$(basename "$sample")
    echo "Processing single-end sample: $sample_name"
    fastp \
        -i "$sample" \
        -o "$outputdir_trimmed/${sample_name}_trimmed.fastq" \
        --adapter_sequence AGATCGGAAGAGC \
        --trim_poly_x \
        --cut_right\
         --cut_right_window_size 4 \
         --cut_right_mean_quality 20 \
         --length_required 50 \
         --dedup \
         --poly_x_min_len 4 \
         --qualified_quality_phred 30
        
done

#### fastqc after trimming ####
output_fastqc_trimmed="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_trimmed/added_dataset"

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
