#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=55G
#SBATCH --time=24:00:00
#SBATCH --job-name=kallisto_paired
#SBATCH --array=1-142%4

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_kallisto_paired/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_kallisto_paired/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"


# get sample paths 
reads_R1=$(cat $base_dir/trimmed_files/sample_R1.paths | head -n $SLURM_ARRAY_TASK_ID | tail -n1)
reads_R2=$(cat $base_dir/trimmed_files/sample_R2.paths | head -n $SLURM_ARRAY_TASK_ID | tail -n1)

# Ensure both R1 and R2 exist
if [ -z "$reads_R1" ] || [ -z "$reads_R2" ]; then
    echo "Error: Missing paired-end reads for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# output file for each sample
outputfile=$base_dir/kallisto_output_paired/$(basename $reads_R1 _R1.fastq)

mkdir -p $( dirname $outputfile)

# Paired-end alignment
kallisto quant -i $base_dir/genome/Arabidopsis_thaliana_transcript/athaliana.idx \
    -o $outputfile $reads_R1 $reads_R2
