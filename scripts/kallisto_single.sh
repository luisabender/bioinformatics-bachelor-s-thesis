#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=55G
#SBATCH --job-name=kallisto_single
#SBATCH --array=1-141%4

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_kallisto_single/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_kallisto_single/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny


source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"

reads=$(cat $base_dir/trimmed_files/sample_single.paths | head -n $SLURM_ARRAY_TASK_ID | tail -n1)

# output file for each sample
outputfile=$base_dir/kallisto_single/$(basename $reads .fastq)

mkdir -p $( dirname $outputfile)

# single end alignment
kallisto quant -i $base_dir/genome/Arabidopsis_thaliana_transcript/athaliana.idx -o $outputfile --single -l 116.41 -s 25.45 $reads
