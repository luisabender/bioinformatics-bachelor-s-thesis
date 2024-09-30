#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=55G
#SBATCH --time=06:00:00
#SBATCH --job-name=STAR_arabidopsis
#SBATCH --array=1-17%4

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=serial
#SBATCH --partition=serial_long

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"
reads=$(cat $base_dir/sequences/sample.paths | head -n $SLURM_ARRAY_TASK_ID | tail -n1)

# single end alignment
kallisto quant -i $base_dir/genome/Arabidopsis_thaliana_transcript/athaliana.idx -o $base_dir/genome/Arabidopsis_thaliana_transcript --single -l 112 -s 20 $reads
