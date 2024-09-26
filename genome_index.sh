#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=55G
#SBATCH --time=01:00:00
#SBATCH --job-name=STAR_arabidopsis

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm/%x.%j.%a.err
 
#This is where you allocate some of cluster specific parameters.
#SBATCH --clusters=serial
#SBATCH --partition=serial_long

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env

base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"



# generate Genome index
STAR --runMode genomeGenerate \
     --genomeDir $base_dir/genome/Arabidospis_thaliana/  \
     --genomeFastaFiles $base_dir/genome/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
     --sjdbGTFfile $base_dir/genome/Arabidopsis_thaliana/Araport11_GTF_genes_transposons.current.gtf.gz \
     --runThreadN $SLURM_CPUS_PER_TASK
