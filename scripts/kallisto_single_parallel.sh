#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=50G
#SBATCH --job-name=kallisto_single
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --get-user-env
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_kallisto_single/%x.%j.%a.out
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm_kallisto_single/%x.%j.%a.err
#SBATCH --mail-user=luisa.bender@tum.de
#SBATCH --mail-type=ALL
 

source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
module load parallel/20220522
 
base_dir="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa"
sample_file="$base_dir/trimmed_files/added_dataset/samples_added.paths"
 
# Check if sample file exists
if [[ ! -f $sample_file ]]; then
    echo "Error: Sample file $sample_file not found!"
    exit 1
fi
 
output_base="$base_dir/kallisto_paired"
mkdir -p "$output_base"
 
# Number of parallel tasks to run (adjust this to control parallelism)
parallel_tasks=8
threads_per_task=$((SLURM_CPUS_PER_TASK/parallel_tasks))
if [[ $threads_per_task -lt 1 ]]; then
    echo "Error: Threads per task cannot be less than 1. Check cpus-per-task and parallel_tasks."
    exit 1
fi
 
# run kallisto quant in parallel
export kallisto_idx="$base_dir/genome/Arabidopsis_thaliana_transcript/athaliana.idx"
export output_base
export threads_per_task
 
echo "Running kallisto in parallel for $(wc -l < "$sample_file") samples..."
 
cat "$sample_file" | parallel --jobs $parallel_tasks --linebuffer \
    'reads={}; \
     output_dir="$output_base/$(basename "$reads" .fastq)"; \
     output_dir="${output_dir%.fq}"; \
     mkdir -p "$output_dir"; \
     kallisto quant -i "$kallisto_idx" -o "$output_dir" --single -l 200 -s 30 -t $threads_per_task "$reads"'
 
echo "All kallisto jobs completed!"



