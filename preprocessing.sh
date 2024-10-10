#!/bin/bash
#unless you are running a parallelized task and know what you're doing, ntasks remain 1, modify the next lines to your task.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=55G
#SBATCH --time=04:00:00
#SBATCH --job-name=preprocessing

 
#make sure you have the folder ~/logs/slurm/ for the next lines.
#SBATCH -o /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm/%x.%j.%a.out 
#SBATCH -e /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/slurm/%x.%j.%a.err
 
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

cd $base_dir/samples/
echo "changed to directory samples"

while IFS= read -r ACC; do

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

###### fastqc #######

#perform fastqc on every extracted .fastq file

output_fastqc="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_results"
source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
# mkdir "$output_fastqc"
echo "Beginning FastQC"

for dir in "$base_dir"/"samples"/*/"fastq"/; do
    echo "Processing directory: $dir"

    # find all fastq files
    for fastq_file in "$dir"*.fastq; do
        if [ -f "$fastq_file" ]; then
            echo "Running FastQC on $fastq_file"
            fastqc "$fastq_file" --outdir="$output_fastqc"
        fi
    done
done

#### adapter trimming ####
echo "Beginning adapter trimming"

for dir in "$base_dir"/"samples"/*/"fastq"/; do
    echo "Processing directory: $dir"

    for fastq_file in "$dir"*.fastq; do
        if [ -f "$fastq_file" ]; then
            
            #running trimmomatic
            echo "Running trimmomatic on $fastq_file"
            cd $base_dir/tools/Trimmomatic-0.39
            java -jar trimmomatic-0.39.jar SE -phred33 $fastq_file $fastq_file-trim ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
            
            #running cutadapt to remove polyA tails
            if [ -f "$fastq_file-trim" ]; then
            source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
            echo "Running cutadapt on $fastq_file-trim"
            cutadapt --poly-a -o $fastq_file-cut $fastq_file-trim
            conda deactivate
            
            fi
        fi
    done
done

#### fastqc after trimming ####
source /dss/dsshome1/0A/ge58rom2/miniconda3/bin/activate bio_env
output_fastqc_trimmed="/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/fastqc_results_trimmed"

#mkdir "$output_fastqc"

for dir in "$base_dir"/"samples"/*/"fastq"/; do
    echo "Processing directory: $dir"

    
    for fastq_file in "$dir"*.fastq-cut; do
        if [ -f "$fastq_file" ]; then
            echo "Running FastQC on $fastq_file"
            fastqc "$fastq_file" --outdir="$output_fastqc_trimmed"
        fi
    done
done

# move trimmed fastq files
for dir in "$base_dir"/"samples"/*/"fastq"/; do
    
    for fastq_file in "$dir"*.fastq-cut; do
        mv $fastq_file "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/sequences"
    done
done

echo "moved fastq files into folder sequences"

# remove -cut from files
for file in *-cut; do mv "$file" "${file%-cut}"; done
