# Bachelor thesis Luisa

## Goals for PA
1. Read in and extract the SRA samples via sratoolkit \n
Dataset from Paper: A primary cell wall cellulose-dependent defense mechanism against vascular pathogens revealed by time-resolved dual transcriptomics 
https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01100-6#Sec13
GEO repository: GSE168919, Accession via PRJNA714597
mRNA-profiles of Arabidopsis thaliana both untreated and infected with pathogen Fusarium oxysporum. 
Downloaded 17 Samples from Day 0 and Day 6 post treatment (untreated and infected) to compare between later. 

2. FASTQC
Used the package fastqc from conda and stored the results in folder fastqc_results.

3. Adapter trimming via fastp
Performed fastp on the .fastq files. 
Note: couldn't remove any adapters, need to adapt the fastp commands.

4. FASTQC again
Note: fastqc after adapter trimming didn't show any big changes on quality score.

TODO:

5. Alignment to an Arabidopsis thaliana genome (TAIR10) with STAR
Download STAR package via bioconda. Download whole Arabidopsis thaliana genome. 
Perform STAR Alignment and use HTSeq-count to count the reads per gene for each sample.


6. Differential gene expression analysis via DESeq2
Load the result files in R and run DESeq.
Enrichment analysis with logfoldchange to identify genes which are most upregulated.
Compare between different samples, infected and control, Day 0 and Day 6.

7. Visualizations