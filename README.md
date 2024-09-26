# Bachelor thesis Luisa

## Goals for PA
1. Read in and extract the SRA samples via **sratoolkit**

Dataset from Paper: A primary cell wall cellulose-dependent defense mechanism against vascular pathogens revealed by time-resolved dual transcriptomics. <br>
https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01100-6#Sec13 <br>
GEO repository: GSE168919, Accession via PRJNA714597<br>
mRNA-profiles of Arabidopsis thaliana both untreated and infected with pathogen Fusarium oxysporum. <br>
Downloaded 17 Samples from Day 0 and Day 6 post treatment (untreated and infected) to compare between later.

2. **FASTQC**

Used the package fastqc from conda and stored the results in folder fastqc_results.

3. **Adapter trimming** via trimmomatic and cutadapt

Performed trimmomatic on the .fastq files for basic adapter trimming. Then run cutadapt on the results to eliminate the polyA adapters.<br>

4. **FASTQC** again and **MultiQC report**

FastQC showed good results on adapter content.
MultiQC report is generated and uploaded.

TODO:

5. **Alignment** to an Arabidopsis thaliana genome (TAIR10) with **STAR**

Downloaded STAR package via conda. 
Downloaded whole Arabidopsis thaliana genome in ENSEMBL: TAIR10.dna.toplevel.fa.qz <br>
Downloaded Arabidopsis GTF annotation file (TAIR11). <br>
Perform STAR genome indexing and single-end alignment in SLURM. <br>
Use HTSeq-count or featureCount to count the reads per gene for each sample.


6. **Differential gene expression analysis** via DESeq2

Load the result files in R and run DESeq.<br>
Enrichment analysis with logfoldchange to identify genes which are most upregulated.<br>
Compare between different samples, infected and control, Day 0 and Day 6.<br>

7. **Visualizations**