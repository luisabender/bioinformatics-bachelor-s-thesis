# Bachelor thesis Luisa

## Goals for PA
1. Read in and extract the SRA samples with **sratoolkit**

Dataset from Paper: A primary cell wall cellulose-dependent defense mechanism against vascular pathogens revealed by time-resolved dual transcriptomics. <br>
https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01100-6#Sec13 <br>
GEO repository: GSE168919, Accession via PRJNA714597<br>
mRNA-profiles of Arabidopsis thaliana both untreated and infected with pathogen Fusarium oxysporum. <br>
Downloaded 17 Samples from Day 0 and Day 6 post treatment (untreated and infected) to compare between later.

2. **FASTQC**

Used the package fastqc from conda and stored the results in folder fastqc_results.

3. **Adapter trimming** with trimmomatic and cutadapt

Performed trimmomatic on the .fastq files for basic adapter trimming. Then run cutadapt on the results to eliminate the polyA adapters.<br>

4. **FASTQC** again and **MultiQC report**

FastQC showed good results on adapter content.
MultiQC report on trimmed reads.

5. **Alignment** to an Arabidopsis thaliana genome (TAIR10) with **Kallisto**

Downloaded Kallisto via conda. 
Downloaded whole Arabidopsis thaliana transcriptome in TAIR: TAIR10.cdna.all.fa <br>
Built Kallisto index and run kallisto quantification algorithm in SLURM. <br>
Output: abundance.tsv file for each of the 17 samples. <br>
Done MultiQC with Kallisto alignment results.

Repeating the preprocessing step with the rest of the samples, altogether 50 samples.
Build raw and tpm count matrix to load into R.

6. **Differential gene expression analysis** with DESeq2

Load the result files in R and make a count and metadata table.<br>
Analyze the count table, make graphs and slides about the metadata too. <br>
PCA with the tpm count matrix to detect potential outliers, check with fastqc result. <br>
Variance stabilization transformation with raw count matrix. <br>
DESeq2 with vst result. <br>


Functional enrichment analysis with logfoldchange to identify genes which are most upregulated.<br>


7. **Visualizations**