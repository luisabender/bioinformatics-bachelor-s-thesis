# Bachelor thesis Luisa

## Goals for PA
### Preprocessing

1. Read in and extract the SRA samples with **sratoolkit**

Bash-Script: preprocessing.sh

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

allignment_kallisto.sh

Downloaded Kallisto via conda. 
Downloaded whole Arabidopsis thaliana transcriptome in TAIR: TAIR10.cdna.all.fa <br>
Built Kallisto index and run kallisto quantification algorithm in SLURM. <br>
Output: abundance.tsv file for each of the 17 samples. <br>
Done MultiQC with Kallisto alignment results.

Repeating the preprocessing step with the rest of the samples, altogether 50 samples.
Build raw and tpm count matrix to load into R. (build_count_matrix.py)

### Differential gene expression analysis in R

DownstreamAnalysis.Rmd

- Loaded count matrices SRA runtable (for metainfo) in R.
- Dropped Samples with very low alignment rate (1.2%-1.6%). 
- Filtered low expressed genes, where the tpm is at least 0.5 for 51% of the genes. 
- Used variance stabilization transformation on the DESeqDataset and created a heatmap of sample-to-sample distances and PCA to search for potential outliers. 

![PCA with vst counts and marked outliers.](plots/PCA.png) 

- DESeq with Interaction term: design: ~treatment + days_post_treatment + treatment:days_post_treatment. 
- Explored the DESeq results and saved the interaction results for time effect in untreated vs treated and the significant genes in each dataset. 
- Created a venn diagram and upset to look for overlaps in significant genes from each day. 

![Venn Diagram with significant genes](plots/upset_sig_genes.png) 

- Loaded Mercator results for Gene level annotations and merged them with the DESeq results. 
- Counted significant genes in each functional category to visualize the counts for each timepoint.  

    | 
:-------------------------:|:-------------------------:
![Barplot with category counts for each day](plots/sig_genes_category.png)  |  ![](lots/sig_genes_category_bubble.png)


- Created a table with the top 30 significant genes for each timepoint and their functional categories. <br>



Next up: Pathway level enrichment