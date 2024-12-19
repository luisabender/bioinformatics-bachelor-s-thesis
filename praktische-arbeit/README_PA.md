## Praktische Arbeit

**Following Dataset was used:**

Paper: A primary cell wall cellulose-dependent defense mechanism against vascular pathogens revealed by time-resolved dual transcriptomics. <br>
https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01100-6#Sec13 <br>
GEO repository: GSE168919, Accession via PRJNA714597<br>

Samples: mRNA-profiles of Arabidopsis thaliana roots both untreated and infected with pathogen Fusarium oxysporum over a time course of six days. <br>

### Preprocessing

Run preprocessing.sh

1. Read in and extract the SRA samples with **sratoolkit**

    Downloaded 48 Samples from Day 0 and Day 6 post treatment (untreated and infected) to compare between later.

2. Quality control with FASTQC

3. Adapter trimming with trimmomatic and cutadapt

    Performed trimmomatic on the .fastq files for basic adapter trimming. Then run cutadapt on the results for removing specifically the polyA adapters.<br>

4. FASTQC again and MultiQC report

    [MultiQC-Report](files/multiqc_report_1.html)

### Alignment

Run index_kallisto.sh and alignment_kallisto.sh

1. Build Kallisto Index with Arabidopsis thaliana transcriptome (from TAIR)
2. Run Kallisto quantification algorithm with the build index 

    [Kallisto MultiQC-Report](files/kallisto_multiqc_report_1.html)

3. Build raw and tpm count matrix with build_count_matrix.py and export them

### Differential gene expression analysis

Run DownstreamAnalysis.Rmd

- Loaded count matrices SRA runtable (for metainfo) in R. Metainfo contains treatment and days post treatment.
- Dropped Samples with very low alignment rate (1.2%-1.6%). 
- Filtered low expressed genes, where the tpm is at least 0.5 for 51% of the genes. 
- Used variance stabilization transformation on the DESeqDataset and created a heatmap of sample-to-sample distances and PCA to search for potential outliers. 

**DESeq2**

- DESeq with Interaction term: design: ~treatment + days_post_treatment + treatment:days_post_treatment. 
- Explored the DESeq results and saved the interaction results for time effect in untreated vs treated and the significant genes in each dataset. 
- Created a venn diagram and upset to look for overlaps in significant genes from each day. 


**Mercator**

- Mercator protein annotation with A. thaliana transcript file from TAIR
- Loaded Mercator results (mapping file and fasta file) for Gene level annotations and merged them with the DESeq results. 
- Counted significant genes in each functional category to visualize the counts for each timepoint.  


- Created a table with the top 30 significant genes for each timepoint and their functional categories. 

GENE_ID     |CATEGORY    |FUNCTION|  P_VALUE|  	LOG_FOLD_CHANGE|	days_post_treatment
:----------:|:----------:|:------:|:-------:|:----------------:|:-----------------------:
AT1G16030.1	|Protein homeostasis|	Protein homeostasis.protein quality control.cytosolic Hsp70 chaperone system.Hsp70 chaperone activities.molecular chaperone *(Hsp70-1/2/3/4/5)|	0|	-11.119594	|3
AT2G22170.1|	Lipid metabolism|	Lipid metabolism.lipid trafficking.endoplasmic reticulum-plasma membrane lipid transfer.lipid trafficking protein *(PLAT)|	0	|5.416741|	3
<br>

- all different expressed annotated genes for each timepoint can be found here:
[annotated DEG](files/)

**MapMan (Mercator4 BIN enrichment analysis)**

- pathway enrichment analysis with Mercator result mapping file
- genes of interest: significant genes for each timepoint, background genes: all genes from annotated fasta file
- result: significant pathways with genes in MapMan category for each timepoint
- focusing on cone category, e.g. cell wall organisation

### Results
All result plots can be found here: [plots](R/plots/).


