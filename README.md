## Bachelor's thesis
##### Author: Luisa Bender

This Github contains the code, plots and files for the bachelor's thesis on WGCNA.

Topic: Unveiling Commonalities and Divergences in Biotrophic Plant-Microbe Interactions using Gene Co-Expression Networks in Arabidopsis thaliana.

### Internship (Praktische Arbeit)
Code and files, as well as another README of the internship are found here: [praktische-arbeit](praktische-arbeit/).

### Data
The data folder contains accession lists, metadata and an Excel file containing information on the downloaded datasets.

### Count tables
This folder contains all raw gene, abundance and vst transformed count tables.

[Raw gene counts](count_tables/counts_genelev_corrected.csv)

[Abundance (TPM) counts](count_tables/counts_tpm_abundance2.csv)

[VST transformed gene counts](count_tables/expr_data.csv)

[Transposed count table as input for WGCNA](count_tables/expr_data_transposed.csv)


### Supplementary
Following supplementary files and figures are referenced in the thesis.

File S1: [RNA-Seq data on Arabidopsis thaliana](data/arabidopsis_datasets_info.xlsx)
File S2: [Metadata](data/metaData_modified.csv)
Figure S1: [MultiQC-Report](quality_control/multiqc_after_trim.html)
Figure S2: [Sample dendrogram with metadata](plots/dendro_metadata_all.png)

Figures for Gene Ontology (GO) enrichment for specific modules are stored here: [GO plots](plots/WGCNA/Gene%20Ontology/)
