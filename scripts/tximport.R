options(repos = c(CRAN = "https://cran.r-project.org"))

if (!require("tximport", quietly = TRUE))
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("tximport")

library(tximport)
library(rhdf5)

base_dir <- '/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/'
samples <- read.table('/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/bachelor-thesis-luisa/data/accession_lists/all.txt', header = TRUE)
files <- file.path(base_dir, "kallisto_paired", samples$sample_ID, "abundance.h5")
names(files) <- samples$sample_ID

# loading tx2gene.tsv file, created with Araport11_GTF_genes_transposons.current.gtf
# awk -F'\t' '$3 == "mRNA" {print $9}' Araport11_GTF_genes_transposons.current.gtf | awk -F';' '{gsub(/ /, "", $0); print $2, $1}' | awk -F'"' '{print $2, $4}' > tx2gene.tsv
tx2gene <- readLines(file.path(base_dir, "genome/Arabidopsis_thaliana/tx2gene.tsv"))
tx2gene <- do.call(rbind, strsplit(tx2gene, "\\s+")) # Splits on whitespace
tx2gene <- data.frame(TXNAME = tx2gene[, 2], GENEID = tx2gene[, 1], stringsAsFactors = FALSE)
tx2gene <- tx2gene[,c("GENEID", "TXNAME")]

# tximport
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write.csv(txi.kallisto$counts, "gene_counts_w_added_dataset.csv")
write.csv(txi.kallisto$abundance, "abundance_counts_w_added_dataset.csv")