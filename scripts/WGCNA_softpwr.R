# WGCNA adjacency with different soft thresholds
library(WGCNA)
# load expression data
base_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/bachelor-thesis-luisa/"
options(stringsAsFactors = FALSE)
expr_data <- read.csv(file.path(base_dir,"count_tables/vst_expression_data.csv"), row.names = 1)

# transpose it for WGCNA
input_mat <- t(expr_data)

print("Running adjacency...")
softPower <- 5
adjacency <- adjacency(input_mat, power = softPower)

# Topological Overlap Matrix
print("Calculating TOM...")
TOM.dissimilarity <- 1-TOMsimilarity(adjacency)


### Hierarchical Clustering Analysis

# creating the dendrogram
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

# start with minimum module size 30 as recommended by the authors of WGCNA
print("Creating Modules...")
minModuleSize = 30
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
ModuleColors <- labels2colors(Modules) # assigns each module number a color

modules.table <- write.csv(table(Modules),file.path(base_dir,"/wgcna_results/modules_table_new.csv"))
moduleColors.table <- write.csv(table(ModuleColors), file.path(base_dir,"/wgcna_results/moduleColors_table_new.csv")) # returns the counts for each color (aka the number of genes within each module)

# plots the gene dendrogram with the module colors
pdf(file.path(base_dir,"/wgcna_results/dendro_new.pdf"))
dendro_plot <- plotDendroAndColors(geneTree, 
                    ModuleColors,
                    "Module",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Module Eigengene Identification
MElist <- moduleEigengenes(input_mat, colors = ModuleColors)
MEs <- MElist$eigengenes
write.csv(MEs,file.path(base_dir,"/wgcna_results/MEs_new.csv"))

# eigengene dissimilarity
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

# construct a cluster tree
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic

pdf(file.path(base_dir,"/wgcna_results/METree_new.pdf"))
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
dev.off()

print("Merging...")
merge <- mergeCloseModules(input_mat, ModuleColors, cutHeight = 0.25)

# merged module colors
mergedColors <- merge$colors
# eigengenes of the new merged modules
mergedMEs <- merge$newMEs
write.csv(mergedMEs, file.path(base_dir,"wgcna_results/mergedMEs_new.csv"))

# plot dendrogram of original and merged module colors
pdf(file.path(base_dir,"/wgcna_results/dendro_merged_new.pdf"))
dendro_merged <- plotDendroAndColors(geneTree, 
                    cbind(ModuleColors, mergedColors),
                    c("original Module", "merged Module"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

pdf(file.path(base_dir, "/wgcna_results/TOMplot_new.pdf"))
TOMplot(TOM.dissimilarity, geneTree, mergedColors)
dev.off()

print("Finished succesfully.")