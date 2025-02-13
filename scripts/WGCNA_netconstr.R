# WGCNA adjacency with different soft thresholds
library(WGCNA)
# load expression data
base_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/bachelor-thesis-luisa/"
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

expr_data <- read.csv(file.path(base_dir,"count_tables/vst_expression_data.csv"), row.names = 1)

# transpose it for WGCNA
input_mat <- t(expr_data)

# change the adjacancy into signed hybrid
print("Running adjacency...")
softPower <- 4
adjacency <- adjacency(input_mat, power = softPower, type = "signed hybrid")

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

# Module Eigengene Identification
MElist <- moduleEigengenes(input_mat, colors = ModuleColors)
MEs <- MElist$eigengenes

# eigengene dissimilarity
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

# construct a cluster tree
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic

pdf(file.path(base_dir,"/wgcna_results/METree_sp4_signed.pdf"))
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
dev.off()

print("Merging...")
merge <- mergeCloseModules(input_mat, ModuleColors, cutHeight = 0.25)

# merged module colors
mergedColors <- merge$colors
# eigengenes of the new merged modules
mergedMEs <- merge$newMEs

# plot dendrogram of original and merged module colors
pdf(file.path(base_dir,"/wgcna_results/dendro_merged_sp4_signed.pdf"))
dendro_merged <- plotDendroAndColors(geneTree, 
                    cbind(ModuleColors, mergedColors),
                    c("original Module", "merged Module"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()


save(mergedMEs, mergedColors, geneTree, file = file.path(base_dir,"/wgcna_results/networkConstruction-pow4_signed.RData"))


print("Finished succesfully")