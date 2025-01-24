# WGCNA adjacency with differen soft thresholds
library(WGCNA)
# load expression data
options(stringsAsFactors = FALSE)
expr_data <- read.csv("~/Studium/bachelor-thesis-luisa/count_tables/vst_expression_data.csv", row.names = 1)

# transpose it for WGCNA
input_mat <- t(expr_data)

softPower <- 8
adjacency <- adjacency(input_mat, power = softPower)

# Topological Overlap Matrix
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM


### Hierarchical Clustering Analysis

# creating the dendrogram
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
# plotting
sizeGrWindow(12,9)

# start with minimum module size 30 as recommended by the authors of WGCNA
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
ModuleColors <- labels2colors(Modules) # assigns each module number a color

table(Modules)
table(ModuleColors) # returns the counts for each color (aka the number of genes within each module)

# plots the gene dendrogram with the module colors
pdf("~/Studium/bachelor-thesis-luisa/plots/dendro_new.pdf")
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

# eigengene dissimilarity
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

# construct a cluster tree
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic

pdf("~/Studium/bachelor-thesis-luisa/plots/METree_new.pdf")
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
dev.off()

merge <- mergeCloseModules(input_mat, ModuleColors, cutHeight = 0.25)

# merged module colors
mergedColors <- merge$colors
# eigengenes of the new merged modules
mergedMEs <- merge$newMEs

# plot dendrogram of original and merged module colors
pdf("~/Studium/bachelor-thesis-luisa/plots/dendro_merged_new.pdf")
dendro_merged <- plotDendroAndColors(geneTree, 
                    cbind(ModuleColors, mergedColors),
                    c("original Module", "merged Module"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()