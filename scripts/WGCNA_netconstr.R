# WGCNA adjacency with different soft thresholds
library(WGCNA)
# load expression data
base_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/bachelor-thesis-luisa"
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

input_mat <- read.csv(file.path(base_dir,"count_tables/expr_data_transposed.csv"), row.names = 1)

powers <- c(4)

# change the adjacancy into signed hybrid
for (power in powers){

print(paste("Starting with power ", power))
print("Running adjacency...")
softPower <- power
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
ME.dissimilarity <- 1-cor(MElist$eigengenes, use="complete") 

#Clustering eigengenes 
METree <- hclust(as.dist(ME.dissimilarity), method = "average") 
par(mar = c(0,4,2,0)) 
par(cex = 0.6);

#pdf(file.path(base_dir, "wgcna_results", "new", paste0("METree_sp", power, "_TOMsigned.pdf")))

plot(METree)
abline(h=.25, col = "red") # a height of .25 corresponds to correlation of .75
#dev.off()

print("Merging...")
merge <- mergeCloseModules(input_mat, ModuleColors, cutHeight = 0.25)

# merged module colors
mergedColors <- merge$colors
# eigengenes of the new merged modules
mergedMEs <- merge$newMEs

# plot dendrogram of original and merged module colors
#pdf(file.path(base_dir, "wgcna_results", "new",paste0("dendro_merged_sp",power,"_TOMsigned.pdf")))
dendro_merged <- plotDendroAndColors(geneTree, 
                    cbind(ModuleColors, mergedColors),
                    c("original Module", "merged Module"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
#dev.off()


save(mergedMEs, ModuleColors, mergedColors, geneTree, file = file.path(base_dir, "wgcna_results", "new",paste0("networkConstruction-pow",power,"_signed_mer.RData")))


print(paste("Finished succesfully with power", power))
}
