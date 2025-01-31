library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

base_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/luisa/bachelor-thesis-luisa/"
expr_data <- read.csv(file.path(base_dir,"count_tables/vst_expression_data.csv"), row.names = 1)

# transpose it for WGCNA
input_mat <- t(expr_data)

print("Calling blockwiseModules...")

net = blockwiseModules(input_mat, power = 5,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "TOM",
verbose = 3)


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf(file.path(base_dir,"wgcna_results/dendro_pow5.pdf"))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = file.path(base_dir,"/wgcna_results/networkConstruction-auto-pow5.RData"))


print("saved data successfully.")