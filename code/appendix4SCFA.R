# association between SCFAs and metabolites module
rm(list = ls())
library(WGCNA)
# load SCFA
scfa <- read.csv("data/histonePTM/SCFA_results.csv", row.names = 1, stringsAsFactors = FALSE)
load("data/derived_data/WGCNA-dataInput.Rdata")
id_sel <- intersect(rownames(scfa), rownames(datExpr))
datTraits <- scfa[id_sel,5:8] # 37 samples left
### Select the soft threshold 
# choose a set of soft-thresholding powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
pickbeta <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = 0.90, verbose = 0)
softpower <- pickbeta$powerEstimate
# Scale-free topology fit index as a function of the soft-thresholding power
plot(pickbeta$fitIndices[,1], -sign(pickbeta$fitIndices[,3])*pickbeta$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",main = paste("Scale independence"));
text(pickbeta$fitIndices[,1], -sign(pickbeta$fitIndices[,3])*pickbeta$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

### Check the scale-free topology 
k <- softConnectivity(datE = datExpr, power = softpower, verbose = 0)
scaleFreePlot(k, main="Check scale free topology\n")


### Network Construction
# calculate the adjacencies, using the soft thresholding power 6
adjacency <- adjacency(datExpr, power = softpower)
# Turn adjacency into topological overlap 
TOM <- TOMsimilarity(adjacency, verbose = 0)
dissTOM <- 1 - TOM
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")
# We like large modules, so we set the minimum module size relatively high
minModuleSize <- 20
# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize, verbose = 0)
# pandoc.table(table(dynamicMods), style = "rmarkdown")
# plot the dendrogram and colors underneath
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Metabolites dendrogram and module colors")

### Merge Module
# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# plot the result
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
# choose a height cut of 0.25, correponding to correlation of 0.75, to merge
MEDissThres <- 0.25
# plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
# call an automatic merging function
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 0)
# merged module colors
mergedColors <- merge$colors
# eigengenes of the nuw merged modules
mergedMEs <- merge$newMEs
# to see what the merging did to our module colors
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# actually, the code above has no effect on the final result
# rename to moduleColors
moduleColors <- mergedColors
# construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLables <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs
# save module eigenvalues
meta.mod <- MEs
rownames(meta.mod) <- rownames(datExpr)
# expand "grey"
idx.grey <- names(datExpr)[moduleColors == "grey"] # 43 unclassified metabolites
meta.mod.expandGrey <- cbind(meta.mod[, -which(names(meta.mod) == "MEgrey")], datExpr[, idx.grey])
# rescale module variable so that they have variance = 1
meta.mod.expandGrey[,1:8] <- apply(meta.mod.expandGrey[, 1:8], 2, function(x) x/sd(x) )
save(meta.mod.expandGrey, file="data/derived_data/WGCNA-meta.mod.expandGrey.Rdata") 

### Identify module associations with the measured phenotypes
# Define numbers of metabolites and samples
nMetas <- ncol(datExpr)
nSamples <- nrow(datExpr)
# # recalculate MEs with color labels (if merged)
# MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs <- orderMEs(MEs0)
MEs <- MEs[id_sel,]
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = 0.7,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Metabolites module-SCFA relationships"))
