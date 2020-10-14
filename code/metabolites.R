# code for analysis related to metabolites
rm(list = ls())
library(openxlsx)
library(WGCNA)
library(ggplot2)
library(pander)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(circlize)
################ Data Cleaning Start ##################
# load phenotypes
pheno <- read.table("data/Phenotype/FatMassBW_LiverTriglycerides_glucose.txt", sep = "\t", 
                    skip = 1, na.strings = c("NA", "#VALUE!"), fill = TRUE, stringsAsFactors = FALSE)
# BW_sac(g)	GDAT_sac(g)	GDAT/BW(%)	TG(nmol/g wet weight)
colnames(pheno) <- c("AnimalID", "WLS", "Diet", "Treatment", "GDAT_BW", "TG_nmol_gWW", "glucose_norm")
pheno <- na.omit(pheno)
# load mouse information
miceInfo <- read.xlsx("data/Metabolomic_data/UWIM-01-19VW CDT sent to the client 06 06 2019.xlsx", sheet = 3, rows = 1:12, cols = 12:61)

# load metabolite
meta <- read.xlsx("data/Metabolomic_data/UWIM-01-19VW CDT sent to the client 06 06 2019.xlsx", sheet = 3, startRow = 13)
metaColNames <- names(meta)
metaColNames[14:61] <- names(miceInfo)[-1]
names(meta) <- metaColNames
# select the necessary columns
meta.val <- meta[, -which(names(meta) %in% 
                            c("PATHWAY_SORTORDER", "SUPER_PATHWAY", "SUB_PATHWAY", "COMP_ID", 
                              "PLATFORM", "CHEMICAL_ID", "RI", "MASS", "PUBCHEM", "CAS", "KEGG", "Group.HMDB"))]
# find the intersect
crossNames <- intersect(names(miceInfo)[-1], pheno$AnimalID)
pheno <- pheno[match(crossNames, pheno$AnimalID), ]
rownames(pheno) <- pheno$AnimalID
pheno$AnimalID <- NULL
miceInfo <- miceInfo[, c("CLIENT_IDENTIFIER", crossNames)]
rownames(miceInfo) <- miceInfo$CLIENT_IDENTIFIER
miceInfo$CLIENT_IDENTIFIER <- NULL

meta.val <- meta.val[, c("BIOCHEMICAL", crossNames)]
rownames(meta.val) <- meta.val$BIOCHEMICAL
meta.val$BIOCHEMICAL <- NULL
meta.val <- as.data.frame(t(meta.val))
# do inverse transformation and test normality
meta.val.INT <- as.data.frame(apply(meta.val, 2, function(x) qnorm((rank(x) - 0.5)/length(x)) ))
rm.idx <- NULL
for (j in 1:ncol(meta.val.INT)) {
  tmp <- try(shapiro.test(meta.val.INT[, j])$p.value)
  if(inherits(tmp, 'try-error')){
    rm.idx <- c(rm.idx, j)
  }else if(tmp < 0.05){
    rm.idx <- c(rm.idx, j)
  }
}
meta.val.filt.INT <- meta.val.INT[, -rm.idx] ## remove 62 non-normal metabolites (774 --> 712)

# compare the histogram 
hist(meta.val$`1-(1-enyl-oleoyl)-GPE (P-18:1)*`, xlab = "1-(1-enyl-oleoyl)-GPE (P-18:1)*", main = "Original data")
hist(meta.val.INT$`1-(1-enyl-oleoyl)-GPE (P-18:1)*`, xlab = "1-(1-enyl-oleoyl)-GPE (P-18:1)*", main = "Transformed data")
rownames(meta.val.filt.INT) <- rownames(meta.val.INT)
names(meta.val.filt.INT) <- names(meta.val.INT)[-rm.idx]

# check for metabolites and samples with too many missing values
gsg <- goodSamplesGenes(meta.val.filt.INT, verbose = 0)
gsg$allOK
# if returns TRUE, all metabolites have passed the cuts.
# if not, we remove the offending metabolites and samples from the data
# if(!gsg$allOK){
#   # Optionally, print the metabolome and sample names that were removed:
#   if (sum(!gsg$goodGenes) > 0)
#     printFlush(paste("Removing genes:", paste(names(meta.val.filt.INT)[!gsg$goodGenes], collapse = ", ")))
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(meta.val.filt.INT)[!gsg$goodSamples], collapse = ", ")))
#   # remove the offending metabolites and samples from the data
#   meta.val.filt.INT <- meta.val.filt.INT[gsg$goodSamples, gsg$goodGenes]
# }
datExpr <- meta.val.filt.INT
datTraits <- pheno[, c("GDAT_BW", "TG_nmol_gWW", "glucose_norm")]
save(datExpr, datTraits, file = "data/derived_data/WGCNA-dataInput.Rdata")  # 712 metabolites
################ Data Cleaning End ##################

################ WGCNA Start ##################
load("data/derived_data/WGCNA-dataInput.Rdata")
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
               main = paste("Metabolites module-phenotype relationships"))

### Intramodular analysis
# define variable GDATdivBW containing the GDAT_BWBW column of datTrait
GDATdivBW <- datTraits[, "GDAT_BW", drop = FALSE]
# define variable LiverTriglycerides containing the TG_nmol_gWW column of datTrait
LiverTriglycerides <- datTraits[, "TG_nmol_gWW", drop = FALSE]
# define varable Glucose containing the glucose_norm column of datTrait
Glucose <- datTraits[,"glucose_norm", drop = FALSE]
# names of the modules
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

# GDATdivBW
geneTraitSignificance <- as.data.frame(cor(datExpr, GDATdivBW, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(GDATdivBW), sep = "")
names(GSPvalue) <- paste("p.GS.", names(GDATdivBW), sep = "")
# LiverTriglycerides
geneTraitSignificance_LT <- as.data.frame(cor(datExpr, LiverTriglycerides, use = "p"))
GSPvalue_LT <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_LT), nSamples))
names(geneTraitSignificance_LT) <- paste("GS.", names(LiverTriglycerides), sep = "")
names(GSPvalue_LT) <- paste("p.GS.", names(LiverTriglycerides), sep = "")
# Glucose
geneTraitSignificance_Glu <- as.data.frame(cor(datExpr, Glucose, use = "p"))
GSPvalue_Glu <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Glu), nSamples))
names(geneTraitSignificance_Glu) <- paste("GS.", names(LiverTriglycerides), sep = "")
names(GSPvalue_Glu) <- paste("p.GS.", names(Glucose), sep = "")
# intramodular analysis: identifying genes with high GS and MM
module <- "turquoise"
column <- match(module, modNames)
moduleGenes <- moduleColors == module

metaInfo <- data.frame(BIOCHEMICAL = names(datExpr),
                       SUPER_PATHWAY = meta[match(names(datExpr), meta$BIOCHEMICAL), "SUPER_PATHWAY"],
                       SUB_PATHWAY = meta[match(names(datExpr), meta$BIOCHEMICAL), "SUB_PATHWAY"],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue,
                       geneTraitSignificance_LT,
                       GSPvalue_LT,
                       geneTraitSignificance_Glu,
                       GSPvalue_Glu)
# order modules by their significance for GDATdivBW
modOrder <- order(-abs(cor(MEs, GDATdivBW, use = "p")))
# add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(metaInfo)
  metaInfo <- data.frame(metaInfo, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]])
  names(metaInfo) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                       paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}
metaOrder <- order(metaInfo$moduleColor)
metaInfo <- metaInfo[metaOrder,]
rownames(metaInfo) <- 1:nrow(metaInfo)
write.csv(metaInfo, file = "result/metaInfo.csv")
save(GDATdivBW, LiverTriglycerides, Glucose,
     miceInfo, pheno, meta.mod.expandGrey, file = "data/derived_data/MedTest-meta-dataInput.Rdata")
################ WGCNA End ##################


################ Mediation Analyis Start ##################
load("data/derived_data/MedTest-meta-dataInput.Rdata")
# # batch effect
X <- t(as.matrix(miceInfo["COHORT", ]))
# exposure
S <- as.matrix(c(rep("147", 23), rep("165", 22)))
rownames(S) <- rownames(X)
# mediator
M <- as.matrix(meta.mod.expandGrey)
diet <- t(as.matrix(miceInfo["DIET", ]))
medTest <- function(X, S, M, Y, diet, target = "Cellulose"){
  p <- ncol(M)
  alpha <- matrix(NA, p, 3)
  colnames(alpha) <- c("Estimate", "P.value.aysm", "P.value.boot")
  beta = gamma = IE <- alpha
  alpha.se <- numeric(length = p)
  gamma.se = beta.se = IE.se <- alpha.se
  idx.diet <- which(diet[,1] == target)
  x <- X[idx.diet,]
  s <- S[idx.diet,]
  y <- Y[idx.diet,]
  fit.y.x <- lm(y ~ x)
  y <- residuals(fit.y.x)
  for (i in 1:p) {
    m <- M[idx.diet, i]
    fit.m.x <- lm(m ~ x)
    m <- residuals(fit.m.x)
    fit.m <- lm(m ~ s)
    fit.y <- lm(y ~ s + m)
    alpha[i, 1:2] <- c(coef(fit.m)[2], summary(fit.m)$coefficients[2, 4]) # alpha_S
    alpha.se[i] <- summary(fit.m)$coefficients[2, 2]
    gamma[i, 1:2] <- c(coef(fit.y)[2], summary(fit.y)$coefficients[2, 4]) # beta_S
    gamma.se[i] <- summary(fit.y)$coefficients[2, 2]
    beta[i, 1:2] <- c(coef(fit.y)[3], summary(fit.y)$coefficients[3, 4]) # beta_M
    beta.se[i] <- summary(fit.y)$coefficients[3, 2]
    
    IE[i, 1] <- alpha[i, 1] * beta[i, 1]
    IE.se[i] <- sqrt((alpha[i, 1] * beta.se[i])^2 + 
                       (alpha.se[i] * beta[i, 1])^2)
    IE[i, 2] <- 2 * pnorm(abs(IE[i, 1]/IE.se[i]), lower.tail = FALSE)
  }
  
  IE.total <- numeric(length = 2)
  DE = TE = L2ME = IE.total.boot <- IE.total
  IE.total[1] <- sum(IE[, 1])
  IE.total.se <- sqrt(sum(IE.se^2))
  IE.total[2] <- 2 * pnorm(abs(IE.total[1]/IE.total.se), lower.tail = FALSE)
  
  ## total effect model
  fit <- lm(y ~ s)
  TE.se <- summary(fit)$coefficients[2, 2]
  # total effect
  TE[1] <- fit$coefficients[2]
  TE[2] <- 2 * pnorm(abs(TE[1]/TE.se), lower.tail = FALSE)
  # direct effect
  DE[1] <- TE[1] - IE.total[1]
  DE.se <- sqrt(IE.total.se^2 + TE.se^2)
  DE[2] <- 2 * pnorm(abs(DE[1]/DE.se), lower.tail = FALSE)
  
  n.draw <- 1000
  bcov.p <- Matrix::bdiag(diag(alpha.se^2), diag(beta.se^2))
  theta.p.mc <- mvtnorm::rmvnorm(n.draw, mean = c(alpha[,1], beta[,1]), sigma = as.matrix(bcov.p))	## sample from the distribution of regression coefficients
  alpha.p.mc <- theta.p.mc[, 1:p]
  beta.p.mc <- theta.p.mc[, p + 1:p]
  ies.p.mc <- beta.p.mc * alpha.p.mc
  ies.p.mc.c <- t(t(ies.p.mc) - apply(ies.p.mc, 2, mean))
  for (j in 1:p){
    IE[j,3] <- 2*mean(ies.p.mc.c[,j] > abs(IE[j,1]))
  }
  total.ide.mc <- apply(ies.p.mc.c, 1, sum)
  total.ide.mc.p <- 2*mean(total.ide.mc > abs(IE.total[1]))
  IE.total.boot[1] <- IE.total[1] 
  IE.total.boot[2] <- total.ide.mc.p
  vie.p.mc <- apply(ies.p.mc.c^2, 1, sum)
  pval.v.p.mc <- mean(vie.p.mc > sum(IE[,1]^2))
  L2ME <- c(sum(IE[,1]^2), pval.v.p.mc)
  
  res <- cbind(TE, DE, IE.total, IE.total.boot, L2ME)
  colnames(res) <- c("TotalEffect", "DirectEffect", "MarginalME.asym", "MarginalME.boot", "L2normME")
  rownames(res) <- c("Estimate", "P-value")
  res <- round(res, 4)
  rownames(IE) = rownames(alpha) = rownames(beta) = rownames(gamma) <- colnames(M)
  
  re <- list(result = res, IE = IE, alpha = alpha, beta = beta, gamma = gamma)
  return(re)
}

res_print <- function(x){
  idx <- which(x$IE[,3] < 0.1)
  if(length(idx) == 0) {
    print("**None**")
  } else {
    if(length(idx) == 1){
      res.sel <- t(x$IE[idx, ])
    }else {
      res.sel <- x$IE[idx, ]
    }
    rownames(res.sel) <- names(idx)
    pandoc.table(res.sel, style = "rmarkdown", split.tables = Inf)
  }
  pandoc.table(x$result, style = "rmarkdown", split.tables = Inf)
}

# GDATdivBW
Y <- as.matrix(scale(GDATdivBW))
res4Cellulose <- medTest(X, S, M, Y, diet, target = "Cellulose")
res_print(res4Cellulose)
res4Inulin <- medTest(X, S, M, Y, diet, target = "Inulin")
res_print(res4Inulin)
res4Pectin <- medTest(X, S, M, Y, diet, target = "Pectin")
res_print(res4Pectin)
res4GDATdivBW <- list(Cellulose = res4Cellulose, Inulin = res4Inulin, Pectin = res4Pectin)

# LiverTriglycerides
Y <- as.matrix(scale(LiverTriglycerides))
res4Cellulose <- medTest(X, S, M, Y, diet, target = "Cellulose")
res_print(res4Cellulose)
res4Inulin <- medTest(X, S, M, Y, diet, target = "Inulin")
res_print(res4Inulin)
res4Pectin <- medTest(X, S, M, Y, diet, target = "Pectin")
res_print(res4Pectin)
res4LiverTriglycerides <- list(Cellulose = res4Cellulose, Inulin = res4Inulin, Pectin = res4Pectin)

# Glucose
Y <- as.matrix(scale(Glucose))
res4Cellulose <- medTest(X, S, M, Y, diet, target = "Cellulose")
res_print(res4Cellulose)
res4Inulin <- medTest(X, S, M, Y, diet, target = "Inulin")
res_print(res4Inulin)
res4Pectin <- medTest(X, S, M, Y, diet, target = "Pectin")
res_print(res4Pectin)
res4Glucose <- list(Cellulose = res4Cellulose, Inulin = res4Inulin, Pectin = res4Pectin)
save(res4GDATdivBW, res4LiverTriglycerides, res4Glucose, file = "result/metaMedRes.Rdata")
################ Mediation Analyis End ##################