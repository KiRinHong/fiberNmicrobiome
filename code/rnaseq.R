# code for analysis related to RNASeq
rm(list = ls())
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(WGCNA)
library(microbiome)
library(phyloseq)
library(dplyr)
library(doParallel)
library(ComplexHeatmap)

################ Data Cleaning Start ##################
### filtering
# read metadata
sampleinfo <- read.csv("data/RNAseq/FIB_WLS_LiverRNAseq_phenotype_data.csv", header = T)
# Reading in the count data
seqdata <- read.delim("data/RNAseq/Ensembl97_STAR_counts.txt", header = T, row.names = "Geneid", stringsAsFactors = FALSE)
# Remove first five columns from seqdata
countdata <- seqdata[,c(-1:-5)]
names(countdata) <- sub("s_", "", sampleinfo[,1])
dim(countdata)
# sum up two different lanes
for (i in seq(1,80,2)) {
  countdata[,i] <- countdata[,i]+countdata[,i+1]
}
countdata <- countdata[,seq(1,80,2)]
# Filtering to remove lowly expressed genes
# Obtain CPMs
myCPM <- cpm(countdata)
# Which values in myCPM are greater than 0.3?
# This produces a logical matrix with TRUEs and FALSEs
thresh <- myCPM > 0.3
# Summary of how many TRUEs there are in each row
table(rowSums(thresh))
# we would like to keep genes that have at least 5 TRUES in each row of thresh
keep <- rowSums(thresh) >= 5
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
# Let's have a look and see whether our threshold of 0.3 does indeed correspond to a count of about 10-15
# We will look at the first sample
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.3)
abline(h=10)
### quality control
dgeObj <- DGEList(counts.keep)
# Library size information is stored in the samples slot
dgeObj$samples
#Library sizes and distribution plots
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
# Add a title to the plot
title("Barplot of library sizes")
# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts
# Get log2 counts per million
logcounts <- cpm(dgeObj,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
# Apply normalisation to DGEList object
dgeObj <- calcNormFactors(dgeObj)
# update the normalisation factors in the DGEList object (their default values are 1)
logcounts1 <- cpm(dgeObj,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts1, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts1),col="blue")
title("Boxplots of logCPMs (normalised)")
# The sample **165-18_S25_L001**(AF) has smallest normalization factor (0.93) and **165J_S40_L002**(Inulin) has largest normalization factor (1.07). 
# Note the plot that shows a few highly expressed genes in the Assorted Fiber sample results in the majority of other genes in the sample having the appearance of being expressed lower.
# Plotting mean difference plots, we could see the composition bias problem. 
# In this step, we use *logcounts*, which have been normalised for library size, but not for composition bias
par(mfrow=c(1,2))
plotMD(logcounts, column = 25)
abline(h=0,col="grey")
plotMD(logcounts,column = 40)
abline(h=0,col="grey")
# As we can see, the composition bias problem has been solved.
sampleinfo <- sampleinfo[seq(1,80,2),]
rownames(sampleinfo) <- colnames(logcounts)
save(dgeObj, sampleinfo, file = "data/derived_data/RNAseq.Rdata")
################ Data Cleaning End ##################

################ WGCNA Start #####################
load("data/derived_data/RNAseq.Rdata")
load("data/derived_data/WGCNA-dataInput.Rdata")
datExpr <- as.data.frame(t(dgeObj$counts))
datConf <- sampleinfo[rownames(datExpr), c("community", "diet", "treatment", "cohort")]
rownames(datExpr) <- unique(sapply(strsplit(rownames(datExpr), "_"),"[",1))
rownames(datConf) <- rownames(datExpr)
# loading clinical trait data
datTraits <- na.omit(datTraits[match(rownames(datExpr), rownames(datTraits)),])
dim(datTraits) # 37 * 4
datExpr <- datExpr[rownames(datTraits),]
datConf <- datConf[rownames(datTraits),]
# As read counts follow a negative binomial distribution, which has a mathematical theory less tractable than that of the normal distribution, RNAseq data was normalised with the voom methodology (Charity W Law et al. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. In: Genome biology 15.2 (Jan. 2014), R29–R29.)
# he voom method estimates the mean-variance of the log-counts and generates a precision weight for each observation. This way, a comparative analysis can be performed with all bioinformatic workflows originally developed for microarray analyses.
datExpr_voom <- voom(t(datExpr))$E
# A large fraction of genes are not differentially expressed between samples. These have to be excluded from WGCNA, as two genes without notable variance in expression between patients will be highly correlated. As a heuristic cutoff, the top 5000 most variant genes have been used in most WGCNA studies. In detail the median absolute devision (MAD) was used as a robust measure of variability.
datExpr <- t(datExpr_voom[order(apply(datExpr_voom,1,mad), decreasing = T)[1:5000],])
rownames(datExpr) <- rownames(datTraits)
sampleTree <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE) # white means low, red means high
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
# Select the soft threshold via visualization
# We considered $R^2=0.90$ as a cut-off to choose soft-thresholding power. As the plot shown, the power 3 is the lowest power for which the scale-free topology fit index reaches 0.90.
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
#   The scale-free topology is confirmed.
k <- softConnectivity(datE = datExpr, power = softpower, verbose = 0)
scaleFreePlot(k, main="Check scale free topology\n")

### Network Construction
# We calculaterd the adjacency by using the soft-thresholding power 6. To minimize the effects of noise and spurious associations, we transformed the adjacency into Topological Overalp Matrix (TOM) and calculated the corresponding dissimilarity. 
# Then we set the minimum module size to 20 and used dynamic tree cut to do module identification. Finally, we plot the dendrogram and colors underneath.
# calculate the adjacencies, using the soft thresholding power 3
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
                    main = "RNAseq dendrogram and module colors")
# The number of each module is listed below: 
pander::pandoc.table(table(dynamicColors), style = "rmarkdown", split.tables = Inf)

### Merge Module
# While the Dynamic Tree Cut may identify modules whose expression profiles are very similar. 
# It may be prudent to merge such modules since their genes are highly co-expressed.   
# To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation. 
# We chose a height cut of 0.2, corresponding to correlation of 0.8, to merge previous dynamic tree. 
# The dendrogram shows that the merged cluster should be consistent with pervious result.
# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# plot the result
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
# choose a height cut of 0.2, correponding to correlation of 0.8, to merge
MEDissThres <- 0.2
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
moduleColors <- mergedColors
# construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLables <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs
# save module eigenvalues
rna.mod <- MEs
rownames(rna.mod) <- rownames(datExpr)
# rescale module variable so that they have variance = 1
rna.mod <- apply(rna.mod, 2, function(x) x/sd(x) )
save(rna.mod, file="data/derived_data/WGCNA-rna.mod.Rdata")
# Identify module associations with the measured phenotypes
# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# # recalculate MEs with color labels (if merged)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
par(mar = c(5.1, 5.1, 4.1, 2.1))
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
               main = paste("RNASeq module-phenotype relationships"))
### Intramodular analysis
# define variable GDATdivBW containing the GDAT_BW column of datTrait
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
module <- "salmon"
column <- match(module, modNames)
moduleGenes <- moduleColors == module
rnaInfo <- data.frame(GENE = colnames(datExpr),
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
  oldNames <- names(rnaInfo)
  rnaInfo <- data.frame(rnaInfo, geneModuleMembership[, modOrder[mod]], 
                        MMPvalue[, modOrder[mod]])
  names(rnaInfo) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                      paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}
rnaOrder <- order(rnaInfo$moduleColor)
rnaInfo <- rnaInfo[rnaOrder,]
rownames(rnaInfo) <- 1:nrow(rnaInfo)
EntrezId <- mapIds(org.Mm.eg.db, keys=as.character(rnaInfo$GENE), keytype = "ENSEMBL", column = "ENTREZID")
Symbol <- mapIds(org.Mm.eg.db, keys=as.character(rnaInfo$GENE), keytype = "ENSEMBL", column = "SYMBOL")
GeneName <- mapIds(org.Mm.eg.db, keys = as.character(rnaInfo$GENE), keytype = "ENSEMBL", column = "GENENAME")
tmp <- data.frame(EntrezId=EntrezId, Symbol=Symbol, GeneName=GeneName)
rnaInfo <- cbind(tmp, rnaInfo)
write.csv(rnaInfo, file = "result/rnaInfo.csv")
################ WGCNA End #####################

################ Association Analysis Start #####################
load("data/derived_data/RNAseq_preprocessing.Rdata")
load("data/derived_data/WGCNA-dataInput.Rdata")
datExpr <- as.data.frame(t(dgeObj$counts))
datConf <- sampleinfo[rownames(datExpr), c("community", "diet", "treatment", "cohort")]
rownames(datExpr) <- unique(sapply(strsplit(rownames(datExpr), "_"),"[",1))
rownames(datConf) <- rownames(datExpr)
# loading clinical trait data
datTraits <- na.omit(datTraits[match(rownames(datExpr), rownames(datTraits)),])
# dimension of datTraits 37 * 4
datExpr <- datExpr[rownames(datTraits),]
datConf <- datConf[rownames(datTraits),]
# As read counts follow a negative binomial distribution, which has a mathematical theory less tractable than that of the normal distribution, RNAseq data was normalised with the voom methodology (Charity W Law et al. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. In: Genome biology 15.2 (Jan. 2014), R29–R29.)
# he voom method estimates the mean-variance of the log-counts and generates a precision weight for each observation. This way, a comparative analysis can be performed with all bioinformatic workflows originally developed for microarray analyses.
datExpr_voom <- voom(t(datExpr))$E
# A large fraction of genes are not differentially expressed between samples. These have to be excluded from WGCNA, as two genes without notable variance in expression between patients will be highly correlated. As a heuristic cutoff, the top 5000 most variant genes have been used in most WGCNA studies. In detail the median absolute devision (MAD) was used as a robust measure of variability.
datExpr <- t(datExpr_voom[order(apply(datExpr_voom,1,mad), decreasing = T)[1:5000],])
rownames(datExpr) <- rownames(datTraits)

tmp <- read.csv("extdata/asv_table_full.csv", row.names = 1)
tmp <- tmp[rownames(datExpr),]
write.csv(tmp, "extdata/asv_table_pairRnaMicro.csv")
tmp <- read.csv("extdata/meta_full.csv", row.names = 1)
tmp <- tmp[rownames(datExpr),]
write.csv(tmp, "extdata/meta_pairRnaMicro.csv")

otu.file <- "extdata/asv_table_pairRnaMicro.csv"
tax.file <- "extdata/taxonomy_table_full.csv"
meta.file <- "extdata/meta_pairRnaMicro.csv"
tree <- read_tree("data/Microbiota16SrRNA/tree.nwk")
pseq <- read_csv2phyloseq(
  otu.file = otu.file,
  taxonomy.file = tax.file,
  metadata.file = meta.file)
pseq@phy_tree <- tree
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(pseq),
               MARGIN = ifelse(taxa_are_rows(pseq), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(pseq),
                    tax_table(pseq))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(pseq) # 1.85
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
pseq <- prune_taxa(keepTaxa, pseq)
# length(get_taxa_unique(pseq, taxonomic.rank = "Genus"))
pseq.genus <- tax_glom(pseq, "Genus", NArm = TRUE)
otu.genus <- abundances(pseq.genus)
meta.genus <- meta(pseq.genus)
tax.genus <- tax_table(pseq.genus)
tree.genus <- pseq.genus@phy_tree
# filter out the taxa with nonzero proportion less than 20%
keepTaxa1 <- rownames(otu.genus)[which(rowMeans(otu.genus!=0)>=0.2)]
pseq.genus<- prune_taxa(keepTaxa1, pseq.genus)
otu.genus <- abundances(pseq.genus)
meta.genus <- meta(pseq.genus)
tax.genus <- tax_table(pseq.genus)
tree.genus <- pseq.genus@phy_tree
otu.genus.comp <- t(otu.genus)
otu.genus.comp[otu.genus.comp == 0] <- 0.5 # use 0.5 to replace 0
otu.genus.comp <- otu.genus.comp/rowSums(otu.genus.comp)
otu.genus.comp.log <- log(otu.genus.comp) # take logrithm
colnames(otu.genus.comp.log) <- tax.genus@.Data[,"Genus"]
# adjust batch effect, diet, community
y <- rna
for (i in 1:ncol(y)) {
  mod <- lm(scale(y[,i]) ~ as.factor(datConf$community) + as.factor(datConf$cohort) + as.factor(datConf$diet))
  y[,i] <- resid(mod)
}

z <- otu.genus.comp.log

dyn.load("code/cvs/cdmm.so")
source("code/cvs/cdmm.R")
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
result_list_gene <- foreach(y=y) %dopar% cdmm_loop(y, z = z)
stopCluster(cl)
list2df <- function(list, dat, colnames = c("metabolites", "genus", "bcv.prob", "stab.prob", "refitted.coef")){
  # convert list to data.frame
  df<- setNames(data.frame(matrix(ncol = 5, nrow = 0)), colnames)
  for (i in 1:length(list)) {
    p <- dim(list[[i]])[1] 
    q <- dim(list[[i]])[2]
    if(p == 0){
      next
    } else{
      tmp <- data.frame(diet = rep(colnames(dat)[i], p))
      colnames(tmp) <- colnames[1]
      tmp <- cbind(tmp, list[[i]])
      df <- rbind(df, tmp)
    }
  }
  return(df)
}
df.gene <- list2df(result_list_gene, dat = y, colnames = c("gene", "genus", "bcv.prob", "stab.prob", "refitted.coef"))
dat.rnainfo <- read.csv("result/rnaInfo.csv", row.names = 1, stringsAsFactors = FALSE)
tmp <- dat.rnainfo[, c("GENE", "GeneName", "moduleColor")]
names(tmp)[1] <- "gene"
res <- merge(df.gene, tmp, by = "gene")
write.csv(res, file = "result/rnaVarSel_full_37samples.csv")
################  Association Analysis End  #####################

################ Mediation Analysis Start #####################
load("data/derived_data/WGCNA-rna.mod.Rdata")
# # batch effect
X <- as.matrix(datConf[, "cohort"])
rownames(X) <- rownames(datConf)
# exposure
S <- as.matrix(datConf[, "community"])
rownames(S) <- rownames(datConf)
# mediator
M <- as.matrix(rna.mod)
diet <- as.matrix(datConf[, "diet"])
rownames(diet) <- rownames(datConf)
# define variable GDATdivBW containing the GDAT_BW column of datTrait
GDATdivBW <- datTraits[, "GDAT_BW", drop = FALSE]
# define variable LiverTriglycerides containing the TG_nmol_gWW column of datTrait
LiverTriglycerides <- datTraits[, "TG_nmol_gWW", drop = FALSE]
# define variable Glucose containing the glocose_norm column of datTrait
Glucose <- datTraits[, "glucose_norm", drop = FALSE]

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
  
  n.draw <- 2000
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
  ##pval.tau.p (MC):		 	p-value for variance component test using Monte Carlo simulation
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
save(res4GDATdivBW, res4LiverTriglycerides,res4Glucose, file = "result/rnaMedRes.Rdata")
################  Mediation Analysis End  #####################