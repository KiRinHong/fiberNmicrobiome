# code for analysis related to Differential Gene Expression

rm(list=ls())
library(gplots)
library(RColorBrewer)
library(edgeR)
library(limma)
library(openxlsx)

### read metadata
fib <- read.csv("FIB_WLS_LiverRNAseq_phenotype_data.csv", header = T)

### read counts table
cnt <- read.delim("Ensembl97_STAR_counts.txt", header = T, row.names = "Geneid")

### delete column from featureCounts
cnt <- cnt[,c(-1:-5)]


## DGEList preparation and annotation
namefib <- fib[,1]
names(cnt)<- namefib

## Create DGList object
x <- DGEList(counts=cnt)


##---------------------------------------------------
## Building the experimental design information  ----
##---------------------------------------------------
group <- as.factor(as.character(fib[,6]))
x$samples$group <- group

lane <- fib[,8]
x$samples$lane <- lane

id <- fib[,9]
x$samples$id <- id

geneid <- rownames(x)

### Check Ids and samples
x$samples


head(x$counts)


##-------------------------------------
## filtering low expressed genes
##----------------------------------

keep <- filterByExpr(x)
x <- x[keep, , keep.lib.sizes=FALSE]

### Addition of Human-readable gene symbols and Entrez identifiers for
##each gene, using the annotation in the org.Mm.eg.db package.

require(org.Mm.eg.db)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

EntrezId <- mapIds(org.Mm.eg.db, keys=rownames(x), keytype="ENSEMBL", column="ENTREZID")
Symbol <- mapIds(org.Mm.eg.db, keys=rownames(x), keytype="ENSEMBL", column="SYMBOL")
GeneName <- mapIds(org.Mm.eg.db, keys = rownames(x), keytype = "ENSEMBL", column = "GENENAME")

x$genes <- data.frame(EntrezId=EntrezId, Symbol=Symbol, GeneName=GeneName)

head(x$genes)

### Some boxplot
nomcol <- x$samples$id
par(mar=c(10,2,1,1))
boxplot(cpm(x, log=T), main="cpms in log scale (not normalized)", names= nomcol, las=2, colours = x$samples$group)


###-----------------------------------
### Normalization
### ----------------------------------

x <- calcNormFactors(x, method = "TMM")
boxplot(cpm(x, log = T), main="cpms in log scale (TMM normalization)", names= nomcol, las=2)

#### PCA's analysis

col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Paired")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-c("black","blue","orange")##si no esta RcolorBrewer
col.lane <- as.character(col.lane)

lcpm <- cpm(x, log = TRUE) 
par(mar=c(2,2,2,2))

prcTMM <- prcomp(lcpm)

### graphics

### TMM normalization
plot(prcTMM$rotation[,1], prcTMM$rotation[,2], col=col.group, main = "PCA (TMM normalization)")
text(prcTMM$rotation[,1],prcTMM$rotation[,2], col= col.group, labels = nomcol)

### PCA 147
prc147 <- prcomp(lcpm[,fib$community=="147"])
plot(prc147$rotation[,1], prc147$rotation[,2], col=col.group, main="PCA 147 community")
text(prc147$rotation[,1],prc147$rotation[,2], col= col.group, labels = nomcol[c(fib$community =="147")])

### PCA 165
prc165 <- prcomp(lcpm[,fib$community=="165"])
plot(prc165$rotation[,1], prc165$rotation[,2], col=col.group, main="PCA 165 community")
text(prc165$rotation[,1],prc165$rotation[,2], col= col.group, labels = nomcol[c(fib$community =="165")])

### PCA Inulin
prcI <- prcomp(lcpm[,fib$diet =="Inulin"])
plot(prcI$rotation[,1], prcI$rotation[,2], col=col.group, main="PCA Inulin diet")
text(prcI$rotation[,1],prcI$rotation[,2], col= col.group, labels = nomcol[c(fib$diet =="Inulin")])

### PCA Pectin
prcP <- prcomp(lcpm[,fib$diet =="Pectin"])
plot(prcP$rotation[,1], prcP$rotation[,2], col=col.group, main="PCA Pectin diet")
text(prcP$rotation[,1],prcP$rotation[,2], col= col.group, labels = nomcol[c(fib$diet =="Pectin")])

### PCA Assorted Fiber
prcAF <- prcomp(lcpm[,fib$diet =="Assorted Fiber"])
plot(prcAF$rotation[,1], prcAF$rotation[,2], col=col.group, main="PCA Assorted Fiber diet")
text(prcAF$rotation[,1],prcAF$rotation[,2], col= col.group, labels = nomcol[c(fib$diet =="Assorted Fiber")])

### PCA Cellulose
prcC <- prcomp(lcpm[,fib$diet =="Cellulose"])
plot(prcC$rotation[,1], prcC$rotation[,2], col=col.group, main="PCA Cellulose diet")
text(prcC$rotation[,1],prcC$rotation[,2], col= col.group, labels = nomcol[c(fib$diet =="Cellulose")])

###-------------------
### building the design
###--------------------

design <- model.matrix(~0+group, data=x$samples)
colnames(design) <- levels(x$samples$group)

###----------------------------------------
### estimating dispersions
###----------------------------------------

x <- estimateDisp(x, design)

x <- estimateGLMTagwiseDisp(x, design)


###-----------------------------------
### Testing for DE genes
###-----------------------------------

fit <- glmQLFit(x, design)

### defining contrast

matrix.contr <- makeContrasts(c147_I.vs.c147_C = c147_I-c147_C,
                              c147_P.vs.c147_C = c147_P-c147_C,
                              c147_AF.vs.c147_C = c147_AF-c147_C,
                              c165_I.vs.c165_C = c165_I-c165_C,
                              c165_P.vs.c165_C = c165_P-c165_C,
                              c165_AF.vs.c165_C = c165_AF-c165_C,
                              c165_I.vs.c147_I = c165_I-c147_I,
                              c165_P.vs.c147_P = c165_P-c147_P,
                              c165_AF.vs.c147_AF = c165_AF-c147_AF,
                              c165_C.vs.c147_C = c165_C-c147_C,
                              levels = design
)

### comparisons

qlf.c165vsc147_I <- glmQLFTest(fit, contrast = matrix.contr[,"c165_I.vs.c147_I"])
qlf.c165vsc147_P <- glmQLFTest(fit, contrast = matrix.contr[,"c165_P.vs.c147_P"])
qlf.c165vsc147_AF <- glmQLFTest(fit, contrast = matrix.contr[, "c165_AF.vs.c147_AF"])
qlf.c165vsc147_C <- glmQLFTest(fit, contrast = matrix.contr[,"c165_C.vs.c147_C"])

## GO Analysis

go.c165vsc147_I <- goana(qlf.c165vsc147_I, species ="Mm", geneid = qlf.c165vsc147_I$genes$EntrezId)
go.c165vsc147_P <- goana(qlf.c165vsc147_P, species ="Mm", geneid = qlf.c165vsc147_P$genes$EntrezId)
go.c165vsc147_C <- goana(qlf.c165vsc147_C, species ="Mm", geneid = qlf.c165vsc147_C$genes$EntrezId)
go.c165vsc147_AF <- goana(qlf.c165vsc147_AF, species ="Mm", geneid = qlf.c165vsc147_AF$genes$EntrezId)

### KEGG Pathways

keg.c165vsc147_I <- kegga(qlf.c165vsc147_I, species ="Mm", geneid = qlf.c165vsc147_I$genes$EntrezId)
keg.c165vsc147_P <- kegga(qlf.c165vsc147_P, species ="Mm", geneid = qlf.c165vsc147_P$genes$EntrezId)
keg.c165vsc147_C <- kegga(qlf.c165vsc147_C, species ="Mm", geneid = qlf.c165vsc147_C$genes$EntrezId)
keg.c165vsc147_AF <- kegga(qlf.c165vsc147_AF, species ="Mm", geneid = qlf.c165vsc147_AF$genes$EntrezId)


### MD Plots
plotMD(qlf.c165vsc147_I, main = "SubB vs SubA: Inulin")
abline(h=c(-1, 1), col="blue")
plotMD(qlf.c165vsc147_P, main = "SubB vs SubA: Pectin")
abline(h=c(-1, 1), col="blue")
plotMD(qlf.c165vsc147_C, main = "SubB vs SubA: Cellulose")
abline(h=c(-1, 1), col="blue")
plotMD(qlf.c165vsc147_AF, main = "SubB vs SubA: A. Fiber")
abline(h=c(-1, 1), col="blue")

### DE genes summary
summary(decideTests(qlf.c165vsc147_I))
summary(decideTests(qlf.c165vsc147_P))
summary(decideTests(qlf.c165vsc147_C))
summary(decideTests(qlf.c165vsc147_AF))

### Write tables

## Contrast 
Inulin <- topTags(qlf.c165vsc147_I, n = Inf)
Pectin <- topTags(qlf.c165vsc147_P, n = Inf)
Cellul <- topTags(qlf.c165vsc147_C, n = Inf)
AFiber <- topTags(qlf.c165vsc147_AF, n = Inf)

## KEGG Pathways
keg.Inulin.UP <- topKEGG(keg.c165vsc147_I, sort="Up", n=100)
keg.Pectin.UP <- topKEGG(keg.c165vsc147_P, sort="Up", n=120)
keg.Cellul.UP <- topKEGG(keg.c165vsc147_C, sort="Up", n=100)
keg.AFiber.UP <- topKEGG(keg.c165vsc147_AF, sort="Up", n=100)

keg.Inulin.DOWN <- topKEGG(keg.c165vsc147_I, sort="Down", n=100)
keg.Pectin.DOWN <- topKEGG(keg.c165vsc147_P, sort="Down", n=100)
keg.Cellul.DOWN <- topKEGG(keg.c165vsc147_C, sort="Down", n=100)
keg.AFiber.DOWN <- topKEGG(keg.c165vsc147_AF, sort="Down", n=100)

## GO Analysis
BP.Inulin.UP <- topGO(go.c165vsc147_I, ont="BP", sort = "Up", n=500)
MF.Inulin.UP <- topGO(go.c165vsc147_I, ont="MF", sort = "Up", n=200)
CC.Inulin.UP <- topGO(go.c165vsc147_I, ont="CC", sort = "Up", n=100)

BP.Pectin.UP <- topGO(go.c165vsc147_P, ont="BP", sort = "Up", n=600)
MF.Pectin.UP <- topGO(go.c165vsc147_P, ont="MF", sort = "Up", n=300)
CC.Pectin.UP <- topGO(go.c165vsc147_P, ont="CC", sort = "Up", n=150)

BP.Cellul.UP <- topGO(go.c165vsc147_C, ont="BP", sort = "Up", n=300)
MF.Cellul.UP <- topGO(go.c165vsc147_C, ont="MF", sort = "Up", n=100)
CC.Cellul.UP <- topGO(go.c165vsc147_C, ont="CC", sort = "Up", n=100)

BP.AFiber.UP <- topGO(go.c165vsc147_AF, ont="BP", sort = "UP", n=700)
MF.AFiber.UP <- topGO(go.c165vsc147_AF, ont="MF", sort = "Up", n=200)
CC.AFiber.UP <- topGO(go.c165vsc147_AF, ont="CC", sort = "Up", n=150)

BP.Inulin.DOWN <- topGO(go.c165vsc147_I, ont="BP", sort = "Down", n=550)
MF.Inulin.DOWN <- topGO(go.c165vsc147_I, ont="MF", sort = "Down", n=200)
CC.Inulin.DOWN <- topGO(go.c165vsc147_I, ont="CC", sort = "Down", n=250)

BP.Pectin.DOWN <- topGO(go.c165vsc147_P, ont="BP", sort = "Down", n=550)
MF.Pectin.DOWN <- topGO(go.c165vsc147_P, ont="MF", sort = "Down", n=250)
CC.Pectin.DOWN <- topGO(go.c165vsc147_P, ont="CC", sort = "Down", n=250)

BP.Cellul.DOWN <- topGO(go.c165vsc147_C, ont="BP", sort = "Down", n=550)
MF.Cellul.DOWN <- topGO(go.c165vsc147_C, ont="MF", sort = "Down", n=250)
CC.Cellul.DOWN <- topGO(go.c165vsc147_C, ont="CC", sort = "Down", n=250)

BP.AFiber.DOWN <- topGO(go.c165vsc147_AF, ont="BP", sort = "Down", n=550)
MF.AFiber.DOWN <- topGO(go.c165vsc147_AF, ont="MF", sort = "Down", n=250)
CC.AFiber.DOWN <- topGO(go.c165vsc147_AF, ont="CC", sort = "Down", n=250)




