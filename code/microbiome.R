# code for analysis related to microbiome
rm(list = ls())
library(microbiome)
library(phyloseq)
library(dplyr)
library(tidyr)
library(reshape2)
library(GUniFrac)
library(vegan)
library(abind)
library(data.table)
library(doParallel)
library(pander)

################ Data Cleaning Start ##################
######## Generate asv_table_subset.csv, taxonomy_table_subset.csv, meta_subset.csv
######## focus on diet cellulose/inulin/pectin, phenotype GDATdivBW, LiverTriglycerides, Glucose
######## The value of phenotype are not scaled
# load phenotypes
pheno <- read.table("data/Phenotype/FatMassBW_LiverTriglycerides_glucose.txt", sep = "\t", 
                    skip = 1, na.strings = c("NA", "#VALUE!"), fill = TRUE, stringsAsFactors = FALSE)
# GDAT/BW(%)	TG(nmol/g wet weight)
colnames(pheno) <- c("AnimalID", "WLS", "Diet", "Treatment", "GDATdivBW", "LiverTriglycerides", "Glucose")
pheno <- na.omit(pheno)
# load mouse information
miceInfo <- read.xlsx("data/Metabolomic_data/UWIM-01-19VW CDT sent to the client 06 06 2019.xlsx", sheet = 3, rows = 1:12, cols = 12:61)
# find the intersect
crossNames <- intersect(names(miceInfo)[-1], pheno$AnimalID)
pheno <- pheno[match(crossNames, pheno$AnimalID), ]
rownames(pheno) <- pheno$AnimalID
pheno$AnimalID <- NULL
miceInfo <- miceInfo[, c("CLIENT_IDENTIFIER", crossNames)]
rownames(miceInfo) <- miceInfo$CLIENT_IDENTIFIER
miceInfo$CLIENT_IDENTIFIER <- NULL
pheno$Cohort <- as.vector(t(miceInfo["COHORT",]))
pheno$WLS <- as.character(pheno$WLS)
pheno <- pheno[, c("WLS", "Diet", "Treatment", "Cohort", "GDATdivBW", "LiverTriglycerides", "Glucose")]
idx.rm <- which(pheno$Diet == "Assorted Fiber")
pheno <- pheno[-idx.rm,]
rm(miceInfo)
write.csv(pheno, "extdata/meta_subset.csv")

asv <- read.table("data/Microbiota16SrRNA/feature-table.txt", skip = 1, stringsAsFactors = FALSE)
colnames(asv) <- asv[1,]
asv <- asv[-1,]
rownames(asv) <- asv$OTU_ID
asv$OTU_ID <- NULL
asv <- asv[,rownames(pheno)]
asv[] <- lapply(asv, as.numeric)
asv <- as.data.frame(t(asv))
write.csv(asv, "extdata/asv_table_subset.csv")

dat <- read.csv("data/Microbiota16SrRNA/taxonomy.csv", sep = "\t", stringsAsFactors = FALSE, row.names = 1)
dat <- dat[match(colnames(asv), rownames(dat)),]
tax <- dat %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") %>% 
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species)

rm(dat)
tax[] <- lapply(tax, function(x) {
  tmp <- stringr::str_replace(x, "[:lower:]__", "")
  ifelse(tmp == "", NA, tmp)
})
tax <- as.data.frame(tax) # class %in% 'tbl_df', 'tbl' cannot use ifelse
# For the taxa that are unclassified at the lower level, their identities at higher levels were used
unclassify <- c("k.", "p.", "c.", "o.", "f.", "g.")
for (i in 2:7) {
  tax[, i] <- ifelse(is.na(tax[, i]) & !stringr::str_detect(tax[, i-1], "\\."),
                     paste0(unclassify[i-1], tax[, i-1]), tax[, i])
}
for (i in 2:7) {
  tax[, i] <- ifelse(is.na(tax[, i]), tax[, i-1], tax[, i])
}
write.csv(tax, "extdata/taxonomy_table_subset.csv")

######## Generate asv_table_full.csv, taxonomy_table_full.csv, meta_full.csv
######## focus on diet cellulose/inulin/pectin/assorted fiber, phenotype GDATdivBW, LiverTriglycerides, Glucose
######## The value of phenotype are not scaled
# load phenotypes
pheno <- read.table("data/Phenotype/FatMassBW_LiverTriglycerides_glucose.txt", sep = "\t", 
                    skip = 1, na.strings = c("NA", "#VALUE!"), fill = TRUE, stringsAsFactors = FALSE)
# BW_sac(g)	GDAT_sac(g)	GDAT/BW(%)	TG(nmol/g wet weight)
colnames(pheno) <- c("AnimalID", "WLS", "Diet", "Treatment", "GDATdivBW", "LiverTriglycerides", "Glucose")
pheno <- na.omit(pheno)
# load mouse information
miceInfo <- read.xlsx("data/Metabolomic_data/UWIM-01-19VW CDT sent to the client 06 06 2019.xlsx", sheet = 3, rows = 1:12, cols = 12:61)
# find the intersect
crossNames <- intersect(names(miceInfo)[-1], pheno$AnimalID)
pheno <- pheno[match(crossNames, pheno$AnimalID), ]
rownames(pheno) <- pheno$AnimalID
pheno$AnimalID <- NULL
miceInfo <- miceInfo[, c("CLIENT_IDENTIFIER", crossNames)]
rownames(miceInfo) <- miceInfo$CLIENT_IDENTIFIER
miceInfo$CLIENT_IDENTIFIER <- NULL
pheno$Cohort <- as.vector(t(miceInfo["COHORT",]))
pheno$WLS <- as.character(pheno$WLS)
pheno <- pheno[, c("WLS", "Diet", "Treatment", "Cohort", "GDATdivBW", "LiverTriglycerides", "Glucose")]
rm(miceInfo)
write.csv(pheno, "extdata/meta_full.csv")

asv <- read.table("data/Microbiota16SrRNA/feature-table.txt", skip = 1, stringsAsFactors = FALSE)
colnames(asv) <- asv[1,]
asv <- asv[-1,]
rownames(asv) <- asv$OTU_ID
asv$OTU_ID <- NULL
asv <- asv[,rownames(pheno)]
asv[] <- lapply(asv, as.numeric)
asv <- as.data.frame(t(asv))
write.csv(asv, "extdata/asv_table_full.csv")

dat <- read.csv("data/Microbiota16SrRNA/taxonomy.csv", sep = "\t", stringsAsFactors = FALSE, row.names = 1)
dat <- dat[match(colnames(asv), rownames(dat)),]
tax <- dat %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") %>% 
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species)

rm(dat)
tax[] <- lapply(tax, function(x) {
  tmp <- stringr::str_replace(x, "[:lower:]__", "")
  ifelse(tmp == "", NA, tmp)
})
tax <- as.data.frame(tax) # class %in% 'tbl_df', 'tbl' cannot use ifelse
# For the taxa that are unclassified at the lower level, their identities at higher levels were used
unclassify <- c("k.", "p.", "c.", "o.", "f.", "g.")
for (i in 2:7) {
  tax[, i] <- ifelse(is.na(tax[, i]) & !stringr::str_detect(tax[, i-1], "\\."),
                     paste0(unclassify[i-1], tax[, i-1]), tax[, i])
}
for (i in 2:7) {
  tax[, i] <- ifelse(is.na(tax[, i]), tax[, i-1], tax[, i])
}
write.csv(tax, "extdata/taxonomy_table_full.csv")


######## Generate asv_table_scfa.csv, taxonomy_table_scfa.csv, meta_scfa.csv
######## The value of phenotype are not scaled
# load phenotypes
scfa <- read.csv("data/histonePTM/SCFA_results.csv", row.names = 1, stringsAsFactors = FALSE)
asv <- read.table("data/Microbiota16SrRNA/feature-table.txt", skip = 1, stringsAsFactors = FALSE)
colnames(asv) <- asv[1,]
asv <- asv[-1,]
rownames(asv) <- asv$OTU_ID
asv$OTU_ID <- NULL
id.keep <- intersect(rownames(scfa), colnames(asv))
asv <- asv[,id.keep]
asv[] <- lapply(asv, as.numeric)
asv <- as.data.frame(t(asv))
write.csv(asv, "extdata/asv_table_scfa.csv")
scfa <- scfa[id.keep,]
write.csv(scfa, "extdata/meta_scfa.csv")

dat <- read.csv("data/Microbiota16SrRNA/taxonomy.csv", sep = "\t", stringsAsFactors = FALSE, row.names = 1)
dat <- dat[match(colnames(asv), rownames(dat)),]
tax <- dat %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") %>% 
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
rm(dat)
tax[] <- lapply(tax, function(x) {
  tmp <- stringr::str_replace(x, "[:lower:]__", "")
  ifelse(tmp == "", NA, tmp)
})
tax <- as.data.frame(tax) # class %in% 'tbl_df', 'tbl' cannot use ifelse
# For the taxa that are unclassified at the lower level, their identities at higher levels were used
unclassify <- c("k.", "p.", "c.", "o.", "f.", "g.")
for (i in 2:7) {
  tax[, i] <- ifelse(is.na(tax[, i]) & !stringr::str_detect(tax[, i-1], "\\."),
                     paste0(unclassify[i-1], tax[, i-1]), tax[, i])
}
for (i in 2:7) {
  tax[, i] <- ifelse(is.na(tax[, i]), tax[, i-1], tax[, i])
}
write.csv(tax, "extdata/taxonomy_table_scfa.csv")
################ Data Cleaning End ##################

################ Association Anslysis Start ##################
### Variable Selection: SCFA ~ Microbiome
otu.file <- "extdata/asv_table_scfa.csv"
tax.file <- "extdata/taxonomy_table_scfa.csv"
meta.file <- "extdata/meta_scfa.csv"
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
# Define phyla to filter
filterPhyla = c("Synergistetes")
# Filter entries with unidentified Phylum.
pseq = subset_taxa(pseq, !Phylum %in% filterPhyla)
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(pseq, "Phylum"))
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(pseq) # 2
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
pseq <- prune_taxa(keepTaxa, pseq)
pseq.genus <- tax_glom(pseq, "Genus", NArm = TRUE)
otu.genus <- abundances(pseq.genus)
# filter out the taxa with nonzero proportion less than 20%
keepTaxa1 <- rownames(otu.genus)[which(rowMeans(otu.genus!=0)>=0.2)]
pseq.genus<- prune_taxa(keepTaxa1, pseq.genus)
# pseq.genus
pseq.genus.rel <- microbiome::transform(pseq.genus, "compositional")
otu.genus <- abundances(pseq.genus)
otu.genus.comp <- abundances(pseq.genus.rel)
meta.genus <- meta(pseq.genus)
tax.genus <- tax_table(pseq.genus)
tree.genus <- pseq.genus@phy_tree
otu.genus.rff <- t(Rarefy(t(otu.genus))$otu.tab.rff)
# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq)
otu.comp <- abundances(pseq.rel)
meta <- meta(pseq)
tree <- pseq@phy_tree
otu.rff <- t(Rarefy(t(otu))$otu.tab.rff)
# calculate the UniFraces
BC <- as.matrix(vegdist(t(otu.rff), method = "bray"))
JAC <- as.matrix(vegdist(t(otu.rff), 'jaccard', binary = TRUE))
unifracs <- GUniFrac(t(otu.rff), tree, alpha = c(0, 0.5, 1))$unifracs # do not need relative abundance
unifracs <- abind(unifracs,BC,JAC)
dimnames(unifracs)[[3]][6:7] <- c("d_BC", "d_JAC")
res4trait <- PermanovaG(unifracs[, , c("d_BC", "d_JAC", "d_1", "d_UW")] ~ Acetate, data = meta, permutations = 2000)
rownames(res4trait$p.tab) <- "p.value"
colnames(res4trait$p.tab) <- c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac", "Omnibus Test")
pandoc.table(round(as.matrix(res4trait$p.tab), 4), style = "rmarkdown", split.tables = Inf, 
             caption = "PERMANOVA significance test for Acetate (Rarefied ASV)")
res4trait <- PermanovaG(unifracs[, , c("d_BC", "d_JAC", "d_1", "d_UW")] ~ Propionate, data = meta, permutations = 2000)
rownames(res4trait$p.tab) <- "p.value"
colnames(res4trait$p.tab) <- c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac", "Omnibus Test")
pandoc.table(round(as.matrix(res4trait$p.tab), 4), style = "rmarkdown", split.tables = Inf, 
             caption = "PERMANOVA significance test for Propionate (Rarefied ASV)")
res4trait <- PermanovaG(unifracs[, , c("d_BC", "d_JAC", "d_1", "d_UW")] ~ Butyrate, data = meta, permutations = 2000)
rownames(res4trait$p.tab) <- "p.value"
colnames(res4trait$p.tab) <- c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac", "Omnibus Test")
pandoc.table(round(as.matrix(res4trait$p.tab), 4), style = "rmarkdown", split.tables = Inf, 
             caption = "PERMANOVA significance test for Butyrate (Rarefied ASV)")
y <- meta.genus[, 5:7]
otu.genus.comp <- t(otu.genus)
otu.genus.comp[otu.genus.comp == 0] <- 0.5 # use 0.5 to replace 0
otu.genus.comp <- otu.genus.comp/rowSums(otu.genus.comp)
otu.genus.comp.log <- log(otu.genus.comp) # take logrithm
colnames(otu.genus.comp.log) <- tax.genus@.Data[,"Genus"]
z <- otu.genus.comp.log
dyn.load("code/cvs/cdmm.so")
source("code/cvs/cdmm.R")
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
result_list_scfa <- foreach(y=y) %dopar% cdmm_loop(y, z = z)
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
df.scfa <- list2df(result_list_scfa, dat = y, colnames = c("SCFA", "genus", "bcv.prob", "stab.prob", "refitted.coef"))


### Variable Selection: Metabolites ~ Microbiome
### data preprocess: filter taxa
otu.file <- "extdata/asv_table_full.csv"
tax.file <- "extdata/taxonomy_table_full.csv"
meta.file <- "extdata/meta_full.csv"
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
# Define phyla to filter
filterPhyla = c("Synergistetes")
# Filter entries with unidentified Phylum.
pseq = subset_taxa(pseq, !Phylum %in% filterPhyla)
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(pseq, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(pseq),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(pseq) # 2.25
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
pseq<- prune_taxa(keepTaxa, pseq)
# length(get_taxa_unique(pseq, taxonomic.rank = "Genus"))
pseq.genus = tax_glom(pseq, "Genus", NArm = TRUE)
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
load("data/derived_data/WGCNA-dataInput.Rdata")
dat.metainfo <- read.csv("result/metaInfo.csv", row.names = 1, stringsAsFactors = FALSE)
# adjust batch effect, diet, community
y <- datExpr[rownames(meta.genus), match(dat.metainfo$BIOCHEMICAL, colnames(datExpr))]
for (i in 1:ncol(y)) {
  mod <- lm(scale(y[,i]) ~ as.factor(meta.genus$WLS) + as.factor(meta.genus$Cohort) + meta.genus$Diet)
  y[,i] <- resid(mod)
}
otu.genus.comp.log <- log(otu.genus.comp) # take logrithm
colnames(otu.genus.comp.log) <- tax.genus@.Data[,"Genus"]
z <- otu.genus.comp.log
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
result_list_metabolite <- foreach(y=y) %dopar% cdmm_loop(y, z = z) # ~ 5hrs 7cores
stopCluster(cl)
df.meta <- list2df(result_list_metabolite, dat = y, colnames = c("metabolites", "genus", "bcv.prob", "stab.prob", "refitted.coef"))
tmp <- dat.metainfo[, c("BIOCHEMICAL", "SUPER_PATHWAY", "SUB_PATHWAY", "moduleColor")]
names(tmp)[1] <- "metabolites"
res <- merge(df.meta, tmp, by = "metabolites")
################ Association Anslysis End ##################

################ Mediation Anslysis Start ##################
### Distance-based Mediation Model
res <- data.frame(Diet = rep(c("Cellulose", "Inulin", "Pectin"), times = 3),
                  Trait = rep(c("GDATdivBW", "LiverTriglycerides", "Glucose"), each = 3))
res[] <- lapply(res, as.character)
res4cmt <- c()
res4perm <- c()
res4beta <- c()
res4alpha <- c()
n.draw <- 2000
for (i in 1:nrow(res)) {
  idx <- which(meta$Diet == res$Diet[i])
  otu.sub <- otu.genus[,idx]
  otu.sub <- t(Rarefy(t(otu.sub))$otu.tab.rff)
  meta.sub <- meta[idx,]
  BC <- as.matrix(vegdist(t(otu.sub), method = "bray"))
  JAC <- as.matrix(vegdist(t(otu.sub), 'jaccard', binary = TRUE))
  unifracs <- GUniFrac(t(otu.sub), tree, alpha = c(0, 0.5, 1))$unifracs
  unifracs <- abind(unifracs, BC, JAC)
  dimnames(unifracs)[[3]][6:7] <- c("d_BC", "d_JAC")
  trait <- scale(meta.sub[,res$Trait[i]])
  maximum_cohort <- max(meta.sub$Cohort)
  x <- meta.sub$Cohort
  x <- ifelse(x == maximum_cohort, 1, 0)
  s <- ifelse(meta.sub$WLS == 165, 1, 0)
  tmp <- PermanovaG(unifracs[, , c("d_BC", "d_JAC", "d_UW", "d_1")] ~ s, permutations = 2000)
  res4cmt <- rbind(res4cmt, round(as.matrix(tmp$p.tab), 4))
  s <- t(t(s))
  fit.y.x <- lm(trait ~ x)
  y <- t(t(resid(fit.y.x)))
  tmp <- PermanovaG(unifracs[, , c("d_BC", "d_JAC", "d_UW", "d_1")] ~ y, permutations = 2000)
  res4perm <- rbind(res4perm, round(as.matrix(tmp$p.tab), 4))
  m.list <- list(BC = vegdist(t(otu.sub), method = "bray"), 
                 JAC = as.matrix(vegdist(t(otu.sub), 'jaccard', binary = TRUE)),
                 UniFrac = unifracs[, , c('d_UW')],
                 WUniFrac = unifracs[, , c('d_1')])
  tmp <- MedOmniTest(s, y, m.list, nperm = 2000)
  res4beta <- rbind(res4beta, round(c(tmp$margPs, tmp$permP), 4))
  
  otu.sub4alpha <- otu_table(otu.sub, taxa_are_rows = TRUE)
  alphadiv <- estimate_richness(otu.sub4alpha, measures = c("Shannon", "InvSimpson"))
  faithpd <- picante::pd(t(otu.sub4alpha), tree = tree)
  alphadiv <- cbind(alphadiv, faithpd)
  alphadiv$SR <- NULL
  alpha.comb <- c()
  for (j in 1:3) {
    m <- alphadiv[,j]
    # asymptotic method/ delta method
    fit.m <- lm(m ~ x + s)
    fit.y <- lm(y ~ x + s + m)
    DE <- coef(fit.y)[3]
    DE.se <- sqrt(vcov(fit.y)[3,3])
    p.DE <- 2 * pnorm(abs(DE/DE.se), lower.tail = FALSE)
    IE <- coef(fit.m)[3] * coef(fit.y)[4]
    IE.se <- sqrt(coef(fit.m)[3]^2 * vcov(fit.y)[4,4] + coef(fit.y)[4]^2 * vcov(fit.m)[3,3])
    p.IE <- 2 * pnorm(abs(IE)/IE.se, lower.tail = FALSE)
    
    de.boot <- numeric(length = n.draw)
    ide.boot <- numeric(length = n.draw)
    for (bb in 1:n.draw){
      set.seed(bb*13)
      id.boot <- sample(1:ncol(otu.sub), replace = TRUE)
      x.boot <- x[id.boot]
      s.boot <- s[id.boot, ]
      m.boot <- m[id.boot]
      y.boot <- y[id.boot, ]
      # fit the mediator model
      fit.m <- lm(m.boot ~ x.boot + s.boot) 
      s.coef <- coef(fit.m)["s.boot"]
      # then fit the outcome model
      fit.y <- lm(y.boot ~ m.boot + s.boot + x.boot)
      m.coef <- coef(fit.y)["m.boot"]
      de.boot[bb] <- coef(fit.y)["s.boot"]
      ide.boot[bb] <- s.coef * m.coef 
    }
    de.boot <- na.omit(de.boot)
    de.boot <- de.boot - mean(de.boot)
    de.boot.p <- mean(abs(de.boot) > abs(DE))
    ide.boot <- na.omit(ide.boot)
    ide.boot <- ide.boot - mean(ide.boot)
    ide.boot.p <- mean(abs(ide.boot) > abs(IE))
    
    tmp <- round(unname(c(DE, p.DE, de.boot.p, IE, p.IE, ide.boot.p)), 4)
    tmp <- c(paste0(tmp[1], " (", tmp[2], ")"), paste0(tmp[1], " (", tmp[3], ")"),
             paste0(tmp[4], " (", tmp[5], ")"), paste0(tmp[4], " (", tmp[6], ")"))
    alpha.comb <- c(alpha.comb, tmp)
  }
  res4alpha <- rbind(res4alpha, alpha.comb)
}
res <- cbind(res, res4perm, res4cmt, res4beta, res4alpha)
names(res)[-(1:2)] <- c("PERMANOVA.Y.pval.BC", "PERMANOVA.Y.pval.JAC", "PERMANOVA.Y.pval.UWUniFrac", "PERMANOVA.Y.pval.WUniFrac", "PERMANOVA.Y.pval.Omnibus", 
                        "PERMANOVA.S.pval.BC", "PERMANOVA.S.pval.JAC", "PERMANOVA.S.pval.UWUniFrac", "PERMANOVA.S.pval.WUniFrac", "PERMANOVA.S.pval.Omnibus", 
                        "Mediation.pval.BC", "Mediation.pval.JAC", "Mediation.pval.UWUniFrac", "Mediation.pval.WUniFrac", "Mediation.pval.Omnibus",
                        "Mediation.DE.Shannon (asym.pval)", "Mediation.DE.Shannon (perm.pval)",
                        "Mediation.IDE.Shannon (asym.pval)", "Mediation.IDE.Shannon (perm.pval)", 
                        "Mediation.DE.InvSimpson (asym.pval)",  "Mediation.DE.InvSimpson (perm.pval)", 
                        "Mediation.IDE.InvSimpson (asym.pval)", "Mediation.IDE.InvSimpson (perm.pval)", 
                        "Mediation.DE.FaithsPD (asym.pval)", "Mediation.DE.FaithsPD (perm.pval)",
                        "Mediation.IDE.FaithsPD (asym.pval)", "Mediation.IDE.FaithsPD (perm.pval)")
rownames(res) <- 1:nrow(res)
pandoc.table(t(res), style = "rmarkdown", split.tables = Inf, caption = "The P-values for PERMANOVA and distance-based mediation tests (Rarefied, Genus level)")


### Causal Mediation Model
# without covariates, continuous outcome, p>>n, no exposure-mediator interaction
medTest <- function(X, S, M, Y, adaptive = FALSE, tsfm = FALSE, n.draw = 1000){
  
  ### X:	matrix for covariates			
  ### S:	matrix for exposure variable		
  ### M:  matrix for high-dimensional correlated mediators 
  ### Y:	matrix for outcome
  ### adaptive: TRUE if selecting the number of transformed mediators up to 80% variability
  ### tsfm: TRUE if doing the spectral decomposition
  ### n.draw: number of Monte-Carlo resampling
  n <- dim(M)[1]; p <- dim(M)[2]; q <- dim(X)[2] + 1 
  # n: sample size, p: # of mediators, q: # of covariates (include intercept)
  if (tsfm){
    # fit the mediator model 
    fit.m <- lm(M ~ S + X) # as.formula(paste0("X[,",1:(q-1),"]", collapse = " + "))
    
    Sigma <- cov(as.matrix(fit.m$residuals)) # do spectral decomposition and select the PCs could explain up to 80%
    Sigma_decomp <- eigen(Sigma)
    if (adaptive){
      n.factor <- sum(cumsum(Sigma_decomp$values)/sum(Sigma_decomp$values) < 0.8) + 1
      U <- Sigma_decomp$vectors[, 1:n.factor]
      p.prime <- n.factor # p.prime: number of independent mediators
    }else{
      U <- Sigma_decomp$vectors
      p.prime <- p
    }
    P <- M %*% U # construct the independent mediators
  }else{
    P <<- M
    p.prime <- p
  }
  # refit the mediator model
  fit.p <- lm(P ~ X + S) # as.formula(paste0("X[,",1:(q-1),"]", collapse = " + "))
  S.var <- vcov(fit.p)[(1:p.prime)*(q+1), (1:p.prime)*(q+1)] # calculate the variance of each exposure
  if(p.prime == 1){
    S.coef <- coef(fit.p)["S"]
    resid.p <<- resid(fit.p) + S * S.coef
  }else{
    S.coef <- coef(fit.p)["S",] # calculate the coefficient of each exposure
    resid.p <<- resid(fit.p) + S %*% S.coef
  }
  
  # then fit the outcome model jointly
  P.var <- NULL; P.coef <- NULL
  if (tsfm){
    for (j in 1:p.prime) {
      Pj <- P[,j]
      fit.y2pj <- lm(Y ~ Pj + S + X) # as.formula(paste0("X[,",1:(q-1),"]", collapse = " + "))
      P.var <- c(P.var, vcov(fit.y2pj)[2, 2]) # calculate the variance of each mediator
      P.coef <- c(P.coef, coef(fit.y2pj)[2]) # calculate the coefficient of each mediator
    }
    P.var <- diag(P.var)
  } else{
    fit.y <- lm(Y ~ P + S + X) # as.formula(paste0("X[,",1:(q-1),"]", collapse = " + "))
    P.var <- vcov(fit.y)[1:p.prime + 1, 1:p.prime + 1]
    P.coef <- coef(fit.y)[1:p.prime + 1]
    P.resid <<- P
    P.coef.resid <<- P.coef
    resid.y <<- resid(fit.y)
  }
  
  mat.var <- Matrix::bdiag(S.var, P.var)
  ide <- S.coef * P.coef # individual indirect effect 
  print(paste0("S.coef: ",S.coef))
  print(paste0("P.coef: ",P.coef))
  ######### p-value for component-wise mediation test using Delta method
  if(p.prime == 1){
    mat.coef <- t(t(c(P.coef, S.coef)))
  } else {
    mat.coef <- rbind(diag(P.coef), diag(S.coef))
  }
  ide.var <- t(mat.coef) %*% as.matrix(mat.var) %*% mat.coef
  stat.ide <- (t(ide) %*% solve(ide.var) %*% ide)[1]
  ide.p <- pchisq(stat.ide, df = p.prime, lower.tail = FALSE)
  
  ######### p-value for marginal mediation test using Delta method
  vec.coef <- c(P.coef, S.coef)
  total.ide <- sum(ide) # marginal mediation(indirect) effect
  total.ide.var <- (t(vec.coef) %*% as.matrix(mat.var) %*% vec.coef)[1]
  total.ide.p <- pchisq(total.ide^2/total.ide.var, df = 1, lower.tail = FALSE)
  # total.ide.var <- sum(ide.var) # variance of marginal effect
  # total.ide.p <- 2 * pnorm(abs(total.ide)/sqrt(total.ide.var), lower.tail = FALSE)
  
  l2.ide <- sum(ide^2) # l2 norm mediation(indirect) effect
  
  ######### p-value for marginal mediation test using Resampling methods [marginal ME & L2 norm ME]
  
  ## [Monte Carlo procedure for approximating the dist. of the marginal ME]
  set.seed(137)
  # sample from the distribution of regression coefficients
  coef.mc <- mvtnorm::rmvnorm(n.draw, mean = c(S.coef, P.coef), sigma = as.matrix(mat.var))
  S.coef.mc <- coef.mc[, 1:p.prime]
  P.coef.mc <- coef.mc[, p.prime + 1:p.prime]
  ide.mc <- t(t(S.coef.mc * P.coef.mc))  # each row denotes each mc / each column denotes the indiv. indirect effect
  
  ide.mc.var <- cov(ide.mc)
  stat.ide.mc <- (t(ide) %*% solve(ide.mc.var) %*% t(t(ide)))[1]
  ide.mc.p <- pchisq(stat.ide.mc, df = p.prime, lower.tail = FALSE)
  
  ide.mc.c <- t(t(ide.mc) - apply(ide.mc, 2, mean)) # same dimension with ide.mc
  l2.ide.mc <- apply(ide.mc.c^2, 1, sum)
  l2.ide.mc.p <- mean(l2.ide.mc > l2.ide)
  
  total.ide.mc <- apply(ide.mc.c, 1, sum)
  total.ide.mc.c <- total.ide.mc - mean(total.ide.mc)
  total.ide.mc.p <- 2*mean(total.ide.mc.c > abs(total.ide))
  
  
  res <- c(total.ide, sqrt(total.ide.var), total.ide.p, total.ide.mc.p, ide.p, ide.mc.p, l2.ide.mc.p)
  names(res) <- c("MarginalME", "se.MarginalME", "pval.Asym.MarginalME", "pval.MC.MarginalME", "pval.Asym.CompME", "pval.MC.CompME", "pval.L2norm.ME")
  return(res)
}

res <- data.frame(Diet = rep(c("Cellulose", "Inulin", "Pectin"), times = 3),
                  Trait = rep(c("GDATdivBW", "LiverTriglycerides", "Glucose"), each = 3))
res[] <- lapply(res, as.character)
### rank 1 ~ rank 2
res4rank.effect <- c()
for (i in 1:nrow(res)) {
  print(paste0("Now processing: ", res$Diet[i], " ~ ", res$Trait[i]))
  idx <- which(meta.genus$Diet == res$Diet[i])
  otu.sub <- otu.genus[,idx]
  meta.sub.i <- meta.genus[idx,]
  W.data <- data.table(data.frame(tax.genus@.Data[,1:6], otu.sub))
  setnames(W.data, 1:6, paste0("Rank", 1:6))
  otucols <- names(W.data)[-(1:6)]
  Rank.low <- paste0("Rank",6)  # change the low rank
  Rank.high <- paste0("Rank",5) # change the high rank
  tt <- W.data[, lapply(.SD, sum, na.rm = TRUE), .SDcols = otucols, by = list(get(Rank.low), get(Rank.high))]
  setnames(tt, 1:2, c(Rank.low, Rank.high))
  W.taxLow <- as.vector(unlist(tt[, Rank.low, with = FALSE]))
  W.taxHigh <- as.vector(unlist(tt[, Rank.high, with = FALSE]))
  W.P <- data.matrix(tt[, otucols, with = FALSE])
  # remove extremly rare genera (at least 90% of the observations are zeros)
  perc.zero.taxa <- rowMeans(W.P == 0)
  rm.idx.taxa <- which(perc.zero.taxa > 0.9)
  if(length(rm.idx.taxa) == 0){
    W2.P <- W.P
    W2.taxLow <- W.taxLow
    W2.taxHigh <- W.taxHigh
  }else{
    W2.P <- W.P[-rm.idx.taxa,]
    W2.taxLow <- W.taxLow[-rm.idx.taxa]
    W2.taxHigh <- W.taxHigh[-rm.idx.taxa]
  }
  # remove zero sequence depth samples
  perc.zero.sample <- colMeans(W2.P != 0)
  rm.idx.sample <- which(perc.zero.sample == 0)
  if(length(rm.idx.sample) == 0){
    W3.P <- W2.P
  }else{
    W3.P <- W2.P[,-rm.idx.sample]
  }
  W3.taxLow <- W2.taxLow
  W3.taxHigh <- W2.taxHigh
  W3.P[W3.P == 0] <- 0.5
  W3.P <- t(W3.P)/rowSums(t(W3.P))
  
  sel <- names(sort(table(W3.taxHigh)))[sort(table(W3.taxHigh)) >= 2]
  sel = "Rikenellaceae"
  effect = c()
  for (j in sel) {
    print(j)
    B.idx <- which(W3.taxHigh == j)
    B.P <- W3.P[,B.idx]
    colnames(B.P) <- W3.taxLow[B.idx]
    m <- B.P/rowSums(B.P)
    K <- which.max(colMeans(m))
    tmp <- colnames(m)[-K]
    print(colnames(m)[K])
    print(tmp)
    m <- log(m[,-K]/m[,K]) # use the most abundant taxon in m as reference
    if(is.vector(m)){
      m <- matrix(m, ncol = 1)
      colnames(m) <- tmp
    }
    rownames(m) <- sub("X","",rownames(B.P))
    rownames(m) <- sub("\\.", "-", rownames(m))
    meta.sub <- meta.sub.i[rownames(m),]
    y <- scale(meta.sub[,res$Trait[i]])
    y <- matrix(y, ncol = 1)
    maximum_cohort <- max(meta.sub$Cohort)
    x <- meta.sub$Cohort
    x <- ifelse(x == maximum_cohort, 1, 0)
    x <- matrix(x, ncol = 1)
    s <- ifelse(meta.sub$WLS == 165, 1, 0)
    s <- matrix(s, ncol = 1)
    if(length(unique(s)) == 1){
      effect <- paste(effect, "Unique exposure", sep = " | ")
      print(paste0(j, " done!"))
      next
    }else if(nrow(y) > ncol(x)+ncol(s)+ncol(m)+1){
      res.vec <- medTest(X = x, S = s, M = m, Y = y, adaptive = FALSE, tsfm = FALSE, n.draw = 1000)
    }else{
      res.vec <- medTest(X = x, S = s, M = m, Y = y, adaptive = TRUE, tsfm = TRUE, n.draw = 1000)
    }
    
    res.vec <- round(res.vec, 4)
    print(paste0(j, " done!"))
    effect <- paste(effect, paste0(j, " [", 
                                   paste(paste0("MarginalME:", res.vec[1]), 
                                         paste0("pval.MC.MarginalME:", res.vec[4]), 
                                         paste0("pval.L2norm.ME:", res.vec[7]),
                                         sep = ", "), 
                                   "]"), sep =  " | ")
  }
  res4rank.effect <- c(res4rank.effect, effect)
}
if(length(res4rank.effect)==0){
  res4rank.effect <- rep("Null", nrow(res))
}
res <- cbind(res, res4rank.effect)
names(res)[-(1:2)] <- c("Kingdom.effect","Phylum.effect", "Class.effect", "Order.effect", "Family.effect")
write.csv(t(res), file = "result/microMedRes_tree.csv")
################ Mediation Anslysis End ##################