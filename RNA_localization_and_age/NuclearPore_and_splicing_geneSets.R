library(readxl)
library(DESeq2)
library("ggplot2")

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

# Nucleoporin genes from http://amigo.geneontology.org/grebe
nups = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/nuclearPoreGenes.xlsx", col_names = FALSE))
nups = read.table("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/nuclear_pore_genes.txt", sep = "\t", header = F)
nupGenes = c(unique(nups$V3), "NUPL1")
nups = geneMap[which(geneMap$Symbol %in% nupGenes),]
nupCounts = geneCounts.down[which(rownames(geneCounts) %in% rownames(nups)),grep("polyA", colnames(geneCounts.down))]

# Splicing genes
splicepos = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splicingPositive.xlsx", col_names = FALSE))
splicing = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splicing.xlsx", col_names = FALSE))
spliceneg = as.data.frame(read_excel("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splicingNegative.xlsx", col_names = FALSE))
spl = rbind(splicepos, splicing, spliceneg)
geneMap[which(rownames(geneMap)=="ENSG00000255663"),6] = "I3L521"
splGenes = c(unique(spl$X2), "C11orf35", "ZNF259") # excludes RBMY1A1
spl = geneMap[which(geneMap$Symbol %in% splGenes),]
splCounts = geneCounts[which(rownames(geneCounts) %in% spl$EnsID),]

# Differential expression of 2 gene sets
Nupdds.int <- DESeqDataSetFromMatrix(countData = nupCounts, colData = pd, design = ~ Library + Zone + Fetal + Fetal:Zone)
Nupdds.int <- DESeq(Nupdds.int)
Nupres.int = results(Nupdds.int)
Nupdds.age <- DESeqDataSetFromMatrix(countData = nupCounts, colData = pd, design = ~ Library + Zone + Fetal)
Nupdds.age <- DESeq(Nupdds.age)
Nupres.age = results(Nupdds.age)
Nupdds.frac <- DESeqDataSetFromMatrix(countData = nupCounts, colData = pd, design = ~ Library + Fetal + Zone)
Nupdds.frac <- DESeq(Nupdds.frac)
Nupres.frac = results(Nupdds.frac)

dat <- Nupres.int[order(rownames(Nupres.int)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Nupres.int <- data.frame(dat, g)
Nupres.int <- Nupres.int[order(Nupres.int$padj),]
write.csv(Nupres.int, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nupres.int.csv")
dat <- Nupres.age[order(rownames(Nupres.age)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Nupres.age <- data.frame(dat, g)
Nupres.age <- Nupres.age[order(Nupres.age$padj),]
write.csv(Nupres.age, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nupres.age.csv")
dat <- Nupres.frac[order(rownames(Nupres.frac)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Nupres.frac <- data.frame(dat, g)
Nupres.frac <- Nupres.frac[order(Nupres.frac$padj),]
write.csv(Nupres.frac, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nupres.frac.csv")

nupres = list(Age = Nupres.age[which(Nupres.age$padj<=0.05),], Fraction = Nupres.frac[which(Nupres.frac$padj<=0.05),], 
              Interaction = Nupres.int[which(Nupres.int$padj<=0.05),])
nupres = lapply(nupres, function(x) x[order(abs(x$log2FoldChange), decreasing =T),])

# Age genes
ENSG00000183486 = plotCounts(Nupdds, "ENSG00000183486", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000183486, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("MX2") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000030066 = plotCounts(Nupdds, "ENSG00000030066", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000030066, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NUP160") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000102900 = plotCounts(Nupdds, "ENSG00000102900", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000102900, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NUP93") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000204764 = plotCounts(Nupdds, "ENSG00000204764", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000204764, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("RANBP17") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000198088 = plotCounts(Nupdds, "ENSG00000198088", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000198088, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NUP62CL") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000132661 = plotCounts(Nupdds, "ENSG00000132661", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000132661, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NPIPA1") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000213024 = plotCounts(Nupdds, "ENSG00000213024", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000213024, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NUP62") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Interaction genes
ENSG00000183426 = plotCounts(Nupdds, "ENSG00000183426", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000183426, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NPIPA1") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000162231 = plotCounts(Nupdds, "ENSG00000162231", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000162231, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NXF1") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000094914 = plotCounts(Nupdds, "ENSG00000094914", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000094914, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("AAAS") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000132341 = plotCounts(Nupdds, "ENSG00000132341", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000132341, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("RAN") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ENSG00000030066 = plotCounts(Nupdds, "ENSG00000030066", intgroup = c("Zone", "Fetal"), returnData=TRUE)
ggplot(ENSG00000030066, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() +
  geom_jitter() + ggtitle("NUP160") + ylab("Count") + 
  xlab("") + scale_y_log10(breaks=c(500,1000,2000,4000,8000)) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Splicing genes
spldds.int <- DESeqDataSetFromMatrix(countData = splCounts, colData = pd, design = ~ Library + Zone + Fetal + Fetal:Zone)
spldds.int <- DESeq(spldds.int)
splres.int = results(spldds.int)
spldds.age <- DESeqDataSetFromMatrix(countData = splCounts, colData = pd, design = ~ Library + Zone + Fetal)
spldds.age <- DESeq(spldds.age)
splres.age = results(spldds.age)
spldds.frac <- DESeqDataSetFromMatrix(countData = splCounts, colData = pd, design = ~ Library + Fetal + Zone)
spldds.frac <- DESeq(spldds.frac)
splres.frac = results(spldds.frac)

dat <- splres.int[order(rownames(splres.int)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
splres.int <- data.frame(dat, g)
splres.int <- splres.int[order(splres.int$padj),]
write.csv(splres.int, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splres.int.csv")
dat <- splres.age[order(rownames(splres.age)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
splres.age <- data.frame(dat, g)
splres.age <- splres.age[order(splres.age$padj),]
write.csv(splres.age, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splres.age.csv")
dat <- splres.frac[order(rownames(splres.frac)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
splres.frac <- data.frame(dat, g)
splres.frac <- splres.frac[order(splres.frac$padj),]
write.csv(splres.frac, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/splres.frac.csv")

splres = list(Age = splres.age[which(splres.age$padj<=0.05),], Fraction = splres.frac[which(splres.frac$padj<=0.05),], 
              Interaction = splres.int[which(splres.int$padj<=0.05),])
splres = lapply(splres, function(x) x[order(abs(x$log2FoldChange), decreasing =T),])
splres

plotMA(Nupres, alpha = 0.05, main = "Nuclear Pore Genes")
plotMA(splres.frac, alpha = 0.05, main = "Splicing Genes")

