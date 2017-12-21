library(DESeq2)
library(pheatmap)
library(RColorBrewer)

load("./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/rawCounts_Amanda_ENCODE_n65.rda")
encCounts = geneCounts
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# load the phenotype table
encpd = read.table("./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/ENCODEpd.txt")
length(unique(encpd$CellType)) # 11 cell types included
encpd = encpd[grep("FastqRd1", encpd$File_Uploaded_to_SRA),]
encpd = encpd[match(colnames(encCounts), encpd$SRA_Run),]
rownames(encpd) = encpd$SRA_Run
encpd$WorkingID = paste0(encpd$CellType,":",encpd$Fraction)

# Make DESeq2 Objects for samples both with and without our samples
combinedCounts = cbind(encCounts[match(rownames(geneCounts), rownames(encCounts)),],geneCounts)
encCounts = encCounts[rowSums(encCounts)>0,] # 51502 genes expressed
dim(combinedCounts)

combPD = data.frame(sampleID = c(as.character(encpd$SRA_Run), rownames(pd)), Fraction = c(as.character(encpd$Fraction), pd$Zone), 
                    CellType = c(as.character(encpd$CellType), rep.int("Brain Tissue",23)))
combPD = combPD[-grep("RiboZero", combPD$sampleID),]
rownames(combPD) = combPD$sampleID
combPD$Age = c(rep.int("Cell Line",65), rep.int("Adult",6), rep.int("Prenatal",6))
combPD$WorkingID = paste0(combPD$CellType, ":",combPD$Fraction,":", combPD$Age)
combinedCounts = combinedCounts[,which(colnames(combinedCounts) %in% combPD$sampleID)]

dds = DESeqDataSetFromMatrix(countData = encCounts, colData = encpd, design = ~ Fraction + CellType)
combdds = DESeqDataSetFromMatrix(countData = combinedCounts, colData = combPD, design = ~ Fraction + CellType)

## Show that relationship between fractions is more different in neuronal lineage than others 
# Euclidean distance heatmap
rlog  = rlog(encCounts)
combRLog = rlog(combinedCounts)
save(rlog,combRLog,file="./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/rlog_transformed_encode_counts.rda")
load("./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/rlog_transformed_encode_counts.rda")

# Euclidean distance between samples heatmap
sampleDists <- dist(t(rlog))
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
match(colnames(encCounts),encpd$SRA_Run)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = encpd$WorkingID
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples")
# encode_euclidean_distance_heatmap.pdf

sampleDists <- dist(t(combRLog))
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
match(colnames(combinedCounts),combPD$sampleID)
colnames(sampleDistMatrix) = rownames(sampleDistMatrix) = combPD$WorkingID
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples")
# encode_andBrain_euclidean_distance_heatmap.pdf


# PCA
load("./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/rlog_transformed_encode_counts.rda")
rlogdds  = rlog(dds)
rlog.comb.dds = rlog(combdds)
save(rlog,combRLog,rlogdds,rlog.comb.dds,file="./Dropbox/sorted_figures/new/github_controlled/other_datasets/data/rlog_transformed_encode_counts.rda")

# PCA Functions:

plotPCA <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

plotPCA2 <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
    geom_point(size = 3) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

plotPCA3 <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

plotPCA4 <- function (x, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PC4 = pca$x[, 4], PC5 = pca$x[, 5], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC4", y = "PC5", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) +
    ylab(paste0("PC5: ", round(percentVar[5] * 100), "% variance")) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
}

# Run PCA
require(genefilter)

plotPCA(rlogdds, intgroup = c("CellType", "Fraction"))
plotPCA2(rlogdds, intgroup = c("CellType", "Fraction"))
plotPCA3(rlogdds, intgroup = c("CellType", "Fraction"))
plotPCA4(rlogdds, intgroup = c("CellType", "Fraction"))

plotPCA(rlog.comb.dds, intgroup = c("CellType", "Fraction"))
plotPCA2(rlog.comb.dds, intgroup = c("CellType", "Fraction"))
plotPCA3(rlog.comb.dds, intgroup = c("CellType", "Fraction"))
plotPCA4(rlog.comb.dds, intgroup = c("CellType", "Fraction"))

plotPCA(rlog.comb.dds, intgroup = c("Fraction","Age"))
plotPCA2(rlog.comb.dds, intgroup = c("Fraction","Age"))
plotPCA3(rlog.comb.dds, intgroup = c("Fraction","Age"))
plotPCA4(rlog.comb.dds, intgroup = c("Fraction","Age"))

#Show that there are more DEGs by fraction in SKNH-SH than in the H1hesc
SKNSHpd = encpd[which(encpd$CellType=="Sknsh"),]
SKNSHCounts = encCounts[,which(colnames(encCounts) %in% SKNSHpd$SRA_Run)]
SKNSHdds = DESeqDataSetFromMatrix(countData = SKNSHCounts, colData = SKNSHpd, design = ~ Fraction)
SKNSHdds = DESeq(SKNSHdds)
SKNSHres = results(SKNSHdds)

H1hescpd = encpd[which(encpd$CellType=="H1hesc"),]
H1hescCounts = encCounts[,which(colnames(encCounts) %in% H1hescpd$SRA_Run)]
H1hescdds = DESeqDataSetFromMatrix(countData = H1hescCounts, colData = H1hescpd, design = ~ Fraction)
H1hescdds = DESeq(H1hescdds)
H1hescres = results(H1hescdds)

dim(SKNSHres[which(SKNSHres$padj<=0.05 & abs(SKNSHres$log2FoldChange)>=1),]) # 14033
dim(H1hescres[which(H1hescres$padj<=0.05 & abs(H1hescres$log2FoldChange)>=1),]) # 7 (!!!)

plotMA(SKNSHres, alpha = 0.05, main="SK-N-SH: Nucleus vs.Cytosol", ylim=c(-6,6)) # SK-N-SH_Nucleus_vs_Cytosol_MAplot.pdf
plotMA(H1hescres, alpha = 0.05, main="H1-HESC: Nucleus vs.Cytosol", ylim=c(-6,6)) # H1hesc_Nucleus_vs_Cytosol_MAplot.pdf

# Check fractionation quality

ens = c("ENSG00000075624","ENSG00000102081","ENSG00000229807","ENSG00000251562")
gencode = geneMap[which(geneMap$ensemblID %in%ens),"gencodeID"] 
encCounts[gencode,]

SKNSHres[gencode,]
H1hescres[gencode,]

