library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(pvclust)
library(rafalib)

load("./Dropbox/sorted_figures/github_controlled/other_datasets/data/rawCounts_Amanda_ENCODE_n65.rda")
encCounts = geneCounts

# load the phenotype table
encpd = read.table("./Dropbox/sorted_figures/github_controlled/other_datasets/data/ENCODEpd.txt")
length(unique(encpd$CellType)) # 11 cell types included
encpd = encpd[grep("FastqRd1", encpd$File_Uploaded_to_SRA),]
encpd = encpd[match(colnames(encCounts), encpd$SRA_Run),]
rownames(encpd) = encpd$SRA_Run
encpd$WorkingID = paste0(encpd$CellType,":",encpd$Fraction)
encpd[encpd$CellType=="Gm12878","Description"] = "Lymphoblasts"
encpd[encpd$CellType=="Huvec","Description"] = "Umbilical Vein Endothelial Cells"
encpd[encpd$CellType=="Hepg2","Description"] = "Hepatocytes"
encpd[encpd$CellType=="H1hesc","Description"] = "hESCs"
encpd[encpd$CellType=="K562","Description"] = "Myelogenous Leukemia"
encpd[encpd$CellType=="Nhek","Description"] = "Keratinocytes"
encpd[encpd$CellType=="Helas3","Description"] = "Cervix Adenocarcinoma"
encpd[encpd$CellType=="Imr90","Description"] = "Myofibroblasts"
encpd[encpd$CellType=="Mcf7","Description"] = "Mammary Epithelium"
encpd[encpd$CellType=="A549","Description"] = "Adenocarcinomic Alveolar Basal Epithelial"
encpd[encpd$CellType=="Sknsh","Description"] = "Neuroblastoma"
head(encpd)


# Make DESeq2 Objects
encCounts = encCounts[rowSums(encCounts)>0,] # 51502 genes expressed
dim(encCounts)

celltypes = as.character(unique(encpd$CellType))
dds = list()
for (i in 1:length(celltypes)) {
  dds[[i]] = DESeqDataSetFromMatrix(countData = encCounts[,which(colnames(encCounts) %in% encpd[which(encpd$CellType==celltypes[i]),"SRA_Run"])], 
                                    colData = encpd[which(encpd$CellType==celltypes[i]),], design = ~ Fraction)
}
names(dds) = celltypes
dds$Fraction = DESeqDataSetFromMatrix(countData = encCounts, colData = encpd, design = ~ CellType + Fraction)

dds = lapply(dds, DESeq)
res = lapply(dds, results)
sigres1 = lapply(res, function(x) list(Cytoplasmic = x[which(x$log2FoldChange<0 & x$padj<=0.05 & abs(x$log2FoldChange)>1),],
                                      Nuclear = x[which(x$log2FoldChange>0 & x$padj<=0.05 & abs(x$log2FoldChange)>1),]))
sigres = lapply(res, function(x) list(Cytoplasmic = x[which(x$log2FoldChange<0 & x$padj<=0.05),],
                                       Nuclear = x[which(x$log2FoldChange>0 & x$padj<=0.05),]))

save(res, dds,sigres, encpd, file = "./Dropbox/sorted_figures/github_controlled/other_datasets/data/ENCODE_DESeq2_output.rda")

load("./Dropbox/sorted_figures/github_controlled/other_datasets/data/ENCODE_DESeq2_output.rda")

x = data.frame(CellType = names(res), DEGs = unlist(lapply(res, function(x) nrow(x[which(x$padj<=0.05 & abs(x$log2FoldChange)>=1),]))), row.names = NULL)
#   CellType  DEGs
#1   Gm12878  7651
#2     Huvec  8931
#3     Hepg2  9502
#4    H1hesc     9
#5      K562  6860
#6      Nhek  8227
#7    Helas3  6180
#8     Imr90  6518
#9      Mcf7 14013
#10     A549  8648
#11    Sknsh 13985
#12 Fraction 15372

pdf("./Dropbox/sorted_figures/github_controlled/other_datasets/figures/MAplots_encode_byFraction.pdf", width= 4,height = 4)
for (i in 1:length(res)) {
  p = plotMA(res[[i]], alpha = 0.05, main=paste0(names(res)[i], ": Nucleus vs.Cytosol"), ylim=c(-8,8))
  print(p)
}
dev.off()

res = Map(cbind, lapply(res, as.data.frame), CellType = as.list(names(res)))
res = do.call(rbind, res[1:11])
res$FDR = "FDR>0.05"
res[which(res$padj<=0.05),"FDR"] = "FDR<0.05"
res$CellType = factor(res$CellType, levels = c("H1hesc","Helas3","Imr90","K562","Gm12878","Nhek","A549","Huvec","Hepg2","Sknsh","Mcf7"))
head(res)

pdf("./Dropbox/sorted_figures/github_controlled/other_datasets/figures/ENCODE_MAplots.pdf", width=6.5, height=5.25)
ggplot(res, aes(x=baseMean/1000, y=log2FoldChange)) + 
  geom_point(aes(colour = factor(FDR))) + scale_colour_manual(values=c("red3","gray47")) +
  facet_wrap(. ~ CellType) + geom_hline(aes(yintercept=0), linetype="dashed") +
  ylab("Log2 Fold Change") + ylim(-15,15) +
  xlab("Mean Normalized Counts/1000") + xlim(0,100) +
  ggtitle("MA Plot: Nucleus Vs. Cytoplasm") + 
  theme(title = element_text(size = 16), text = element_text(size = 16)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"), legend.position = c(.87, .15))
dev.off()


# Cluster the LFC between fractions for each cell type and see how it looks

lapply(res, head)
mat = do.call(cbind, lapply(res, function(x) x$log2FoldChange))
rownames(mat) = rownames(res[[1]])
head(mat)
mat = mat[complete.cases(mat),1:11]
sigmat = mat[which(rownames(mat) %in% unlist(lapply(res[1:11], function(x) rownames(x[which(x$padj<=0.05 & abs(x$log2FoldChange)>=1),])))),]
dim(sigmat) # 18803

hc = hclust(dist(t(mat[,1:11])), method="ward")
hc_cut = cutree(hc, k= 9)
shc = hclust(dist(t(sigmat[,1:11])), method="ward")
shc_cut = cutree(shc, k= 9)

pdf("./Dropbox/sorted_figures/github_controlled/other_datasets/figures/hclust_encode_Fraction_LFCs.pdf", h = 4, w = 4)
palette(brewer.pal(12,"Paired"))
myplclust(hc, lab.col=hc_cut, xlab="", hang=0.05, cex=1.1, main = "Cluster by Fraction LFC")
myplclust(shc, lab.col=shc_cut, xlab="", hang=0.05, cex=1.1, main = "Cluster by Fraction LFC")
dev.off()


# Euclidean distance between samples heatmap

pdf("./Dropbox/sorted_figures/github_controlled/other_datasets/figures/heatmap_encode_Fraction_LFCs.pdf", h = 4, w = 4)
sampleDists <- dist(t(mat))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples\nby Fraction LFC: All Genes")
sampleDists <- dist(t(sigmat))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=colors, main="Euclidean Distance Between Samples\nby Fraction LFC: FDR<0.05")
dev.off()


## look for enrichment in nuclear-enriched genes (FDR<0.05) for all cell types and gene sets

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
pgc2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]


## Enrichment in genes differentially expressed by fraction in ENCODE samples

load("./Dropbox/sorted_figures/github_controlled/other_datasets/data/ENCODE_DESeq2_output.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"Symbol"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$Gene.Symbol) %in% geneuniverse ), ] # drop genes that are not present in the test set
aej_sets_expressed$Gene.Symbol = as.character(aej_sets_expressed$Gene.Symbol)
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
splitSets$PGC2 = data.frame(Gene.Symbol = unique(pgc2$Symbol))
splitSets = splitSets[which(names(splitSets)!="SCZ PGC GWAS")]

inGroup = unlist(lapply(sigres, function(x) lapply(x, function(y) 
                          as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(y)),"Symbol"]))))), recursive = F)
outGroup = unlist(lapply(sigres, function(x) lapply(x, function(y) 
                    geneuniverse[!(geneuniverse %in% as.character(geneMap[which(geneMap$gencodeID %in% rownames(y)),"Symbol"]))])), recursive = F)

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$Gene.Symbol),sum(!(inG %in% x$Gene.Symbol)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$Gene.Symbol), sum(!(outG %in% x$Gene.Symbol)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)

enrich = lapply(enrich, data.frame)
enrich = do.call(rbind, Map(cbind, CellType = as.list(gsub("\\..*","",names(enrich))), Compartment = as.list(gsub("^.*\\.","",names(enrich))), lapply(enrich, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")
enrich$Description = encpd[match(enrich$CellType, encpd$CellType),"Description"] 

write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/github_controlled/other_datasets/data/Birnbaum_geneSet_enrichment_ENCODE_FractionDEGs.csv")
bb = read.csv("./Dropbox/sorted_figures/github_controlled/other_datasets/data/Birnbaum_geneSet_enrichment_ENCODE_FractionDEGs.csv")


bb[bb$CellType=="Fraction" & bb$FDR<=0.05,]
#        X CellType Compartment      GeneSet      P.Value OddsRatio          FDR Description
#  222 222 Fraction Cytoplasmic ASD.DATABASE 8.490516e-04 0.5170464 0.0135848257        <NA>
#  231 231 Fraction     Nuclear      ASD.CNV 2.885493e-03 1.6313473 0.0329770656        <NA>
#  232 232 Fraction     Nuclear ASD.DATABASE 1.781139e-03 1.5253235 0.0237485258        <NA>
#  233 233 Fraction     Nuclear    BPAD.GWAS 6.068541e-03 1.6729367 0.0485483286        <NA>
#  240 240 Fraction     Nuclear         PGC2 7.861570e-06 1.6088334 0.0003493339        <NA>


bb[bb$CellType=="Fraction" & bb$Compartment=="Nuclear" & bb$FDR<=0.05,"GeneSet"]
# ASD.CNV      BPAD.GWAS    PGC2 in nuclear
# ASD.DATABASE in both compartments


bb[bb$CellType=="Fraction" & bb$FDR>0.05,]
#        X CellType Compartment           GeneSet     P.Value OddsRatio        FDR Description
#  221 221 Fraction Cytoplasmic           ASD.CNV 0.185233606 1.3006196 0.51693099        <NA>
#  223 223 Fraction Cytoplasmic         BPAD.GWAS 0.558069172 1.1354921 0.88699736        <NA>
#  224 224 Fraction Cytoplasmic                ID 0.786822993 1.0741156 1.00000000        <NA>
#  225 225 Fraction Cytoplasmic               NDD 0.495322855 0.6422010 0.84912489        <NA>
#  226 226 Fraction Cytoplasmic Neurodegenerative 0.134184808 1.6467325 0.42374150        <NA>
#  227 227 Fraction Cytoplasmic           SCZ.CNV 0.028278599 1.6741088 0.13850742        <NA>
#  228 228 Fraction Cytoplasmic SCZ.Meta.analysis 0.289077372 0.5384979 0.65546974        <NA>
#  229 229 Fraction Cytoplasmic           SCZ.SNV 0.173583839 0.7487096 0.49011908        <NA>
#  230 230 Fraction Cytoplasmic              PGC2 0.692002009 1.0517921 0.96558420        <NA>
#  234 234 Fraction     Nuclear                ID 0.367592766 1.2267108 0.71725418        <NA>
#  235 235 Fraction     Nuclear               NDD 0.011474319 2.5489231 0.07396458        <NA>
#  236 236 Fraction     Nuclear Neurodegenerative 0.877463963 1.0384367 1.00000000        <NA>
#  237 237 Fraction     Nuclear           SCZ.CNV 0.024649841 1.5893570 0.12587153        <NA>
#  238 238 Fraction     Nuclear SCZ.Meta.analysis 0.859134726 0.8922439 1.00000000        <NA>
#  239 239 Fraction     Nuclear           SCZ.SNV 0.008026742 1.4822411 0.06020057        <NA>
  

bb[bb$CellType!="Fraction" & bb$FDR<=0.05,]
#      X CellType Compartment           GeneSet      P.Value OddsRatio          FDR                               Description
#12   12  Gm12878     Nuclear      ASD.DATABASE 2.697138e-03 1.6858121 3.261696e-02                              Lymphoblasts
#22   22    Huvec Cytoplasmic      ASD.DATABASE 5.088443e-03 0.4593873 4.854833e-02          Umbilical Vein Endothelial Cells
#23   23    Huvec Cytoplasmic         BPAD.GWAS 5.806018e-03 1.9951368 4.854833e-02          Umbilical Vein Endothelial Cells
#27   27    Huvec Cytoplasmic           SCZ.CNV 4.043824e-03 2.2089899 4.072896e-02          Umbilical Vein Endothelial Cells
#112 112     Nhek     Nuclear      ASD.DATABASE 5.813664e-05 1.9250985 1.395279e-03                             Keratinocytes
#115 115     Nhek     Nuclear               NDD 1.463389e-03 3.6286198 2.065961e-02                             Keratinocytes
#119 119     Nhek     Nuclear           SCZ.SNV 6.484416e-08 2.5071082 7.781299e-06                             Keratinocytes
#149 149    Imr90 Cytoplasmic           SCZ.SNV 4.072896e-03 0.2832474 4.072896e-02                            Myofibroblasts
#160 160    Imr90     Nuclear              PGC2 5.933291e-03 1.5188819 4.854833e-02                            Myofibroblasts
#166 166     Mcf7 Cytoplasmic Neurodegenerative 1.768734e-05 3.5981732 5.306202e-04                        Mammary Epithelium
#170 170     Mcf7 Cytoplasmic              PGC2 5.320262e-03 1.3612497 4.854833e-02                        Mammary Epithelium
#171 171     Mcf7     Nuclear           ASD.CNV 5.780917e-03 1.6309852 4.854833e-02                        Mammary Epithelium
#173 173     Mcf7     Nuclear         BPAD.GWAS 1.388786e-03 1.8954313 2.065961e-02                        Mammary Epithelium
#177 177     Mcf7     Nuclear           SCZ.CNV 3.045959e-04 2.1565118 6.091918e-03                        Mammary Epithelium
#191 191     A549     Nuclear           ASD.CNV 2.598726e-08 2.7662965 6.236942e-06 Adenocarcinomic Alveolar Basal Epithelial
#192 192     A549     Nuclear      ASD.DATABASE 8.113216e-04 1.7174944 1.358483e-02 Adenocarcinomic Alveolar Basal Epithelial
#193 193     A549     Nuclear         BPAD.GWAS 7.968518e-04 2.1037381 1.358483e-02 Adenocarcinomic Alveolar Basal Epithelial
#197 197     A549     Nuclear           SCZ.CNV 4.060575e-07 3.1146066 3.248460e-05 Adenocarcinomic Alveolar Basal Epithelial
#199 199     A549     Nuclear           SCZ.SNV 7.183570e-05 1.9752339 1.567324e-03 Adenocarcinomic Alveolar Basal Epithelial
#200 200     A549     Nuclear              PGC2 1.554133e-05 1.7369899 5.306202e-04 Adenocarcinomic Alveolar Basal Epithelial
#206 206    Sknsh Cytoplasmic Neurodegenerative 2.718080e-03 2.5028067 3.261696e-02                             Neuroblastoma
#211 211    Sknsh     Nuclear           ASD.CNV 4.733912e-06 2.2330831 2.840347e-04                             Neuroblastoma
#213 213    Sknsh     Nuclear         BPAD.GWAS 3.567218e-03 1.8393558 3.891510e-02                             Neuroblastoma
#217 217    Sknsh     Nuclear           SCZ.CNV 2.114809e-05 2.5216408 5.639490e-04                             Neuroblastoma
#220 220    Sknsh     Nuclear              PGC2 8.733346e-06 1.6941233 3.493339e-04                             Neuroblastoma

dim(bb[bb$CellType!="Fraction" & bb$FDR<=0.05,]) # 25
dim(bb[bb$CellType!="Fraction" & bb$FDR<=0.05 & bb$OddsRatio<1,]) # 2
bb[bb$CellType!="Fraction" & bb$FDR<=0.05 & bb$OddsRatio<1,]
#22   22    Huvec Cytoplasmic ASD.DATABASE 0.0050884432 0.4593873 0.04854833 Umbilical Vein Endothelial Cells
#149 149    Imr90 Cytoplasmic      SCZ.SNV 0.0040728962 0.2832474 0.04072896                   Myofibroblasts

  
x = split(bb[bb$FDR<=0.05,], bb[bb$FDR<=0.05,"GeneSet"])
x = lapply(x, function(y) split(y, y$CellType))
do.call(rbind, Map(cbind, GeneSet = as.list(names(lapply(x, elementNROWS))), lapply(x, function(y) data.frame(rbind(elementNROWS(y))))))
#             GeneSet A549 Fraction Gm12878 H1hesc Helas3 Hepg2 Huvec Imr90 K562 Mcf7 Nhek Sknsh
#1            ASD.CNV    1        1       0      0      0     0     0     0    0    1    0     1
#2       ASD.DATABASE    1        2       1      0      0     0     1     0    0    0    1     0
#3          BPAD.GWAS    1        1       0      0      0     0     1     0    0    1    0     1
#4                 ID    0        0       0      0      0     0     0     0    0    0    0     0
#5                NDD    0        0       0      0      0     0     0     0    0    0    1     0
#6  Neurodegenerative    0        0       0      0      0     0     0     0    0    1    0     1
#7               PGC2    1        1       0      0      0     0     0     1    0    1    0     1
#8            SCZ.CNV    1        0       0      0      0     0     1     0    0    1    0     1
#9  SCZ.Meta.analysis    0        0       0      0      0     0     0     0    0    0    0     0
#10           SCZ.SNV    1        0       0      0      0     0     0     1    0    0    1     0


x = split(bb[bb$FDR<=0.05 & bb$Compartment=="Nuclear",], bb[bb$FDR<=0.05 & bb$Compartment=="Nuclear","GeneSet"])
x = lapply(x, function(y) split(y, y$CellType))
do.call(rbind, Map(cbind, GeneSet = as.list(names(lapply(x, elementNROWS))), lapply(x, function(y) data.frame(rbind(elementNROWS(y))))))
#             GeneSet A549 Fraction Gm12878 H1hesc Helas3 Hepg2 Huvec Imr90 K562 Mcf7 Nhek Sknsh
#1            ASD.CNV    1        1       0      0      0     0     0     0    0    1    0     1
#2       ASD.DATABASE    1        1       1      0      0     0     0     0    0    0    1     0
#3          BPAD.GWAS    1        1       0      0      0     0     0     0    0    1    0     1
#4                 ID    0        0       0      0      0     0     0     0    0    0    0     0
#5                NDD    0        0       0      0      0     0     0     0    0    0    1     0
#6  Neurodegenerative    0        0       0      0      0     0     0     0    0    0    0     0
#7               PGC2    1        1       0      0      0     0     0     1    0    0    0     1
#8            SCZ.CNV    1        0       0      0      0     0     0     0    0    1    0     1
#9  SCZ.Meta.analysis    0        0       0      0      0     0     0     0    0    0    0     0
#10           SCZ.SNV    1        0       0      0      0     0     0     0    0    0    1     0

c("ASD.CNV","ASD.DATABASE","BPAD.GWAS","ID","NDD","Neurodegenerative","PGC2","SCZ.CNV","SCZ.Meta.analysis","SCZ.SNV")

x = split(bb[bb$FDR<=0.05 & bb$Compartment=="Cytoplasmic",], bb[bb$FDR<=0.05 & bb$Compartment=="Cytoplasmic","GeneSet"])
x = lapply(x, function(y) split(y, y$CellType))
do.call(rbind, Map(cbind, GeneSet = as.list(names(lapply(x, elementNROWS))), lapply(x, function(y) data.frame(rbind(elementNROWS(y))))))

#GeneSet A549 Fraction Gm12878 H1hesc Helas3 Hepg2 Huvec Imr90 K562 Mcf7 Nhek Sknsh
#1            ASD.CNV    0        0       0      0      0     0     0     0    0    0    0     0
#2       ASD.DATABASE    0        1       0      0      0     0     1     0    0    0    0     0
#3          BPAD.GWAS    0        0       0      0      0     0     1     0    0    0    0     0
#4                 ID    0        0       0      0      0     0     0     0    0    0    0     0
#5                NDD    0        0       0      0      0     0     0     0    0    0    0     0
#6  Neurodegenerative    0        0       0      0      0     0     0     0    0    1    0     1
#7               PGC2    0        0       0      0      0     0     0     0    0    1    0     0
#8            SCZ.CNV    0        0       0      0      0     0     1     0    0    0    0     0
#9  SCZ.Meta.analysis    0        0       0      0      0     0     0     0    0    0    0     0
#10           SCZ.SNV    0        0       0      0      0     0     0     1    0    0    0     0
                
x = bb[bb$CellType!="Fraction" & bb$FDR<=0.05,]
x$GeneSet = gsub("ASD.CNV", "ASD\n(CNV)", x$GeneSet)
x$GeneSet = gsub("ASD.DATABASE", "ASD\n(Database)", x$GeneSet)
x$GeneSet = gsub("BPAD.GWAS", "BPAD\n(GWAS)", x$GeneSet)
x$GeneSet = gsub("ID", "Intellectual\nDisability", x$GeneSet)
x$GeneSet = gsub("NDD", "Neuro-\ndevel.", x$GeneSet)
x$GeneSet = gsub("Neurodegenerative", "Neuro-\ndegen.", x$GeneSet)
x$GeneSet = gsub("SCZ.Meta.analysis", "SCZ\n(Meta\nanalysis)", x$GeneSet)
x$GeneSet = gsub("SCZ.SNV", "SCZ\n(SNV)", x$GeneSet)
x$GeneSet = gsub("PGC2", "SCZ\n(PGC2)", x$GeneSet)
x$GeneSet = gsub("SCZ.CNV", "SCZ\n(CNV)", x$GeneSet)
x$cat = ifelse(x$GeneSet %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)"), "Nuclear\nin Both", "Not\nNuclear")
x[which(x$GeneSet %in% c("BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")),"cat"] = "Nuclear in\nAdult Only"
x$cat = factor(x$cat, levels = c("Nuclear\nin Both","Nuclear in\nAdult Only","Not\nNuclear"))
x$Compartment = factor(x$Compartment, levels = c("Nuclear","Cytoplasmic"))
x

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/ENCODE_OR_plot.pdf", width = 7.5, height = 3.25)
ggplot(x, aes(cat, OddsRatio, fill = cat)) +
  geom_boxplot(outlier.shape = NA) + ylab("Odds Ratio") +  
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  xlab("") + facet_grid(. ~ Compartment) + geom_jitter() +
  geom_hline(yintercept=1, linetype="dotted") + ggtitle("ENCODE Cell Line Enrichment\nIn Disease-Associated Gene Sets") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.title=element_blank(),
        legend.position = "none", legend.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


## look for enrichment in nuclear-enriched genes (FDR<0.05 AND abs(LFC)>1) for all cell types and gene sets

sigres1 = lapply(sigres, function(x) lapply(x, function(y) y[which(abs(y$log2FoldChange)>=1),]))


inGroup = unlist(lapply(sigres1, function(x) lapply(x, function(y) as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(y)),"Symbol"]))))), recursive = F)
outGroup = unlist(lapply(sigres1, function(x) lapply(x, function(y) geneuniverse[!(geneuniverse %in% as.character(geneMap[which(geneMap$gencodeID %in% rownames(y)),"Symbol"]))])), recursive = F)

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$Gene.Symbol),sum(!(inG %in% x$Gene.Symbol)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$Gene.Symbol), sum(!(outG %in% x$Gene.Symbol)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)


enrich = lapply(enrich, data.frame)
enrich = do.call(rbind, Map(cbind, CellType = as.list(gsub("\\..*","",names(enrich))), Compartment = as.list(gsub("^.*\\.","",names(enrich))), lapply(enrich, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/github_controlled/other_datasets/data/Birnbaum_geneSet_enrichment_ENCODE_FractionDEGs_LFC1.csv")
bb1 = read.csv("./Dropbox/sorted_figures/github_controlled/other_datasets/data/Birnbaum_geneSet_enrichment_ENCODE_FractionDEGs_LFC1.csv")

bb1[bb1$CellType=="Fraction" & bb1$FDR<=0.05,]
#      X CellType Compartment GeneSet      P.Value OddsRatio         FDR
#240 240 Fraction     Nuclear    PGC2 0.0003558155  1.501109 0.009012085


bb1[bb1$CellType!="Fraction" & bb1$FDR<=0.05,]
#      X CellType Compartment           GeneSet      P.Value OddsRatio          FDR
#7     7  Gm12878 Cytoplasmic           SCZ.CNV 1.138842e-03 2.4488785 2.102477e-02
#22   22    Huvec Cytoplasmic      ASD.DATABASE 3.233504e-03 0.4323234 4.850256e-02
#27   27    Huvec Cytoplasmic           SCZ.CNV 2.128097e-03 2.2790003 3.404954e-02
#112 112     Nhek     Nuclear      ASD.DATABASE 9.451364e-05 1.9166235 4.536655e-03
#115 115     Nhek     Nuclear               NDD 1.065794e-03 3.7915148 2.102477e-02
#119 119     Nhek     Nuclear           SCZ.SNV 6.273450e-08 2.5549336 1.505628e-05
#166 166     Mcf7 Cytoplasmic Neurodegenerative 4.135252e-04 2.9636353 9.022369e-03
#174 174     Mcf7     Nuclear                ID 3.755036e-04 0.1822970 9.012085e-03
#191 191     A549     Nuclear           ASD.CNV 2.489870e-05 2.3634899 1.991896e-03
#197 197     A549     Nuclear           SCZ.CNV 1.877208e-04 2.5981471 6.436141e-03
#200 200     A549     Nuclear              PGC2 2.360451e-04 1.6984844 7.081353e-03
#211 211    Sknsh     Nuclear           ASD.CNV 9.279279e-05 2.0380465 4.536655e-03
#213 213    Sknsh     Nuclear         BPAD.GWAS 1.540505e-03 1.9542603 2.640866e-02
#217 217    Sknsh     Nuclear           SCZ.CNV 1.764940e-04 2.3251975 6.436141e-03
#220 220    Sknsh     Nuclear              PGC2 1.181839e-06 1.8051601 1.418207e-04  

b = bb[bb$CellType!="Fraction" & bb$FDR<=0.05,"X"]
b1 = bb1[bb1$CellType!="Fraction" & bb1$FDR<=0.05,"X"]
b1[-which(b1 %in% b)] # 174
nob = b[-which(b %in% b1)]

bb1[bb1$X==174,]
#      X CellType Compartment GeneSet      P.Value OddsRatio         FDR
#174 174     Mcf7     Nuclear      ID 0.0003755036  0.182297 0.009012085
bb[bb$X==174,]
#      X CellType Compartment GeneSet  P.Value OddsRatio FDR        Description
#174 174     Mcf7     Nuclear      ID 0.900371  1.019166   1 Mammary Epithelium

bb1[bb1$X %in% nob,]
#      X CellType Compartment           GeneSet     P.Value OddsRatio        FDR
#12   12  Gm12878     Nuclear      ASD.DATABASE 0.022864452 1.5138991 0.15242968
#23   23    Huvec Cytoplasmic         BPAD.GWAS 0.004943481 2.0584853 0.05741093
#149 149    Imr90 Cytoplasmic           SCZ.SNV 0.004072896 0.2832474 0.05144711
#160 160    Imr90     Nuclear              PGC2 0.005933291 1.5188819 0.06191261
#170 170     Mcf7 Cytoplasmic              PGC2 0.122364090 1.2078957 0.43831913
#171 171     Mcf7     Nuclear           ASD.CNV 0.022051726 1.5590826 0.15121183
#173 173     Mcf7     Nuclear         BPAD.GWAS 0.078484352 1.4832901 0.33636151
#177 177     Mcf7     Nuclear           SCZ.CNV 0.098573044 1.5024110 0.39429218
#192 192     A549     Nuclear      ASD.DATABASE 0.083488606 1.3699560 0.35153097
#193 193     A549     Nuclear         BPAD.GWAS 0.005612454 1.9555974 0.06122677
#199 199     A549     Nuclear           SCZ.SNV 0.007107099 1.6960846 0.06560399
#206 206    Sknsh Cytoplasmic Neurodegenerative 0.008580452 2.3342404 0.07230832

bb[bb$X %in% nob,]
#      X CellType Compartment           GeneSet      P.Value OddsRatio         FDR                               Description
#12   12  Gm12878     Nuclear      ASD.DATABASE 0.0026971384 1.6858121 0.032616962                              Lymphoblasts
#23   23    Huvec Cytoplasmic         BPAD.GWAS 0.0058060182 1.9951368 0.048548329          Umbilical Vein Endothelial Cells
#149 149    Imr90 Cytoplasmic           SCZ.SNV 0.0040728962 0.2832474 0.040728962                            Myofibroblasts
#160 160    Imr90     Nuclear              PGC2 0.0059332914 1.5188819 0.048548329                            Myofibroblasts
#170 170     Mcf7 Cytoplasmic              PGC2 0.0053202617 1.3612497 0.048548329                        Mammary Epithelium
#171 171     Mcf7     Nuclear           ASD.CNV 0.0057809174 1.6309852 0.048548329                        Mammary Epithelium
#173 173     Mcf7     Nuclear         BPAD.GWAS 0.0013887859 1.8954313 0.020659614                        Mammary Epithelium
#177 177     Mcf7     Nuclear           SCZ.CNV 0.0003045959 2.1565118 0.006091918                        Mammary Epithelium
#192 192     A549     Nuclear      ASD.DATABASE 0.0008113216 1.7174944 0.013584826 Adenocarcinomic Alveolar Basal Epithelial
#193 193     A549     Nuclear         BPAD.GWAS 0.0007968518 2.1037381 0.013584826 Adenocarcinomic Alveolar Basal Epithelial
#199 199     A549     Nuclear           SCZ.SNV 0.0000718357 1.9752339 0.001567324 Adenocarcinomic Alveolar Basal Epithelial
#206 206    Sknsh Cytoplasmic Neurodegenerative 0.0027180802 2.5028067 0.032616962                             Neuroblastoma




