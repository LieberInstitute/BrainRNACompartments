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

pdf("./Dropbox/sorted_figures/github_controlled/other_datasets/figures/ENCODE_MAplots.pdf", width=22, height=3)
ggplot(res, aes(x=baseMean/1000, y=log2FoldChange)) + 
  geom_point(aes(colour = factor(FDR))) + scale_colour_manual(values=c("red3","gray47")) +
  facet_grid(. ~ CellType) + geom_hline(aes(yintercept=0), linetype="dashed") +
  ylab("Log2 Fold Change") + ylim(-15,15) +
  xlab("Mean Normalized Counts/1000") + xlim(0,100) +
  ggtitle(paste0("MA Plot: Nucleus Vs. Cytoplasm")) + 
  theme(title = element_text(size = 16)) +
  theme(text = element_text(size = 16)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
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


# look for enrichment in nuclear-enriched genes (FDR<0.05) for all cell types and gene sets

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(res[[1]])),"Symbol"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$Gene.Symbol) %in% geneuniverse ), ] # drop genes that are not present in the test set
aej_sets_expressed$Gene.Symbol = as.character(aej_sets_expressed$Gene.Symbol)
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
inGroup = unlist(lapply(sigres, function(x) lapply(x, function(y) as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(y)),"Symbol"]))))), recursive = F)
outGroup = unlist(lapply(sigres, function(x) lapply(x, function(y) geneuniverse[!(geneuniverse %in% as.character(geneMap[which(geneMap$gencodeID %in% rownames(y)),"Symbol"]))])), recursive = F)


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
write.csv(enrich[,which(colnames(enrich)!="X")], quote=F, file="./Dropbox/sorted_figures/github_controlled/other_datasets/data/Birnbaum_geneSet_enrichment_ENCODE_FractionDEGs.csv")
bb = read.csv("./Dropbox/sorted_figures/github_controlled/other_datasets/data/Birnbaum_geneSet_enrichment_ENCODE_FractionDEGs.csv")


bb[bb$CellType=="Fraction" & bb$FDR<=0.05,]
#      X CellType Compartment      GeneSet      P.Value OddsRatio          FDR
#227 227 Fraction Cytoplasmic      SCZ.CNV 6.935834e-03  1.906897 0.0361869599
#229 229 Fraction Cytoplasmic SCZ.PGC.GWAS 2.053568e-03  2.044055 0.0154017633
#231 231 Fraction     Nuclear      ASD.CNV 1.168678e-04  1.886624 0.0017530173
#232 232 Fraction     Nuclear ASD.DATABASE 1.299575e-05  1.792615 0.0003465534
#233 233 Fraction     Nuclear    BPAD.GWAS 4.012418e-04  1.948994 0.0046795189
#235 235 Fraction     Nuclear          NDD 4.274167e-03  3.015421 0.0256450034
#237 237 Fraction     Nuclear      SCZ.CNV 2.946567e-03  1.846487 0.0186098970
#239 239 Fraction     Nuclear SCZ.PGC.GWAS 8.984007e-03  1.710841 0.0422776785
#240 240 Fraction     Nuclear      SCZ.SNV 3.655184e-04  1.694176 0.0046170744

bb[bb$CellType=="Fraction" & bb$Compartment=="Nuclear" & bb$FDR<=0.05,"GeneSet"]
# ASD.CNV, ASD.DATABASE, BPAD.GWAS, NDD, SCZ.SNV in nuclear
# SCZ.CNV, SCZ.PGC.GWAS in both compartments

range(bb[bb$CellType=="Fraction" & bb$Compartment=="Nuclear" & bb$GeneSet!="SCZ.CNV" & bb$GeneSet!="SCZ.PGC.GWAS" 
         & bb$FDR<=0.05,"OddsRatio"])
# 1.694176 3.015421
range(bb[bb$CellType=="Fraction" & bb$Compartment=="Nuclear" & bb$GeneSet!="SCZ.CNV" & bb$GeneSet!="SCZ.PGC.GWAS" 
         & bb$FDR<=0.05,"FDR"])
# 0.0003465534 0.0256450034
range(bb[bb$CellType=="Fraction" & bb$GeneSet %in% c("SCZ.CNV","SCZ.PGC.GWAS") & bb$FDR<=0.05,"OddsRatio"])
# 1.710841 2.044055
range(bb[bb$CellType=="Fraction" & bb$GeneSet %in% c("SCZ.CNV","SCZ.PGC.GWAS") & bb$FDR<=0.05,"FDR"])
# 0.01540176 0.04227768

bb[bb$CellType=="Fraction" & bb$FDR>0.05,]
#      X CellType Compartment           GeneSet    P.Value OddsRatio        FDR
#221 221 Fraction Cytoplasmic           ASD.CNV 0.04388209 1.4782612 0.14833384
#222 222 Fraction Cytoplasmic      ASD.DATABASE 0.01164650 0.5950604 0.05176222
#223 223 Fraction Cytoplasmic         BPAD.GWAS 0.26983949 1.2979531 0.50993290
#224 224 Fraction Cytoplasmic                ID 0.39686144 1.2410368 0.66049601
#225 225 Fraction Cytoplasmic               NDD 0.80832069 0.7417979 1.00000000
#226 226 Fraction Cytoplasmic Neurodegenerative 0.05145949 1.9021453 0.16440723
#228 228 Fraction Cytoplasmic SCZ.Meta.analysis 0.50504213 0.6025676 0.76846105
#230 230 Fraction Cytoplasmic           SCZ.SNV 0.45282288 0.8452208 0.71971848
#234 234 Fraction     Nuclear                ID 0.10307203 1.4517861 0.28764287
#236 236 Fraction     Nuclear Neurodegenerative 0.52191313 1.2286688 0.76846105
#238 238 Fraction     Nuclear SCZ.Meta.analysis 1.00000000 1.0133675 1.00000000

bb[bb$CellType!="Fraction" & bb$FDR<=0.05,]
#      X CellType Compartment           GeneSet      P.Value OddsRatio          FDR                               Description
#3     3  Gm12878 Cytoplasmic         BPAD.GWAS 2.268665e-03 2.1912290 1.555656e-02                              Lymphoblasts
#7     7  Gm12878 Cytoplasmic           SCZ.CNV 2.855340e-03 2.2699173 1.852112e-02                              Lymphoblasts
#11   11  Gm12878     Nuclear           ASD.CNV 6.896549e-03 1.8177740 3.618696e-02                              Lymphoblasts
#12   12  Gm12878     Nuclear      ASD.DATABASE 2.358840e-04 1.9101130 3.330127e-03                              Lymphoblasts
#23   23    Huvec Cytoplasmic         BPAD.GWAS 1.291297e-03 2.2527276 1.068660e-02          Umbilical Vein Endothelial Cells
#27   27    Huvec Cytoplasmic           SCZ.CNV 7.386960e-04 2.4885013 7.091482e-03          Umbilical Vein Endothelial Cells
#29   29    Huvec Cytoplasmic      SCZ.PGC.GWAS 8.540930e-03 2.1298471 4.099647e-02          Umbilical Vein Endothelial Cells
#51   51    Hepg2     Nuclear           ASD.CNV 7.864402e-03 1.7498413 3.932201e-02                               Hepatocytes
#87   87     K562 Cytoplasmic           SCZ.CNV 4.889415e-03 2.3107631 2.793951e-02                      Myelogenous Leukemia
#109 109     Nhek Cytoplasmic      SCZ.PGC.GWAS 1.597850e-03 2.3812111 1.278280e-02                             Keratinocytes
#112 112     Nhek     Nuclear      ASD.DATABASE 2.047846e-06 2.1927572 7.108360e-05                             Keratinocytes
#114 114     Nhek     Nuclear                ID 5.202552e-03 2.1109897 2.903750e-02                             Keratinocytes
#115 115     Nhek     Nuclear               NDD 5.405785e-04 4.1530196 5.640819e-03                             Keratinocytes
#120 120     Nhek     Nuclear           SCZ.SNV 7.672816e-10 2.8620782 9.207380e-08                             Keratinocytes
#149 149    Imr90 Cytoplasmic      SCZ.PGC.GWAS 3.449386e-04 2.9544361 4.599181e-03                            Myofibroblasts
#150 150    Imr90 Cytoplasmic           SCZ.SNV 1.078405e-02 0.3152731 4.977255e-02                            Myofibroblasts
#153 153    Imr90     Nuclear         BPAD.GWAS 6.043496e-03 2.0360327 3.296452e-02                            Myofibroblasts
#159 159    Imr90     Nuclear      SCZ.PGC.GWAS 7.986704e-04 2.4648729 7.372342e-03                            Myofibroblasts
#162 162     Mcf7 Cytoplasmic      ASD.DATABASE 7.862166e-03 1.4512202 3.932201e-02                        Mammary Epithelium
#163 163     Mcf7 Cytoplasmic         BPAD.GWAS 4.710757e-03 1.7399317 2.757516e-02                        Mammary Epithelium
#164 164     Mcf7 Cytoplasmic                ID 1.026951e-03 2.0698543 8.802433e-03                        Mammary Epithelium
#166 166     Mcf7 Cytoplasmic Neurodegenerative 1.909261e-06 4.2371583 7.108360e-05                        Mammary Epithelium
#169 169     Mcf7 Cytoplasmic      SCZ.PGC.GWAS 6.491040e-05 2.2795566 1.298208e-03                        Mammary Epithelium
#170 170     Mcf7 Cytoplasmic           SCZ.SNV 8.352313e-03 1.4881046 4.090929e-02                        Mammary Epithelium
#171 171     Mcf7     Nuclear           ASD.CNV 5.369691e-04 1.8570119 5.640819e-03                        Mammary Epithelium
#173 173     Mcf7     Nuclear         BPAD.GWAS 9.846163e-05 2.1701057 1.687914e-03                        Mammary Epithelium
#177 177     Mcf7     Nuclear           SCZ.CNV 3.077517e-05 2.4602293 6.714582e-04                        Mammary Epithelium
#180 180     Mcf7     Nuclear           SCZ.SNV 1.663127e-03 1.6494239 1.287582e-02                        Mammary Epithelium
#189 189     A549 Cytoplasmic      SCZ.PGC.GWAS 2.507873e-03 1.9791729 1.671915e-02 Adenocarcinomic Alveolar Basal Epithelial
#191 191     A549     Nuclear           ASD.CNV 4.610656e-10 3.1161131 9.207380e-08 Adenocarcinomic Alveolar Basal Epithelial
#192 192     A549     Nuclear      ASD.DATABASE 2.941906e-05 1.9621266 6.714582e-04 Adenocarcinomic Alveolar Basal Epithelial
#193 193     A549     Nuclear         BPAD.GWAS 8.197559e-05 2.3868926 1.513396e-03 Adenocarcinomic Alveolar Basal Epithelial
#197 197     A549     Nuclear           SCZ.CNV 2.408413e-08 3.5189242 1.926730e-06 Adenocarcinomic Alveolar Basal Epithelial
#199 199     A549     Nuclear      SCZ.PGC.GWAS 2.182599e-03 2.1286607 1.555656e-02 Adenocarcinomic Alveolar Basal Epithelial
#200 200     A549     Nuclear           SCZ.SNV 3.462616e-06 2.2070325 1.038785e-04 Adenocarcinomic Alveolar Basal Epithelial
#202 202    Sknsh Cytoplasmic      ASD.DATABASE 3.750254e-03 1.5353344 2.307849e-02                             Neuroblastoma
#204 204    Sknsh Cytoplasmic                ID 2.220974e-03 2.0603960 1.555656e-02                             Neuroblastoma
#206 206    Sknsh Cytoplasmic Neurodegenerative 7.364629e-04 2.9127229 7.091482e-03                             Neuroblastoma
#209 209    Sknsh Cytoplasmic      SCZ.PGC.GWAS 1.150541e-04 2.2870572 1.753017e-03                             Neuroblastoma
#211 211    Sknsh     Nuclear           ASD.CNV 1.819176e-07 2.5266929 1.091505e-05                             Neuroblastoma
#213 213    Sknsh     Nuclear         BPAD.GWAS 4.094579e-04 2.0965476 4.679519e-03                             Neuroblastoma
#217 217    Sknsh     Nuclear           SCZ.CNV 2.073272e-06 2.8619126 7.108360e-05                             Neuroblastoma
#219 219    Sknsh     Nuclear      SCZ.PGC.GWAS 9.749406e-04 2.1275973 8.666139e-03                             Neuroblastoma

dim(bb[bb$FDR<=0.05,]) # 52
dim(bb[bb$FDR<=0.05 & bb$OddsRatio<1,]) # 1
bb[bb$FDR<=0.05 & bb$OddsRatio<1,]
#150 150    Imr90 Cytoplasmic SCZ.SNV 0.01078405 0.3152731 0.04977255 Myofibroblasts

x = split(bb[bb$FDR<=0.05,], bb[bb$FDR<=0.05,"GeneSet"])
x = lapply(x, function(y) split(y, y$CellType))
do.call(rbind, Map(cbind, GeneSet = as.list(names(lapply(x, elementNROWS))), lapply(x, function(y) data.frame(rbind(elementNROWS(y))))))
#             GeneSet A549 Fraction Gm12878 H1hesc Helas3 Hepg2 Huvec Imr90 K562 Mcf7 Nhek Sknsh
#1            ASD.CNV    1        1       1      0      0     1     0     0    0    1    0     1
#2       ASD.DATABASE    1        1       1      0      0     0     0     0    0    1    1     1
#3          BPAD.GWAS    1        1       1      0      0     0     1     1    0    2    0     1
#4                 ID    0        0       0      0      0     0     0     0    0    1    1     1
#5                NDD    0        1       0      0      0     0     0     0    0    0    1     0
#6  Neurodegenerative    0        0       0      0      0     0     0     0    0    1    0     1
#7            SCZ.CNV    1        2       1      0      0     0     1     0    1    1    0     1
#8  SCZ.Meta.analysis    0        0       0      0      0     0     0     0    0    0    0     0
#9       SCZ.PGC.GWAS    2        2       0      0      0     0     1     2    0    1    1     2
#10           SCZ.SNV    1        1       0      0      0     0     0     1    0    2    1     0

x = split(bb[bb$FDR<=0.05 & bb$Compartment=="Nuclear",], bb[bb$FDR<=0.05 & bb$Compartment=="Nuclear","GeneSet"])
x = lapply(x, function(y) split(y, y$CellType))
do.call(rbind, Map(cbind, GeneSet = as.list(names(lapply(x, elementNROWS))), lapply(x, function(y) data.frame(rbind(elementNROWS(y))))))
#             GeneSet A549 Fraction Gm12878 H1hesc Helas3 Hepg2 Huvec Imr90 K562 Mcf7 Nhek Sknsh
#1            ASD.CNV    1        1       1      0      0     1     0     0    0    1    0     1
#2       ASD.DATABASE    1        1       1      0      0     0     0     0    0    0    1     0
#3          BPAD.GWAS    1        1       0      0      0     0     0     1    0    1    0     1
#4                 ID    0        0       0      0      0     0     0     0    0    0    1     0
#5                NDD    0        1       0      0      0     0     0     0    0    0    1     0
#6  Neurodegenerative    0        0       0      0      0     0     0     0    0    0    0     0
#7            SCZ.CNV    1        1       0      0      0     0     0     0    0    1    0     1
#8  SCZ.Meta.analysis    0        0       0      0      0     0     0     0    0    0    0     0
#9       SCZ.PGC.GWAS    1        1       0      0      0     0     0     1    0    0    0     1
#10           SCZ.SNV    1        1       0      0      0     0     0     0    0    1    1     0

c("ASD.CNV","ASD.DATABASE","BPAD.GWAS","ID","NDD","Neurodegenerative","SCZ.CNV","SCZ.Meta.analysis","SCZ.PGC.GWAS","SCZ.SNV")
range(bb[bb$FDR<=0.05 & bb$CellType=="A549" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE","BPAD.GWAS","SCZ.CNV","SCZ.Meta.analysis","SCZ.PGC.GWAS","SCZ.SNV"),"OddsRatio"])
# 1.962127 3.518924
range(bb[bb$FDR<=0.05 & bb$CellType=="A549" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE","BPAD.GWAS","SCZ.CNV","SCZ.Meta.analysis","SCZ.PGC.GWAS","SCZ.SNV"),"FDR"])
# 9.207380e-08 1.555656e-02
range(bb[bb$FDR<=0.05 & bb$CellType %in% c("Gm12878","Hepg2") & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE"),"OddsRatio"])
# 1.749841 1.910113
range(bb[bb$FDR<=0.05 & bb$CellType %in% c("Gm12878","Hepg2") & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE"),"FDR"])
# 0.003330127 0.039322008
bb[bb$FDR<=0.05 & bb$CellType=="Imr90" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("BPAD.GWAS"),]

range(bb[bb$FDR<=0.05 & bb$CellType %in% c("Mcf7","Nhek") & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE","SCZ.CNV","SCZ.Meta.analysis","SCZ.PGC.GWAS","SCZ.SNV"),"OddsRatio"])
# 1.649424 2.862078
range(bb[bb$FDR<=0.05 & bb$CellType %in% c("Mcf7","Nhek") & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE","SCZ.CNV","SCZ.Meta.analysis","SCZ.PGC.GWAS","SCZ.SNV"),"FDR"])
# 9.207380e-08 1.287582e-02
range(bb[bb$FDR<=0.05 & bb$CellType=="Sknsh" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE","BPAD.GWAS","SCZ.CNV","SCZ.Meta.analysis","SCZ.SNV"),"OddsRatio"])
# 2.096548 2.861913
range(bb[bb$FDR<=0.05 & bb$CellType=="Sknsh" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ASD.CNV","ASD.DATABASE","BPAD.GWAS","SCZ.CNV","SCZ.Meta.analysis","SCZ.SNV"),"FDR"])
# 1.091505e-05 4.679519e-03
range(bb[bb$FDR<=0.05 & bb$CellType=="Nhek" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ID","NDD","Neurodegenerative"),"OddsRatio"])
# 2.11099 4.15302
range(bb[bb$FDR<=0.05 & bb$CellType=="Nhek" & bb$Compartment=="Nuclear" & bb$GeneSet %in% c("ID","NDD","Neurodegenerative"),"FDR"])
# 0.005640819 0.029037500

range(bb[bb$FDR<=0.05 & bb$CellType %in% c("Mcf7","Sknsh") & bb$Compartment=="Cytoplasmic" & bb$GeneSet %in% c("ID","NDD","Neurodegenerative"),"OddsRatio"])
# 2.060396 4.237158
range(bb[bb$FDR<=0.05 & bb$CellType %in% c("Mcf7","Sknsh") & bb$Compartment=="Cytoplasmic" & bb$GeneSet %in% c("ID","NDD","Neurodegenerative"),"FDR"])
# 0.0000710836 0.0155565621

                
# look for enrichment in nuclear-enriched genes (FDR<0.05 AND abs(LFC)>1) for all cell types and gene sets

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
# none

bb1[bb1$CellType!="Fraction" & bb1$FDR<=0.05,]
#      X CellType Compartment           GeneSet      P.Value OddsRatio          FDR
#3     3  Gm12878 Cytoplasmic         BPAD.GWAS 1.639030e-03 2.3656619 1.966836e-02
#7     7  Gm12878 Cytoplasmic           SCZ.CNV 2.778852e-04 2.7535543 7.410272e-03
#11   11  Gm12878     Nuclear           ASD.CNV 3.252017e-03 1.9211447 3.121937e-02
#12   12  Gm12878     Nuclear      ASD.DATABASE 4.498024e-03 1.7131525 3.811649e-02
#23   23    Huvec Cytoplasmic         BPAD.GWAS 9.833657e-04 2.3234472 1.573385e-02
#27   27    Huvec Cytoplasmic           SCZ.CNV 5.583409e-04 2.5664733 1.030783e-02
#29   29    Huvec Cytoplasmic      SCZ.PGC.GWAS 4.605743e-03 2.1964418 3.811649e-02
#87   87     K562 Cytoplasmic           SCZ.CNV 2.448267e-03 2.3847960 2.554713e-02
#109 109     Nhek Cytoplasmic      SCZ.PGC.GWAS 4.159897e-03 2.2326750 3.811649e-02
#112 112     Nhek     Nuclear      ASD.DATABASE 3.610057e-06 2.1816086 2.888045e-04
#115 115     Nhek     Nuclear               NDD 3.889998e-04 4.3357383 8.487269e-03
#120 120     Nhek     Nuclear           SCZ.SNV 7.479022e-10 2.9159805 1.794965e-07
#149 149    Imr90 Cytoplasmic      SCZ.PGC.GWAS 3.449386e-04 2.9544361 8.278526e-03
#153 153    Imr90     Nuclear         BPAD.GWAS 6.043496e-03 2.0360327 4.834797e-02
#159 159    Imr90     Nuclear      SCZ.PGC.GWAS 7.986704e-04 2.4648729 1.369149e-02
#166 166     Mcf7 Cytoplasmic Neurodegenerative 7.398176e-05 3.4274677 2.536517e-03
#169 169     Mcf7 Cytoplasmic      SCZ.PGC.GWAS 6.299519e-03 1.8826591 4.877047e-02
#171 171     Mcf7     Nuclear           ASD.CNV 4.323802e-03 1.7583979 3.811649e-02
#174 174     Mcf7     Nuclear                ID 1.222969e-03 0.2090022 1.655182e-02
#189 189     A549 Cytoplasmic      SCZ.PGC.GWAS 2.108961e-03 2.2287057 2.300685e-02
#191 191     A549     Nuclear           ASD.CNV 2.302351e-06 2.6532092 2.762821e-04
#193 193     A549     Nuclear         BPAD.GWAS 1.241386e-03 2.2081535 1.655182e-02
#197 197     A549     Nuclear           SCZ.CNV 2.050530e-05 2.9249160 9.065772e-04
#200 200     A549     Nuclear           SCZ.SNV 1.081210e-03 1.8892447 1.621815e-02
#204 204    Sknsh Cytoplasmic                ID 3.127997e-03 2.0402252 3.121937e-02
#206 206    Sknsh Cytoplasmic Neurodegenerative 1.787327e-03 2.7003548 2.042659e-02
#209 209    Sknsh Cytoplasmic      SCZ.PGC.GWAS 1.593387e-03 2.0421020 1.966836e-02
#211 211    Sknsh     Nuclear           ASD.CNV 5.771018e-06 2.3012873 3.462611e-04
#213 213    Sknsh     Nuclear         BPAD.GWAS 2.058806e-04 2.2201534 6.176418e-03
#217 217    Sknsh     Nuclear           SCZ.CNV 2.266443e-05 2.6331005 9.065772e-04
#219 219    Sknsh     Nuclear      SCZ.PGC.GWAS 5.486993e-04 2.2233295 1.030783e-02


b = bb[bb$CellType!="Fraction" & bb$FDR<=0.05,"X"]
b1 = bb1[bb1$CellType!="Fraction" & bb1$FDR<=0.05,"X"]
b1[-which(b1 %in% b)] # 174
nob = b[-which(b %in% b1)]

bb1[bb1$X==174,]
#      X CellType Compartment GeneSet     P.Value OddsRatio        FDR
#174 174     Mcf7     Nuclear      ID 0.001222969 0.2090022 0.01655182
bb[bb$X==174,]
#      X CellType Compartment GeneSet   P.Value OddsRatio       FDR
#174 174     Mcf7     Nuclear      ID 0.5119175  1.182585 0.7684611

bb1[bb1$X %in% nob,]
#      X CellType Compartment GeneSet     P.Value OddsRatio        FDR
#51   51    Hepg2     Nuclear      ASD.CNV 0.01058025 1.7193358 0.07612272
#114 114     Nhek     Nuclear           ID 0.01897759 1.9259030 0.10081739
#150 150    Imr90 Cytoplasmic      SCZ.SNV 0.01078405 0.3152731 0.07612272
#162 162     Mcf7 Cytoplasmic ASD.DATABASE 0.55442901 1.1061383 0.88121167
#163 163     Mcf7 Cytoplasmic    BPAD.GWAS 0.15163668 1.3846428 0.43324766
#164 164     Mcf7 Cytoplasmic           ID 0.78343267 1.0643393 1.00000000
#170 170     Mcf7 Cytoplasmic      SCZ.SNV 0.41110567 0.8405873 0.73841194
#173 173     Mcf7     Nuclear    BPAD.GWAS 0.02469570 1.6820900 0.11621504
#177 177     Mcf7     Nuclear      SCZ.CNV 0.03153466 1.6999787 0.14015404
#180 180     Mcf7     Nuclear      SCZ.SNV 0.91987877 0.9607012 1.00000000
#192 192     A549     Nuclear ASD.DATABASE 0.01821790 1.5566359 0.10081739
#199 199     A549     Nuclear SCZ.PGC.GWAS 0.02108687 1.8726343 0.10635032
#202 202    Sknsh Cytoplasmic ASD.DATABASE 0.01334259 1.4690118 0.08654652

bb[bb$X %in% nob,]
#      X CellType Compartment      GeneSet      P.Value OddsRatio          FDR
#51   51    Hepg2     Nuclear      ASD.CNV 7.864402e-03 1.7498413 0.0393220081
#114 114     Nhek     Nuclear           ID 5.202552e-03 2.1109897 0.0290375001
#150 150    Imr90 Cytoplasmic      SCZ.SNV 1.078405e-02 0.3152731 0.0497725455
#162 162     Mcf7 Cytoplasmic ASD.DATABASE 7.862166e-03 1.4512202 0.0393220081
#163 163     Mcf7 Cytoplasmic    BPAD.GWAS 4.710757e-03 1.7399317 0.0275751650
#164 164     Mcf7 Cytoplasmic           ID 1.026951e-03 2.0698543 0.0088024331
#170 170     Mcf7 Cytoplasmic      SCZ.SNV 8.352313e-03 1.4881046 0.0409092904
#173 173     Mcf7     Nuclear    BPAD.GWAS 9.846163e-05 2.1701057 0.0016879137
#177 177     Mcf7     Nuclear      SCZ.CNV 3.077517e-05 2.4602293 0.0006714582
#180 180     Mcf7     Nuclear      SCZ.SNV 1.663127e-03 1.6494239 0.0128758182
#192 192     A549     Nuclear ASD.DATABASE 2.941906e-05 1.9621266 0.0006714582
#199 199     A549     Nuclear SCZ.PGC.GWAS 2.182599e-03 2.1286607 0.0155565621
#202 202    Sknsh Cytoplasmic ASD.DATABASE 3.750254e-03 1.5353344 0.0230784862





