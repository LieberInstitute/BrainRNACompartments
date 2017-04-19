library(DESeq2)
library(VennDiagram)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

## Does the pattern of fewer DE genes by fraction in Prenatal hold when measuring the last exon only?
# extract the last exon from each transcript

head(exonMap)
exonPlus = exonMap[which(exonMap$Strand=="+"),]
exonMinus = exonMap[which(exonMap$Strand=="-"),]
exonPlus$Num = seq_len(nrow(exonPlus))
exonPlus = exonPlus[order(exonPlus$Num, decreasing = T),]
firstPlus = exonPlus[match(unique(exonPlus$gencodeTx), exonPlus$gencodeTx),]
exonMinus$Num = seq_len(nrow(exonMinus))
exonMinus = exonMinus[order(exonMinus$Num, decreasing = T),]
firstMinus = exonMinus[match(unique(exonMinus$Geneid), exonMinus$Geneid),]
LastExonMap = rbind(firstPlus, firstMinus)
LastExonMap = LastExonMap[which(LastExonMap$gencodeID %in% rownames(geneCounts)),]
lastexonCounts = exonCounts[which(rownames(exonCounts) %in% rownames(LastExonMap) & rowSums(exonCounts)>0),]
LastExonMap = LastExonMap[which(rownames(LastExonMap) %in% rownames(lastexonCounts)),]

LastExonMap.down = rbind(firstPlus, firstMinus)
LastExonMap.down = LastExonMap.down[which(LastExonMap.down$gencodeID %in% rownames(geneCounts.down)),]
lastexonCounts.down = exonCounts.down[which(rownames(exonCounts.down) %in% rownames(LastExonMap.down) & rowSums(exonCounts.down)>0),]
LastExonMap.down = LastExonMap.down[which(rownames(LastExonMap.down) %in% rownames(lastexonCounts.down)),]

# DESeq2 by fraction in four groups
Ap.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Adult" & pd$Library=="polyA"),"SampleID"])],
                                   colData = pd[which(pd$Fetal=="Adult" & pd$Library=="polyA"),], design = ~ Zone)
Ap.ddsLE <- DESeq(Ap.ddsLE)
ApresLE <- results(Ap.ddsLE)
Fp.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),"SampleID"])],
                                   colData = pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),], design = ~ Zone)
Fp.ddsLE <- DESeq(Fp.ddsLE)
FpresLE <- results(Fp.ddsLE)

Fp.ddsLE.down <- DESeqDataSetFromMatrix(countData = lastexonCounts.down[,which(colnames(lastexonCounts.down) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),"SampleID"])],
                                   colData = pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),], design = ~ Zone)
Fp.ddsLE.down <- DESeq(Fp.ddsLE.down)
FpresLE.down <- results(Fp.ddsLE.down)

Ar.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Adult" & pd$Library=="RiboZero"),"SampleID"])],
                                   colData = pd[which(pd$Fetal=="Adult" & pd$Library=="RiboZero"),], design = ~ Zone)
Ar.ddsLE <- DESeq(Ar.ddsLE)
ArresLE <- results(Ar.ddsLE)
Fr.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="RiboZero"),"SampleID"])],
                                   colData = pd[which(pd$Fetal=="Prenatal" & pd$Library=="RiboZero"),], design = ~ Zone)
Fr.ddsLE <- DESeq(Fr.ddsLE)
FrresLE <- results(Fr.ddsLE)
save(ApresLE, FpresLE,FpresLE.down,ArresLE,FrresLE, 
     file = "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/last.exon.DESeq2.results.rda")

# MA Plots
plotMA(ApresLE, alpha=0.05, main="Adult PolyA Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(FpresLE, alpha=0.05, main="Prenatal PolyA Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(FpresLE.down, alpha=0.05, main="Prenatal PolyA Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(ArresLE, alpha=0.05, main="Adult Ribozero Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(FrresLE, alpha=0.05, main="Prenatal Ribozero Samples (Last Exon Only)", ylim=c(-3,3))

## Are similar numbers of genes expressed in the four groups? 
  # 134806 tested transcripts
dim(lastexonCounts[which(rowMeans(lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Adult" & pd$Library=="polyA"),"SampleID"])]) > 0),
                   which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Adult" & pd$Library=="polyA"),"SampleID"])]) 
# 121779
dim(lastexonCounts[which(rowMeans(lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),"SampleID"])]) > 0),
                   which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),"SampleID"])]) 
# 123404
dim(lastexonCounts.down[which(rowMeans(lastexonCounts.down[,which(colnames(lastexonCounts.down) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),"SampleID"])]) > 0),
                   which(colnames(lastexonCounts.down) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="polyA"),"SampleID"])]) 
# 119845 
dim(lastexonCounts[which(rowMeans(lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Adult" & pd$Library=="RiboZero"),"SampleID"])]) > 0),
                   which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Adult" & pd$Library=="RiboZero"),"SampleID"])]) 
# 122112
dim(lastexonCounts[which(rowMeans(lastexonCounts[,which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="RiboZero"),"SampleID"])]) > 0),
                   which(colnames(lastexonCounts) %in% pd[which(pd$Fetal=="Prenatal" & pd$Library=="RiboZero"),"SampleID"])]) 
# 119977

## What are the DE genes?
# Make list of significant genes
LEResults = list(AP = as.data.frame(ApresLE), FP = as.data.frame(FpresLE),
                 AR = as.data.frame(ArresLE), FR = as.data.frame(FrresLE), 
                 FP.down = as.data.frame(FpresLE.down))
sigLE = lapply(LEResults, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(sigLE)
    #AP      FP      AR      FR FP.down 
    #3812     142    3932     287      79
Sign = lapply(sigLE, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigLE = Map(cbind, sigLE, Sign = Sign)
sigLE = lapply(sigLE, function(x) split(x, x$Sign))
sigLE = unlist(sigLE, recursive = F)
sigLE = Map(cbind, sigLE, lapply(sigLE, function(x) exonMap[match(rownames(x), rownames(exonMap)),]))
sigGenesLE = lapply(sigLE, function(x) unique(x$gencodeID))
elementNROWS(sigGenesLE)

# Compare to   
FracList = list(Apres = as.data.frame(Apres), Fpres = as.data.frame(Fpres),
                Arres = as.data.frame(Arres),Frres = as.data.frame(Frres), 
                Fpres.down = as.data.frame(Fpres.down))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
    #Apres      Fpres      Arres      Frres Fpres.down 
    #1894         52       1892         30         40 
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
SigList = Map(cbind, SigFracList, Sign = Sign)
SigList = lapply(SigList, function(x) split(x, x$Sign))
SigList = unlist(SigList, recursive=F)
SigList = lapply(SigList, function(x) rownames(x))
elementNROWS(SigList)
    #Apres.DownNuc        Apres.UpNuc      Fpres.DownNuc        Fpres.UpNuc      Arres.DownNuc        Arres.UpNuc 
    #938                956                  1                 51               1024                868 
    #Frres.DownNuc        Frres.UpNuc Fpres.down.DownNuc   Fpres.down.UpNuc 
    #7                 23                  1                 39 

# Create Venn Diagrams of the overlap of the genes DE in the whole gene and by last exon (LFC≥1, FDR<0.05)
LE.GeneP.UP <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.UpNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.UpNuc"]], 
                     "Adult:PolyA\nGene"=SigList[["Apres.UpNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.UpNuc"]])
LE.GeneR.UP <- list("Adult:Ribozero\nLast Exon"=sigGenesLE[["AR.UpNuc"]], "Prenatal:Ribozero\nLast Exon"=sigGenesLE[["FR.UpNuc"]], 
                    "Adult:Ribozero\nGene"=SigList[["Arres.UpNuc"]], "Prenatal:Ribozero\nGene"=SigList[["Frres.UpNuc"]])
LE.GeneP.D <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.DownNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.DownNuc"]], 
                    "Adult:PolyA\nGene"=SigList[["Apres.DownNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.DownNuc"]])
LE.GeneR.D <- list("Adult:Ribozero\nLast Exon"=sigGenesLE[["AR.DownNuc"]], "Prenatal:Ribozero\nLast Exon"=sigGenesLE[["FR.DownNuc"]], 
                    "Adult:Ribozero\nGene"=SigList[["Arres.DownNuc"]], "Prenatal:Ribozero\nGene"=SigList[["Frres.DownNuc"]])
venn.LEvsGeneP.UP <- venn.diagram(LE.GeneP.UP, 
                                  "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Nuclear_polya.jpeg", 
                                  main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                            col = "transparent",
                            fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                            alpha = 0.50,
                            label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                          "white", "white", "white", "white", "palevioletred4", "white",
                                          "white", "white", "white", "darkblue", "white"),
                            fontfamily = "Arial",
                            fontface = "bold",
                            cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                            cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneP.D <- venn.diagram(LE.GeneP.D, 
                                 "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Cytosolic_polya.jpeg", 
                                 main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                           col = "transparent",
                           fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                           alpha = 0.50,
                           label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                         "white", "white", "white", "white", "palevioletred4", "white",
                                         "white", "white", "white", "darkblue", "white"),
                           fontfamily = "Arial",
                           fontface = "bold",
                           cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                           cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneR.UP <- venn.diagram(LE.GeneR.UP, 
                                  "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Nuclear_ribo.jpeg", 
                                  main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                  col = "transparent",
                                  fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                                  alpha = 0.50,
                                  label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                                "white", "white", "white", "white", "palevioletred4", "white",
                                                "white", "white", "white", "darkblue", "white"),
                                  fontfamily = "Arial",
                                  fontface = "bold",
                                  cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                  cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneR.D <- venn.diagram(LE.GeneR.D, 
                                 "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Cytosolic_ribo.jpeg", 
                                 main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                 col = "transparent",
                                 fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                                 alpha = 0.50,
                                 label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                               "white", "white", "white", "white", "palevioletred4", "white",
                                               "white", "white", "white", "darkblue", "white"),
                                 fontfamily = "Arial",
                                 fontface = "bold",
                                 cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                 cat.fontfamily = "Arial", margin=0.2)

LE.GeneP.UP.down <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.UpNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.down.UpNuc"]], 
                    "Adult:PolyA\nGene"=SigList[["Apres.UpNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.UpNuc"]])
LE.GeneP.D.down <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.DownNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.down.DownNuc"]], 
                   "Adult:PolyA\nGene"=SigList[["Apres.DownNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.DownNuc"]])
venn.LEvsGeneP.UP <- venn.diagram(LE.GeneP.UP.down, 
                                  "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Nuclear_polya_downsampled.jpeg", 
                                  main="Differentially Expressed Genes\nUp-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                  col = "transparent",fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
                                  label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                                "white", "white", "white", "white", "palevioletred4", "white",
                                                "white", "white", "white", "darkblue", "white"),
                                  fontfamily = "Arial",fontface = "bold",cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                  cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneP.D <- venn.diagram(LE.GeneP.D.down, 
                                 "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Cytosolic_polya_downsampled.jpeg", 
                                 main="Differentially expressed Genes\nDown-regulated in Nucleus (LFC≥1, FDR<0.05)",
                                 col = "transparent",fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
                                 label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                               "white", "white", "white", "white", "palevioletred4", "white",
                                               "white", "white", "white", "darkblue", "white"),
                                 fontfamily = "Arial",fontface = "bold",cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                 cat.fontfamily = "Arial", margin=0.2)

# Find overlaps just for P<0.05, no LFC criteria for Gene-level expression (FDR<0.05)
sigLE = lapply(LEResults, function(x) x[which(x$padj<=0.05),])
elementNROWS(sigLE)
Sign = lapply(sigLE, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigLE = Map(cbind, sigLE, Sign = Sign)
sigLE = lapply(sigLE, function(x) split(x, x$Sign))
sigLE = unlist(sigLE, recursive = F)
sigLE = Map(cbind, sigLE, lapply(sigLE, function(x) exonMap[match(rownames(x), rownames(exonMap)),]))
sigGenesLE = lapply(sigLE, function(x) unique(x$gencodeID))
elementNROWS(sigGenesLE)

# Compare to   
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
SigList = Map(cbind, SigFracList, Sign = Sign)
SigList = lapply(SigList, function(x) split(x, x$Sign))
SigList = unlist(SigList, recursive=F)
SigList = lapply(SigList, function(x) rownames(x))
elementNROWS(SigList)

LE.GeneP.UP <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.UpNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.UpNuc"]], 
                    "Adult:PolyA\nGene"=SigList[["Apres.UpNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.UpNuc"]])
LE.GeneR.UP <- list("Adult:Ribozero\nLast Exon"=sigGenesLE[["AR.UpNuc"]], "Prenatal:Ribozero\nLast Exon"=sigGenesLE[["FR.UpNuc"]], 
                    "Adult:Ribozero\nGene"=SigList[["Arres.UpNuc"]], "Prenatal:Ribozero\nGene"=SigList[["Frres.UpNuc"]])
LE.GeneP.D <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.DownNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.DownNuc"]], 
                   "Adult:PolyA\nGene"=SigList[["Apres.DownNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.DownNuc"]])
LE.GeneR.D <- list("Adult:Ribozero\nLast Exon"=sigGenesLE[["AR.DownNuc"]], "Prenatal:Ribozero\nLast Exon"=sigGenesLE[["FR.DownNuc"]], 
                   "Adult:Ribozero\nGene"=SigList[["Arres.DownNuc"]], "Prenatal:Ribozero\nGene"=SigList[["Frres.DownNuc"]])
venn.LEvsGeneP.UP <- venn.diagram(LE.GeneP.UP, 
                                  "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Nuclear_polya_noLFC.jpeg", 
                                  main="Differentially Expressed Genes\nUp-regulated in Nucleus (FDR<0.05)",
                                  col = "transparent",
                                  fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                                  alpha = 0.50,
                                  label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                                "white", "white", "white", "white", "palevioletred4", "white",
                                                "white", "white", "white", "darkblue", "white"),
                                  fontfamily = "Arial",
                                  fontface = "bold",
                                  cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                  cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneP.D <- venn.diagram(LE.GeneP.D, 
                                 "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Cytosolic_polya_noLFC.jpeg", 
                                 main="Differentially expressed Genes\nDown-regulated in Nucleus (FDR<0.05)",
                                 col = "transparent",
                                 fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                                 alpha = 0.50,
                                 label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                               "white", "white", "white", "white", "palevioletred4", "white",
                                               "white", "white", "white", "darkblue", "white"),
                                 fontfamily = "Arial",
                                 fontface = "bold",
                                 cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                 cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneR.UP <- venn.diagram(LE.GeneR.UP, 
                                  "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Nuclear_ribo_noLFC.jpeg", 
                                  main="Differentially Expressed Genes\nUp-regulated in Nucleus (FDR<0.05)",
                                  col = "transparent",
                                  fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                                  alpha = 0.50,
                                  label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                                "white", "white", "white", "white", "palevioletred4", "white",
                                                "white", "white", "white", "darkblue", "white"),
                                  fontfamily = "Arial",
                                  fontface = "bold",
                                  cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                  cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneR.D <- venn.diagram(LE.GeneR.D, 
                                 "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Cytosolic_ribo_noLFC.jpeg", 
                                 main="Differentially expressed Genes\nDown-regulated in Nucleus (FDR<0.05)",
                                 col = "transparent",
                                 fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                                 alpha = 0.50,
                                 label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                               "white", "white", "white", "white", "palevioletred4", "white",
                                               "white", "white", "white", "darkblue", "white"),
                                 fontfamily = "Arial",
                                 fontface = "bold",
                                 cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                 cat.fontfamily = "Arial", margin=0.2)

LE.GeneP.UP.down <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.UpNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.down.UpNuc"]], 
                         "Adult:PolyA\nGene"=SigList[["Apres.UpNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.UpNuc"]])
LE.GeneP.D.down <- list("Adult:PolyA\nLast Exon"=sigGenesLE[["AP.DownNuc"]], "Prenatal:PolyA\nLast Exon"=sigGenesLE[["FP.down.DownNuc"]], 
                        "Adult:PolyA\nGene"=SigList[["Apres.DownNuc"]], "Prenatal:PolyA\nGene"=SigList[["Fpres.DownNuc"]])
venn.LEvsGeneP.UP <- venn.diagram(LE.GeneP.UP.down, 
                                  "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Nuclear_polya_noLFC_downsampled.jpeg", 
                                  main="Differentially Expressed Genes\nUp-regulated in Nucleus (FDR<0.05)",
                                  col = "transparent",fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
                                  label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                                "white", "white", "white", "white", "palevioletred4", "white",
                                                "white", "white", "white", "darkblue", "white"),
                                  fontfamily = "Arial",fontface = "bold",cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                  cat.fontfamily = "Arial", margin=0.2)
venn.LEvsGeneP.D <- venn.diagram(LE.GeneP.D.down, 
                                 "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/venn.LE_vs_gene_Cytosolic_polya_noLFC_downsampled.jpeg", 
                                 main="Differentially expressed Genes\nDown-regulated in Nucleus (FDR<0.05)",
                                 col = "transparent",fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
                                 label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                               "white", "white", "white", "white", "palevioletred4", "white",
                                               "white", "white", "white", "darkblue", "white"),
                                 fontfamily = "Arial",fontface = "bold",cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                                 cat.fontfamily = "Arial", margin=0.2)