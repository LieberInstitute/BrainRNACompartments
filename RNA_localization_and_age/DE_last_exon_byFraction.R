library('DESeq2')
library(VennDiagram)

## Does the pattern of fewer DE genes by fraction in fetal hold when measuring the last exon only?
# extract the last exon from each gene

head(exonMap)
exonPlus = exonMap[which(exonMap$Strand=="+"),]
exonMinus = exonMap[which(exonMap$Strand=="-"),]
exonPlus$Num = seq_len(nrow(exonPlus))
exonPlus = exonPlus[order(exonPlus$Num, decreasing = T),]
firstPlus = exonPlus[match(unique(exonPlus$Geneid), exonPlus$Geneid),]
exonMinus$Num = seq_len(nrow(exonMinus))
exonMinus = exonMinus[order(exonMinus$Num, decreasing = T),]
firstMinus = exonMinus[match(unique(exonMinus$Geneid), exonMinus$Geneid),]
LastExonMap = rbind(firstPlus, firstMinus)
LastExonMap = LastExonMap[which(LastExonMap$Geneid %in% rownames(geneCounts)),]
lastexonCounts = exonCounts[which(rownames(exonCounts)%in%rownames(LastExonMap)),]

# DESeq2 by fraction in four groups

Ap.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts)%in%Adult.polyA$SampleID)],
                                   colData = Adult.polyA, design = ~ Zone)
Ap.ddsLE <- DESeq(Ap.ddsLE)
ApresLE <- results(Ap.ddsLE)
Fp.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts)%in%Fetal.polyA$SampleID)],
                                   colData = Fetal.polyA, design = ~ Zone)
Fp.ddsLE <- DESeq(Fp.ddsLE)
FpresLE <- results(Fp.ddsLE)
Ar.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts)%in%Adult.Ribo$SampleID)],
                                   colData = Adult.Ribo, design = ~ Zone)
Ar.ddsLE <- DESeq(Ar.ddsLE)
ArresLE <- results(Ar.ddsLE)
Fr.ddsLE <- DESeqDataSetFromMatrix(countData = lastexonCounts[,which(colnames(lastexonCounts)%in%Fetal.Ribo$SampleID)], 
                                   colData = Fetal.Ribo, design = ~ Zone)
Fr.ddsLE <- DESeq(Fr.ddsLE)
FrresLE <- results(Fr.ddsLE)
plotMA(ApresLE, alpha=0.05, main="Adult PolyA Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(FpresLE, alpha=0.05, main="Fetal PolyA Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(ArresLE, alpha=0.05, main="Adult Ribozero Samples (Last Exon Only)", ylim=c(-3,3))
plotMA(FrresLE, alpha=0.05, main="Fetal Ribozero Samples (Last Exon Only)", ylim=c(-3,3))

## Are similar numbers of genes expressed in the four groups? 
  # 44441 tested genes

dim(Adult.polyA.counts[which(rowMeans(Adult.polyA.counts) > 0),]) # 31389
dim(Fetal.polyA.counts[which(rowMeans(Fetal.polyA.counts) > 0),]) # 31627
dim(Adult.Ribo.counts[which(rowMeans(Adult.Ribo.counts) > 0),]) # 38401 
dim(Fetal.Ribo.counts[which(rowMeans(Fetal.Ribo.counts) > 0),]) # 38119

## What are the DE genes?

# Make list of significant genes
LEResults = list(AP = ApresLE, FP = FpresLE, AR = ArresLE, FR = FrresLE)
sigLE = lapply(LEResults, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
    elementLengths(sigLE)
    #AP   FP   AR   FR 
    #841   58 1709  616
sigLE = lapply(sigLE, function(x) as.data.frame(x))
Sign = lapply(sigLE, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigLE = Map(cbind, sigLE, Sign = Sign)
sigLE = lapply(sigLE, function(x) split(x, x$Sign))
sigLE = list("Adult\nPolyA\nUp" = sigLE[["AP"]][["UpNuc"]], "Adult\nPolyA\nDown" = sigLE[["AP"]][["DownNuc"]],
             "Fetal\nPolyA\nUp" = sigLE[["FP"]][["UpNuc"]], "Fetal\nPolyA\nDown" = sigLE[["FP"]][["DownNuc"]],
             "Adult\nRibozero\nUp" = sigLE[["AR"]][["UpNuc"]], "Adult\nRibozero\nDown" = sigLE[["AR"]][["DownNuc"]],
             "Fetal\nRibozero\nUp" = sigLE[["FR"]][["UpNuc"]], "Fetal\nRibozero\nDown" = sigLE[["FR"]][["DownNuc"]])
sigLE = lapply(sigLE, function(x) rownames(x))
sigLE = lapply(sigLE, function(x) exonMap$Geneid[which(rownames(exonMap)%in%x)])
  elementLengths(sigLE)
    #Adult\nPolyA\nUp    Adult\nPolyA\nDown      Fetal\nPolyA\nUp    Fetal\nPolyA\nDown   Adult\nRibozero\nUp Adult\nRibozero\nDown 
    #273                   568                    57                     1                   901                   808 
    #Fetal\nRibozero\nUp Fetal\nRibozero\nDown 
    #598                    18 
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
  elementLengths(SigFracList)
    #Apres Fpres Arres Frres 
    #1648    75  2869   349 
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracList = Map(cbind, SigFracList, Sign = Sign)
sigFracList = lapply(sigFracList, function(x) split(x, x$Sign))
SigList = list("Adult\nPolyA\nUp" = sigFracList[["Apres"]][["UpNuc"]], "Adult\nPolyA\nDown" = sigFracList[["Apres"]][["DownNuc"]],
               "Fetal\nPolyA\nUp" = sigFracList[["Fpres"]][["UpNuc"]], "Fetal\nPolyA\nDown" = sigFracList[["Fpres"]][["DownNuc"]],
               "Adult\nRibozero\nUp" = sigFracList[["Arres"]][["UpNuc"]], "Adult\nRibozero\nDown" = sigFracList[["Arres"]][["DownNuc"]],
               "Fetal\nRibozero\nUp" = sigFracList[["Frres"]][["UpNuc"]], "Fetal\nRibozero\nDown" = sigFracList[["Frres"]][["DownNuc"]]) 
SigList = lapply(SigList, function(x) as.character(x$X))
  elementLengths(SigList)
    #Adult\nPolyA\nUp    Adult\nPolyA\nDown      Fetal\nPolyA\nUp    Fetal\nPolyA\nDown   Adult\nRibozero\nUp Adult\nRibozero\nDown 
    #954                   694                    74                     1                  1507                  1362 
    #Fetal\nRibozero\nUp Fetal\nRibozero\nDown 
    #339                    10 
save(LEResults, LastExonMap, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Last.Exon>DEResults.rda")

# Create Venn Diagrams of the overlap of the genes DE in the whole gene and by last exon (LFC≥1, FDR<0.05)

LE.GeneP.UP <- list("Adult:PolyA\nLast Exon"=sigLE[["Adult\nPolyA\nUp"]], "Fetal:PolyA\nLast Exon"=sigLE[["Fetal\nPolyA\nUp"]], 
                     "Adult:PolyA\nGene"=SigList[["Adult\nRibozero\nUp"]], "Fetal:PolyA\nGene"=SigList[["Fetal\nRibozero\nUp"]])
LE.GeneR.UP <- list("Adult:Ribozero\nLast Exon"=sigLE[["Adult\nRibozero\nUp"]], "Fetal:Ribozero\nLast Exon"=sigLE[["Fetal\nRibozero\nUp"]], 
                    "Adult:Ribozero\nGene"=SigList[["Adult\nRibozero\nUp"]], "Fetal:Ribozero\nGene"=SigList[["Fetal\nRibozero\nUp"]])

LE.GeneP.D <- list("Adult:PolyA\nLast Exon"=sigLE[["Adult\nPolyA\nDown"]], "Fetal:PolyA\nLast Exon"=sigLE[["Fetal\nPolyA\nDown"]], 
                    "Adult:PolyA\nGene"=SigList[["Adult\nRibozero\nDown"]], "Fetal:PolyA\nGene"=SigList[["Fetal\nRibozero\nDown"]])
LE.GeneR.D <- list("Adult:Ribozero\nLast Exon"=sigLE[["Adult\nRibozero\nDown"]], "Fetal:Ribozero\nLast Exon"=sigLE[["Fetal\nRibozero\nDown"]], 
                    "Adult:Ribozero\nGene"=SigList[["Adult\nRibozero\nDown"]], "Fetal:Ribozero\nGene"=SigList[["Fetal\nRibozero\nDown"]])

venn.LEvsGeneP.UP <- venn.diagram(LE.GeneP.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GeneP.UP.jpeg", 
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


venn.LEvsGeneP.D <- venn.diagram(LE.GeneP.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GeneP.D.jpeg", 
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

venn.LEvsGeneR.UP <- venn.diagram(LE.GeneR.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GeneR.UP.jpeg", 
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


venn.LEvsGeneR.D <- venn.diagram(LE.GeneR.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GeneR.D.jpeg", 
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

# Find overlaps just for P<0.05, no LFC criteria for Gene-level expression (FDR<0.05)

SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05),])
  elementLengths(SigFracList)
    #Apres Fpres Arres Frres 
    #4417   297  5605  1515
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracList = Map(cbind, SigFracList, Sign = Sign)
sigFracList = lapply(sigFracList, function(x) split(x, x$Sign))
SigList = list(Apres.Up = sigFracList[["Apres"]][["UpNuc"]], Apres.Down = sigFracList[["Apres"]][["DownNuc"]],
               Fpres.Up = sigFracList[["Fpres"]][["UpNuc"]], Fpres.Down = sigFracList[["Fpres"]][["DownNuc"]],
               Arres.Up = sigFracList[["Arres"]][["UpNuc"]], Arres.Down = sigFracList[["Arres"]][["DownNuc"]],
               Frres.Up = sigFracList[["Frres"]][["UpNuc"]], Frres.Down = sigFracList[["Frres"]][["DownNuc"]]) 
SigList = lapply(SigList, function(x) as.character(x$X))
names(SigList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", "Fetal\nRibozero\nDown")
  elementLengths(SigList)
    #Adult\nPolyA\nUp    Adult\nPolyA\nDown      Fetal\nPolyA\nUp    Fetal\nPolyA\nDown   Adult\nRibozero\nUp Adult\nRibozero\nDown 
    #2115                  2302                   274                    23                  2261                  3344 
    #Fetal\nRibozero\nUp Fetal\nRibozero\nDown 
    #930                   585 

LE.GeneP.UP <- list("Adult:PolyA\nLast Exon"=sigLE[["Adult\nPolyA\nUp"]], "Fetal:PolyA\nLast Exon"=sigLE[["Fetal\nPolyA\nUp"]], 
                    "Adult:PolyA\nGene"=SigList[["Adult\nRibozero\nUp"]], "Fetal:PolyA\nGene"=SigList[["Fetal\nRibozero\nUp"]])
LE.GeneR.UP <- list("Adult:Ribozero\nLast Exon"=sigLE[["Adult\nRibozero\nUp"]], "Fetal:Ribozero\nLast Exon"=sigLE[["Fetal\nRibozero\nUp"]], 
                    "Adult:Ribozero\nGene"=SigList[["Adult\nRibozero\nUp"]], "Fetal:Ribozero\nGene"=SigList[["Fetal\nRibozero\nUp"]])

LE.GeneP.D <- list("Adult:PolyA\nLast Exon"=sigLE[["Adult\nPolyA\nDown"]], "Fetal:PolyA\nLast Exon"=sigLE[["Fetal\nPolyA\nDown"]], 
                   "Adult:PolyA\nGene"=SigList[["Adult\nRibozero\nDown"]], "Fetal:PolyA\nGene"=SigList[["Fetal\nRibozero\nDown"]])
LE.GeneR.D <- list("Adult:Ribozero\nLast Exon"=sigLE[["Adult\nRibozero\nDown"]], "Fetal:Ribozero\nLast Exon"=sigLE[["Fetal\nRibozero\nDown"]], 
                   "Adult:Ribozero\nGene"=SigList[["Adult\nRibozero\nDown"]], "Fetal:Ribozero\nGene"=SigList[["Fetal\nRibozero\nDown"]])

venn.LEvsGenePnoLFC.UP <- venn.diagram(LE.GeneP.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GenePnoLFC.UP.jpeg", 
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


venn.LEvsGenePnoLFC.D <- venn.diagram(LE.GeneP.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GenePnoLFC.D.jpeg", 
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

venn.LEvsGeneRnoLFC.UP <- venn.diagram(LE.GeneR.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GeneRnoLFC.UP.jpeg", 
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


venn.LEvsGeneRnoLFC.D <- venn.diagram(LE.GeneR.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/LE.GeneRnoLFC.D.jpeg", 
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

