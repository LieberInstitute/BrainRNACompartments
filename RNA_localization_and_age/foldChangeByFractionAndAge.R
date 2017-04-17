library(GenomicRanges)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

# Prepare significant genes in a list
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres), 
                Arres = data.frame(Arres), Frres = data.frame(Frres),
                Fpres_down = data.frame(Fpres.down))
DirList = Map(cbind, FracList, 
              GeneID = lapply(FracList, function(x) rownames(x)),
              Sign = lapply(FracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc")),
              Comparison = list("Adult:PolyA","Prenatal:PolyA",
                                "Adult:RiboZero","Prenatal:RiboZero","Prenatal:PolyA"),
              Sig=lapply(FracList, function(x) ifelse(x$padj<=0.05,"YES", "NO")))

# in all genes
DirList = lapply(DirList, function(x) x[order(x$GeneID),])
polya = data.frame(gene.A = DirList[["Apres"]][,7], gene.F = DirList[["Fpres"]][,7],
                   LFC.A = DirList[["Apres"]][,2], LFC.F = DirList[["Fpres"]][,2],
                   Sign.A = DirList[["Apres"]][,8], Sign.F = DirList[["Fpres"]][,8],
                   Sig.A = DirList[["Apres"]][,10], Sig.F = DirList[["Fpres"]][,10])
ribo = data.frame(gene.A = DirList[["Arres"]][,7], gene.F = DirList[["Frres"]][,7],
                   LFC.A = DirList[["Arres"]][,2], LFC.F = DirList[["Frres"]][,2],
                   Sign.A = DirList[["Arres"]][,8], Sign.F = DirList[["Frres"]][,8],
                   Sig.A = DirList[["Arres"]][,10], Sig.F = DirList[["Frres"]][,10])
x = DirList[["Apres"]]
x = x[match(as.character(DirList[["Fpres_down"]][,7]),as.character(x$GeneID)),]
polya.down = data.frame(gene.A = x[,7], gene.F = DirList[["Fpres_down"]][,7],
                   LFC.A = x[,2], LFC.F = DirList[["Fpres_down"]][,2],
                   Sign.A = x[,8], Sign.F = DirList[["Fpres_down"]][,8],
                   Sig.A = x[,10], Sig.F = DirList[["Fpres_down"]][,10])
identical(as.character(polya.down$gene.A),as.character(polya.down$gene.F))
x = ifelse(as.character(polya.down$gene.A)==as.character(polya.down$gene.F), "YES","NO")
x[which(x=="NO")]
identical(as.character(polya$gene.A),as.character(polya$gene.F))
identical(as.character(ribo$gene.A),as.character(ribo$gene.F))
polya$quad = ifelse(polya$Sign.A==polya$Sign.F, "same","diff")
polya.down$quad = ifelse(polya.down$Sign.A==polya.down$Sign.F, "same","diff")
ribo$quad = ifelse(ribo$Sign.A==ribo$Sign.F, "same","diff")

ggplot(polya[which(polya$quad!="NA"),], aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad))) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(polya.down[which(polya.down$quad!="NA"),], 
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad))) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo[which(ribo$quad!="NA"),], aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad))) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# in significant genes
ggplot(polya[which((polya$Sig.A=="YES" | polya$Sig.F=="YES") & polya$quad!="NA"),], 
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad))) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(polya.down[which((polya.down$Sig.A=="YES" | polya.down$Sig.F=="YES") & polya.down$quad!="NA"),], 
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad))) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo[which((ribo$Sig.A=="YES" | ribo$Sig.F=="YES") & ribo$quad!="NA"),],
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad))) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# 


names(DirList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", 
                   "Fetal\nRibozero\nDown")
APu = DirList[["Adult\nPolyA\nUp"]]
FPu = DirList[["Fetal\nPolyA\nUp"]]
ARu = DirList[["Adult\nRibozero\nUp"]]
FRu = DirList[["Fetal\nRibozero\nUp"]]
APd = DirList[["Adult\nPolyA\nDown"]]
FPd = DirList[["Fetal\nPolyA\nDown"]]
ARd = DirList[["Adult\nRibozero\nDown"]]
FRd = DirList[["Fetal\nRibozero\nDown"]]
upP = unique(c(as.character(APu$X), as.character(FPu$X)))
upR = unique(c(as.character(ARu$X), as.character(FRu$X)))
dP = unique(c(as.character(APd$X), as.character(FPd$X)))
dR = unique(c(as.character(APd$X), as.character(FPd$X)))
polyA=polyA[which(polyA$quad=="same" | polyA$quad=="diff"),]
ribo=ribo[which(ribo$quad=="same" | ribo$quad=="diff"),]
venn.polyA.UP = list("Different Direction" = polyA[which(polyA$quad=="diff" & polyA$Gene %in% upP),],
                     "Same Direction" = polyA[which(polyA$quad=="same" & polyA$Gene %in% upP),],
                  "Adult:Nucleus"=DirList[["Adult\nPolyA\nUp"]], "Fetal:Nucleus"=DirList[["Fetal\nPolyA\nUp"]])
venn.polyA.D = list("Different Direction"=polyA[which(polyA$quad=="diff" & polyA$Gene %in% dP),], 
                    "Same Direction"=polyA[which(polyA$quad=="same" & polyA$Gene %in% dP),],
                    "Adult:Cytosol"=DirList[["Adult\nPolyA\nDown"]], "Fetal:Cytosol"=DirList[["Fetal\nPolyA\nDown"]])
venn.polyA.UP = lapply(venn.polyA.UP, function(x) as.character(x[,1]))
venn.polyA.D = lapply(venn.polyA.D, function(x) as.character(x[,1]))
venn.ribo.UP = list("Different Direction" = ribo[which(ribo$quad=="diff" & ribo$Gene %in% upR),],
                     "Same Direction" = ribo[which(ribo$quad=="same" & ribo$Gene %in% upR),],
                     "Adult:Nucleus"=DirList[["Adult\nRibozero\nUp"]], "Fetal:Nucleus"=DirList[["Fetal\nRibozero\nUp"]])
venn.ribo.D = list("Different Direction"=ribo[which(ribo$quad=="diff" & ribo$Gene %in% dR),], 
                    "Same Direction"=ribo[which(ribo$quad=="same" & ribo$Gene %in% dR),],
                    "Adult:Cytosol"=DirList[["Adult\nRibozero\nDown"]], "Fetal:Cytosol"=DirList[["Fetal\nRibozero\nDown"]])
venn.ribo.UP = lapply(venn.ribo.UP, function(x) as.character(x[,1]))
venn.ribo.D = lapply(venn.ribo.D, function(x) as.character(x[,1]))

venn.p.UP <- venn.diagram(venn.polyA.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.polyA_UP_LFC_and_DEGs.jpeg", 
                          main="Overlap of Nuclear DEGs and LCF Direction (polyA)",
                          col = "transparent",
                          fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                          alpha = 0.50,
                          label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                        "white", "white", "white", "white", "palevioletred4", "white",
                                        "white", "white", "white", "darkblue", "white"),
                          fontfamily = "Arial",
                          fontface = "bold",
                          cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                          cat.fontfamily = "Arial", margin=0.2, cat.dist=0.3)
venn.p.D <- venn.diagram(venn.polyA.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.polyA_D_LFC_and_DEGs.jpeg", 
                         main="Overlap of Cytosolic DEGs and LCF Direction (polyA)",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2, cat.dist=0.3)


venn.r.UP <- venn.diagram(venn.ribo.UP, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.ribo_UP_LFC_and_DEGs.jpeg", 
                           main="Overlap of Nuclear DEGs and LCF Direction (Ribozero)",
                           col = "transparent",
                           fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                           alpha = 0.50,
                           label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                         "white", "white", "white", "white", "palevioletred4", "white",
                                         "white", "white", "white", "darkblue", "white"),
                           fontfamily = "Arial",
                           fontface = "bold",
                           cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                           cat.fontfamily = "Arial", margin=0.2, cat.dist=0.3)
venn.r.D <- venn.diagram(venn.ribo.D, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/venn.ribo_D_LFC_and_DEGs.jpeg", 
                          main="Overlap of Cytosolic DEGs and LCF Direction (Ribozero)",
                          col = "transparent",
                          fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                          alpha = 0.50,
                          label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                        "white", "white", "white", "white", "palevioletred4", "white",
                                        "white", "white", "white", "darkblue", "white"),
                          fontfamily = "Arial",
                          fontface = "bold",
                          cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                          cat.fontfamily = "Arial", margin=0.2, cat.dist=0.3)


### By Age
# Prepare significant genes in a list
AgeList = list(Cpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Cpres.csv"),
                Npres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Npres.csv"),
                Crres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Crres.csv"),
                Nrres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nrres.csv"))

SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05),])
elementLengths(SigAgeList)
CP = SigAgeList[[1]]
NP = SigAgeList[[2]]
CR = SigAgeList[[3]]
NR = SigAgeList[[4]]
sigGenesP = unique(c(as.character(CP[,1]),as.character(NP[,1])))
sigGenesR = unique(c(as.character(CR[,1]),as.character(NR[,1])))
SigP = lapply(AgeList[1:2], function(x) x[which(x$X %in% sigGenesP),])
SigR = lapply(AgeList[3:4], function(x) x[which(x$X %in% sigGenesR),])

# in all genes
x = lapply(AgeList, function(x) x[order(x$X),])
lfc = list()
for (i in 1:length(AgeList)){
  tmp = x[[i]]
  lfc[[i]] = data.frame(Gene=tmp$X, baseMean = tmp$baseMean, log2FoldChange = tmp$log2FoldChange, padj = tmp$padj, comparison = names[i])
  tmp2 = lfc[[i]]
  tmp2$Sig= ifelse(tmp2$padj<=0.05,"YES", "NO")
  lfc[[i]] = tmp2
}
CP = lfc[[1]]
NP = lfc[[2]]
CR = lfc[[3]]
NR = lfc[[4]]

polyA = data.frame(Gene = CP$Gene, Adult.LFC = CP$log2FoldChange, Fetal.LFC = NP$log2FoldChange)
ribo = data.frame(Gene = CR$Gene, Adult.LFC = CR$log2FoldChange, Fetal.LFC = NR$log2FoldChange)
polyA$Cs=ifelse(polyA$Adult.LFC>0, "pos","neg")
polyA$Ns=ifelse(polyA$Fetal.LFC>0, "pos","neg")
polyA$quad=ifelse(polyA$Cs==polyA$Ns, "same","diff")
ribo$Cs=ifelse(ribo$Adult.LFC>0, "pos","neg")
ribo$Ns=ifelse(ribo$Fetal.LFC>0, "pos","neg")
ribo$quad=ifelse(ribo$Cs==ribo$Ns, "same","diff")

ggplot(polyA, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Nuclear") + 
  xlab("Cytosolic") +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Nuclear") + 
  xlab("Cytosolic") +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# in significant genes
SigP = lapply(SigP, function(x) x[order(x$X),])
SigR = lapply(SigR, function(x) x[order(x$X),])
lfcP = lfcR = list()
for (i in 1:length(SigP)){
  tmpP = SigP[[i]]
  tmpR = SigR[[i]]
  lfcP[[i]] = data.frame(Gene=tmpP$X, baseMean = tmpP$baseMean, log2FoldChange = tmpP$log2FoldChange, padj = tmpP$padj)
  lfcR[[i]] = data.frame(Gene=tmpR$X, baseMean = tmpR$baseMean, log2FoldChange = tmpR$log2FoldChange, padj = tmpR$padj)
}
CP = lfcP[[1]]
NP = lfcP[[2]]
CR = lfcR[[1]]
NR = lfcR[[2]]

polyA = data.frame(Gene = CP$Gene, Adult.LFC = CP$log2FoldChange, Fetal.LFC = NP$log2FoldChange)
ribo = data.frame(Gene = CR$Gene, Adult.LFC = CR$log2FoldChange, Fetal.LFC = NR$log2FoldChange)
polyA$Cs=ifelse(polyA$Adult.LFC>0, "pos","neg")
polyA$Ns=ifelse(polyA$Fetal.LFC>0, "pos","neg")
polyA$quad=ifelse(polyA$Cs==polyA$Ns, "same","diff")
ribo$Cs=ifelse(ribo$Adult.LFC>0, "pos","neg")
ribo$Ns=ifelse(ribo$Fetal.LFC>0, "pos","neg")
ribo$quad=ifelse(ribo$Cs==ribo$Ns, "same","diff")

ggplot(polyA, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Nuclear") + 
  xlab("Cytosolic") +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Nuclear") + 
  xlab("Cytosolic") +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
