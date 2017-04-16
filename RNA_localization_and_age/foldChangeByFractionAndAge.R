library("GenomicRanges")
library("ggplot2")

# Prepare significant genes in a list
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))

SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05),])
elementLengths(SigFracList)
AP = SigFracList[[1]]
FP = SigFracList[[2]]
AR = SigFracList[[3]]
FR = SigFracList[[4]]
sigGenesP = unique(c(as.character(AP[,1]),as.character(FP[,1])))
sigGenesR = unique(c(as.character(AR[,1]),as.character(FR[,1])))
SigP = lapply(FracList[1:2], function(x) x[which(x$X %in% sigGenesP),])
SigR = lapply(FracList[3:4], function(x) x[which(x$X %in% sigGenesR),])

# in all genes
x = lapply(FracList, function(x) x[order(x$X),])
lfc = list()
for (i in 1:length(FracList)){
  tmp = x[[i]]
  lfc[[i]] = data.frame(Gene=tmp$X, baseMean = tmp$baseMean, log2FoldChange = tmp$log2FoldChange, padj = tmp$padj, comparison = names[i])
  tmp2 = lfc[[i]]
  tmp2$Sig= ifelse(tmp2$padj<=0.05,"YES", "NO")
  lfc[[i]] = tmp2
}
AP = lfc[[1]]
FP = lfc[[2]]
AR = lfc[[3]]
FR = lfc[[4]]

polyA = data.frame(Gene = AP$Gene, Adult.LFC = AP$log2FoldChange, Fetal.LFC = FP$log2FoldChange)
ribo = data.frame(Gene = AR$Gene, Adult.LFC = AR$log2FoldChange, Fetal.LFC = FR$log2FoldChange)
polyA$As=ifelse(polyA$Adult.LFC>0, "pos","neg")
polyA$Fs=ifelse(polyA$Fetal.LFC>0, "pos","neg")
polyA$quad=ifelse(polyA$As==polyA$Fs, "same","diff")
ribo$As=ifelse(ribo$Adult.LFC>0, "pos","neg")
ribo$Fs=ifelse(ribo$Fetal.LFC>0, "pos","neg")
ribo$quad=ifelse(ribo$As==ribo$Fs, "same","diff")

ggplot(polyA, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Fetal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Fetal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (RiboZero)") + 
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
AP = lfcP[[1]]
FP = lfcP[[2]]
AR = lfcR[[1]]
FR = lfcR[[2]]

polyA = data.frame(Gene = AP$Gene, Adult.LFC = AP$log2FoldChange, Fetal.LFC = FP$log2FoldChange)
ribo = data.frame(Gene = AR$Gene, Adult.LFC = AR$log2FoldChange, Fetal.LFC = FR$log2FoldChange)
polyA$As=ifelse(polyA$Adult.LFC>0, "pos","neg")
polyA$Fs=ifelse(polyA$Fetal.LFC>0, "pos","neg")
polyA$quad=ifelse(polyA$As==polyA$Fs, "same","diff")
ribo$As=ifelse(ribo$Adult.LFC>0, "pos","neg")
ribo$Fs=ifelse(ribo$Fetal.LFC>0, "pos","neg")
ribo$quad=ifelse(ribo$As==ribo$Fs, "same","diff")

dim(polyA)
#[1] 4474    6
dim(polyA[which(polyA$quad=="same"),])
#[1] 4034    6
dim(ribo)
#[1] 6063    6
dim(ribo[which(ribo$quad=="same"),])
#[1] 5475    6

ggplot(polyA, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Fetal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo, aes(Adult.LFC, Fetal.LFC)) + geom_point(aes(colour = factor(quad))) +
  ylab("Fetal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))


elementLengths(DirList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]]) 
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