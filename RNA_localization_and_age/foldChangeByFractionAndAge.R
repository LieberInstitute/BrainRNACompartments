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

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/sigLFC_byFractions_allGenes.pdf", width=7, height=5)
ggplot(polya[which(polya$quad!="NA"),], aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  scale_colour_manual(values=c("dimgray","black")) +
  ylab("Prenatal") + 
  xlab("Adult") +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(polya.down[which(polya.down$quad!="NA"),], 
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + 
  xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo[which(ribo$quad!="NA"),], aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + 
  xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# in significant genes
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/sigLFC_byFractions_sigGenes.pdf", width=7, height=5)
ggplot(polya[which((polya$Sig.A=="YES" | polya$Sig.F=="YES") & polya$quad!="NA"),], 
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + 
  xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(polya.down[which((polya.down$Sig.A=="YES" | polya.down$Sig.F=="YES") & polya.down$quad!="NA"),], 
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + 
  xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo[which((ribo$Sig.A=="YES" | ribo$Sig.F=="YES") & ribo$quad!="NA"),],
       aes(LFC.A, LFC.F)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + 
  xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBetween Fractions (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# Calculate the correlation
cor.test(x = polya.down$LFC.A, y = polya.down$LFC.F)
#Pearson's product-moment correlation
#data:  polya.down$LFC.A and polya.down$LFC.F
#t = 125.91, df = 28961, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#0.5872882 0.6021739
#sample estimates:
#cor 
#0.594782 


### By Age
# Prepare significant genes in a list
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres), 
                Crres = data.frame(Crres), Nrres = data.frame(Nrres),
                Cpres_down = data.frame(Cpres.down))
DirList = Map(cbind, AgeList, 
              GeneID = lapply(AgeList, function(x) rownames(x)),
              Sign = lapply(AgeList, function(x) ifelse(x$log2FoldChange > 0,"UpPrenatal", "DownPrenatal")),
              Comparison = list("Cytosol:PolyA","Nucleus:PolyA",
                                "Cytosol:RiboZero","Nucleus:RiboZero","Nucleus:PolyA"),
              Sig=lapply(AgeList, function(x) ifelse(x$padj<=0.05,"YES", "NO")))
elementNROWS(DirList)

# in all genes
DirList = lapply(DirList, function(x) x[order(x$GeneID),])
polya = data.frame(gene.C = DirList[["Cpres"]][,7], gene.F = DirList[["Npres"]][,7],
                   LFC.C = DirList[["Cpres"]][,2], LFC.N = DirList[["Npres"]][,2],
                   Sign.C = DirList[["Cpres"]][,8], Sign.N = DirList[["Npres"]][,8],
                   Sig.C = DirList[["Cpres"]][,10], Sig.N = DirList[["Npres"]][,10])
ribo = data.frame(gene.C = DirList[["Crres"]][,7], gene.N = DirList[["Nrres"]][,7],
                  LFC.C = DirList[["Crres"]][,2], LFC.N = DirList[["Nrres"]][,2],
                  Sign.C = DirList[["Crres"]][,8], Sign.N = DirList[["Nrres"]][,8],
                  Sig.C = DirList[["Crres"]][,10], Sig.N = DirList[["Nrres"]][,10])
x = DirList[["Npres"]]
x = x[match(as.character(DirList[["Cpres_down"]][,7]),as.character(x$GeneID)),]
polya.down = data.frame(gene.N = x[,7], gene.C = DirList[["Cpres_down"]][,7],
                        LFC.N = x[,2], LFC.C = DirList[["Cpres_down"]][,2],
                        Sign.N = x[,8], Sign.C = DirList[["Cpres_down"]][,8],
                        Sig.N = x[,10], Sig.C = DirList[["Cpres_down"]][,10])
identical(as.character(polya.down$gene.C),as.character(polya.down$gene.N))
x = ifelse(as.character(polya.down$gene.C)==as.character(polya.down$gene.N), "YES","NO")
x[which(x=="NO")]
identical(as.character(polya$gene.C),as.character(polya$gene.N))
identical(as.character(ribo$gene.C),as.character(ribo$gene.N))
polya$quad = ifelse(polya$Sign.C==polya$Sign.N, "same","diff")
polya.down$quad = ifelse(polya.down$Sign.C==polya.down$Sign.N, "same","diff")
ribo$quad = ifelse(ribo$Sign.C==ribo$Sign.N, "same","diff")

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/sigLFC_byAge_allGenes.pdf", width=7, height=5)
ggplot(polya[which(polya$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo[which(ribo$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(polya.down[which(polya.down$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

# in significant genes
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/sigLFC_byAge_sigGenes.pdf", width=7, height=5)
ggplot(polya[which((polya$Sig.C=="YES" | polya$Sig.N=="YES") & polya$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(ribo[which((ribo$Sig.C=="YES" | ribo$Sig.N=="YES") & ribo$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (RiboZero)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(polya.down[which((polya.down$Sig.C=="YES" | polya.down$Sig.N=="YES") & polya.down$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBetween Ages (PolyA)") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# calculate correlation

cor.test(x = polya.down$LFC.C, y = polya.down$LFC.N)
#	Pearson's product-moment correlation
#data:  polya.down$LFC.C and polya.down$LFC.N
#t = 335.78, df = 30560, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.8845739 0.8893560
#sample estimates:
#  cor 
#0.8869887 
