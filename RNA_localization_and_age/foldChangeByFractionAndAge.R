library(GenomicRanges)
library(ggplot2)

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

# Prepare significant genes in a list
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres.down), 
                Arres = data.frame(Arres), Frres = data.frame(Frres))
DirList = Map(cbind, FracList, GeneID = lapply(FracList, function(x) rownames(x)),
              Sign = lapply(FracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc")),
              Comparison = list("Adult:PolyA","Prenatal:PolyA", "Adult:RiboZero","Prenatal:RiboZero"))

# in all genes
polya = data.frame(geneID = as.character(DirList$Apres$GeneID[which(DirList$Apres$GeneID %in% DirList$Fpres$GeneID)]))
polya = cbind(polya, data.frame(LFC.A = DirList[["Apres"]][match(polya$geneID, DirList[["Apres"]][,"GeneID"]),"log2FoldChange"],
                                LFC.P = DirList[["Fpres"]][match(polya$geneID, DirList[["Fpres"]][,"GeneID"]),"log2FoldChange"],
                                Sign.A = DirList[["Apres"]][match(polya$geneID, DirList[["Apres"]][,"GeneID"]),"Sign"],
                                Sign.P = DirList[["Fpres"]][match(polya$geneID, DirList[["Fpres"]][,"GeneID"]),"Sign"],
                                Sig.A = DirList[["Apres"]][match(polya$geneID, DirList[["Apres"]][,"GeneID"]),"padj"],
                                Sig.P = DirList[["Fpres"]][match(polya$geneID, DirList[["Fpres"]][,"GeneID"]),"padj"]))
ribo = data.frame(geneID = DirList$Arres$GeneID[which(DirList$Arres$GeneID %in% DirList$Frres$GeneID)])
ribo = data.frame(geneID = DirList$Arres$GeneID[which(DirList$Arres$GeneID %in% DirList$Frres$GeneID)],
                   LFC.A = DirList[["Arres"]][match(ribo$geneID, DirList[["Arres"]][,"GeneID"]),"log2FoldChange"],
                   LFC.P = DirList[["Frres"]][match(ribo$geneID, DirList[["Frres"]][,"GeneID"]),"log2FoldChange"],
                   Sign.A = DirList[["Arres"]][match(ribo$geneID, DirList[["Arres"]][,"GeneID"]),"Sign"],
                   Sign.P = DirList[["Frres"]][match(ribo$geneID, DirList[["Frres"]][,"GeneID"]),"Sign"],
                   Sig.A = DirList[["Arres"]][match(ribo$geneID, DirList[["Arres"]][,"GeneID"]),"padj"],
                   Sig.P = DirList[["Frres"]][match(ribo$geneID, DirList[["Frres"]][,"GeneID"]),"padj"])

polya$quad = ifelse(polya$Sign.A==polya$Sign.P, "same","diff")
ribo$quad = ifelse(ribo$Sign.A==ribo$Sign.P, "same","diff")
polya = polya[which(polya$quad!="NA"),]
ribo = ribo[which(ribo$quad!="NA"),]

head(polya)
polya$sig = "NO"
polya[which(as.character(polya$geneID) %in% c(as.character(polya[which(polya$Sig.A<=0.05),"geneID"]), 
                                              as.character(polya[which(polya$Sig.P<=0.05),"geneID"]))),"sig"] = "YES"
ribo$sig = "NO"
ribo[which(as.character(ribo$geneID) %in% c(as.character(ribo[which(ribo$Sig.A<=0.05),"geneID"]), 
                                              as.character(ribo[which(ribo$Sig.P<=0.05),"geneID"]))),"sig"] = "YES"
ribo$sig = ifelse((ribo$Sig.A<=0.05 | ribo$Sig.P<=0.05), "YES", "NO")
polya$sig = factor(polya$sig, levels = c("YES","NO"))
ribo$sig = factor(ribo$sig, levels = c("YES","NO"))


pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/sigLFC_byFractions_allGenes.pdf", width=3.5, height=4.25)
ggplot(polya, aes(LFC.A, LFC.P)) + geom_point(aes(colour = sig), alpha = 1/2) +
  scale_colour_manual(values=c("red3","gray32")) +
  ylab("In Prenatal") + xlab("In Adult") + ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBy Fraction (PolyA)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
ggplot(ribo[which(ribo$quad!="NA"),], aes(LFC.A, LFC.P)) + geom_point(aes(colour = factor(sig)), alpha = 1/2) +
  ylab("In Prenatal") + xlab("In Adult") + ylim(-3,3) + xlim(-3,3) +
  scale_colour_manual(values=c("red3","gray32")) +
  ggtitle("Log2 Fold Change\nBy Fraction (RiboZero)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
dev.off()


# in significant genes
pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/sigLFC_byFractions_sigGenes.pdf", width=7, height=5)
ggplot(polya[which((polya$Sig.A<=0.05 | polya$Sig.P<=0.05) & polya$quad!="NA"),], 
       aes(LFC.A, LFC.P)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBy Fractions (PolyA)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
ggplot(ribo[which((ribo$Sig.A<=0.05 | ribo$Sig.P<=0.05) & ribo$quad!="NA"),],
       aes(LFC.A, LFC.P)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Prenatal") + xlab("Adult") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-3,3) + xlim(-3,3) +
  ggtitle("Log2 Fold Change\nBy Fractions (RiboZero)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
dev.off()


# Calculate the correlation
cor.test(x = polya$LFC.A, y = polya$LFC.P)
#Pearson's product-moment correlation
#data:  polya.down$LFC.A and polya.down$LFC.P
#t = 125.91, df = 28961, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#0.5872882 0.6021739
#sample estimates:
#cor 
#0.594782 


### By Age
# Prepare significant genes in a list
AgeList = list(Cpres = data.frame(Cpres.down), Npres = data.frame(Npres), 
                Crres = data.frame(Crres), Nrres = data.frame(Nrres))
DirList = Map(cbind, AgeList, GeneID = lapply(AgeList, rownames),
              Sign = lapply(AgeList, function(x) ifelse(x$log2FoldChange > 0,"UpPrenatal", "DownPrenatal")),
              Comparison = list("Cytosol:PolyA","Nucleus:PolyA","Cytosol:RiboZero","Nucleus:RiboZero"))
elementNROWS(DirList)

# in all genes
polya = data.frame(geneID = as.character(DirList$Cpres$GeneID[which(DirList$Cpres$GeneID %in% DirList$Npres$GeneID)]))
polya = cbind(polya, data.frame(LFC.C = DirList[["Cpres"]][match(polya$geneID, DirList[["Cpres"]][,"GeneID"]),"log2FoldChange"],
                                LFC.N = DirList[["Npres"]][match(polya$geneID, DirList[["Npres"]][,"GeneID"]),"log2FoldChange"],
                                Sign.C = DirList[["Cpres"]][match(polya$geneID, DirList[["Cpres"]][,"GeneID"]),"Sign"],
                                Sign.N = DirList[["Npres"]][match(polya$geneID, DirList[["Npres"]][,"GeneID"]),"Sign"],
                                Sig.C = DirList[["Cpres"]][match(polya$geneID, DirList[["Cpres"]][,"GeneID"]),"padj"],
                                Sig.N = DirList[["Npres"]][match(polya$geneID, DirList[["Npres"]][,"GeneID"]),"padj"]))
ribo = data.frame(geneID = DirList$Crres$GeneID[which(DirList$Crres$GeneID %in% DirList$Nrres$GeneID)])
ribo = data.frame(geneID = DirList$Crres$GeneID[which(DirList$Crres$GeneID %in% DirList$Nrres$GeneID)],
                  LFC.C = DirList[["Crres"]][match(ribo$geneID, DirList[["Crres"]][,"GeneID"]),"log2FoldChange"],
                  LFC.N = DirList[["Nrres"]][match(ribo$geneID, DirList[["Nrres"]][,"GeneID"]),"log2FoldChange"],
                  Sign.C = DirList[["Crres"]][match(ribo$geneID, DirList[["Crres"]][,"GeneID"]),"Sign"],
                  Sign.N = DirList[["Nrres"]][match(ribo$geneID, DirList[["Nrres"]][,"GeneID"]),"Sign"],
                  Sig.C = DirList[["Crres"]][match(ribo$geneID, DirList[["Crres"]][,"GeneID"]),"padj"],
                  Sig.N = DirList[["Nrres"]][match(ribo$geneID, DirList[["Nrres"]][,"GeneID"]),"padj"])

polya$quad = ifelse(polya$Sign.C==polya$Sign.N, "same","diff")
ribo$quad = ifelse(ribo$Sign.C==ribo$Sign.N, "same","diff")
polya = polya[which(polya$quad!="NA"),]
ribo = ribo[which(ribo$quad!="NA"),]

head(polya)
polya$sig = "NO"
polya[which(as.character(polya$geneID) %in% c(as.character(polya[which(polya$Sig.C<=0.05),"geneID"]), 
                                              as.character(polya[which(polya$Sig.N<=0.05),"geneID"]))),"sig"] = "YES"
ribo$sig = "NO"
ribo[which(as.character(ribo$geneID) %in% c(as.character(ribo[which(ribo$Sig.C<=0.05),"geneID"]), 
                                            as.character(ribo[which(ribo$Sig.N<=0.05),"geneID"]))),"sig"] = "YES"
ribo$sig = ifelse((ribo$Sig.C<=0.05 | ribo$Sig.N<=0.05), "YES", "NO")
polya$sig = factor(polya$sig, levels = c("YES","NO"))
ribo$sig = factor(ribo$sig, levels = c("YES","NO"))


pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/sigLFC_byAge_allGenes.pdf", width=3.5, height=4.25)
ggplot(polya, aes(x=LFC.C, y=LFC.N)) + geom_point(aes(colour = sig), alpha = 1/2) +
  ylab("Nucleus") + xlab("Cytoplasm") +
  scale_colour_manual(values=c("red3","gray32")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBy Ages (PolyA)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
ggplot(ribo,aes(x=LFC.C, y=LFC.N)) + geom_point(aes(colour = sig), alpha = 1/2) +
  ylab("Nucleus") + 
  xlab("Cytoplasm") +
  scale_colour_manual(values=c("red3","gray32")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBy Ages (RiboZero)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
dev.off()

# in significant genes
pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/sigLFC_byAge_sigGenes.pdf", width=4, height=5)
ggplot(polya[which((polya$Sig.C<=0.05 | polya$Sig.N<=0.05) & polya$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBy Ages (PolyA)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
ggplot(ribo[which((ribo$Sig.C<=0.05 | ribo$Sig.N<=0.05) & ribo$quad!="NA"),],
       aes(LFC.C, LFC.N)) + geom_point(aes(colour = factor(quad)), alpha = 1/2) +
  ylab("Nucleus") + xlab("Cytoplasm") +
  scale_colour_manual(values=c("dimgray","black")) +
  ylim(-8,8) + xlim(-8,8) +
  ggtitle("Log2 Fold Change\nBy Ages (RiboZero)") + 
  theme(title = element_text(size = 20), legend.position="none", text = element_text(size = 20))
dev.off()


# calculate correlation

cor.test(x = polya$LFC.C, y = polya$LFC.N)
#	Pearson's product-moment correlation
#data:  polya.down$LFC.C and polya.down$LFC.N
#t = 335.78, df = 30560, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.8845739 0.8893560
#sample estimates:
#  cor 
#0.8869887 
