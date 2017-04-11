library("ggplot2")
library(DESeq2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
vsd = varianceStabilizingTransformation(geneCounts)
testgenes = rbind(geneRpkm["ENSG00000075624",], geneRpkm["ENSG00000102081",], geneRpkm["ENSG00000229807",], geneRpkm["ENSG00000251562",])
rownames(testgenes)=c("ACTB", "FMR1", "XIST", "MALAT1")

FracList = list(Apres = read.csv("/Users/amandaprice/Downloads/Apres.csv"),
                Fpres = read.csv("/Users/amandaprice/Downloads/Fpres.csv"),
                Arres = read.csv("/Users/amandaprice/Downloads/Arres.csv"),
                Frres = read.csv("/Users/amandaprice/Downloads/Frres.csv"))
genes = c("ENSG00000075624","ENSG00000251562")
lfc = lapply(FracList, function(x) x[which(x$X %in% genes),])
x = do.call(rbind, lfc)
write.table(x, "./Dropbox/sorted_figures/new/MALAT1.ACTB.txt", quote=F)

# Exclude the unpaired sample
tg = t(testgenes[,which(colnames(testgenes)!="Br1113C1_RiboZero")])
tg = data.frame(tg)
Cytosol <- pd[which(pd$Zone=="Cytosol"),]
Nucleus <- pd[which(pd$Zone=="Nucleus"),]
testnuc = tg[which(rownames(tg)%in%Nucleus$SampleID),]
testcyt = tg[which(rownames(tg)%in%Cytosol$SampleID),]

# one-tailed paired t test: ACTB
t.test(testnuc$ACTB, testcyt$ACTB, paired=TRUE, alternative = "less")

#data (using VSD):  testnuc$ACTB and testcyt$ACTB
#t = -5.7857, df = 10, p-value = 8.818e-05
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.6481196
#sample estimates:
#  mean of the differences 
#-0.9437691

#data (using RPKM):  testnuc$ACTB and testcyt$ACTB
#t = -4.5757, df = 10, p-value = 0.0005086
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -128.8618
#sample estimates:
#  mean of the differences 
#-213.3835

# one-tailed paired t test: FMR1
t.test(testnuc$FMR1, testcyt$FMR1, paired=TRUE, alternative = "less")

#data (using VSD):  testnuc$FMR1 and testcyt$FMR1
#t = -0.9437, df = 10, p-value = 0.1838
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.05699625
#sample estimates:
#  mean of the differences 
#-0.06191235 

#data (using RPKM):  testnuc$FMR1 and testcyt$FMR1
#t = -2.9087, df = 10, p-value = 0.0078
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.7006774
#sample estimates:
#  mean of the differences 
#-1.859164

# one-tailed paired t test in PolyA samples: FMR1
FMR1pN = rbind(testnuc[1,],  testnuc[2,], testnuc[3,], testnuc[4,], testnuc[5,], testnuc[6,])
FMR1pC = rbind(testcyt[1,],  testcyt[2,], testcyt[3,], testcyt[4,], testcyt[5,], testcyt[6,])
t.test(FMR1pN$FMR1, FMR1pC$FMR1, paired=TRUE, alternative = "less")

#data (usng VSD):  FMR1pN$FMR1 and FMR1pC$FMR1
#t = -3.1722, df = 5, p-value = 0.01238
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.05447295
#sample estimates:
#  mean of the differences 
#-0.1493327

#data (using RPKM):  FMR1pN$FMR1 and FMR1pC$FMR1
#t = -2.0609, df = 5, p-value = 0.04717
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.02070892
#sample estimates:
#  mean of the differences 
#-0.9316126 

# one-tailed paired t test in females: XIST
females <- pd[which(pd$Sex=="F"),]
XISTnuc = rbind(testnuc[1,],  testnuc[4,], testnuc[9,])
XISTcyt = rbind(testcyt[1,],  testcyt[4,], testcyt[9,])
t.test(XISTnuc$XIST, XISTcyt$XIST, paired=TRUE, alternative = "greater")

#data (using VSD):  XISTnuc$XIST and XISTcyt$XIST
#t = 2.2905, df = 2, p-value = 0.07456
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.6471344        Inf
#sample estimates:
#  mean of the differences 
#2.354768 

#data (using RPKM):  XISTnuc$XIST and XISTcyt$XIST
#t = 1.9355, df = 2, p-value = 0.09628
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -2.763293       Inf
#sample estimates:
#  mean of the differences 
#5.432986

# one-tailed paired t test: MALAT1
t.test(testnuc$MALAT1, testcyt$MALAT1, paired=TRUE, alternative = "greater")

#data (using VSD):  testnuc$MALAT1 and testcyt$MALAT1
#t = 6.1758, df = 10, p-value = 5.236e-05
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.9017213       Inf
#sample estimates:
#  mean of the differences 
#1.276278

#data (using RPKM):  testnuc$MALAT1 and testcyt$MALAT1
#t = 2.8484, df = 10, p-value = 0.008649
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  36.94485      Inf
#sample estimates:
#  mean of the differences 
#101.5832 


# Plot the Log2 Fold Change and SE by library for these genes
Zpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Zpres.csv")
Zrres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Zrres.csv")
genes = rbind(Zpres[which(Zpres$X=="ENSG00000075624"),], Zrres[which(Zrres$X=="ENSG00000075624"),], 
              Zpres[which(Zpres$X=="ENSG00000102081"),], Zrres[which(Zrres$X=="ENSG00000102081"),],
             Zpres[which(Zpres$X=="ENSG00000229807"),], Zrres[which(Zrres$X=="ENSG00000229807"),], 
             Zpres[which(Zpres$X=="ENSG00000251562"),], Zrres[which(Zrres$X=="ENSG00000251562"),])
genes = data.frame(genes, "gene"=c("ACTB", "ACTB", "FMR1", "FMR1", "XIST", "XIST", "MALAT1", "MALAT1"), 
                   Library = c("PolyA", "RiboZero", "PolyA", "RiboZero", "PolyA", "RiboZero", "PolyA", "RiboZero"))

dodge <- position_dodge(width=0.9)
limits <- aes(ymax = log2FoldChange + lfcSE, ymin=log2FoldChange - lfcSE)
ggplot(genes[which(genes$gene!="XIST" & genes$gene!="FMR1"),], aes(x=gene, y=log2FoldChange, fill=Library), color=Library) + 
  stat_summary(position=position_dodge(),geom="bar") +
  geom_errorbar(mapping = limits, position = dodge, width=0.25) +
  ylim(-2,2) +
  ylab("Log2 Fold Change (±SE)") + 
  xlab("") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(.83, 0.3)) + 
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))