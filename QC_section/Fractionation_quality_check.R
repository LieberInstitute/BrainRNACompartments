library(ggplot2)
library(DESeq2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

ens = c("ENSG00000075624","ENSG00000102081","ENSG00000229807","ENSG00000251562")
map = data.frame(ids = rownames(geneMap[which(geneMap$ensemblID %in% ens),]), 
                 ens = geneMap[which(geneMap$ensemblID %in% ens),"ensemblID"],
                 sym = geneMap[which(geneMap$ensemblID %in% ens),"Symbol"])
testgenes = geneRpkm[which(rownames(geneRpkm) %in% map$ids),]
rownames(testgenes)=geneMap[match(map$ids,rownames(geneMap)), "Symbol"]
FracList = list(full = list(Apres = Apres, Fpres = Fpres, Arres = Arres, Frres = Frres),
                downsampled = list(Apres= Apres.down, Fpres = Fpres.down, 
                                   Arres = Arres.down, Frres = Frres.down))
FracList = lapply(FracList, function(x) lapply(x, function(y) data.frame(y)))
lfc = lapply(FracList, function(x) lapply(x, function(y) 
  y[which(rownames(y) %in% map[which(map$sym=="ACTB" | map$sym=="MALAT1"),"ids"]),]))
x = list()
for (i in 1:2){x[[i]] = do.call(rbind, lfc[[i]])}
x[[1]]["Comparison"] = x[[2]]["Comparison"] = c("Adult:PolyA","Adult:PolyA","Prenatal:PolyA","Prenatal:PolyA",
                   "Adult:RiboZero","Adult:RiboZero","Prenatal:RiboZero","Prenatal:RiboZero")
x[[1]]["Gene"] = x[[2]]["Gene"] = rep.int(c("ACTB","MALAT1"),4)
write.table(x[[1]], "./Dropbox/sorted_figures/new/github_controlled/QC_section/data/MALAT1.ACTB.txt", 
            row.names=F, quote=F, sep="\t")
write.table(x[[2]], "./Dropbox/sorted_figures/new/github_controlled/QC_section/data/MALAT1.ACTB.downsampled.txt", 
            row.names=F, quote=F, sep="\t")

# Exclude the unpaired sample
tg = data.frame(t(testgenes[,which(colnames(testgenes)!="Br1113C1_RiboZero")]))
Cytosol <- pd[which(pd$Zone=="Cytosol"),]
Nucleus <- pd[,]
testnuc = tg[which(rownames(tg) %in% pd[which(pd$Zone=="Nucleus"),"SampleID"]),]
testcyt = tg[which(rownames(tg) %in% pd[which(pd$Zone=="Cytosol"),"SampleID"]),]

# one-tailed paired t test: ACTB
t.test(log(testnuc$ACTB), log(testcyt$ACTB), paired=TRUE, alternative = "less")
#data:  log(testnuc$ACTB) and log(testcyt$ACTB)
#t = -5.8346, df = 10, p-value = 8.251e-05
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.5340181
#sample estimates:
#  mean of the differences 
#-0.7746596

# one-tailed paired t test: FMR1
t.test(log(testnuc$FMR1), log(testcyt$FMR1), paired=TRUE, alternative = "less")
#data:  log(testnuc$FMR1) and log(testcyt$FMR1)
#t = -2.4166, df = 10, p-value = 0.01814
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.0338547
#sample estimates:
#  mean of the differences 
#-0.1354167 

# one-tailed paired t test in PolyA samples: FMR1
t.test(testnuc[c(1:2,4,6,8,10),"FMR1"], testcyt[c(1:2,4,6,8,10),"FMR1"], paired=TRUE, alternative = "less")
#data:  testnuc[c(1:2, 4, 6, 8, 10), "FMR1"] and testcyt[c(1:2, 4, 6, 8, 10), "FMR1"]
#t = -2.5587, df = 5, p-value = 0.02536
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.1752863
#sample estimates:
#  mean of the differences 
#-0.824947 

# one-tailed paired t test in females: XIST
t.test(testnuc[which(rownames(testnuc) %in% rownames(pd[which(pd$Sex=="F"),])),"XIST"], 
       testcyt[which(rownames(testcyt) %in% rownames(pd[which(pd$Sex=="F"),])),"XIST"], 
       paired=TRUE, alternative = "greater")
#data:  testnuc[which(rownames(testnuc) %in% rownames(pd[which(pd$Sex ==  and testcyt[which(rownames(testcyt) %in% rownames(pd[which(pd$Sex ==     "F"), ])), "XIST"] and     "F"), ])), "XIST"]
#t = 1.2936, df = 2, p-value = 0.1625
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -14.1975      Inf
#sample estimates:
#  mean of the differences 
#11.29229 

# one-tailed paired t test: MALAT1
t.test(testnuc$MALAT1, testcyt$MALAT1, paired=TRUE, alternative = "greater")
#data:  testnuc$MALAT1 and testcyt$MALAT1
#t = 2.6881, df = 10, p-value = 0.01139
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  33.81784      Inf
#sample estimates:
#  mean of the differences 
#103.8188

# Plot the Log2 Fold Change and SE by library for these genes
zone.res = list(list(Zpres = data.frame(Zpres), Zrres = data.frame(Zrres)), 
                list(Zpres.down = data.frame(Zpres.down), Zrres.down = data.frame(Zrres.down)))
zone.res = lapply(zone.res, function(x) lapply(x, function(y) y[which(rownames(y) %in% map$ids),]))
zone.res[[1]][[1]]["Gene"] = zone.res[[2]][[1]]["Gene"] = 
  zone.res[[1]][[2]]["Gene"] = zone.res[[2]][[2]]["Gene"] = map$sym
zone.res = lapply(zone.res, function(x) do.call(rbind,x))
zone.res[[1]]["Library"] = zone.res[[2]]["Library"] = c(rep.int("PolyA", 4), rep.int("RiboZero",4))

dodge <- position_dodge(width=0.9)
limits <- aes(ymax = log2FoldChange + lfcSE, ymin=log2FoldChange - lfcSE)
lapply(zone.res, function(x) {ggplot(x[which(x$Gene!="XIST" & x$Gene!="FMR1"),], 
       aes(x=Gene, y=log2FoldChange, fill=Library), color=Library) + 
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
        legend.key = element_rect(fill = "transparent", color = "transparent"))})