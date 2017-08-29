library(GenomicRanges)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

## read in Differential IR results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
nonconst = list()
for (i in 1:length(comps)){
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
}
names(nonconst) = comps
elementNROWS(nonconst)
string = lapply(nonconst, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE))
IR.diff = lapply(nonconst, function(y) y$A.IRratio - y$B.IRratio)
dIR = Map(cbind, nonconst, ensID = lapply(string, function(y) y[grep("ENSG", y)]), 
          comments = lapply(string, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])), 
          IR.diff = IR.diff, Sign = lapply(IR.diff, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")), 
          intronID = lapply(nonconst, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End)))
dIRclean = lapply(dIR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" &
                                             y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                                             y$comments=="clean"),])
sigdIR = lapply(dIRclean, function(x) x[which(x$p.diff<=0.05),])
sigdIR = lapply(sigdIR, function(x) split(x, x$Sign))
sigdIR = unlist(sigdIR, recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                  "Cytosol-Increased", "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)]))
introns = c("Introns (Fraction)" = list(do.call(rbind, lapply(dIRclean[1:2], function(x) x[,c(1:4,6,8,30,32,33,34)]))), 
            "Introns (Age)" = list(do.call(rbind, lapply(dIRclean[3:4], function(x) x[,c(1:4,6,8,30,32,33,34)]))), sigdIR)
introns = lapply(introns[c(1:10,12:14)], function(x) data.frame(Chr = paste0("chr",x$Chr), x[2:10]))

### How many genes are represented by how many retained introns?
lapply(introns, head)
introns = c(introns, list(Cytosolic = rbind(introns[["Adult:Cytosol-Increased"]],introns[["Prenatal:Cytosol-Increased"]]),
                          Nuclear = rbind(introns[["Adult:Nucleus-Increased"]],introns[["Prenatal:Nucleus-Increased"]]),
                          Adult = rbind(introns[["Cytosol:Adult-Increased"]],introns[["Nucleus:Adult-Increased"]]),
                          Prenatal = rbind(introns[["Cytosol:Prenatal-Increased"]],introns[["Nucleus:Prenatal-Increased"]])))
genes = data.frame(numgenes = unlist(lapply(introns, function(x) length(unique(x$ensID)))),
                   numintrons = unlist(lapply(introns, function(x) length(unique(x$intronID)))))
genes$group = c("All", "All", "Cytosolic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytosolic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                "Adult Retention\nIn Cytosol","Prenatal Retention\nIn Cytosol","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus",
                "Nuclear Retention", "Adult Retention", "Prenatal Retention", "Pooled Cytosolic\nRetention","Pooled Nuclear\nRetention",
                "Pooled Adult\nRetention","Pooled Prenatal\nRetention")
genes$retention = c("NA","NA","Cytosolic","Nuclear","Cytosolic","Nuclear","Adult","Prenatal","Adult","Prenatal","Nuclear","Adult", "Prenatal","Cytosolic","Nuclear","Adult","Prenatal")

## by Fraction: cytosolic (is there a relationship between the number of introns found dIR and the number of genes?)
fisher.test(data.frame(c(3,3),c(1164-3,1407-3)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1616479 9.0475613
#sample estimates:
#  odds ratio 
#1.209209
## by Fraction: nuclear
fisher.test(data.frame(c(218,250),c(1164-218,1407-250)))
#p-value = 0.5381
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8675996 1.3103089
#sample estimates:
#  odds ratio 
#1.066454
## by Fraction: cytosolic vs nuclear
fisher.test(data.frame(c(3,3),c(218,250)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1519708 8.6501791
#sample estimates:
#  odds ratio 
#1.146455
## by Fraction in adults: cytosolic
fisher.test(data.frame(c(2,2),c(1164-2,1407-2)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.08753648 16.69911644
#sample estimates:
#  odds ratio 
#1.209028
## by Fraction in adults: nuclear
fisher.test(data.frame(c(143,160),c(1164-143,1407-160)))
#p-value = 0.4993
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.851767 1.397904
#sample estimates:
#  odds ratio 
#1.091545
## by Fraction in adults: cytosolic vs nuclear
fisher.test(data.frame(c(2,2),c(143,160)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.08013867 15.61042299
#sample estimates:
#  odds ratio 
#1.11847
## by Fraction in prenatal: cytosolic
fisher.test(data.frame(c(1,1),c(1164-1,1407-1)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0153953 94.8975538
#sample estimates:
#  odds ratio 
#1.208849
## by Fraction in prenatal: nuclear
fisher.test(data.frame(c(85,92),c(1164-85,1407-92)))
#p-value = 0.4815
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8189755 1.5465049
#sample estimates:
#  odds ratio 
#1.125939
## by Fraction in prenatal: cytosolic vs nuclear
fisher.test(data.frame(c(1,1),c(85,92)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.01363748 85.81742953
#sample estimates:
#  odds ratio 
#1.081855
## by Age: adult
fisher.test(data.frame(c(110,116),c(973-110,1145-116)))
#p-value = 0.3973
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8492254 1.5045466
#sample estimates:
#  odds ratio 
#1.130612
## by Age: prenatal
fisher.test(data.frame(c(120,131),c(973-120,1145-131)))
#p-value = 0.5442
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8285683 1.4300551
#sample estimates:
#  odds ratio 
#1.088883
## by Age: adult vs prenatal
fisher.test(data.frame(c(110,116),c(120,131)))
#p-value = 0.855
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.710850 1.507415
#sample estimates:
#  odds ratio 
#1.035145
## by Age in cytosol: adult
fisher.test(data.frame(c(17,17),c(973-17,1145-17)))
#p-value = 0.7293
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.562748 2.473472
#sample estimates:
#  odds ratio 
#1.179821
## by Age in cytosol: prenatal
fisher.test(data.frame(c(60,63),c(973-60,1145-63)))
#p-value = 0.5159
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7701461 1.6525699
#sample estimates:
#  odds ratio 
#1.128604 
## by Age in cytosol: adult vs prenatal
fisher.test(data.frame(c(17,17),c(60,63)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4569037 2.4118642
#sample estimates:
#  odds ratio 
#1.049699
## by Age in nucleus: adult
fisher.test(data.frame(c(98,102),c(973-98,1145-102)))
#p-value = 0.3717
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8461161 1.5493104
#sample estimates:
#  odds ratio 
#1.145179
## by Age in nucleus: prenatal
fisher.test(data.frame(c(72,77),c(973-72,1145-77)))
#p-value = 0.5518
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7823529 1.5686686
#sample estimates:
#  odds ratio 
#1.108324
## by Age in nucleus: adult vs prenatal
fisher.test(data.frame(c(72,77),c(98,102)))
#p-value = 0.9142
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.622518 1.521118
#sample estimates:
#  odds ratio 
#0.973297

## How many introns per gene represented?
pergeneIR = lapply(introns, function(x) data.frame(numIR=count(x$ensID)))
stat = data.frame(mean = unlist(lapply(pergeneIR, function(x) mean(x$numIR.freq))), sd = unlist(lapply(pergeneIR,function(x) mean(x$numIR.freq))))
morethan1 = lapply(pergeneIR, function(x) x[which(x$numIR.freq>1),])
elementNROWS(morethan1)
elementNROWS(pergeneIR)

## number of genes with more than one retained intron by Fraction in adults: cytosolic
fisher.test(data.frame(c(2,0),c(1164-2,208-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.03350136        Inf
#sample estimates:
#  odds ratio 
#Inf
## by Fraction in adults: nuclear
fisher.test(data.frame(c(143-14,14),c(1164-208-143+14,208-14)))
#p-value = 0.006989
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.208367 4.153066
#sample estimates:
#  odds ratio 
#2.160386
## by Fraction in adults: cytosolic vs nuclear
fisher.test(data.frame(c(0,2),c(14,143-14)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00000 51.17491
#sample estimates:
#  odds ratio 
#0
## by Fraction in prenatal: cytosolic
fisher.test(data.frame(c(1,0),c(1164-208-1,208-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.005594338         Inf
#sample estimates:
#  odds ratio 
#Inf
## by Fraction in prenatal: nuclear
fisher.test(data.frame(c(85-6,6),c(1164-208-85+6,208-6)))
#p-value = 0.004889
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.306040 8.626995
#sample estimates:
#  odds ratio 
#3.030696 
## by Fraction in prenatal: cytosolic vs nuclear
fisher.test(data.frame(c(0,1),c(6,85-6)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.000 515.817
#sample estimates:
#  odds ratio 
#0
## number of genes with more than one retained intron by Age in cytosol: adult
fisher.test(data.frame(c(17,0),c(973-179-17,179-0)))
#p-value = 0.05415
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9392915       Inf
#sample estimates:
#  odds ratio 
#Inf
## by Age in cytosol: prenatal
fisher.test(data.frame(c(60-3,3),c(973-179-60+3,179-3)))
#p-value = 0.003167
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.447555 22.889346
#sample estimates:
#  odds ratio 
#4.532908
## by Age in cytosol: adult vs prenatal
fisher.test(data.frame(c(17,0),c(60-3,3)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.11463     Inf
#sample estimates:
#  odds ratio 
#Inf
## by Age in nucleus: adult
fisher.test(data.frame(c(98-2,2),c(973-98-179+2,179-2)))
#p-value = 3.463e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.217669 102.846112
#sample estimates:
#  odds ratio 
#12.15652
## by Age in nucleus: prenatal
fisher.test(data.frame(c(72-5,5),c(973-72-179+5,179-5)))
#p-value = 0.00681
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.27900 10.34251
#sample estimates:
#  odds ratio 
#3.204523
## by Age in nucleus: adult vs prenatal
fisher.test(data.frame(c(98-2,2),c(72-5,5)))
#p-value = 0.1352
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5618737 38.3894865
#sample estimates:
#  odds ratio 
#3.555643

### Are genes with more than one retained intron more likely to be differentially expressed by fraction or age?
lapply(morethan1, head)
fracList = list(Ares = as.data.frame(Ares), Fres = as.data.frame(Fres), Cres = as.data.frame(Cres), Nres = as.data.frame(Nres))
fracList = Map(cbind, fracList, ensID = lapply(fracList, function(x) geneMap[match(rownames(x), geneMap$gencodeID),"ensemblID"]))
lapply(fracList,head)
length(morethan1)
length(fracList)
res = list(list(),list(),list(),list())
for (i in 1:length(fracList)){
  fr = fracList[[i]]
  for (j in 1:length(morethan1)){
    more = morethan1[[j]]
    res[[i]][[j]] = fr[which(more$numIR.x %in% fr$ensID),]
  }
}
names(res) = names(fracList)
names(res[[1]]) = names(res[[2]]) = names(res[[3]]) = names(res[[4]]) = names(morethan1)
elementNROWS(morethan1)
lapply(res, elementNROWS)
lapply(lapply(res, function(x) lapply(x, function(y) y[which(y$padj<=0.05),])), elementNROWS)

## In genes with more than one intron, is there a relationship between having an intron be dIR and the gene being DE? 
#by Fraction in adults: cytosolic
fisher.test(data.frame(c(67-0,0),c(207-67-0,0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
## by Fraction in adults: nuclear
fisher.test(data.frame(c(67-2,2),c(207-67-14+2,14-2)))
#p-value = 0.2346
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6461105 28.7194481
#sample estimates:
#  odds ratio 
#3.033725
## by Fraction in adults: cytosolic vs nuclear
fisher.test(data.frame(c(2,0),c(14-2,0-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
## by Fraction in prenatal: cytosolic
fisher.test(data.frame(c(24-0,0),c(207-24-0,0-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
## by Fraction in prenatal: nuclear
fisher.test(data.frame(c(24-1,1),c(207-24-6+1,6-1)))
#p-value = 0.5272
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.06798533 31.87023517
#sample estimates:
#  odds ratio 
#0.6476653
## by Fraction in prenatal: cytosolic vs nuclear
fisher.test(data.frame(c(1,0),c(6-1,0-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
## number of genes with more than one retained intron by Age in cytosol: adult
fisher.test(data.frame(c(81-0,0),c(178-81-0,0-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
## by Age in cytosol: prenatal
fisher.test(data.frame(c(81-0,0),c(178-81-3,3-0)))
#p-value = 0.2519
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3466785       Inf
#sample estimates:
#  odds ratio 
#Inf 
## by Age in cytosol: adult vs prenatal
fisher.test(data.frame(c(0,0),c(0-0,3-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
## by Age in nucleus: adult
fisher.test(data.frame(c(70-0,0),c(178-70-2,2-0)))
#p-value = 0.5201
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1216475       Inf
#sample estimates:
#  odds ratio 
#Inf
## by Age in nucleus: prenatal
fisher.test(data.frame(c(70-0,0),c(178-70-5,5-0)))
#p-value = 0.1582
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6016422       Inf
#sample estimates:
#  odds ratio 
#Inf
## by Age in nucleus: adult vs prenatal
fisher.test(data.frame(c(0,0),c(2-0,5-0)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0