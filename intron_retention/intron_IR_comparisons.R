library(GenomicRanges)
library(data.table)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

### Differences by intron
## Differential IR by group
# read in results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
IRcomp = rat.1 = nonconst.66warn = nonconst.nowarn = nonconst = list()
for (i in 1:length(comps)){
  IRcomp[[i]] = read.table(paste0(path, comps[i], ".tab"), header = TRUE, comment.char="#")
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
  rat.1[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_0.1.tab"), header = TRUE, comment.char="#")
  nonconst.nowarn[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.nowarn.tab"), header = TRUE, comment.char="#")
  nonconst.66warn[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.66warn.tab"), header = TRUE, comment.char="#")
}
names(IRcomp) = names(nonconst) = names(rat.1) = names(nonconst.nowarn) = names(nonconst.66warn) = comps
dIR = list(IRcomp=IRcomp, nonconst=nonconst, rat.1=rat.1, nonconst.nowarn=nonconst.nowarn, nonconst.66warn=nonconst.66warn)

elementNROWS(IRcomp)
elementNROWS(nonconst)
elementNROWS(rat.1)
elementNROWS(nonconst.nowarn)
elementNROWS(nonconst.66warn)

string = lapply(dIR, function(x) lapply(x, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE)))
genes = lapply(string, function(x) lapply(x, function(y) y[grep("ENSG", y)]))
comments = lapply(string, function(x) lapply(x, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])))
IR.diff = lapply(dIR, function(x) lapply(x, function(y) y$A.IRratio - y$B.IRratio))
Sign = lapply(IR.diff, function(x) lapply(x, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")))
IR = list()
for (i in 1:length(dIR)){
  IR[[i]] = Map(cbind, dIR[[i]], ensID = genes[[i]], comments = comments[[i]], IR.diff = IR.diff[[i]], Sign = Sign[[i]])
}
names(IR) = names(dIR)
IRclean = lapply(IR, function(x) lapply(x, function(y) 
  y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & 
            y$A.warnings!="NonUniformIntronCover" & y$B.warnings!="LowCover" & 
            y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
            y$comments=="clean"),]))

### begin already gone through results from IR_analysis.R ###
# Explore the results
lapply(IRclean, function(x) elementNROWS(x))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1097              1001               689              1214               136               168 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1097              1001               689              1214               136               168 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#785               661               421               877               131               155 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#184               330                94               237                 5                 3 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#646               710               409               694                56                57
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1087               975               576               595               136               104 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1087               975               576               595               136               104 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#776               641               334               380               131                95 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#182               322                80               137                 5                 3 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#638               698               344               355                56                38 
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$p.diff<=0.05),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#314               213               171               351                99                95 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#314               213               171               351                99                95 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#268               192               141               323                97                93 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#57                30                18                48                 3                 2 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#203               134                83               191                42                36
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$p.diff<=0.05 & y$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#312               212               140               165                99                62 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#312               212               140               165                99                62 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#266               191               112               145                97                60 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#56                30                14                23                 3                 2 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#201               134                67                91                42                25
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 3                16                 3                 7 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 3                16                 3                 7 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 3                16                 3                 7 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#0                 0                 0                 1                 0                 0 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2                 1                 0                 9                 0                 1      
lapply(IRclean, function(x) elementNROWS(lapply(lapply(x, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),]), 
                                                function(z) z[which(z$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 1                 6                 3                 5 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 1                 6                 3                 5 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 1                 6                 3                 5 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#0                 0                 0                 1                 0                 0 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2                 1                 0                 4                 0                 1

#moving forward with the nonconstitutively spliced results because they capture more introns, 
#and coverage issues are filtered after differential retention calculation
nonconst = IRclean[["nonconst"]]
elementNROWS(lapply(nonconst, function(x) x[which(x$p.diff<=0.05),]))

# Comparison of significantly vs nonsignificantly retained introns by zone/age
fisher.test(data.frame(c(312,2),c(775,8))) 
# adult zone
#p-value = 0.7332
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3188043 15.6435734
#sample estimates:
#  odds ratio 
#1.609729
fisher.test(data.frame(c(212,1),c(763,25)))
# fetal zone
#p-value = 0.02666
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.122308 286.209119
#sample estimates:
#  odds ratio 
#6.938855
fisher.test(data.frame(c(141,30),c(435,83)))
# cytosol age
#p-value = 0.6352
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5572585 1.4726430
#sample estimates:
#  odds ratio 
#0.8969333 
fisher.test(data.frame(c(165,186),c(430,433)))
# nucleus age
#p-value = 0.3762
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6912077 1.1541152
#sample estimates:
#  odds ratio 
#0.8933697

# Get the DEG pval and LFC sign by Fraction for differentially retained introns
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down), Ares = data.frame(Ares), Fres.down = data.frame(Fres.down))
FracList = Map(cbind, FracList, lapply(FracList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
nonconst = Map(cbind, nonconst, 
               AP.sig = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"padj"]),
               AP.LFC = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"log2FoldChange"]),
               FP.sig = lapply(nonconst, function(x) FracList[["Fpres"]][match(x$ensID, FracList[["Fpres"]][,"ensemblID"]),"padj"]),
               FP.LFC = lapply(nonconst, function(x) FracList[["Fpres"]][match(x$ensID, FracList[["Fpres"]][,"ensemblID"]),"log2FoldChange"]),
               AR.sig = lapply(nonconst, function(x) FracList[["Arres"]][match(x$ensID, FracList[["Arres"]][,"ensemblID"]),"padj"]),
               AR.LFC = lapply(nonconst, function(x) FracList[["Arres"]][match(x$ensID, FracList[["Arres"]][,"ensemblID"]),"log2FoldChange"]),
               FR.sig = lapply(nonconst, function(x) FracList[["Frres"]][match(x$ensID, FracList[["Frres"]][,"ensemblID"]),"padj"]),
               FR.LFC = lapply(nonconst, function(x) FracList[["Frres"]][match(x$ensID, FracList[["Frres"]][,"ensemblID"]),"log2FoldChange"]),
               A.sig = lapply(nonconst, function(x) FracList[["Ares"]][match(x$ensID, FracList[["Ares"]][,"ensemblID"]),"padj"]),
               A.LFC = lapply(nonconst, function(x) FracList[["Ares"]][match(x$ensID, FracList[["Ares"]][,"ensemblID"]),"log2FoldChange"]),
               F.down.sig = lapply(nonconst, function(x) FracList[["Fres.down"]][match(x$ensID, FracList[["Fres.down"]][,"ensemblID"]),"padj"]),
               F.down.LFC = lapply(nonconst, function(x) FracList[["Fres.down"]][match(x$ensID, FracList[["Fres.down"]][,"ensemblID"]),"log2FoldChange"]))
elementNROWS(nonconst)
lapply(nonconst, head)

# Are genes with significantly differentially retained introns more likely to be significantly DEG by fraction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]>=0.05),])
fisher.test(data.frame(c(164,148), c(286,489))) #adult polya
#data:  data.frame(c(164, 148), c(286, 489))
#p-value = 2.539e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.439242 2.493470
#sample estimates:
#  odds ratio 
#1.893514
fisher.test(data.frame(c(42,167), c(47,717))) #prenatal polya
#data:  data.frame(c(42, 167), c(47, 717))
#p-value = 1.507e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.379715 6.151949
#sample estimates:
#  odds ratio 
#3.829787
fisher.test(data.frame(c(194,117), c(368,408))) #adult both libraries
#data:  data.frame(c(194, 117), c(368, 408))
#p-value = 8.927e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.392601 2.430641
#sample estimates:
#  odds ratio 
#1.837345
fisher.test(data.frame(c(42,167), c(47,717))) #prenatal both libraries
#data:  data.frame(c(42, 167), c(47, 717))
#p-value = 1.507e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.379715 6.151949
#sample estimates:
#  odds ratio 
#3.829787

# Are genes with significantly differentially retained introns more likely to have LFC in one direction by fraction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]<0),])
fisher.test(data.frame(c(147,65), c(447,335))) #prenatal polya
#data:  data.frame(c(147, 65), c(447, 335))
#p-value = 0.001547
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.212726 2.384151
#sample estimates:
#  odds ratio 
#1.694006
fisher.test(data.frame(c(151,162), c(439,341))) #adult polya
#data:  data.frame(c(151, 162), c(439, 341))
#p-value = 0.01873
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5517484 0.9501778
#sample estimates:
#  odds ratio 
#0.7242496
fisher.test(data.frame(c(149,164), c(453,327))) #adult in both libraries
#data:  data.frame(c(149, 164), c(453, 327))
#p-value = 0.001954
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4994930 0.8611575
#sample estimates:
#  odds ratio 
#0.6561047
fisher.test(data.frame(c(134,78), c(488,295))) #prenatal in both libraries
#data:  data.frame(c(134, 78), c(488, 295))
#p-value = 0.873
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7504722 1.4433322
#sample estimates:
#  odds ratio 
#1.038481

# Are genes with significantly differentially retained introns more likely to have a higher IR ratio in nuclear RNA?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #312
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #2
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #775
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #8
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #212
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #1
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #763
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #25
fisher.test(data.frame(c(312,2), c(775,8))) #adult
#data:  data.frame(c(312, 2), c(775, 8))
#p-value = 0.7332
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3188043 15.6435734
#sample estimates:
#  odds ratio 
#1.609729
fisher.test(data.frame(c(212,1), c(763,25))) #prenatal
#data:  data.frame(c(212, 1), c(763, 25))
#p-value = 0.02666
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.122308 286.209119
#sample estimates:
#  odds ratio 
#6.938855

# Get the DEG Age p-value and LFC sign for differentially retained introns
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres),
                Crres = data.frame(Crres), Nrres = data.frame(Nrres), 
                Cpres.down = data.frame(Cpres.down), Nres = data.frame(Nres), Cres.down = data.frame(Cres.down))
AgeList = Map(cbind, AgeList, lapply(AgeList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
nonconst = Map(cbind, nonconst, 
               CP.sig = lapply(nonconst, function(x) AgeList[["Cpres"]][match(x$ensID, AgeList[["Cpres"]][,"ensemblID"]),"padj"]),
               CP.LFC = lapply(nonconst, function(x) AgeList[["Cpres"]][match(x$ensID, AgeList[["Cpres"]][,"ensemblID"]),"log2FoldChange"]),
               NP.sig = lapply(nonconst, function(x) AgeList[["Npres"]][match(x$ensID, AgeList[["Npres"]][,"ensemblID"]),"padj"]),
               NP.LFC = lapply(nonconst, function(x) AgeList[["Npres"]][match(x$ensID, AgeList[["Npres"]][,"ensemblID"]),"log2FoldChange"]),
               CR.sig = lapply(nonconst, function(x) AgeList[["Crres"]][match(x$ensID, AgeList[["Crres"]][,"ensemblID"]),"padj"]),
               CR.LFC = lapply(nonconst, function(x) AgeList[["Crres"]][match(x$ensID, AgeList[["Crres"]][,"ensemblID"]),"log2FoldChange"]),
               NR.sig = lapply(nonconst, function(x) AgeList[["Nrres"]][match(x$ensID, AgeList[["Nrres"]][,"ensemblID"]),"padj"]),
               NR.LFC = lapply(nonconst, function(x) AgeList[["Nrres"]][match(x$ensID, AgeList[["Nrres"]][,"ensemblID"]),"log2FoldChange"]),
               CP.down.sig = lapply(nonconst, function(x) AgeList[["Cpres.down"]][match(x$ensID, AgeList[["Cpres.down"]][,"ensemblID"]),"padj"]),
               CP.down.LFC = lapply(nonconst, function(x) AgeList[["Cpres.down"]][match(x$ensID, AgeList[["Cpres.down"]][,"ensemblID"]),"log2FoldChange"]),
               N.sig = lapply(nonconst, function(x) AgeList[["Nres"]][match(x$ensID, AgeList[["Nres"]][,"ensemblID"]),"padj"]),
               N.LFC = lapply(nonconst, function(x) AgeList[["Nres"]][match(x$ensID, AgeList[["Nres"]][,"ensemblID"]),"log2FoldChange"]),
               C.down.sig = lapply(nonconst, function(x) AgeList[["Cres.down"]][match(x$ensID, AgeList[["Cres.down"]][,"ensemblID"]),"padj"]),
               C.down.LFC = lapply(nonconst, function(x) AgeList[["Cres.down"]][match(x$ensID, AgeList[["Cres.down"]][,"ensemblID"]),"log2FoldChange"]))
elementNROWS(nonconst)
lapply(nonconst, head)

# Are genes with significantly differentially retained introns more likely to be significantly DEG by age?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]>=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]>=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.sig"]>=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.sig"]>=0.05),])
fisher.test(data.frame(c(118,50), c(302,207))) #cytosol polya
#data:  data.frame(c(118, 50), c(302, 207))
#p-value = 0.0132
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.096525 2.406682
#sample estimates:
#  odds ratio 
#1.61649
fisher.test(data.frame(c(235,112), c(513,339))) #nucleus polya
#data:  data.frame(c(235, 112), c(513, 339))
#p-value = 0.01509
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.057376 1.823272
#sample estimates:
#  odds ratio 
#1.386161
fisher.test(data.frame(c(127,40), c(356,153))) #cytosol in both libraries
#data:  data.frame(c(127, 40), c(356, 153))
#p-value = 0.1393
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8994604 2.0982478
#sample estimates:
#  odds ratio 
#1.363929
fisher.test(data.frame(c(244,103), c(559,293))) #nucleus in both libraries
#data:  data.frame(c(244, 103), c(559, 293))
#p-value = 0.1199
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9405281 1.6451632
#sample estimates:
#  odds ratio 
#1.241454

# Are genes with significantly differentially retained introns more likely to have LFC in one direction by age?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]>0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]<0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]<05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]>0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.LFC"]>0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"C.down.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.LFC"]<0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.LFC"]<05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"N.LFC"]>0),])
fisher.test(data.frame(c(109,59), c(248,261))) #cytosol polya
#data:  data.frame(c(109, 59), c(248, 261))
#p-value = 0.0003479
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.335597 2.844027
#sample estimates:
#  odds ratio 
#1.942424
fisher.test(data.frame(c(200,148), c(853,392))) #nucleus polya
#data:  data.frame(c(200, 148), c(853, 392))
#p-value = 0.0001536
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4833510 0.7992836
#sample estimates:
#  odds ratio 
#0.6211912
fisher.test(data.frame(c(105,63), c(262,247))) #cytosol in both libraries
#data:  data.frame(c(105, 63), c(262, 247))
#p-value = 0.01577
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.083506 2.287870
#sample estimates:
#  odds ratio 
#1.570188
fisher.test(data.frame(c(206,142), c(853,354))) #nucleus in both libraries
#data:  data.frame(c(206, 142), c(853, 354))
#p-value = 6.697e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4670682 0.7776257
#sample estimates:
#  odds ratio 
#0.6022291

# Are genes with significantly differentially retained introns more likely to have a higher IR ratio in prenatal samples?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #140
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #31
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #436
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #82
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #165
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #186
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #430
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #433
fisher.test(data.frame(c(140,31), c(436,82))) #cytosol
#data:  data.frame(c(140, 31), c(436, 82))
#p-value = 0.4768
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5296991 1.3881400
#sample estimates:
#  odds ratio 
#0.8495898
fisher.test(data.frame(c(165,186), c(430,433))) #nucleus
#data:  data.frame(c(165, 186), c(430, 433))
#p-value = 0.3762
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6912077 1.1541152
#sample estimates:
#  odds ratio 
#0.8933697

# IR regulated genes 
DirList = lapply(nonconst, function(x) split(x, x$Sign))
IRlist = unlist(DirList, recursive = F)
names(IRlist) = names = c("Adult:Cytosol-Enriched","Adult:Nuclear-Enriched","Prenatal:Cytosol-Enriched","Prenatal:Nuclear-Enriched",
                  "Cytosol:Adult-Enriched","Cytosol:Prenatal-Enriched","Nucleus:Adult-Enriched","Nucleus:Prenatal-Enriched",
                  "Cytosol-Enriched","Nuclear-Enriched","Adult-Enriched","Prenatal-Enriched")
degs = list(Ares = as.data.frame(Ares), Fres.down = as.data.frame(Fres.down),Cres.down = as.data.frame(Cres.down), Nres = as.data.frame(Nres))
degs = lapply(degs, function(x) x[,1:6])
degs = Map(cbind, degs, Comparison = list("Adult","Prenatal","Cytosol","Nucleus"), Collapsed.Comparison = list("Fraction","Fraction","Age","Age"), 
           Sign = lapply(degs, function(x) ifelse(x$log2FoldChange>0, "Pos","Neg")),
           FDR = lapply(degs, function(x) ifelse(x$padj<=0.05, "FDR<0.05", "FDR>0.05")),
           ensID = lapply(degs, function(x) geneMap[match(as.character(rownames(x)),as.character(rownames(geneMap))),"ensemblID"]),
           Symbol = lapply(degs, function(x) geneMap[match(as.character(rownames(x)),as.character(rownames(geneMap))),"Symbol"]),
           EntrezID = lapply(degs, function(x) geneMap[match(as.character(rownames(x)),geneMap$gencodeID),"EntrezID"]),
           Type = lapply(degs, function(x) geneMap[match(as.character(rownames(x)),geneMap$gencodeID),"gene_type"]))
degs = do.call(rbind, degs)
degs = degs[which(degs$padj!="NA"),]
IRlist = lapply(IRlist, function(x) x[which(x$p.diff <=0.05),])
IRlist = lapply(IRlist, function(x) degs[which(degs$ensID %in% x$ensID),])
elementNROWS(IRlist)

# Plot DEG by fraction LFC and FDR for genes that contain introns that are differentially retained

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRgenes_byFraction.pdf", width=8, height=8)
for (i in 1:length(IRlist)){
  g = ggplot(IRlist[[i]][which(IRlist[[i]][,"Collapsed.Comparison"]=="Fraction"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("RNA Localization by IR Ratio:\n",names[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
}
dev.off()

# Plot DEG by age LFC and FDR for genes that contain introns that are differentially retained

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/IRgenes_byAge.pdf", width=8, height=8)
for (i in 1:length(IRlist)){
  g = ggplot(IRlist[[i]][which(IRlist[[i]][,"Collapsed.Comparison"]=="Age"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("Age Expression Changes by IR Ratio:\n",names[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
}
dev.off()

## Are the Fraction LFC values greater in genes containing a retained intron differentially?

t.test(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"], 
       IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"], alternative = "two.sided")
#data:  Adult:Cytosol-Enriched and Adult:Nuclear-Enriched
#t = 2.5258, df = 1.0291, p-value = 0.2343
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3.340238  5.142398
#sample estimates:
#  mean of x   mean of y 
#0.82133680 -0.07974323 
t.test(IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"], 
       IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"], alternative = "two.sided")
#data:  Prenatal:Cytosol-Enriched and Prenatal:Nuclear-Enriched
#not enough 'x' observations
t.test(c(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
         IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"]), 
       c(IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
         IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"]), alternative = "two.sided")
#data:  Combined  Cytosol-Enriched and Nuclear-Enriched
#t = 1.3287, df = 2.0234, p-value = 0.314
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.093380  2.086214
#sample estimates:
#  mean of x  mean of y 
#0.50996684 0.01354981
t.test(IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"], 
       IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"], alternative = "two.sided")
#data: Cytosol:Adult-Enriched and Cytosol:Prenatal-Enriched
#t = 5.8578, df = 39.413, p-value = 7.874e-07
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.9365837 1.9240148
#sample estimates:
#  mean of x  mean of y 
#0.6834974 -0.7468019
t.test(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"], 
       IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"], alternative = "two.sided")
#data:  Nucleus:Adult-Enriched and Nucleus:Prenatal-Enriched
#t = 9.7837, df = 270.56, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.9818406 1.4765372
#sample estimates:
#  mean of x  mean of y 
#0.2691135 -0.9600754
t.test(c(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
         IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"]), 
       c(IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
         IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"]), alternative = "two.sided")
#data:  Combined Adult-Enriched and Combined Prenatal-Enriched
#t = 11.608, df = 457.07, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.986928 1.389202
#sample estimates:
#  mean of x  mean of y 
#0.3267855 -0.8612796 

## Are the Fraction pvalues values more significant in genes containing a retained intron preferentially?

t.test(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"padj"], 
       IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"padj"], alternative = "two.sided")
#data:  Adult:Cytosol-Enriched and Adult:Nuclear-Enriched
#t = -10.085, df = 265.99, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2104949 -0.1417319
#sample estimates:
#  mean of x   mean of y 
#0.001126789 0.177240174 
t.test(IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"padj"], 
       IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"padj"], alternative = "two.sided")
#data:  Prenatal:Cytosol-Enriched and Prenatal:Nuclear-Enriched
#not enough 'x' observations
t.test(c(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"padj"],
         IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"padj"]), 
       c(IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"padj"],
         IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"padj"]), alternative = "two.sided")
#data:  Combined  Cytosol-Enriched and Nuclear-Enriched
#t = 0.16058, df = 2.0111, p-value = 0.8871
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.117272  1.204380
#sample estimates:
#  mean of x mean of y 
#0.2719766 0.2284224
t.test(IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"padj"], 
       IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"padj"], alternative = "two.sided")
#data: Cytosol:Adult-Enriched and Cytosol:Prenatal-Enriched
#t = 0.37388, df = 36.487, p-value = 0.7107
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.08262087  0.11998935
#sample estimates:
#  mean of x  mean of y 
#0.11864635 0.09996211 
t.test(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"padj"], 
       IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"padj"], alternative = "two.sided")
#data:  Nucleus:Adult-Enriched and Nucleus:Prenatal-Enriched
#t = 2.3936, df = 310.23, p-value = 0.01728
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.01360126 0.13926969
#sample estimates:
#  mean of x mean of y 
#0.1881114 0.1116759
t.test(c(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
         IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"padj"]), 
       c(IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
         IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"padj"]), alternative = "two.sided")
#data:  Combined Adult-Enriched and Combined Prenatal-Enriched
#t = 2.7719, df = 357.99, p-value = 0.005865
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.02097326 0.12341456
#sample estimates:
#  mean of x mean of y 
#0.1784436 0.1062497

## Are dIR introns by age more likely to be dIR introns by fraction?
nonconst = Map(cbind, IRclean[["nonconst"]], intronID = lapply(IRclean[["nonconst"]], function(x) paste0(x$Chr, ":", x$Start, "-", x$"End")))
sigIR.intron = lapply(nonconst, function(x) as.character(x[which(x$p.diff<=0.05),"intronID"])) 
nonsigIR.intron = lapply(nonconst, function(x) as.character(x[which(x$p.diff>0.05),"intronID"]))
sigIR.intron.combined = list("Significantly IR\nBy Fraction" = c(sigIR.intron[["Adult_PolyA_Zone"]], sigIR.intron[["Fetal_PolyA_Zone"]]),
                             "Significantly IR\nBy Age" = c(sigIR.intron[["Cytosol_PolyA_Age"]], sigIR.intron[["Nuclear_PolyA_Age"]]))
nonsigIR.intron.combined = list("Non-Significantly IR\nBy Fraction" = c(nonsigIR.intron[["Adult_PolyA_Zone"]], nonsigIR.intron[["Fetal_PolyA_Zone"]]),
                             "Non-Significantly IR\nBy Age" = c(nonsigIR.intron[["Cytosol_PolyA_Age"]], nonsigIR.intron[["Nuclear_PolyA_Age"]]))

venn.diagram(c(sigIR.intron[c(1,3)], nonsigIR.intron[c(1,3)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_adult_cytosol.jpeg", 
             main="dIR_FractionbyAge_adult_cytosol", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.intron[c(1,4)], nonsigIR.intron[c(1,4)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_adult_nucleus.jpeg", 
             main="dIR_FractionbyAge_adult_nucleus", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.intron[c(1:2)], nonsigIR.intron[c(1:2)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_Fraction_adult_fetal.jpeg", 
             main="dIR_Fraction_adult_fetal", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.intron[c(2,3)], nonsigIR.intron[c(2,3)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_prenatal_cytosol.jpeg", 
             main="dIR_FractionbyAge_prenatal_cytosol", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.intron[c(2,4)], nonsigIR.intron[c(2,4)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_prenatal_nucleus.jpeg", 
             main="dIR_FractionbyAge_prenatal_nucleus", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.intron[c(3,4)], nonsigIR.intron[c(3,4)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_Age_cytosol_nucleus.jpeg", 
             main="dIR_Age_cytosol_nucleus", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.intron.combined, nonsigIR.intron.combined), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_combined.jpeg", 
             main="dIR_FractionbyAge_combined", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
# all 7 saved together at ./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_overlap_Fraction_Age.pdf

fisher.test(data.frame(c(1+4+10+85,5+15+78+285),c(6+20+51+322,322+927+1088))) 
#data: dIR_FractionbyAge_combined.jpeg  
#  p-value = 0.001026
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.184734 1.961597
#sample estimates:
#  odds ratio 
#1.529056
fisher.test(data.frame(c(35,11+125),c(19+260,52+447+720)))
#data: dIR_FractionbyAge_adult_cytosol.jpeg
#p-value = 0.5371
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7354046 1.6826860
#sample estimates:
#  odds ratio 
#1.12436
fisher.test(data.frame(c(5,19+147),c(6+202,26+486+743)))
#data: dIR_FractionbyAge_prenatal_cytosol.jpeg
#p-value = 3.313e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.05758012 0.44079056
#sample estimates:
#  odds ratio 
#0.1818474
fisher.test(data.frame(c(31,34+286),c(21+262,111+731+638)))
#data: dIR_FractionbyAge_adult_nucleus.jpeg
#p-value = 0.0003822
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3313704 0.7524600
#sample estimates:
#  odds ratio 
#0.5067602
fisher.test(data.frame(c(41,49+261),c(19+153,170+674+569)))
#data: dIR_FractionbyAge_prenatal_nucleus.jpeg
#p-value = 0.6382
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7367658 1.5731319
#sample estimates:
#  odds ratio 
#1.086481
fisher.test(data.frame(c(28,14+171),c(23+263,41+724+728)))
#data: dIR_Fraction_adult_fetal.jpeg
#p-value = 0.3194
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5009596 1.2076840
#sample estimates:
#  odds ratio 
#0.790188
fisher.test(data.frame(c(39,20+292),c(5+127,50+808+448)))
#data: dIR_Age_cytosol_nucleus.jpeg
#p-value = 0.2665
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8242523 1.8220293
#sample estimates:
#  odds ratio 
#1.236588

# What about only the introns reported in both comparisons?
fisher.test(data.frame(c(1+4+10+85,5+78),c(6+51,322))) 
#data: dIR_FractionbyAge_combined.jpeg  
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  4.449856 10.415955
#sample estimates:
#  odds ratio 
#6.777491
fisher.test(data.frame(c(35,11),c(19,52)))
#data: dIR_FractionbyAge_adult_cytosol.jpeg
#p-value = 1.851e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.431975 22.659230
#sample estimates:
#  odds ratio 
#8.51932
fisher.test(data.frame(c(5,19),c(6,26)))
#data: dIR_FractionbyAge_prenatal_cytosol.jpeg
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2367084 5.2428228
#sample estimates:
#  odds ratio 
#1.137662
fisher.test(data.frame(c(31,34),c(21,111)))
#data: dIR_FractionbyAge_adult_nucleus.jpeg
#p-value = 4.731e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.327024 10.006823
#sample estimates:
#  odds ratio 
#4.77432
fisher.test(data.frame(c(41,49),c(19,170)))
#data: dIR_FractionbyAge_prenatal_nucleus.jpeg
#p-value = 7.144e-11
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.820863 14.858858
#sample estimates:
#  odds ratio 
#7.419343
fisher.test(data.frame(c(28,14),c(23,41)))
#data: dIR_Fraction_adult_fetal.jpeg
#p-value = 0.00277
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.459885 8.825248
#sample estimates:
#  odds ratio 
#3.519946
fisher.test(data.frame(c(39,20),c(5,50)))
#data: dIR_Age_cytosol_nucleus.jpeg
#p-value = 1.683e-10
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  6.250532 70.547304
#sample estimates:
#  odds ratio 
#18.87679


## Are genes that have a dIR intron by age more likely to have a dIR intron by fraction?
sigIR.gene = lapply(nonconst, function(x) unique(as.character(x[which(x$p.diff<=0.05),"ensID"])))
names(sigIR.gene) = paste0(names(sigIR.gene),"\nsig")
nonsigIR.gene = lapply(nonconst, function(x) unique(as.character(x[which(x$p.diff>0.05),"ensID"])))
names(nonsigIR.gene) = paste0(names(nonsigIR.gene),"\nnonsig")
sigIR.gene.combined = list("Significantly IR\nBy Fraction" = c(sigIR.gene[["Adult_PolyA_Zone\nsig"]], sigIR.gene[["Fetal_PolyA_Zone\nsig"]]),
                             "Significantly IR\nBy Age" = c(sigIR.gene[["Cytosol_PolyA_Age\nsig"]], sigIR.gene[["Nuclear_PolyA_Age\nsig"]]))
nonsigIR.gene.combined = list("Non-Significantly IR\nBy Fraction" = c(nonsigIR.gene[["Adult_PolyA_Zone\nnonsig"]], nonsigIR.gene[["Fetal_PolyA_Zone\nnonsig"]]),
                                "Non-Significantly IR\nBy Age" = c(nonsigIR.gene[["Cytosol_PolyA_Age\nnonsig"]], nonsigIR.gene[["Nuclear_PolyA_Age\nnonsig"]]))

venn.diagram(c(sigIR.gene[c(1,3)], nonsigIR.gene[c(1,3)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_FractionbyAge_adult_cytosol.jpeg", 
             main="dIR_FractionbyAge_adult_cytosol", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.gene[c(1,4)], nonsigIR.gene[c(1,4)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_FractionbyAge_adult_nucleus.jpeg", 
             main="dIR_FractionbyAge_adult_nucleus", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.gene[c(1:2)], nonsigIR.gene[c(1:2)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_Fraction_adult_fetal.jpeg", 
             main="dIR_Fraction_adult_fetal", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.gene[c(2,3)], nonsigIR.gene[c(2,3)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_FractionbyAge_prenatal_cytosol.jpeg", 
             main="dIR_FractionbyAge_prenatal_cytosol", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.gene[c(2,4)], nonsigIR.gene[c(2,4)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_FractionbyAge_prenatal_nucleus.jpeg", 
             main="dIR_FractionbyAge_prenatal_nucleus", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.gene[c(3,4)], nonsigIR.gene[c(3,4)]), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_Age_cytosol_nucleus.jpeg", 
             main="dIR_Age_cytosol_nucleus", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(c(sigIR.gene.combined, nonsigIR.gene.combined), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_gene_FractionbyAge_combined.jpeg", 
             main="dIR_FractionbyAge_combined", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
# all 7 saved together at ./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_byGene_overlap_Fraction_Age.pdf


fisher.test(data.frame(c(21+17+29+53,35+74+46+142),c(51+39+50+153,368+663+529))) 
#data: dIR_gene_FractionbyAge_combined.jpeg  
#p-value = 3.879e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.665630 2.767285
#sample estimates:
#  odds ratio 
#2.150367
fisher.test(data.frame(c(5+3+4+25,19+4+12+85),c(13+30+28+163,101+299+521)))
#data: dIR_gene_FractionbyAge_adult_cytosol.jpeg
#p-value = 0.3448
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7933247 1.8224939
#sample estimates:
#  odds ratio 
#1.213404
fisher.test(data.frame(c(0+1+1+11,24+8+15+97),c(6+13+24+124,64+360+578)))
#data: dIR_gene_FractionbyAge_prenatal_cytosol.jpeg
#p-value = 0.04613
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2751708 0.9850721
#sample estimates:
#  odds ratio 
#0.5418933
fisher.test(data.frame(c(6+2+10+36,46+15+34+164),c(13+29+21+154,153+504+431)))
#data: dIR_gene_FractionbyAge_adult_nucleus.jpeg
#p-value = 0.8005
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7384812 1.4615272
#sample estimates:
#  odds ratio 
#1.045328
fisher.test(data.frame(c(6+4+8+32,8+22+9+91),c(29+18+53+163,182+487+421)))
#data: dIR_gene_FractionbyAge_prenatal_nucleus.jpeg
#p-value = 0.01352
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.095534 2.291578
#sample estimates:
#  odds ratio 
#1.593495
fisher.test(data.frame(c(3+2+8+25,27+7+19+89),c(16+34+23+160,95+529+516)))
#data: dIR_gene_Fraction_adult_fetal.jpeg
#p-value = 0.1748
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8659325 1.9426236
#sample estimates:
#  odds ratio 
#1.309072
fisher.test(data.frame(c(4+7+6+36,32+7+39+182),c(10+3+11+80,103+583+301)))
#data: dIR_gene_Age_cytosol_nucleus.jpeg
#p-value = 0.0004928
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.323990 2.800034
#sample estimates:
#  odds ratio 
#1.933604

# What about only the genes reported in both comparisons?
fisher.test(data.frame(c(21+17+29+53,35+74+46+142),c(51+39+50+153,368+663+529))) 
#data: dIR_gene_FractionbyAge_combined.jpeg  
#p-value = 3.879e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.665630 2.767285
#sample estimates:
#  odds ratio 
#2.150367
fisher.test(data.frame(c(5+3+4+25,19+4+12+85),c(13+30+28+163,101+299+521)))
#data: dIR_gene_FractionbyAge_adult_cytosol.jpeg
#p-value = 0.3448
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7933247 1.8224939
#sample estimates:
#  odds ratio 
#1.213404
fisher.test(data.frame(c(0+1+1+11,24+8+15+97),c(6+13+24+124,64+360+578)))
#data: dIR_gene_FractionbyAge_prenatal_cytosol.jpeg
#p-value = 0.04613
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2751708 0.9850721
#sample estimates:
#  odds ratio 
#0.5418933
fisher.test(data.frame(c(6+2+10+36,46+15+34+164),c(13+29+21+154,153+504+431)))
#data: dIR_gene_FractionbyAge_adult_nucleus.jpeg
#p-value = 0.8005
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7384812 1.4615272
#sample estimates:
#  odds ratio 
#1.045328
fisher.test(data.frame(c(6+4+8+32,8+22+9+91),c(29+18+53+163,182+487+421)))
#data: dIR_gene_FractionbyAge_prenatal_nucleus.jpeg
#p-value = 0.01352
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.095534 2.291578
#sample estimates:
#  odds ratio 
#1.593495
fisher.test(data.frame(c(3+2+8+25,27+7+19+89),c(16+34+23+160,95+529+516)))
#data: dIR_gene_Fraction_adult_fetal.jpeg
#p-value = 0.1748
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8659325 1.9426236
#sample estimates:
#  odds ratio 
#1.309072
fisher.test(data.frame(c(4+7+6+36,32+7),c(10+3,103)))
#data: dIR_gene_Age_cytosol_nucleus.jpeg
#p-value = 7.898e-13
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.059531 23.696676
#sample estimates:
#  odds ratio 
#10.6203

# What about only the genes reported in both comparisons?
fisher.test(data.frame(c(21+17+29+53,35+74),c(51+39,368))) 
#data: dIR_gene_FractionbyAge_combined.jpeg  
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.133920 6.463821
#sample estimates:
#  odds ratio 
#4.490156
fisher.test(data.frame(c(5+3+4+25,19+4),c(13+30,101)))
#data: dIR_gene_FractionbyAge_adult_cytosol.jpeg
#p-value = 3.532e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.919544 7.474978
#sample estimates:
#  odds ratio 
#3.751376
fisher.test(data.frame(c(0+1+1+11,24+8),c(6+13,64)))
#data: dIR_gene_FractionbyAge_prenatal_cytosol.jpeg
#p-value = 0.523
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5457205 3.3482308
#sample estimates:
#  odds ratio 
#1.364944
fisher.test(data.frame(c(6+2+10+36,46+15),c(13+29,153)))
#data: dIR_gene_FractionbyAge_adult_nucleus.jpeg
#p-value = 4.187e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.894974 5.487542
#sample estimates:
#  odds ratio 
#3.21149
fisher.test(data.frame(c(6+4+8+32,8+22),c(18+53,182)))
#data: dIR_gene_FractionbyAge_prenatal_nucleus.jpeg
#p-value = 5.295e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.436011 7.527893
#sample estimates:
#  odds ratio 
#4.250734
fisher.test(data.frame(c(3+2+8+25,27+7),c(16+34,95)))
#data: dIR_gene_Fraction_adult_fetal.jpeg
#p-value = 0.1748
#p-value = 0.01247
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.146241 3.931282
#sample estimates:
#  odds ratio 
#2.115891
fisher.test(data.frame(c(4+7+6+36,32+7),c(10+3,103)))
#data: dIR_gene_Age_cytosol_nucleus.jpeg
#p-value = 7.898e-13
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.059531 23.696676
#sample estimates:
#  odds ratio 
#10.6203