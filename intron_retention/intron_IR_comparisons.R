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
            y$B.warnings!="NonUniformIntronCover" & y$B.warnings!="LowCover" & 
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

### pick up here ###

fracDevel = lapply(sigIR, function(x) total[which(total$ensID %in% x$ensID),])
pdf("/Users/amandaprice/Dropbox/sorted_figures/new/FracDevel_plots.pdf", width=8, height=8)
names = c("Both Nuclear","Both Cytosolic", "Nuclear in Fetal Only",
          "Nuclear in Adult Only", "Nuclear in Adult/Cytosolic in Fetal",
          "Nuclear in Fetal/Cytosolic in Adult", "Interaction Effect")
for (i in 1:length(fracDevel)){
  g = ggplot(fracDevel[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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

# IR regulated genes 
DirList = lapply(IRclean, function(x) split(x, x$Sign))
IRlist = list("Adult:Nucleus"=DirList[["Adult_PolyA_Zone"]][["MoreIRInNuc.Fetal"]],
              "Adult:Cytosol"=DirList[["Adult_PolyA_Zone"]][["MoreIRInCyt.Adult"]],
              "Fetal:Nucleus"=DirList[["Fetal_PolyA_Zone"]][["MoreIRInNuc.Fetal"]],
              "Fetal:Cytosol"=DirList[["Fetal_PolyA_Zone"]][["MoreIRInCyt.Adult"]],
              "Cytosol:Fetal"=DirList[["Cytosol_PolyA_Age"]][["MoreIRInNuc.Fetal"]],
              "Cytosol:Adult"=DirList[["Cytosol_PolyA_Age"]][["MoreIRInCyt.Adult"]],
              "Nucleus:Fetal"=DirList[["Nuclear_PolyA_Age"]][["MoreIRInNuc.Fetal"]],
              "Nucleus:Adult"=DirList[["Nuclear_PolyA_Age"]][["MoreIRInCyt.Adult"]],
              "Nuclear-Enriched"=DirList[["PolyA_Zone"]][["MoreIRInNuc.Fetal"]],
              "Fetal-Enriched"=DirList[["PolyA_Age"]][["MoreIRInNuc.Fetal"]],
              "Adult-Enriched"=DirList[["PolyA_Age"]][["MoreIRInCyt.Adult"]])

totalF = rbind(xA,xF)
totalF$FDR = ifelse(totalF$padj<=0.05, "FDR<0.05", "FDR>0.05")
totalF$IR = factor(ifelse(totalF$IRratio>=0.5, ">0.5", "<0.5"))
totalF = totalF[which(totalF$padj!="NA"),]
combined.IRlist = list(nuc.A.F = rbind(IRlist[["Adult:Nucleus"]], IRlist[["Fetal:Nucleus"]]),
                       cyt.A.F = rbind(IRlist[["Adult:Cytosol"]], IRlist[["Fetal:Cytosol"]]),
                       adult.C.N = rbind(IRlist[["Cytosol:Adult"]], IRlist[["Nucleus:Adult"]]),
                       fetal.C.N = rbind(IRlist[["Cytosol:Fetal"]], IRlist[["Nucleus:Fetal"]]))
IRfrac = lapply(IRlist, function(x) totalF[which(totalF$EnsID %in% x$ensID),])
IRfrac2 = lapply(combined.IRlist, function(x) totalF[which(totalF$EnsID %in% x$ensID),])

pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byFrac_plots.pdf", width=8, height=8)
names = names(IRfrac)
for (i in 1:length(IRfrac)){
  g = ggplot(IRfrac[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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
pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byFrac_combined.pdf", width=8, height=8)
names = names(IRfrac2)
for (i in 1:length(IRfrac2)){
  g = ggplot(IRfrac2[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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


# Make the same plots for Age differences
totalA = rbind(xC,xN)
totalA$FDR = ifelse(totalA$padj<=0.05, "FDR<0.05", "FDR>0.05")
totalA$IR = factor(ifelse(totalA$IRratio>=0.5, ">0.5", "<0.5"))
totalA = totalA[which(totalA$padj!="NA"),]
IRAge = lapply(IRlist, function(x) totalA[which(totalA$ensID %in% x$ensID),])
IRAge2 = lapply(combined.IRlist, function(x) totalA[which(totalA$EnsID %in% x$ensID),])

pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byAge_plots.pdf", width=8, height=8)
names = names(IRlist)
for (i in 1:length(IRAge)){
  g = ggplot(IRAge[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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
pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byAge_combined.pdf", width=8, height=8)
names = names(IRAge2)
for (i in 1:length(IRAge2)){
  g = ggplot(IRAge2[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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

