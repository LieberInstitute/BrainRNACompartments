library(data.table)
library(GenomicRanges)
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

## Do DEG by fraction contain introns with higher IR ratios than non-DEG by fraction?
# Get the IR ratio for Fraction genes
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down))
FracList = Map(cbind, FracList, lapply(FracList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
SigFracList = Map(cbind, SigFracList, Sign = Sign)
SigFracList = lapply(SigFracList, function(x) split(x, x$Sign))
SigList = unlist(SigFracList, recursive = F) 
lapply(SigList, head)
allgenes = as.character(rownames(FracList[[1]]))
sigFracIRratio = list(list(), list(), list(),list(),list(), list(),list(), list(),list(), list()) 
for (i in 1:length(SigList)){
  for (j in 1:length(IRfiltered2)){
    filt = IRfiltered2[[j]][,c(22,23,20)]
    sigFracIRratio[[i]][[j]] = filt[which(as.character(SigList[[i]][,"ensemblID"]) %in% as.character(filt$genes)),]
  }}
names(sigFracIRratio) = names(SigList)
for (i in 1:length(sigFracIRratio)){names(sigFracIRratio[[i]]) = names(IRfiltered2)}
lapply(sigFracIRratio, elementNROWS)
sigFracIRratio = lapply(sigFracIRratio, function(x) lapply(x, function(y) data.table(y, key="genes")))
sigFracIRratio = lapply(sigFracIRratio, function(x) lapply(x, function(y) data.frame(y[, list(IRratio=max(IRratio)), by="genes"]))) # limit to intron with highest IR ratio per gene per sample
sigFracIRratio = lapply(sigFracIRratio, function(x) do.call(rbind,x))
elementNROWS(sigFracIRratio)

# Compare the intron with the greatest retention per gene per sample in different sets of DEGs
t.test(sigFracIRratio[["Apres.DownNuc"]][,"IRratio"], sigFracIRratio[["Apres.UpNuc"]][,"IRratio"], alternative = "less")
#data:  sigFracIRratio[["Apres.DownNuc"]][, "IRratio"] and sigFracIRratio[["Apres.UpNuc"]][, "IRratio"]
#t = -0.77008, df = 2768, p-value = 0.2207
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.00141993
#sample estimates:
#  mean of x  mean of y 
#0.01539956 0.01664877 
t.test(sigFracIRratio[["Arres.DownNuc"]][,"IRratio"], sigFracIRratio[["Arres.UpNuc"]][,"IRratio"], alternative = "two.sided")
#data:  sigFracIRratio[["Arres.DownNuc"]][, "IRratio"] and sigFracIRratio[["Arres.UpNuc"]][, "IRratio"]
#t = 3.6311, df = 1971.5, p-value = 0.0002894
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.002722987 0.009118719
#sample estimates:
#  mean of x  mean of y 
#0.01730208 0.01138122  

# Get the DEG Fraction p-value and LFC sign for differentially retained introns
nonconst = IRclean[["nonconst"]]
elementNROWS(lapply(nonconst, function(x) x[which(x$p.diff<=0.05),]))
nonconst = Map(cbind, nonconst, 
               AP.sig = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"padj"]),
               AP.LFC = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"log2FoldChange"]),
               FP.sig = lapply(nonconst, function(x) FracList[["Fpres"]][match(x$ensID, FracList[["Fpres"]][,"ensemblID"]),"padj"]),
               FP.LFC = lapply(nonconst, function(x) FracList[["Fpres"]][match(x$ensID, FracList[["Fpres"]][,"ensemblID"]),"log2FoldChange"]),
               AR.sig = lapply(nonconst, function(x) FracList[["Arres"]][match(x$ensID, FracList[["Arres"]][,"ensemblID"]),"padj"]),
               AR.LFC = lapply(nonconst, function(x) FracList[["Arres"]][match(x$ensID, FracList[["Arres"]][,"ensemblID"]),"log2FoldChange"]),
               FR.sig = lapply(nonconst, function(x) FracList[["Frres"]][match(x$ensID, FracList[["Frres"]][,"ensemblID"]),"padj"]),
               FR.LFC = lapply(nonconst, function(x) FracList[["Frres"]][match(x$ensID, FracList[["Frres"]][,"ensemblID"]),"log2FoldChange"]),
               FP.down.sig = lapply(nonconst, function(x) FracList[["Fpres.down"]][match(x$ensID, FracList[["Fpres.down"]][,"ensemblID"]),"padj"]),
               FP.down.LFC = lapply(nonconst, function(x) FracList[["Fpres.down"]][match(x$ensID, FracList[["Fpres.down"]][,"ensemblID"]),"log2FoldChange"]))
elementNROWS(nonconst)
lapply(nonconst, head)

# Are significantly dIR genes more likely to be significantly DEG by fraction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]>=0.05),])
fisher.test(data.frame(c(164,148), c(286,489))) #adult
#data:  data.frame(c(164, 148), c(286, 489))
#p-value = 2.539e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.439242 2.493470
#sample estimates:
#  odds ratio 
#1.893514
fisher.test(data.frame(c(42,167), c(47,717))) #prenatal
#data:  data.frame(c(42, 167), c(47, 717))
#p-value = 1.507e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.379715 6.151949
#sample estimates:
#  odds ratio 
#3.829787

# Of dIR Fraction genes, does the LFC go in the same direction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]<0),])
fisher.test(data.frame(c(143,69), c(447,337))) #prenatal
#data:  data.frame(c(143, 69), c(447, 337))
#p-value = 0.007306
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.122700 2.186824
#sample estimates:
#  odds ratio 
#1.561769
fisher.test(data.frame(c(151,162), c(439,341))) #adult
#data:  data.frame(c(151, 162), c(439, 341))
#p-value = 0.01873
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5517484 0.9501778
#sample estimates:
#  odds ratio 
#0.7242496

# Of dIR Fraction genes, does the direction of dIR more in nucleus?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
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


## Do DEG by age contain introns with higher IR ratios than non-DEG by age?
# Get the IR ratio for developmental genes
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres),
               Crres = data.frame(Crres),Nrres = data.frame(Nrres), 
               Cpres.down = data.frame(Cpres.down))
AgeList = Map(cbind, AgeList, lapply(AgeList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigAgeList)
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpPrenatal", "DownPrenatal"))
SigAgeList = Map(cbind, SigAgeList, Sign = Sign)
SigAgeList = lapply(SigAgeList, function(x) split(x, x$Sign))
SigList = unlist(SigAgeList, recursive = F) 
lapply(SigList, head)
allgenes = as.character(rownames(AgeList[[1]]))
sigAgeIRratio = list(list(), list(), list(),list(),list(), list(),list(), list(),list(), list()) 
for (i in 1:length(SigList)){
  for (j in 1:length(IRfiltered2)){
    filt = IRfiltered2[[j]][,c(22,23,20)]
    sigAgeIRratio[[i]][[j]] = filt[which(as.character(SigList[[i]][,"ensemblID"]) %in% as.character(filt$genes)),]
  }}
names(sigAgeIRratio) = names(SigList)
for (i in 1:length(sigAgeIRratio)){names(sigAgeIRratio[[i]]) = names(IRfiltered2)}
lapply(sigAgeIRratio, elementNROWS)
sigAgeIRratio = lapply(sigAgeIRratio, function(x) lapply(x, function(y) data.table(y, key="genes")))
sigAgeIRratio = lapply(sigAgeIRratio, function(x) lapply(x, function(y) data.frame(y[, list(IRratio=max(IRratio)), by="genes"]))) # limit to intron with highest IR ratio per gene per sample
sigAgeIRratio = lapply(sigAgeIRratio, function(x) do.call(rbind,x))
elementNROWS(sigAgeIRratio)

# Compare the intron with the greatest retention per gene per sample in different sets of DEGs
t.test(sigAgeIRratio[["Cpres.DownPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.UpPrenatal"]][,"IRratio"], alternative = "less")
#data:  sigAgeIRratio[["Cpres.DownPrenatal"]][, "IRratio"] and sigAgeIRratio[["Cpres.UpPrenatal"]][, "IRratio"]
#t = -4.6249, df = 11983, p-value = 1.893e-06
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.001788791
#sample estimates:
#  mean of x  mean of y 
#0.01217757 0.01495380 
t.test(sigAgeIRratio[["Npres.DownPrenatal"]][,"IRratio"], sigAgeIRratio[["Npres.UpPrenatal"]][,"IRratio"], alternative = "less")
#data:  sigAgeIRratio[["Npres.DownPrenatal"]][, "IRratio"] and sigAgeIRratio[["Npres.UpPrenatal"]][, "IRratio"]
#t = -5.3515, df = 10065, p-value = 4.457e-08
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.002540159
#sample estimates:
#  mean of x  mean of y 
#0.01248776 0.01615528
t.test(sigAgeIRratio[["Cpres.down.DownPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.down.UpPrenatal"]][,"IRratio"], alternative = "less")
#data:  sigAgeIRratio[["Cpres.down.DownPrenatal"]][, "IRratio"] and sigAgeIRratio[["Cpres.down.UpPrenatal"]][, "IRratio"]
#t = -5.3569, df = 10730, p-value = 4.319e-08
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.002534245
#sample estimates:
#  mean of x  mean of y 
#0.01280146 0.01645878

# Get the DEG Age p-value and LFC sign for differentially retained introns
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
               CP.down.LFC = lapply(nonconst, function(x) AgeList[["Cpres.down"]][match(x$ensID, AgeList[["Cpres.down"]][,"ensemblID"]),"log2FoldChange"]))
elementNROWS(nonconst)
lapply(nonconst, head)

# Are significantly dIR genes more likely to be significantly DEG by age?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]>=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]>=0.05),])
fisher.test(data.frame(c(118,50), c(302,207))) #cytosol
#data:  data.frame(c(164, 148), c(286, 489))
#data:  data.frame(c(118, 50), c(302, 207))
#p-value = 0.0132
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.096525 2.406682
#sample estimates:
#  odds ratio 
#1.61649
fisher.test(data.frame(c(235,112), c(513,339))) #nucleus
#data:  data.frame(c(235, 112), c(513, 339))
#p-value = 0.01509
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.057376 1.823272
#sample estimates:
#  odds ratio 
#1.386161

# Of dIR Age genes, does the LFC go in the same direction?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]>0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]<0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]<05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]>0),])
fisher.test(data.frame(c(109,59), c(248,261))) #cytosol
#data:  data.frame(c(109, 59), c(248, 261))
#p-value = 0.0003479
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.335597 2.844027
#sample estimates:
#  odds ratio 
#1.942424
fisher.test(data.frame(c(200,148), c(853,392))) #nucleus
#data:  data.frame(c(200, 148), c(853, 392))
#p-value = 0.0001536
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4833510 0.7992836
#sample estimates:
#  odds ratio 
#0.6211912

# Of dIR Age genes, does the direction of dIR more in adult?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
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

