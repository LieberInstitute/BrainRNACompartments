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


# Explore the results
lapply(IRclean, function(x) elementNROWS(x))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#792               649               461               740                55                85 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#792               649               461               740                55                85 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#566               413               274               517                52                76 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#145               253                74               181                 3                 2 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#465               471               270               430                19                31
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#783               629               391               390                55                59 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#783               629               391               390                55                59 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#558               399               220               231                52                50 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#144               246                62               113                 3                 2 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#458               462               227               239                19                20 
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$p.diff<=0.05),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#162                93                80               179                32                37 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#162                93                80               179                32                37 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#143                84                69               169                32                35 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#45                21                12                32                 1                 1 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#117                66                42               105                10                17
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$p.diff<=0.05 & y$Sign=="MoreIRInNuc.Fetal"),])))
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#160                92                63                77                32                25 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#160                92                63                77                32                25 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#141                83                53                72                32                23 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#44                21                 9                17                 1                 1 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#115                66                30                45                10                11
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 5                 1                12                 1                 3 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 5                 1                12                 1                 3 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 5                 1                12                 1                 3 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#0                 0                 0                 1                 0                 0 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2                 1                 0                 7                 0                 1 
lapply(IRclean, function(x) elementNROWS(lapply(lapply(x, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),]), 
                                                function(z) z[which(z$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
# Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 5                 1                 4                 1                 3 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 5                 1                 4                 1                 3 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#4                 5                 1                 4                 1                 3 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#0                 0                 0                 1                 0                 0 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2                 1                 0                 2                 0                 1 


#moving forward with the nonconstitutively spliced results because they capture more introns, 
#and coverage issues are filtered after differential retention calculation
nonconst = IRclean[["nonconst"]]
elementNROWS(nonconst)
elementNROWS(lapply(nonconst, function(y) y[which(y$Sign=="MoreIRInNuc.Fetal"),]))
elementNROWS(lapply(nonconst, function(x) x[which(x$p.diff<=0.05),]))
elementNROWS(lapply(nonconst, function(y) y[which(y$p.diff<=0.05 & y$Sign=="MoreIRInNuc.Fetal"),]))
elementNROWS(lapply(nonconst, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),]))
elementNROWS(lapply(nonconst, function(y) y[which((y$A.IRratio>=0.5 | y$B.IRratio>=0.5) & y$Sign=="MoreIRInNuc.Fetal"),]))

# Comparison of significantly vs nonsignificantly retained introns by zone/age
fisher.test(data.frame(c(160,2),c(623,7))) 
# adult zone
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1689921 8.9521084
#sample estimates:
#  odds ratio 
#0.8989833
fisher.test(data.frame(c(92,1),c(537,19)))
# fetal zone
#p-value = 0.3372
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5043948 136.6375340
#sample estimates:
#  odds ratio 
#3.251195
fisher.test(data.frame(c(63,17),c(391-63,461-391-17)))
# cytosol age
#p-value = 0.1214
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3166124 1.1787385
#sample estimates:
#  odds ratio 
#0.5995641
fisher.test(data.frame(c(77,102),c(390-77,740-390-102)))
# nucleus age
#p-value = 0.003416
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4194177 0.8517435
#sample estimates:
#  odds ratio 
#0.5985671

# Get the DEG pval and LFC sign by Fraction for differentially retained introns
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down), Ares = data.frame(Ares), Fres.down = data.frame(Fres.down))
FracList = Map(cbind, FracList, lapply(FracList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
nonconst = Map(cbind, nonconst, 
               AP.sig = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"padj"]),
               AP.LFC = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"log2FoldChange"]),
               FP.sig = lapply(nonconst, function(x) FracList[["Fpres.down"]][match(x$ensID, FracList[["Fpres.down"]][,"ensemblID"]),"padj"]),
               FP.LFC = lapply(nonconst, function(x) FracList[["Fpres.down"]][match(x$ensID, FracList[["Fpres.down"]][,"ensemblID"]),"log2FoldChange"]),
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
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.sig"]>=0.05),])
fisher.test(data.frame(c(85,77), c(230,394))) #adult polya
#p-value = 0.0004274
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.313810 2.721389
#sample estimates:
#  odds ratio 
#1.889471
fisher.test(data.frame(c(22,68), c(36,504))) #prenatal polya
#p-value = 1.949e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.380947 8.425692
#sample estimates:
#  odds ratio 
#4.513507
fisher.test(data.frame(c(97,65), c(291,333))) #adult both libraries
#p-value = 0.002739
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.184965 2.469225
#sample estimates:
#  odds ratio 
#1.706497
fisher.test(data.frame(c(30,62), c(111,442))) #prenatal both libraries
#p-value = 0.009472
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.143450 3.189854
#sample estimates:
#  odds ratio 
#1.9246

# Are genes with significantly differentially retained introns more likely to have LFC in one direction by fraction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"A.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"F.down.LFC"]<0),])
fisher.test(data.frame(c(85,77), c(371,256))) #adult polya
#p-value = 0.1299
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5305547 1.0952479
#sample estimates:
#  odds ratio 
#0.7619954
fisher.test(data.frame(c(70,22), c(334,219))) #prenatal polya
#p-value = 0.003621
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.231986 3.644240
#sample estimates:
#  odds ratio 
#2.084049
fisher.test(data.frame(c(84,78), c(378,249))) #adult in both libraries
#p-value = 0.06017
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4939812 1.0201661
#sample estimates:
#  odds ratio 
#0.7097328
fisher.test(data.frame(c(66,26), c(367,187))) #prenatal in both libraries
#p-value = 0.3389
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7797678 2.1955444
#sample estimates:
#  odds ratio 
#1.292938

# Are significantly differentially retained introns by fraction more likely to have a higher IR ratio in nuclear RNA?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #160
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #2
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #623
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #7
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #92
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #1
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #537
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #19
fisher.test(data.frame(c(160,2), c(623,7))) #adult
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1689921 8.9521084
#sample estimates:
#  odds ratio 
#0.8989833
fisher.test(data.frame(c(92,1), c(537,19))) #prenatal
#p-value = 0.3372
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5043948 136.6375340
#sample estimates:
#  odds ratio 
#3.251195

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
fisher.test(data.frame(c(56,23), c(217,159))) #cytosol polya
#p-value = 0.03199
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.028505 3.169073
#sample estimates:
#  odds ratio 
#1.78183
fisher.test(data.frame(c(116,60), c(326,229))) #nucleus polya
#p-value = 0.0935
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9402847 1.9730512
#sample estimates:
#  odds ratio 
#1.357514
fisher.test(data.frame(c(63,16), c(261,115))) #cytosol in both libraries
#p-value = 0.07552
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9397578 3.3574674
#sample estimates:
#  odds ratio 
#1.732976
fisher.test(data.frame(c(121,55), c(356,199))) #nucleus in both libraries
#p-value = 0.2768
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.844410 1.804591
#sample estimates:
#  odds ratio 
#1.229432

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
fisher.test(data.frame(c(43,36), c(162,214))) #cytosol polya
#p-value = 0.08126
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9413047 2.6522660
#sample estimates:
#  odds ratio 
#1.576209
fisher.test(data.frame(c(89,88), c(555,296))) #nucleus polya
#p-value = 0.0002336
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3839657 0.7584748
#sample estimates:
#  odds ratio 
#0.539719
fisher.test(data.frame(c(43,36), c(175,201))) #cytosol in both libraries
#p-value = 0.2169
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8192296 2.3048746
#sample estimates:
#  odds ratio 
#1.370935
fisher.test(data.frame(c(92,85), c(555,268))) #nucleus in both libraries
#p-value = 0.0001298
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3712786 0.7369890
#sample estimates:
#  odds ratio 
#0.5229967

# Are significantly differentially retained introns more likely to have a higher IR ratio in prenatal samples?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #63
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #17
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #328
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #53
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #77
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #102
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),]) #313
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),]) #248
fisher.test(data.frame(c(63,17), c(328,53))) #cytosol
#p-value = 0.1214
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3166124 1.1787385
#sample estimates:
#  odds ratio 
#0.5995641
fisher.test(data.frame(c(77,102), c(313,248))) #nucleus
#p-value = 0.003416
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4194177 0.8517435
#sample estimates:
#  odds ratio 
#0.5985671

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
#t = 2.2342, df = 1.0548, p-value = 0.2575
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3.235357  4.839308
#sample estimates:
#  mean of x  mean of y 
#0.82133680 0.01936121 
t.test(IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"], 
       IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"], alternative = "two.sided")
#data:  Prenatal:Cytosol-Enriched and Prenatal:Nuclear-Enriched
#not enough 'x' observations
t.test(c(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
         IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"]), 
       c(IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
         IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"]), alternative = "two.sided")
#data:  Combined  Cytosol-Enriched and Nuclear-Enriched
#t = 1.0935, df = 2.0468, p-value = 0.3861
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.167567  1.986998
#sample estimates:
#  mean of x mean of y 
#0.5099668 0.1002516
t.test(IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"], 
       IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"], alternative = "two.sided")
#data: Cytosol:Adult-Enriched and Cytosol:Prenatal-Enriched
#t = 5.2669, df = 24.108, p-value = 2.089e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.9850287 2.2539872
#sample estimates:
#  mean of x  mean of y 
#0.9928103 -0.6266976
t.test(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"], 
       IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"], alternative = "two.sided")
#data:  Nucleus:Adult-Enriched and Nucleus:Prenatal-Enriched
#t = 6.9574, df = 135.18, p-value = 1.356e-10
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.8695872 1.5602946
#sample estimates:
#  mean of x  mean of y 
#0.3771470 -0.8377939 
t.test(c(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
         IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"]), 
       c(IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
         IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"]), alternative = "two.sided")
#data:  Combined Adult-Enriched and Combined Prenatal-Enriched
#t = 8.5545, df = 240.5, p-value = 1.405e-15
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.9272796 1.4820912
#sample estimates:
#  mean of x  mean of y 
#0.4643206 -0.7403648  

## Are the Fraction pvalues more significant in genes containing a retained intron preferentially?

t.test(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"padj"], 
       IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"padj"], alternative = "two.sided")
#data:  Adult:Cytosol-Enriched and Adult:Nuclear-Enriched
#t = -8.0333, df = 142.5, p-value = 3.254e-13
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2483328 -0.1502528
#sample estimates:
#  mean of x   mean of y 
#0.001126789 0.200419598  
t.test(IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"padj"], 
       IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"padj"], alternative = "two.sided")
#data:  Prenatal:Cytosol-Enriched and Prenatal:Nuclear-Enriched
#not enough 'x' observations
t.test(c(IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"padj"],
         IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"padj"]), 
       c(IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"padj"],
         IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"padj"]), alternative = "two.sided")
#data:  Combined  Cytosol-Enriched and Nuclear-Enriched
#t = 0.12266, df = 2.0225, p-value = 0.9135
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.12296  1.18959
#sample estimates:
#  mean of x mean of y 
#0.2719766 0.2386617
t.test(IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"padj"], 
       IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"padj"], alternative = "two.sided")
#data: Cytosol:Adult-Enriched and Cytosol:Prenatal-Enriched
#t = 0.53028, df = 20.638, p-value = 0.6016
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1188533  0.2000960
#sample estimates:
#  mean of x  mean of y 
#0.13665037 0.09602901 
t.test(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"padj"], 
       IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"padj"], alternative = "two.sided")
#data:  Nucleus:Adult-Enriched and Nucleus:Prenatal-Enriched
#t = 1.7406, df = 161.27, p-value = 0.08367
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.01103201  0.17498677
#sample estimates:
#  mean of x mean of y 
#0.2127748 0.1307974
t.test(c(IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
         IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"padj"]), 
       c(IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
         IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"padj"]), alternative = "two.sided")
#data:  Combined Adult-Enriched and Combined Prenatal-Enriched
#t = 2.3053, df = 212.35, p-value = 0.02212
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.01264375 0.16184756
#sample estimates:
#  mean of x mean of y 
#0.2019961 0.1147505

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

fisher.test(data.frame(c(1+3+4+42,7+2+46+142),c(5+2+33+163,873+233+630))) 
#data: dIR_FractionbyAge_combined  
#p-value = 2.972e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.506956 3.083273
#sample estimates:
#  odds ratio 
#2.169531
fisher.test(data.frame(c(15,7+58),c(14+133,581+42+325)))
#data: dIR_FractionbyAge_adult_cytosol
#p-value = 0.1805
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7669728 2.7235930
#sample estimates:
#  odds ratio 
#1.487655
fisher.test(data.frame(c(1,6+73),c(3+89,22+528+356)))
#data: dIR_FractionbyAge_prenatal_cytosol
#p-value = 0.01107
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.003087761 0.734638044
#sample estimates:
#  odds ratio 
#0.1247735
fisher.test(data.frame(c(21,29+129),c(15+126,97+504+449)))
#data: dIR_FractionbyAge_adult_nucleus
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5765156 1.6293415
#sample estimates:
#  odds ratio 
#0.9897722
fisher.test(data.frame(c(16,17+146),c(9+68,93+446+459)))
#data: dIR_FractionbyAge_prenatal_nucleus
#p-value = 0.4399
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6753004 2.2670649
#sample estimates:
#  odds ratio 
#1.271986
fisher.test(data.frame(c(2,5+86),c(6+154,21+529+604)))
#data: dIR_Fraction_adult_fetal
#p-value = 0.001206
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0187539 0.6002508
#sample estimates:
#  odds ratio 
#0.1586279
fisher.test(data.frame(c(12,10+157),c(4+64,30+341+527)))
#data: dIR_Age_cytosol_nucleus
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4572037 1.8166853
#sample estimates:
#  odds ratio 
#0.94895

# What about only the introns reported in both comparisons?
fisher.test(data.frame(c(1+3+4+42,2+46),c(2+33,233))) 
#data: dIR_FractionbyAge_combined  
#p-value = 4.229e-13
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.933238 12.221562
#sample estimates:
#  odds ratio 
#6.885668
fisher.test(data.frame(c(15,7),c(14,42)))
#data: dIR_FractionbyAge_adult_cytosol
#p-value = 0.0006281
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.938127 22.224856
#sample estimates:
#  odds ratio 
#6.249445
fisher.test(data.frame(c(1,6),c(3,22)))
#data: dIR_FractionbyAge_prenatal_cytosol
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.01998962 18.71582978
#sample estimates:
#  odds ratio 
#1.214275
fisher.test(data.frame(c(21,29),c(15,97)))
#data: dIR_FractionbyAge_adult_nucleus
#p-value = 0.0001511
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.993149 11.047539
#sample estimates:
#  odds ratio 
#4.629837
fisher.test(data.frame(c(16,17),c(9,93)))
#data: dIR_FractionbyAge_prenatal_nucleus
#p-value = 2.7e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.340701 28.882652
#sample estimates:
#  odds ratio 
#9.492299
fisher.test(data.frame(c(2,5),c(6,21)))
#data: dIR_Fraction_adult_fetal
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1064076 11.6437411
#sample estimates:
#  odds ratio 
#1.385473
fisher.test(data.frame(c(12,10),c(4,30)))
#data: dIR_Age_cytosol_nucleus
#p-value = 0.0008489
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.030847 45.265290
#sample estimates:
#  odds ratio 
#8.58202



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


fisher.test(data.frame(c(30+9+11+9,44+22+19+84),c(35+18+22+87,606+271+421))) 
#data: dIR_gene_FractionbyAge_combined
#p-value = 1.842e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.955509 3.962748
#sample estimates:
#  odds ratio 
#2.795009
fisher.test(data.frame(c(1+1+2+12,13+4+4+40),c(6+17+89+17,454+75+243)))
#data: dIR_gene_FractionbyAge_adult_cytosol
#p-value = 0.1328
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8180197 2.8602972
#sample estimates:
#  odds ratio 
#1.568876
fisher.test(data.frame(c(2,8+5+5+57),c(3+3+9+69,42+439+293)))
#data: dIR_gene_FractionbyAge_prenatal_cytosol
#p-value = 0.03734
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.02873558 0.95037313
#sample estimates:
#  odds ratio 
#0.2459412
fisher.test(data.frame(c(20+1+4+3,35+7+14+84),c(4+18+15+80,122+382+348)))
#data: dIR_gene_FractionbyAge_adult_nucleus
#p-value = 0.104
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8931127 2.3128402
#sample estimates:
#  odds ratio 
#1.455879
fisher.test(data.frame(c(18+2+2+1,20+7+15+103),c(5+11+4+43,113+354+363)))
#data: dIR_gene_FractionbyAge_prenatal_nucleus
#p-value = 0.007971
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.196070 3.542889
#sample estimates:
#  odds ratio 
#2.088009
fisher.test(data.frame(c(7+3,12+4+8+52),c(10+6+17+102,65+465+413)))
#data: dIR_gene_Fraction_adult_fetal
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4134921 1.8414453
#sample estimates:
#  odds ratio 
#0.919139
fisher.test(data.frame(c(12+2+2+1,20+3+19+109),c(5+2+5+48,58+260+427)))
#data: dIR_gene_Age_cytosol_nucleus
#p-value = 0.2703
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7428577 2.5088940
#sample estimates:
#  odds ratio 
#1.397376


# What about only the genes reported in both comparisons?
fisher.test(data.frame(c(30+9+11+9,44+22),c(35+18,271))) 
#data: dIR_gene_FractionbyAge_combined
#p-value = 7.768e-11
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.813228 7.410751
#sample estimates:
#  odds ratio 
#4.552273
fisher.test(data.frame(c(1+1+2+12,13+4),c(6+17,75)))
#data: dIR_gene_FractionbyAge_adult_cytosol
#p-value = 0.008738
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.228627 7.583875
#sample estimates:
#  odds ratio 
#3.039395
fisher.test(data.frame(c(2,8+5),c(3+3,42)))
#data: dIR_gene_FractionbyAge_prenatal_cytosol
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.09517173 7.03631018
#sample estimates:
#  odds ratio 
#1.075645
fisher.test(data.frame(c(20+1+4+3,35+7),c(4+18,122)))
#data: dIR_gene_FractionbyAge_adult_nucleus
#p-value = 0.0001195
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.812570 7.542576
#sample estimates:
#  odds ratio 
#3.67088
fisher.test(data.frame(c(18+2+2+1,20+7),c(5+11,113)))
#data: dIR_gene_FractionbyAge_prenatal_nucleus
#p-value = 3.752e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.615232 13.878318
#sample estimates:
#  odds ratio 
#5.940231
fisher.test(data.frame(c(7+3,12+4),c(10+6,65)))
#data: dIR_gene_Fraction_adult_fetal
#p-value = 0.06736
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.852056 7.293437
#sample estimates:
#  odds ratio 
#2.514344
fisher.test(data.frame(c(12+2+2+1,20+3),c(5+2,58)))
#data: dIR_gene_Age_cytosol_nucleus
#p-value = 0.0002721
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.043956 19.550674
#sample estimates:
#  odds ratio 
#6.00427