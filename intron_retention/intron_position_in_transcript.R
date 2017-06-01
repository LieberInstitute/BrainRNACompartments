library(GenomicRanges)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

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


## Map introns to transcripts
exonMap$exonID = rownames(exonMap)
x = strsplit(exonMap$gencodeTx, ";", fixed=T)
names(x) = paste0(rownames(exonMap),".")
x = unlist(x, recursive =F)
x = data.frame(exonID = names(x), TxID = x)
x$exonID = gsub("\\..*","",x$exonID)
exonMap = cbind(x, exonMap[match(x$exonID, exonMap$exonID),])
lastinTx = ddply(exonMap, .(TxID), summarise, LastStart=max(Start), FirstStart=min(Start), TxEnd=max(End))
lastinTx$matchID = paste0(lastinTx$TxID, ":", lastinTx$LastStart)
exonMap$matchID = paste0(exonMap$TxID, ":", exonMap$Start)
exonMap$lastinTx = ifelse(exonMap$matchID %in% lastinTx$matchID, "Last", "NotLast")
exonMap = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
introns = lapply(introns, function(x) makeGRangesFromDataFrame(na.omit(x), 
                                                               seqnames.field="Chr",
                                                               start.field="Start",
                                                               end.field="End",
                                                               strand.field="Direction",
                                                               keep.extra.columns = T))
overlaps = lapply(introns, function(x) findOverlaps(exonMap, x))
introns = lapply(introns, function(x) as.data.frame(x))
exonMap = as.data.frame(exonMap)
intronMap = intronTx = Tx = list()
for (i in 1:length(introns)){
  ov = overlaps[[i]]
  int = introns[[i]]
  intronMap[[i]] = cbind(int[subjectHits(ov),],exonMap[queryHits(ov),])
  map = intronMap[[i]]
  map$afterExon = ifelse((map[,2]==map[,14]), "YES","NO")
  map$beforeExon = ifelse((map[,3]==map[,13]), "YES","NO")
  intronTx[[i]] = lastinTx[match(map$TxID, lastinTx$TxID),]
  tx = intronTx[[i]]
  Tx[[i]] = cbind(tx, map[match(tx$TxID, map$TxID),])
  intronMap[[i]] = map
}
names(intronMap) = names(intronTx) = names(Tx) = names(introns)
Tx = lapply(Tx, function(x) x[,c(1:4, 6:16)])


## Check intron position in the transcript
# calculate the distance of midpoint of intron from the start of the last exon divided by the length of the transcript
position = lapply(Tx, function(x) (x$TxEnd-(x$start + x$width/2))/(x$TxEnd-x$FirstStart))
lapply(position, head)

# summary statistics on position
stats = data.frame(mean = unlist(lapply(position, mean)), median = unlist(lapply(position, median)), SD = unlist(lapply(position, sd)))
stats$group = c("All", "All", "Cytosolic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytosolic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                "Adult Retention\nIn Cytosol","Prenatal Retention\nIn Cytosol","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus",
                "Nuclear Retention", "Adult Retention", "Prenatal Retention")
stats$retention = c("NA","NA","Cytosolic","Nuclear","Cytosolic","Nuclear","Adult","Prenatal","Adult","Prenatal","Nuclear","Adult", "Prenatal")


## Is there a difference by fraction
t.test(x=c(position[["Adult:Cytosol-Increased"]], position[["Prenatal:Cytosol-Increased"]]), 
       y=c(position[["Adult:Nucleus-Increased"]], position[["Prenatal:Nucleus-Increased"]]), alternative = "greater")
#t = 5.6563, df = 16.21, p-value = 1.702e-05
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.1933941       Inf
#sample estimates:
#  mean of x mean of y 
#0.7626270 0.4829873

## Is there a difference by fraction in adults
t.test(x=position[["Adult:Cytosol-Increased"]], y=position[["Adult:Nucleus-Increased"]], alternative = "greater")
#t = 4.5057, df = 11.924, p-value = 0.0003655
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.1648034       Inf
#sample estimates:
#  mean of x mean of y 
#0.7170641 0.4443119

## Is there a difference by fraction in prenatal
t.test(x=position[["Prenatal:Cytosol-Increased"]], y=position[["Prenatal:Nucleus-Increased"]], alternative = "greater")
#t = 22.141, df = 436.39, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.3219331       Inf
#sample estimates:
#  mean of x mean of y 
#0.8993157 0.5514879

## Is there a difference by age
t.test(x=c(position[["Cytosol:Adult-Increased"]], position[["Nucleus:Adult-Increased"]]), 
       y=c(position[["Cytosol:Prenatal-Increased"]], position[["Nucleus:Prenatal-Increased"]]), alternative = "greater")
#t = -0.17523, df = 1113.5, p-value = 0.5695
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.03952203         Inf
#sample estimates:
#  mean of x mean of y 
#0.4972683 0.5010705

## Is there a difference by age in cytosol
t.test(x=position[["Cytosol:Adult-Increased"]], y=position[["Cytosol:Prenatal-Increased"]], alternative = "greater")
#t = -0.15552, df = 145.44, p-value = 0.5617
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.07083251         Inf
#sample estimates:
#  mean of x mean of y 
#0.5318243 0.5379075

## Is there a difference by age in nucleus
t.test(x=position[["Nucleus:Adult-Increased"]], y=position[["Nucleus:Prenatal-Increased"]], alternative = "greater")
#t = 0.58646, df = 917.22, p-value = 0.2789
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.02718578         Inf
#sample estimates:
#  mean of x mean of y 
#0.4916112 0.4765712
