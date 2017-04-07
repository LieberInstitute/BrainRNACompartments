library(ggplot2)
library(reshape2)
library(VennDiagram)
library(data.table)
library(stringr)

# Load IRFinder Results
names = scan("./Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = c(unique(gsub( "_.*$", "", names)),"APC_Pooled","APN_Pooled","FPC_Pooled","FPN_Pooled",
                   "AP_Pooled","FP_Pooled","NP_Pooled","CP_Pooled")
path = "./Dropbox/sorted_figures/IRfinder/PolyA/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
}
names(IRres) = shortenedNames
intron = lapply(IRres, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End))
IRres = lapply(IRres, function(x) data.frame(x, name = intron[[1]]))

### Filter introns
dim(IRres[[1]]) # 233535
clean = lapply(IRres, function(x) grep("clean", x$GeneIntronDetails, fixed=T))
    # the options are "clean","anti-near","anti-over","known-exon+anti-near","known-exon","known-exon+anti-near+anti-over"
head(clean[[2]])
IRcleaned = lapply(IRres, function(x) x[clean[[1]],])
dim(IRcleaned[[1]]) # 175475

# Filter by IR ratio
ratio = data.frame(lapply(IRcleaned, function(x) x$IRratio))
head(ratio)
rownames(ratio) = IRcleaned[[1]][,"name"]
ratio.nonconst = ratio[which(rowSums(ratio[,1:12])>0),]
dim(ratio.nonconst) # 125946
IR.nonconst = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(ratio.nonconst)),])
ratio.0.1 = ratio[which(rowSums(ratio[,1:12])>0.1),]
dim(ratio.0.1) # 46046
IR.0.1 = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(ratio.0.1)),])

IR.nonconst.fetal = ratio[which(rowSums(ratio[,7:12])>0),]
dim(IR.nonconst.fetal) # 109413
IR.nonconst.fetal = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.nonconst.fetal)),])
IR.0.1.fetal = ratio[which(rowSums(ratio[,7:12])>0.1),]
dim(IR.0.1.fetal) # 25243
IR.0.1.fetal = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.0.1.fetal)),])

IR.nonconst.adult = ratio[which(rowSums(ratio[,1:6])>0),]
dim(IR.nonconst.adult) # 106159
IR.nonconst.adult = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.nonconst.adult)),])
IR.0.1.adult = ratio[which(rowSums(ratio[,1:6])>0.1),]
dim(IR.0.1.adult) # 29787
IR.0.1.adult = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.0.1.adult)),])

IR.nonconst.nuc = ratio[which(rowSums(ratio[,c(2,4,6,8,10,12)])>0),]
dim(IR.nonconst.nuc) # 118820
IR.nonconst.nuc = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.nonconst.nuc)),])
IR.0.1.nuc = ratio[which(rowSums(ratio[,c(2,4,6,8,10,12)])>0.1),]
dim(IR.0.1.nuc) # 35092
IR.0.1.nuc = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.0.1.nuc)),])

IR.nonconst.cyt = ratio[which(rowSums(ratio[,c(1,3,5,7,9,11)])>0),]
dim(IR.nonconst.cyt) # 101849
IR.nonconst.cyt = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.nonconst.cyt)),])
IR.0.1.cyt = ratio[which(rowSums(ratio[,c(1,3,5,7,9,11)])>0.1),]
dim(IR.0.1.cyt) # 20567
IR.0.1.cyt = lapply(IRcleaned, function(x) x[which(x$name %in% rownames(IR.0.1.cyt)),])

for (i in 1:length(shortenedNames)){
  write.table(IR.nonconst[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_all_nonconst.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.0.1[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_all_0.1.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.fetal[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_fetal_nonconst.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.0.1.fetal[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_fetal_0.1.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.adult[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_adult_nonconst.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.0.1.adult[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_adult_0.1.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.nuc[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_nuc_nonconst.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.0.1.nuc[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_nuc_0.1.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.cyt[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_cyt_nonconst.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.0.1.cyt[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_cyt_0.1.txt"),
              quote =F, row.names = F, sep ="\t")
}

# Filtering by coverage is a bit too severe
SplicesExact = data.frame(lapply(IR.nonconst, function(x) x$SplicesExact), row.names = rownames(ratio.nonconst))
SplicesExact = SplicesExact[which(rowMeans(SplicesExact)>=5),]
dim(SplicesExact) # 68296
ExonToIntronReadsLeft = data.frame(lapply(IR.nonconst, function(x) x$ExonToIntronReadsLeft),
                                   row.names = rownames(ratio.nonconst))
ExonToIntronReadsLeft = ExonToIntronReadsLeft[which(rowMeans(ExonToIntronReadsLeft)>=5),]
dim(ExonToIntronReadsLeft) # 7072
ExonToIntronReadsRight = data.frame(lapply(IR.nonconst, function(x) x$ExonToIntronReadsRight),
                                    row.names = rownames(ratio.nonconst))
ExonToIntronReadsRight = ExonToIntronReadsRight[which(rowMeans(ExonToIntronReadsRight)>=5),]
dim(ExonToIntronReadsRight) # 2189

# Further Filter by warnings
    # Warning Levels: - LowCover LowSplicing MinorIsoform NonUniformIntronCover
nowarn = lapply(IR.nonconst[1:12], function(x) x[which(x$Warnings=="-"),])
nowarnnames = as.character(unlist(lapply(nowarn, function(x) x$name)))
uniquenames = unique(nowarnnames)
length(uniquenames)

nowarnA = lapply(IR.nonconst.adult[1:6], function(x) x[which(x$Warnings=="-"),])
nowarnnamesA = as.character(unlist(lapply(nowarnA, function(x) x$name)))
uniquenamesA = unique(nowarnnamesA)
length(uniquenamesA)
nowarnN = lapply(IR.nonconst.nuc[c(2,4,6,8,10,12)], function(x) x[which(x$Warnings=="-"),])
nowarnnamesN = as.character(unlist(lapply(nowarnN, function(x) x$name)))
uniquenamesN = unique(nowarnnamesN)
length(uniquenamesN)
nowarnC = lapply(IR.nonconst.cyt[c(1,3,5,7,9,11)], function(x) x[which(x$Warnings=="-"),])
nowarnnamesC = as.character(unlist(lapply(nowarnC, function(x) x$name)))
uniquenamesC = unique(nowarnnamesC)
length(uniquenamesC)
nowarnF = lapply(IR.nonconst.fetal[7:12], function(x) x[which(x$Warnings=="-"),])
nowarnnamesF = as.character(unlist(lapply(nowarnF, function(x) x$name)))
uniquenamesF = unique(nowarnnamesF)
length(uniquenamesF)

countsA = countsN = countsC = countsF = counts_all = c()
for (i in 1:length(uniquenamesC)){
  countsC[i] = sum(str_count(nowarnnamesC, uniquenamesC[i]))}
freqC = data.frame(names = uniquenamesC, freq = countsC)
save(freqC, file="./Dropbox/cyt_freq.rda")
for (i in 1:length(uniquenamesF)){
  countsF[i] = sum(str_count(nowarnnamesF, uniquenamesF[i]))}
freqF = data.frame(names = uniquenamesF, freq = countsF)
save(freqF, file="./Dropbox/fetal_freq.rda")
for (i in 1:length(uniquenames)){
  counts_all[i] = sum(str_count(nowarnnames, uniquenames[i]))}
freq = data.frame(names = uniquenames, freq = counts_all)
save(freq, file="./Dropbox/all_freq.rda")
for (i in 1:length(uniquenamesA)){
  countsA[i] = sum(str_count(nowarnnamesA, uniquenamesA[i]))}
freqA = data.frame(names = uniquenamesA, freq = countsA)
save(freqA, file="./Dropbox/adult_freq.rda")
for (i in 1:length(uniquenamesN)){
  countsN[i] = sum(str_count(nowarnnamesN, uniquenamesN[i]))}
freqN = data.frame(names = uniquenamesN, freq = countsN)
save(freqN, file="./Dropbox/nuc_freq.rda")

# Explore warning results
load("./Dropbox/sorted_figures/new/cyt_freq.rda")
load("./Dropbox/sorted_figures/new/nuc_freq.rda")
load("./Dropbox/sorted_figures/new/adult_freq.rda")
load("./Dropbox/sorted_figures/new/fetal_freq.rda")
load("./Dropbox/sorted_figures/new/all_freq.rda")
dim(IR.nonconst[[1]])
dim(freq[which(freq$freq>=12),]) # 3367
dim(freqF[which(freqF$freq>=6),]) # 18851
dim(freqA[which(freqA$freq>=6),]) # 5228
dim(freqC[which(freqC$freq>=6),]) # 4966
dim(freqN[which(freqN$freq>=6),]) # 9076

dim(freq[which(freq$freq>=8),]) # 22738
dim(freqF[which(freqF$freq>=4),]) # 34531
dim(freqA[which(freqA$freq>=4),]) # 19829
dim(freqC[which(freqC$freq>=4),]) # 20846
dim(freqN[which(freqN$freq>=4),]) # 24425

IR.nonconst.nowarn = lapply(IRcleaned, function(x) x[which(x$name %in% freq[which(freq$freq>=12),"names"]),])
IR.nonconst.nowarnF = lapply(IRcleaned, function(x) x[which(x$name %in% freqF[which(freqF$freq>=6),"names"]),])
IR.nonconst.nowarnA = lapply(IRcleaned, function(x) x[which(x$name %in% freqA[which(freqA$freq>=6),"names"]),])
IR.nonconst.nowarnC = lapply(IRcleaned, function(x) x[which(x$name %in% freqC[which(freqC$freq>=6),"names"]),])
IR.nonconst.nowarnN = lapply(IRcleaned, function(x) x[which(x$name %in% freqN[which(freqN$freq>=6),"names"]),])

IR.nonconst.66warn = lapply(IRcleaned, function(x) x[which(x$name %in% freq[which(freq$freq>=8),"names"]),])
IR.nonconst.66warnF = lapply(IRcleaned, function(x) x[which(x$name %in% freqF[which(freqF$freq>=4),"names"]),])
IR.nonconst.66warnA = lapply(IRcleaned, function(x) x[which(x$name %in% freqA[which(freqA$freq>=4),"names"]),])
IR.nonconst.66warnC = lapply(IRcleaned, function(x) x[which(x$name %in% freqC[which(freqC$freq>=4),"names"]),])
IR.nonconst.66warnN = lapply(IRcleaned, function(x) x[which(x$name %in% freqN[which(freqN$freq>=4),"names"]),])

for (i in 1:length(shortenedNames)){
  write.table(IR.nonconst.nowarn[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_all_nonconst.nowarn.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.66warn[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_all_nonconst.66warn.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.nowarnA[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_adult_nonconst.nowarnA.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.66warnA[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_adult_nonconst.66warnA.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.nowarnF[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_fetal_nonconst.nowarnF.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.66warnF[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_fetal_nonconst.66warnF.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.nowarnC[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_cyt_nonconst.nowarnC.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.66warnC[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_cyt_nonconst.66warnC.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.nowarnN[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_nuc_nonconst.nowarnN.txt"),
              quote =F, row.names = F, sep ="\t")
  write.table(IR.nonconst.66warnN[[i]][,1:21], file=paste0(path,shortenedNames[i],"/IRFinder-IR-nondir_nuc_nonconst.66warnN.txt"),
              quote =F, row.names = F, sep ="\t")
}

#cat IRFinder-IR-nondir.txt | grep -P “0\.[1-9]\d+\s+-$” | awk 'BEGIN {OFS="\t"}{if ($9>5 && $19>5 && $20 >0.1)print $0}'
#to get IR events with at least 5 reads depth in the intron and across spliced exons and at least 0.1 IR ratio. 
#In fact for my first pass I will often filter very aggressively and look at these manually in an IGV viewer:
#  cat IRFinder-IR-nondir.txt | grep -P “0\.[1-9]\d+\s+-$” | awk ‘BEGIN {OFS=“\t"}{if ($9>10 && $19>10 && $20 >0.3)print $0}' 

#1) convert the 9th column (intron depth) from IRFinder-IR-nondir.txt to integers and treat them as “raw counts” when you build your DESeqDataSet object.
#2) calculate the sum of the 9th and 19th column (intron depth and spliced reads) from IRFinder-IR-nondir.txt for each sample. Integrate them as the library sizes for your samples into the DESeqDataSet object.
#3) I recommend to apply the filters which William suggested. That is going to remove introns with low sequencing/statistical confidence out of the analysis
#4) You can run DESeq2 as normal.
