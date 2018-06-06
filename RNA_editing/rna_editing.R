### ON SERVER 4 ###
# Create R object of RNA editing pipeline output
ids = scan("../../editingsamples.txt", what="character")
ids = ids[c(grep("L001", ids),25:26)]
editingres = list()
for (i in 1:length(ids)){
  editingres[[i]] = read.table(ids[i])
  colnames(editingres[[i]]) = c("chromosome", "start", "end", "ref", "alt", "depth", "valdepth", "ref.count", "alt.count")
}
names(editingres) = gsub("\\_.*","",ids)
names(editingres) = c(names(editingres)[1:12], "Br5339C1.downsampled", "Br5340C1.downsampled")
save(editingres, file="/media/DATA/Amanda/rna_editing_results.rda")

### ON LOCAL COMPUTER ###

library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(VennDiagram)
library(plyr)
library(clusterProfiler)
require("org.Hs.eg.db")

load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/rna_editing/data/rna_editing_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Explore results
lapply(editingres, head)
editingres = editingres[c(1:6,8,10:14)]
names(editingres) = gsub(".downsampled", "", names(editingres))
editingres = Map(cbind, editingres, conversion = lapply(editingres, function(x) paste0(x$ref, ":", x$alt)),
                 sampleID = lapply(names(editingres), function(x) x), rnum = lapply(editingres, function(x) 1:nrow(x)),
                 editingID = lapply(editingres, function(x) paste0(x$chromosome,":",x$start,"-",x$end,":",x$ref,":",x$alt)))
for (i in 1:length(editingres)){
  tmp = editingres[[i]]
  tmp$Fraction = ifelse((tmp$rnum %in% grep("N", tmp$sampleID)), "Nucleus", "Cytosol")
  tmp$Age = ifelse((tmp$rnum %in% grep("53", tmp$sampleID)), "Prenatal", "Adult")
  tmp$Group = paste0(tmp$Age, ":", tmp$Fraction)
  tmp$collapsedconversion = NA
  tmp[which(tmp$conversion=="A:C" | tmp$conversion=="T:G"), "collapsedconversion"] = "A:C / T:G"
  tmp[which(tmp$conversion=="A:G" | tmp$conversion=="T:C"), "collapsedconversion"] = "A:G / T:C"
  tmp[which(tmp$conversion=="A:T" | tmp$conversion=="T:A"), "collapsedconversion"] = "A:T / T:A"
  tmp[which(tmp$conversion=="C:A" | tmp$conversion=="G:T"), "collapsedconversion"] = "C:A / G:T"
  tmp[which(tmp$conversion=="C:G" | tmp$conversion=="G:C"), "collapsedconversion"] = "C:G / G:C"
  tmp[which(tmp$conversion=="C:T" | tmp$conversion=="G:A"), "collapsedconversion"] = "C:T / G:A"
  editingres[[i]] = tmp
}
editingres_df = do.call(rbind, editingres)
editingres_df$rate = editingres_df$alt.count / editingres_df$valdepth

# how many editing events identified, and what type?
min(elementNROWS(editingres)) # 3064
max(elementNROWS(editingres)) # 5840
dim(editingres_df) # 52498
length(unique(editingres_df$editingID)) # 25051

editingres_df = data.table(editingres_df)
editingres_df[collapsedconversion=="A:G / T:C",length(unique(editingID)),] # 18907
editingres_df[, length(unique(editingID)), by = "collapsedconversion"]
#collapsedconversion    V1
#1:           A:G / T:C 18907
#2:           C:T / G:A  3259
#3:           A:T / T:A   663
#4:           C:G / G:C   802
#5:           C:A / G:T   775
#6:           A:C / T:G   645
editingres_df[, length(unique(editingID)), by = c("Fraction", "Age")]
#Fraction      Age    V1
#1:  Cytosol    Adult  8140
#2:  Nucleus    Adult 12749
#3:  Nucleus Prenatal  9926
#4:  Cytosol Prenatal  8116

editingtypecounts = editingres_df[, length(unique(editingID)), by = c("collapsedconversion", "Fraction", "Age")]
editingtypecounts = as.data.frame(editingtypecounts)
editingtypecounts$collapsedconversion = factor(editingtypecounts$collapsedconversion,
                                               levels = c("A:G / T:C","A:C / T:G","A:T / T:A","C:A / G:T","C:G / G:C","C:T / G:A"))

ggplot(editingtypecounts, aes(x = Fraction, y = V1, fill = collapsedconversion)) + geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Unique RNA Editing Sites by Group") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Unique_RNA_Editing_Sites_byGroup.pdf

# What is the range of editing rate by group?
editingres_df[, list(min = min(na.omit(rate)),max = max(na.omit(rate)), median = median(na.omit(rate)), 
                     mean = mean(na.omit(rate)), sd = sd(na.omit(rate))), by = "sampleID"]
#sampleID min max median      mean        sd
#1: Br1113C1   0   1    0.5 0.5543089 0.2820572
#2: Br1113N1   0   1    0.5 0.5096131 0.2460976
#3:  Br2046C   0   1    0.5 0.5325636 0.2534987
#4:  Br2046N   0   1    0.5 0.5250878 0.2456670
#5:  Br2074C   0   1    0.5 0.5324752 0.2562208
#6:  Br2074N   0   1    0.5 0.5286539 0.2455392
#7: Br5339N1   0   1    0.5 0.5322611 0.2559207
#8: Br5340N1   0   1    0.5 0.5358024 0.2545051
#9: Br5341C1   0   1    0.5 0.5450855 0.2782518
#10: Br5341N1   0   1    0.5 0.5506776 0.2700437
#11: Br5339C1   0   1    0.5 0.5448914 0.2719962
#12: Br5340C1   0   1    0.5 0.5321418 0.2599204

# What is the validated coverage range per sample at edited sites?
editingres_df[, list(min = min(valdepth),max = max(valdepth), median = median(valdepth), 
                     mean = mean(valdepth), sd = sd(valdepth)), by = c("sampleID")]
#    sampleID min max median     mean       sd
#1: Br1113C1   0 109     12 16.91928 14.99633
#2: Br1113N1   2 108     12 16.10279 12.73881
#3:  Br2046C   1 108     12 16.33267 14.73171
#4:  Br2046N   1 109     12 15.61522 13.31537
#5:  Br2074C   1 109     12 16.79319 15.28954
#6:  Br2074N   0 107     12 15.86815 14.13114
#7: Br5339N1   1 107     11 14.41077 11.70882
#8: Br5340N1   0 109     11 14.75442 12.21129
#9: Br5341C1   2 107     12 16.35597 14.74964
#10: Br5341N1   0 107     11 15.53969 14.19272
#11: Br5339C1   0 108     11 15.18437 12.81806
#12: Br5340C1   0 107     11 15.21979 12.80784

ggplot(editingres_df, aes(x = sampleID, y = valdepth, fill = collapsedconversion)) + geom_boxplot() +
  facet_grid(. ~ Group, scales = "free") +
  labs(fill="") +
  ylab("Read Depth\n(Filtered By Read Quality)") + 
  xlab("") +
  ggtitle("Validated Coverage Range By Sample At Edited Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# valDepth_byEditingType_bySampleID_byGroup.pdf

ggplot(editingres_df[which(editingres_df$collapsedconversion=="A:G / T:C"),], 
       aes(x = sampleID, y = valdepth)) + geom_boxplot() +
  facet_grid(. ~ Group, scales = "free") +
  labs(fill="") +
  ylab("Read Depth\n(Filtered By Read Quality)") + 
  xlab("") +
  ggtitle("Validated Coverage Range By Sample At Edited Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# valDepth_bySampleID_byGroup_AtoG_only.pdf

# Annotate editing sites to features in the genome
txdb = loadDb("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
features = list(genes = genes(txdb), CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
lapply(features, head)

editingres_df$strand = ifelse(editingres_df$conversion=="A:G","+","-")
editingres_df[editingres_df$collapsedconversion!="A:G / T:C","strand"] = "*"
grediting = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
editedgenes = findOverlaps(grediting, features$genes)
genes = as.data.frame(features$genes)
x = cbind(editingres_df[queryHits(editedgenes),], overlappingGene = genes[subjectHits(editedgenes),"gene_id"])
editingres_df = rbind(x, data.frame(editingres_df[-unique(queryHits(editedgenes)),], overlappingGene = "NA"))
grediting = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grediting, y))

grediting$rnum = 1:length(grediting)
grediting$cds = ifelse(grediting$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grediting$intron = ifelse(grediting$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grediting$UTR5 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grediting$UTR3 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grediting$anno = paste0(grediting$cds,":",grediting$intron, ":", grediting$UTR5, ":", grediting$UTR3)

editing = as.data.frame(grediting)
editing[which(editing$anno == "NA:NA:NA:NA"),"annotation"] = "Intergenic" 
editing[grep("CDS", editing$cds),"annotation"] = "CDS"
editing[which(is.na(editing$annotation) & editing$UTR3 == "UTR3"),"annotation"] = "UTR3"
editing[which(is.na(editing$annotation) & editing$UTR5 == "UTR5"),"annotation"] = "UTR5"
editing[which(is.na(editing$annotation) & editing$intron == "Intron"),"annotation"] = "Intron"

# Mapping editing sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grediting, geneMapGR)
editing$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
editing$nearestID = names(geneMapGR)[subjectHits(dA)]
editing$distToGene = mcols(dA)$distance
editing$EntrezID = ifelse(editing$overlappingGene!="NA", geneMap[match(editing$overlappingGene, geneMap$gencodeID),"EntrezID"], geneMapGR$EntrezID[subjectHits(dA)])
editing_anno = data.table(editing)
editing_anno$annotation = gsub("UTR3","3'UTR",editing_anno$annotation)
editing_anno$annotation = gsub("UTR5","5'UTR", editing_anno$annotation)
editing_anno$annotation = factor(editing_anno$annotation, levels = c("CDS","Intron","3'UTR","5'UTR","Intergenic"))

# What is the distribution of features edited across groups?

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Genomic_features_editing_allSites.pdf", width = 7,height = 4)
ggplot(editing_anno[,length(unique(editingID)), by = c("annotation", "Fraction", "Age")], 
       aes(x = Fraction, y = V1, fill = annotation)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Number of RNA Editing Sites\n by Feature and Group") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Genomic_features_editing_AtoGOnly.pdf", width = 7,height = 4)
ggplot(editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("annotation", "Fraction", "Age")], 
       aes(x = Fraction, y = V1, fill = annotation)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Number of RNA Editing Sites\nby Feature and Group") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

intron = editing_anno[collapsedconversion=="A:G / T:C" & annotation=="Intron",length(unique(editingID)), by = c("annotation", "Fraction", "Age")]
utr3 = editing_anno[collapsedconversion=="A:G / T:C" & annotation=="3'UTR",length(unique(editingID)), by = c("annotation", "Fraction", "Age")]
total = editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction", "Age")]
(intron$V1/total$V1 * 100)[order(intron$V1/total$V1 * 100)] # 21.67809 33.73761 33.75515 26.02888
(utr3$V1/total$V1 * 100)[order(utr3$V1/total$V1 * 100)] # 37.56576 40.39850 43.84477 50.83252


### In editing sites shared between all four groups:
#How many sites are present in all groups?
intersecting = Reduce(intersect, list(as.character(editing_anno[Group=="Adult:Cytosol", editingID,]),
                                      as.character(editing_anno[Group=="Adult:Nucleus",editingID]), 
                                      as.character(editing_anno[Group=="Prenatal:Cytosol",editingID]),
                                      as.character(editing_anno[Group=="Prenatal:Nucleus",editingID])))
length(intersecting) # 1253
intersectres = editing_anno[editingID %in% intersecting,,]
intersectres = intersectres[collapsedconversion=="A:G / T:C",,]
length(unique(intersectres[,editingID,])) # 1025 of 1253 are A to I editing events
cov = intersectres[, list(min = min(valdepth),max = max(valdepth), mean = mean(valdepth), sd = sd(valdepth)), by = c("editingID","Group")]

# Identify differentially edited sites
subset = intersectres[,list(editingID,sampleID,rate),]
propmat = matrix(NA, nrow=length(unique(subset[,editingID,])), ncol=12)
rownames(propmat) = unique(subset[,editingID,])
colnames(propmat) = unique(subset[,sampleID,])

for (i in 1:nrow(intersectres)) {
  rowIdx=which(intersectres$editingID[i]==rownames(propmat))
  colIdx=which(intersectres$sampleID[i]==colnames(propmat))
  propmat[rowIdx,colIdx]=intersectres$rate[i]
}

result = data.frame(id=rownames(propmat),fstat=NA,df=NA,Age=NA,Fraction=NA,Interaction=NA)
for (i in 1:nrow(propmat)) {
  print(i)
  age = c(rep.int("Adult",6), rep.int("Prenatal",6))
  fraction = rep.int(c("Cytosol","Nucleus"),6)
  y=log(propmat[i,]/(1-propmat[i,])) # logit transformation of rate
  idx=which(!is.nan(y) & is.finite(y))
  if (length(idx)<5) {next}
  else {
    y=y[idx]; age=age[idx]; fraction=fraction[idx];
    if (length(unique(age))>1 & length(unique(fraction))>1) {
      fit=lm(y~age+fraction+age*fraction)
      result[i,2:3]=summary(fit)$fstatistic[1:2]
      result[i,4:6]=summary(fit)$coefficients[,4][2:4]
    }
  }
}
result$fdr.age = p.adjust(result$Age, method = "fdr")
result$fdr.frac = p.adjust(result$Fraction, method = "fdr")
result$fdr.int = p.adjust(result$Interaction, method = "fdr")
head(na.omit(result[order(result$fdr.int),]))
head(na.omit(result[order(result$fdr.frac),]))
head(na.omit(result[order(result$fdr.age),]))
result = melt(result, id.vars = c("id","fstat","df","fdr.age","fdr.frac","fdr.int"))
result$variable = factor(result$variable, levels = c("Age","Fraction","Interaction"))


## Plot the pvalue distributions

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Pvalue_distributions_RNAediting.pdf",height = 4,width = 9)
ggplot(data=result, aes(result$value)) + geom_histogram(bins=20) + facet_grid(. ~ variable) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value Distribution") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

editing.sig = list(byAge = result[which(result$fdr.age<=0.05),"id"],
                   byFrac = result[which(result$fdr.frac<=0.05),"id"],
                   Interaction = result[which(result$fdr.int<=0.05),"id"],
                   Int.Unadjusted = result[which(result$pval.int<=0.05),"id"]) 
elementNROWS(editing.sig)
editing.sig = lapply(editing.sig, function(x) editing_anno[editing_anno$editingID %in% x])

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/differential_editing_byAge.pdf")
dat = as.data.frame(editing.sig$byAge)
for (i in 1:length(unique(dat$editingID))){
  ids = unique(dat$editingID)
  g = ggplot(dat[which(dat$editingID==ids[i]),], 
             aes(x = Age, y = rate, fill = Fraction)) + geom_boxplot() + geom_jitter() +
             ylab("Editing Rate") + ylim(0,1) +
    scale_fill_brewer(palette="Dark2") +
             xlab("") +
             ggtitle(ids[i]) +
             theme(title = element_text(size = 20)) +
             theme(text = element_text(size = 20))
  print(g)
}
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/differential_editing_byFrac.pdf")
dat = as.data.frame(editing.sig$byFrac)
for (i in 1:length(unique(dat$editingID))){
  ids = unique(dat$editingID)
  g = ggplot(dat[which(dat$editingID==ids[i]),], 
             aes(x = Age, y = rate, fill = Fraction)) + geom_boxplot() + geom_jitter() +
    ylab("Editing Rate") + ylim(0,1) +
    scale_fill_brewer(palette="Dark2") +
    xlab("") +
    ggtitle(ids[i]) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)
}
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/differential_editing_Interaction.pdf")
dat = as.data.frame(editing.sig$Interaction)
for (i in 1:length(unique(dat$editingID))){
  ids = as.character(unique(dat$editingID))
  g = ggplot(dat[which(dat$editingID==ids[i]),], 
             aes(x = Age, y = rate, fill = Fraction)) + geom_boxplot() + geom_jitter() +
    ylab("Editing Rate") + ylim(0,1) +
    scale_fill_brewer(palette="Dark2") +
    xlab("") +
    ggtitle(ids[i]) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)
}
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/differential_editing_Interaction_unadjusted.pdf")
dat = as.data.frame(editing.sig$Int.Unadjusted)
for (i in 1:length(unique(dat$editingID))){
  ids = unique(dat$editingID)
  g = ggplot(dat[which(dat$editingID==ids[i]),], 
             aes(x = Age, y = rate, fill = Fraction)) + geom_boxplot() + geom_jitter() +
    ylab("Editing Rate") + ylim(0,1) +
    xlab("") + scale_fill_brewer(palette="Dark2") +
    ggtitle(ids[i]) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)
}
dev.off()


## How about split by age?
# In Adult
result.adult = data.frame(id=rownames(propmat[,1:6]),fstat=NA,df=NA,pval.frac=NA)
for (i in 1:nrow(propmat[,1:6])) {
  print(i)
  age = "Adult"
  fraction = rep.int(c("Cytosol","Nucleus"),3)
  y=log(propmat[i,1:6]/(1-propmat[i,1:6])) # logit transformation of rate
  idx=which(!is.nan(y) & is.finite(y))
  if (length(idx)<5) {next}
  else {
    y=y[idx]; age=age[idx]; fraction=fraction[idx];
    if (length(unique(age))>1 & length(unique(fraction))>1) {
      fit=lm(y~fraction)
      result.adult[i,2:3]=summary(fit)$fstatistic[1:2]
      result.adult[i,4]=summary(fit)$coefficients[,4][2]
    }
  }
}
result.adult$fdr.frac = p.adjust(result.adult$pval.frac, method = "fdr")

# In Prenatal
result.prenatal = data.frame(id=rownames(propmat[,7:12]),fstat=NA,df=NA,pval.frac=NA)
for (i in 1:nrow(propmat[,7:12])) {
  print(i)
  age = "prenatal"
  fraction = rep.int(c("Cytosol","Nucleus"),3)
  y=log(propmat[i,7:12]/(1-propmat[i,7:12])) # logit transformation of rate
  idx=which(!is.nan(y) & is.finite(y))
  if (length(idx)<5) {next}
  else {
    y=y[idx]; age=age[idx]; fraction=fraction[idx];
    if (length(unique(age))>1 & length(unique(fraction))>1) {
      fit=lm(y~fraction)
      result.prenatal[i,2:3]=summary(fit)$fstatistic[1:2]
      result.prenatal[i,4]=summary(fit)$coefficients[,4][2]
    }
  }
}
result.prenatal$fdr.frac = p.adjust(result.prenatal$pval.frac, method = "fdr")

# Plot 
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Pvalue_ByFrac_inAdultOnly_inPrenatalOnly_RNAediting.pdf", height = 6, width = 8)
ggplot(data=result.adult, aes(result.adult$pval.frac)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Fraction in Adult") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(data=result.prenatal, aes(result.prenatal$pval.frac)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Fraction in Prenatal") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

AandP.editing.sig = list(Adult.byFrac = result.adult[which(result.adult$fdr.frac<=0.05),"id"],
                         Prenatal.byFrac = result.prenatal[which(result.prenatal$fdr.frac<=0.05),"id"]) 
elementNROWS(AandP.editing.sig) # None meet the FDR threshold


## Editing sites present in one group but not another
editingID = editing_anno[collapsedconversion=="A:G / T:C",list(editingID), by = "Group"]
editingID = split(editingID, f = editingID$Group)
editingID = lapply(editingID, function(x) as.character(x$editingID))
lapply(editingID, head)

venn.diagram(editingID, "./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/editing_site_overlap.jpeg", 
             main="Overlap of Identified Editing Sites By Group", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)

AGonly = editing_anno[collapsedconversion=="A:G / T:C",,]
cyt = AGonly[Fraction=="Cytosol",,]
nuc = AGonly[Fraction=="Nucleus",,]
ad = AGonly[Age=="Adult",,]
pren = AGonly[Age=="Prenatal",,]
AC = AGonly[Group=="Adult:Cytosol",,]
AN = AGonly[Group=="Adult:Nucleus",,]
PC = AGonly[Group=="Prenatal:Cytosol",,]
PN = AGonly[Group=="Prenatal:Nucleus",,]

unique = list(cytosolOnly = cyt[!(editingID %in% nuc$editingID),,], nucleusOnly = nuc[!(editingID %in% cyt$editingID),,], 
              adultOnly = ad[!(editingID %in% pren$editingID),,], prenatalOnly = pren[!(editingID %in% ad$editingID),,], 
              ANnotAC = AN[!(editingID %in% AC$editingID),,], ACnotAN = AC[!(editingID %in% AN$editingID),,], 
              ANnotPN = AN[!(editingID %in% PN$editingID),,], PNnotAN = PN[!(editingID %in% AN$editingID),,],
              ACnotPC = AC[!(editingID %in% PC$editingID),,], PCnotAC = PC[!(editingID %in% AC$editingID),,], 
              PCnotPN = PC[!(editingID %in% PN$editingID),,], PNnotPC = PN[!(editingID %in% PC$editingID),,])
elementNROWS(unique)
all = list(cytosolAll = cyt, nucleusAll = nuc, adultAll = ad, prenatalAll = pren, allAC = AC, allAN = AN, allPC = PC, allPN = PN)

## Proportion of edited sites identified in x that are unique to x

prop = unlist(lapply(unique, function(x) x[,length(unique(editingID)),]))
total = list(cytosol = AGonly[Fraction=="Cytosol",length(unique(editingID)),], nucleus = AGonly[Fraction=="Nucleus",length(unique(editingID)),],
             adult = AGonly[Age=="Adult",length(unique(editingID)),], prenatal = AGonly[Age=="Prenatal",length(unique(editingID)),], 
             AC = AGonly[Group=="Adult:Cytosol",length(unique(editingID)),], AN = AGonly[Group=="Adult:Nucleus",length(unique(editingID)),], 
             PN = AGonly[Group=="Prenatal:Nucleus",length(unique(editingID)),], PC = AGonly[Group=="Prenatal:Cytosol",length(unique(editingID)),])
prop = cbind(comparison = rownames(prop), total = c(total$cytosol, total$nucleus, total$adult, total$prenatal,
                                                    total$AN, total$AC, total$AN, total$PN, total$AC, total$PC, total$PC, total$PN), 
             unique = prop, proportion = c(prop["cytosolOnly"]/total$cytosol, prop["nucleusOnly"]/total$nucleus,
                                           prop["adultOnly"]/total$adult,prop["prenatalOnly"]/total$prenatal,
                                           prop["ANnotAC"]/total$AN, prop["ACnotAN"]/total$AC,
                                           prop["ANnotPN"]/total$AN, prop["PNnotAN"]/total$PN,
                                           prop["ACnotPC"]/total$AC, prop["PCnotAC"]/total$PC,
                                           prop["PCnotPN"]/total$PC, prop["PNnotPC"]/total$PN))
prop = as.data.frame(prop)
prop[order(prop[,3]),]
prop$Group = factor(rownames(prop), levels = c("cytosolOnly", "nucleusOnly", "adultOnly", "prenatalOnly", "ANnotAC",
                                                 "ACnotAN", "ANnotPN", "PNnotAN", "ACnotPC", "PCnotAC", "PCnotPN", "PNnotPC"))

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/unique_editing_sites_inEach_category.pdf", width = 8, height = 6)
ggplot(prop, aes(x = Group, y = unique)) + geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(100*proportion,1),"%")), vjust = -.5) +
  labs(y = "Percent", fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Unique Editing Sites in Each Category") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()

# Check coverage by unique editing group
for (i in 1:length(unique)){  unique[[i]] = cbind(unique[[i]], Status = names(unique)[i])  }
unique_dt = do.call(rbind, unique)
unique_dt$Status = factor(unique_dt$Status, levels = c("cytosolOnly", "nucleusOnly", "adultOnly", "prenatalOnly", "ANnotAC",
                                                 "ACnotAN", "ANnotPN", "PNnotAN", "ACnotPC", "PCnotAC", "PCnotPN", "PNnotPC"))

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/valdepth_unique_editing_sites_byGroup_bySample.pdf", width = 24, height = 6)
ggplot(unique_dt, aes(x = sampleID, y = valdepth, fill = Group)) + geom_boxplot() +
  facet_grid(. ~ Status, scales = "free") +
  labs(fill="") +
  ylab("High Quality Read Depth)") + ylim(0,30) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  ggtitle("Validated Coverage Range By Group At Uniquely Edited Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Are editing sites in a specific annotation more likely to be unique to a group?

x = lapply(unique, function(x) as.data.frame(cbind(x[,length(unique(editingID)),by="annotation"], total = length(unique(x$editingID)))))
x = Map(cbind, x, diff = lapply(x, function(y) y$total - y$V1))
comps = list(byFraction = c("cytosolOnly","nucleusOnly"), byAge = c("adultOnly","prenatalOnly"), byFracInAdult = c("ACnotAN","ANnotAC"),
             byFracinPrenatal = c("PCnotPN","PNnotPC"), byAgeinCyt = c("ACnotPC","PCnotAC"), byAgeinNuc = c("ANnotPN", "PNnotAN"))
anno = c("CDS","Intron","Intergenic","3'UTR","5'UTR")

## Are  more likely to be in cyt over nuc or vice versa?
anno.comps = list(list(),list(),list(),list(),list(),list())
for (i in (1:length(comps))){
  for (j in (1:length(anno))){
    anno.comps[[i]][[j]] = data.frame(Cyt.Adult = c(x[[comps[[i]][1]]][which(x[[comps[[i]][1]]][,"annotation"]==anno[j]),"V1"],
                                                    x[[comps[[i]][1]]][which(x[[comps[[i]][1]]][,"annotation"]==anno[j]),"diff"]),
                                      Nuc.Prenatal = c(x[[comps[[i]][2]]][which(x[[comps[[i]][1]]][,"annotation"]==anno[j]),"V1"],
                                                       x[[comps[[i]][2]]][which(x[[comps[[i]][1]]][,"annotation"]==anno[j]),"diff"]), row.names = c("inAnno","notAnno"))
  }
  names(anno.comps[[i]]) = anno
}
names(anno.comps) = names(comps)
anno.comps = lapply(anno.comps, lapply, fisher.test)

x = do.call(rbind, Map(cbind, Group = as.list(names(anno.comps)), lapply(anno.comps, function(a) 
            do.call(rbind, Map(cbind, Annotation = as.list(names(a)), lapply(a, function(z) data.frame(pval = z$p.value, OddsRatio = z$estimate, row.names = NULL)))))))
x$FDR = p.adjust(x$pval, method = "fdr")
write.csv(x,quote=F, file= "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_annotationEnrichment_inUniqueSites_byGroup.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_annotationEnrichment_inUniqueSites_byGroup.csv")
df[df$FDR<=0.05 & df$OddsRatio>1,colnames(df)!="X"]
#              Group Annotation          pval OddsRatio           FDR
#3        byFraction Intergenic  2.730282e-06  1.223474  3.561237e-06
#4        byFraction      3'UTR  0.000000e+00 82.370592  0.000000e+00
#8             byAge Intergenic  0.000000e+00 27.361875  0.000000e+00
#9             byAge      3'UTR  1.315546e-20  1.380692  1.973319e-20
#12    byFracInAdult     Intron 1.117287e-256 34.364053 2.578354e-256
#13    byFracInAdult Intergenic  1.526726e-05  1.259822  1.908407e-05
#14    byFracInAdult      3'UTR  1.983926e-10  1.353012  2.834180e-10
#18 byFracinPrenatal Intergenic 6.300473e-270 50.507131 1.575118e-269
#19 byFracinPrenatal      3'UTR  3.522626e-34  1.972863  5.871044e-34
#23       byAgeinCyt Intergenic 2.545671e-280 22.697174 6.942740e-280
#24       byAgeinCyt      3'UTR  1.250038e-63  2.130625  2.343820e-63
#27       byAgeinNuc     Intron  2.230578e-07  1.225238  3.041697e-07
#28       byAgeinNuc Intergenic  0.000000e+00 24.117045  0.000000e+00
#29       byAgeinNuc      3'UTR  0.000000e+00 61.294598  0.000000e+00
df[df$FDR<=0.05 & df$OddsRatio<1,colnames(df)!="X"]
#              Group Annotation          pval  OddsRatio           FDR
#1        byFraction        CDS  0.000000e+00 0.02523414  0.000000e+00
#2        byFraction     Intron  6.411616e-75 0.46490368  1.282323e-74
#6             byAge        CDS  0.000000e+00 0.02737029  0.000000e+00
#11    byFracInAdult        CDS 5.323224e-311 0.02267203 1.774408e-310
#16 byFracinPrenatal        CDS 2.785537e-240 0.02228103 5.969007e-240
#17 byFracinPrenatal     Intron  9.482181e-28 0.54250043  1.497186e-27
#21       byAgeinCyt        CDS 8.794660e-309 0.04676909 2.638398e-308
#22       byAgeinCyt     Intron  3.371831e-57 0.46696500  5.950290e-57
#25       byAgeinCyt      5'UTR  2.643376e-02 0.55454641  3.172052e-02
#26       byAgeinNuc        CDS  0.000000e+00 0.02697268  0.000000e+00
#30       byAgeinNuc      5'UTR  0.000000e+00 0.01657080  0.000000e+00



## Gene ontology of different groups of editing sites

entrez.editing = lapply(unique, function(x) as.character(unique(na.omit(x$EntrezID))))

split.anno = lapply(unique, function(x) split(x, f = x$annotation))
cytosolOnly = split.anno[["cytosolOnly"]]
nucleusOnly = split.anno[["nucleusOnly"]] 
adultOnly = split.anno[["adultOnly"]]   
prenatalOnly = split.anno[["prenatalOnly"]]
ANnotAC = split.anno[["ANnotAC"]]   
ACnotAN = split.anno[["ACnotAN"]]     
ANnotPN = split.anno[["ANnotPN"]]     
PNnotAN = split.anno[["PNnotAN"]]     
ACnotPC = split.anno[["ACnotPC"]]     
PCnotAC = split.anno[["PCnotAC"]]    
PCnotPN = split.anno[["PCnotPN"]]
PNnotPC = split.anno[["PNnotPC"]]

cytosolOnly = lapply(cytosolOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
nucleusOnly = lapply(nucleusOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
adultOnly = lapply(adultOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
prenatalOnly = lapply(prenatalOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
ANnotAC = lapply(ANnotAC, function(x) unique(as.character(na.omit(x$EntrezID))))
ACnotAN = lapply(ACnotAN, function(x) unique(as.character(na.omit(x$EntrezID))))
ANnotPN = lapply(ANnotPN, function(x) unique(as.character(na.omit(x$EntrezID))))
PNnotAN = lapply(PNnotAN, function(x) unique(as.character(na.omit(x$EntrezID))))
ACnotPC = lapply(ACnotPC, function(x) unique(as.character(na.omit(x$EntrezID))))
PCnotAC = lapply(PCnotAC, function(x) unique(as.character(na.omit(x$EntrezID))))
PCnotPN = lapply(PCnotPN, function(x) unique(as.character(na.omit(x$EntrezID))))
PNnotPC = lapply(PNnotPC, function(x) unique(as.character(na.omit(x$EntrezID))))

## Compare the enriched terms between the unsplit groups
# KEGG
compareKegg = compareCluster(entrez.editing, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrez.editing, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrez.editing, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrez.editing, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrez.editing, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/interaction.kegg.GO.DO.objects.RNAediting.rda")
# Plot
pdf(file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/interaction.kegg.GO.DO_unsplit_by_annotation.pdf", width=22, height=14)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO,colorBy="p.adjust",  showCategory = 45, title= "Disease Ontology Enrichment")
dev.off()

## Compare the enriched terms between the split groups
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ANnotAC","ACnotAN","ANnotPN","PNnotAN","ACnotPC","PCnotAC","PCnotPN","PNnotPC")
# KEGG
Kegg.cytosolOnly = compareCluster(cytosolOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.nucleusOnly = compareCluster(nucleusOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.adultOnly = compareCluster(adultOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.prenatalOnly = compareCluster(prenatalOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ANnotAC = compareCluster(ANnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ACnotAN = compareCluster(ACnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ANnotPN = compareCluster(ANnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PNnotAN = compareCluster(PNnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ACnotPC = compareCluster(ACnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PCnotAC = compareCluster(PCnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PCnotPN = compareCluster(PCnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PNnotPC = compareCluster(PNnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(Kegg.cytosolOnly,Kegg.nucleusOnly,Kegg.adultOnly,Kegg.prenatalOnly,Kegg.ANnotAC,Kegg.ACnotAN,Kegg.ANnotPN,Kegg.PNnotAN,Kegg.ACnotPC,Kegg.PCnotAC,Kegg.PCnotPN,Kegg.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/kegg.objects.RNAediting.SplitByAnnotation.rda")

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Kegg_split_by_annotation.pdf", width = 10)
plot(Kegg.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: cytosolOnly")
plot(Kegg.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: adultOnly")
plot(Kegg.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: prenatalOnly")
plot(Kegg.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotAC")
plot(Kegg.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotAN")
plot(Kegg.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotPN")
plot(Kegg.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotAN")
plot(Kegg.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotPC")
plot(Kegg.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotAC")
plot(Kegg.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotPN")
dev.off()

# Biological Process
BP.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(BP.cytosolOnly,BP.nucleusOnly,BP.adultOnly,BP.ANnotAC,BP.ACnotAN,BP.ANnotPN,BP.ACnotPC,BP.PCnotAC,BP.PCnotPN,BP.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/BP.objects.RNAediting.SplitByAnnotation.rda")

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/BP_split_by_annotation.pdf", height = 6, width = 10.5)
plot(BP.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: cytosolOnly")
plot(BP.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: nucleusOnly")
plot(BP.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: adultOnly")
plot(BP.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotAC")
plot(BP.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ACnotAN")
plot(BP.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotPN")
plot(BP.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotAC")
plot(BP.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotPN")
plot(BP.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PNnotPC")
dev.off()

# Molecular Function
MF.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(MF.cytosolOnly,MF.nucleusOnly,MF.adultOnly,MF.ANnotAC,MF.ACnotAN,MF.ANnotPN,MF.ACnotPC,MF.PCnotPN,MF.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/MF.objects.RNAediting.SplitByAnnotation.rda")

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/MF_split_by_annotation.pdf", height = 8, width = 12)
plot(MF.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: cytosolOnly")
plot(MF.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: nucleusOnly")
plot(MF.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: adultOnly")
plot(MF.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotAC")
plot(MF.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotAN")
plot(MF.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotPN")
plot(MF.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotPC")
plot(MF.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PCnotPN")
plot(MF.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PNnotPC")
dev.off()

# Cellular Component
CC.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(CC.cytosolOnly,CC.nucleusOnly,CC.adultOnly,CC.prenatalOnly,CC.ANnotAC,CC.ACnotAN,CC.ANnotPN,CC.PNnotAN,CC.ACnotPC,CC.PCnotAC,CC.PCnotPN,CC.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/CC.objects.RNAediting.SplitByAnnotation.rda")

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/CC_split_by_annotation.pdf", height = 8, width = 10)
plot(CC.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: cytosolOnly")
plot(CC.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: nucleusOnly")
plot(CC.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: adultOnly")
plot(CC.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: prenatalOnly")
plot(CC.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ANnotAC")
plot(CC.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotAN")
plot(CC.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ANnotPN")
plot(CC.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotAN")
plot(CC.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotPC")
plot(CC.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PCnotAC")
plot(CC.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PCnotPN")
plot(CC.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotPC")
dev.off()

# Disease Ontology
DO.cytosolOnly = compareCluster(cytosolOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.nucleusOnly = compareCluster(nucleusOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.adultOnly = compareCluster(adultOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.prenatalOnly = compareCluster(prenatalOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ANnotAC = compareCluster(ANnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ACnotAN = compareCluster(ACnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ANnotPN = compareCluster(ANnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PNnotAN = compareCluster(PNnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ACnotPC = compareCluster(ACnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PCnotAC = compareCluster(PCnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PCnotPN = compareCluster(PCnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PNnotPC = compareCluster(PNnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(DO.adultOnly,DO.prenatalOnly,DO.ACnotAN,DO.ANnotPN,DO.PNnotAN,DO.ACnotPC,DO.PCnotAC,DO.PCnotPN, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/DO.objects.RNAediting.SplitByAnnotation.rda")

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/DO_split_by_annotation.pdf")
plot(DO.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: adultOnly")
plot(DO.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: prenatalOnly")
plot(DO.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ACnotAN")
plot(DO.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotPN")
plot(DO.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PNnotAN")
plot(DO.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ACnotPC")
plot(DO.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotAC")
plot(DO.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotPN")
dev.off()