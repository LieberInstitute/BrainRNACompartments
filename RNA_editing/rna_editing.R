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
#sampleID min max median     mean       sd
#1: Br1113C1   0 109     12 16.99608 15.02656
#2: Br1113N1   2 108     12 16.10694 12.74313
#3:  Br2046C   1 108     12 16.24075 14.81553
#4:  Br2046N   1 109     12 15.50497 13.21504
#5:  Br2074C   1 109     12 16.79001 15.35052
#6:  Br2074N   0 107     12 15.77555 14.04571
#7: Br5339N1   1 107     11 14.41094 11.83242
#8: Br5340N1   0 109     11 14.80892 12.41942
#9: Br5341C1   2 107     12 16.62589 15.18498
#10: Br5341N1   0 107     11 15.73896 14.59073
#11: Br5339C1   0 108     12 15.32501 13.00267
#12: Br5340C1   0 107     11 15.29666 13.01752

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
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
lapply(features, head)

grediting = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grediting, y))

grediting$rnum = 1:length(grediting)
grediting$cds = ifelse(grediting$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grediting$intron = ifelse(grediting$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grediting$UTR5 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grediting$UTR3 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grediting$anno = paste0(grediting$cds,":",grediting$intron, ":", grediting$UTR5, ":", grediting$UTR3)

editing = as.data.frame(grediting)
editing[which(editing$anno == "NA:NA:NA:NA"),"annotation"] = "Other" 
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
editing$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
editing_anno = data.table(editing)

# What is the distribution of features edited across groups?

ggplot(editing_anno[,length(unique(editingID)), by = c("annotation", "Fraction", "Age")], 
       aes(x = Fraction, y = V1, fill = annotation)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Number of RNA Editing Sites\n by Feature and Group") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Genomic_features_editing_allSites.pdf

ggplot(editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("annotation", "Fraction", "Age")], 
       aes(x = Fraction, y = V1, fill = annotation)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Number of RNA Editing Sites\nby Feature and Group") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Genomic_features_editing_AtoGOnly.pdf

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

result = data.frame(id=rownames(propmat),fstat=NA,df=NA,pval.age=NA,pval.frac=NA,pval.int=NA)
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
result$fdr.age = p.adjust(result$pval.age, method = "fdr")
result$fdr.frac = p.adjust(result$pval.frac, method = "fdr")
result$fdr.int = p.adjust(result$pval.int, method = "fdr")
head(na.omit(result[order(result$fdr.int),]))
head(na.omit(result[order(result$fdr.frac),]))
head(na.omit(result[order(result$fdr.age),]))

## Plot the pvalue distributions

ggplot(data=result, aes(result$pval.age)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Pvalue_ByAge_RNAediting.pdf

ggplot(data=result, aes(result$pval.frac)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Pvalue_ByFrac_RNAediting.pdf

ggplot(data=result, aes(result$pval.int)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Pvalue_Interaction_RNAediting.pdf

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
    xlab("") +
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
ggplot(data=result.adult, aes(result.adult$pval.frac)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Fraction in Adult") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Pvalue_ByFrac_AdultOnly_RNAediting.pdf

ggplot(data=result.prenatal, aes(result.prenatal$pval.frac)) + geom_histogram(bins=20) +
  ylab("Count") + 
  xlab("P-value") +
  ggtitle("P-value By Fraction in Prenatal") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# Pvalue_ByFrac_PrenatalOnly_RNAediting.pdf

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
prop

# Check coverage by unique editing group
for (i in 1:length(unique)){  unique[[i]] = cbind(unique[[i]], Status = names(unique)[i])  }
unique_dt = do.call(rbind, unique)
unique_dt$Status = factor(unique$Status, levels = c("cytosolOnly", "nucleusOnly", "adultOnly", "prenatalOnly", "ANnotAC",
                                                 "ACnotAN", "ANnotPN", "PNnotAN", "ACnotPC", "PCnotAC", "PCnotPN", "PNnotPC"))

ggplot(unique_dt, aes(x = sampleID, y = valdepth, fill = Group)) + geom_boxplot() +
  facet_grid(. ~ Status, scales = "free") +
  labs(fill="") +
  ylab("High Quality Read Depth)") + ylim(0,30) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  ggtitle("Validated Coverage Range By Group At Uniquely Edited Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# valdepth_unique_editing_sites_byGroup_bySample.pdf


## Are editing sites in a specific annotation more likely to be unique to a group?

elementNROWS(lapply(unique[1:4], function(x) unique(x$editingID)))
#cytosolOnly  nucleusOnly    adultOnly prenatalOnly
#   3830         8809         9539         6049
x = lapply(unique, function(x) x[,length(unique(editingID)),by="annotation"])
x = lapply(x, function(x) as.data.frame(x))
x = lapply(x, function(x) x[order(x$annotation),])
x = data.frame(annotation = x[[1]][,1], cyt = x[[1]][,2], nuc = x[[2]][,2], adult = x[[3]][,2], prenatal = x[[4]][,2])
#annotation  cyt  nuc adult prenatal
#1        CDS   53  116   173      132
#2     Intron 1255 4237  3635     2635
#3      Other  777 1617  1799     1138
#4       UTR3 1710 2748  3854     2069
#5       UTR5   35   91    78       75

## Are  more likely to be in cyt over nuc or vice versa?
fisher.test(data.frame(c(53,3830-53),c(116,8809-116))) # CDS
#p-value = 0.8005
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7433125 1.4710638
#sample estimates:
#  odds ratio 
#1.051571
fisher.test(data.frame(c(1255,3830-1255),c(4237,8809-4237))) # Intron
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4853220 0.5697537
#sample estimates:
#  odds ratio 
#0.5259342
fisher.test(data.frame(c(777,3830-777),c(1617,8809-1617))) # Other
#p-value = 0.01175
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.027431 1.246496
#sample estimates:
#  odds ratio 
#1.131948
fisher.test(data.frame(c(1710,3830-1710),c(2748,8809-2748))) # 3'UTR
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.644115 1.924855
#sample estimates:
#  odds ratio 
#1.779005
fisher.test(data.frame(c(35,3830-35),c(91,8809-91))) # 5'UTR
#p-value = 0.5607
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5793647 1.3208958
#sample estimates:
#  odds ratio 
#0.8835633

## Are  more likely to be in adult over prenatal or vice versa?
fisher.test(data.frame(c(173,9539-173),c(132,6049-132))) # CDS
#p-value = 0.1093
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6547078 1.0490527
#sample estimates:
#  odds ratio 
#0.8279966
fisher.test(data.frame(c(3635,9539-3635),c(2635,6049-2635))) # Intron
#p-value = 1.421e-11
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7467143 0.8521855
#sample estimates:
#  odds ratio 
#0.7977237
fisher.test(data.frame(c(1799,9539-1799),c(1138,6049-1138))) # Other
#p-value = 0.9497
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9229963 1.0902577
#sample estimates:
#  odds ratio 
#1.003039
fisher.test(data.frame(c(3854,9539-3854),c(2069,6049-2069))) # 3'UTR
#p-value = 7.539e-15
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.218885 1.395338
#sample estimates:
#  odds ratio 
#1.304049
fisher.test(data.frame(c(78,9539-78),c(75,6049-75))) # 5'UTR
#p-value = 0.01212
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4712914 0.9155842
#sample estimates:
#  odds ratio 
#0.6567312


## Break the above down by age
elementNROWS(lapply(unique[5:12], function(x) unique(x$editingID)))
#ANnotAC ACnotAN ANnotPN PNnotAN ACnotPC PCnotAC PCnotPN PNnotPC 
#6732    2469    8044    4688    4558    3972    2335    3828
x = lapply(unique, function(x) x[,length(unique(editingID)),by="annotation"])
x = lapply(x, function(x) as.data.frame(x))
x = lapply(x, function(x) x[order(x$annotation),])
y = data.frame(x[[5]][,1], x[[5]][,2], x[[6]][,2], x[[7]][,2], x[[8]][,2],x[[9]][,2],x[[10]][,2],x[[11]][,2],x[[12]][,2])
colnames(y) = c("annotation", names(x)[5:12])
#annotation ANnotAC ACnotAN ANnotPN PNnotAN ACnotPC PCnotAC PCnotPN PNnotPC
#1        CDS      84      35     150     114     107     100      30      50
#2     Intron    3068     690    3214    2070    1244    1466     834    1846
#3      Other    1200     507    1460     868     834     721     438     652
#4       UTR3    2323    1218    3155    1577    2345    1644    1011    1239
#5       UTR5      57      19      65      59      28      41      22      41

# In Adult: Are  more likely to be in a cyt over nuc or vice versa?
fisher.test(data.frame(c(84,6732-84),c(35,2469-35))) # CDS
#p-value = 0.5325
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5840017 1.3472645
#sample estimates:
#  odds ratio 
#0.8787258
fisher.test(data.frame(c(3068,6732-3068),c(690,2469-690))) # Intron
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.951227 2.390122
#sample estimates:
#  odds ratio 
#2.158708
fisher.test(data.frame(c(1200,6732-1200),c(507,2469-507))) # Other
#p-value = 0.003324
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7466130 0.9446533
#sample estimates:
#  odds ratio 
#0.8394753
fisher.test(data.frame(c(2323,6732-2323),c(1218,2469-1218))) # 3'UTR
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4922896 0.5949168
#sample estimates:
#  odds ratio 
#0.5411784
fisher.test(data.frame(c(57,6732-57),c(19,2469-19))) # 5'UTR
#p-value = 0.7956
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6437852 1.9643228
#sample estimates:
#  odds ratio 
#1.101111

# In Prenatal: Are  more likely to be in a cyt over nuc or vice versa?
fisher.test(data.frame(c(30,2335-30),c(50,3828-50))) # CDS
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.601889 1.582234
#sample estimates:
#  odds ratio 
#0.9834295
fisher.test(data.frame(c(834,2335-834),c(1846,3828-1846))) # Intron
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5358590 0.6640308
#sample estimates:
#  odds ratio 
#0.5966263
fisher.test(data.frame(c(438,2335-438),c(652,3828-652))) # Other
#p-value = 0.08547
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9812938 1.2882884
#sample estimates:
#  odds ratio 
#1.124683 
fisher.test(data.frame(c(1011,2335-1011),c(1239,3828-1239))) # 3'UTR
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.432615 1.776922
#sample estimates:
#  odds ratio 
#1.595525
fisher.test(data.frame(c(22,2335-22),c(41,3828-41))) # 5'UTR
#p-value = 0.6962
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4970463 1.5141666
#sample estimates:
#  odds ratio 
#0.8785313 

# In Cytosol: Are  more likely to be in a adult over prenatal or vice versa?
fisher.test(data.frame(c(107,4558-107),c(100,3972-100))) # CDS
#p-value = 0.6221
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6995324 1.2394529
#sample estimates:
#  odds ratio 
#0.9308348
fisher.test(data.frame(c(1244,4558-1244),c(1466,3972-1466))) # Intron
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5848471 0.7040080
#sample estimates:
#  odds ratio 
#0.6417252
fisher.test(data.frame(c(834,4558-834),c(721,3972-721))) # Other
#p-value = 0.8661
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9030578 1.1293151
#sample estimates:
#  odds ratio 
#1.009806
fisher.test(data.frame(c(2345,4558-2345),c(1644,3972-1644))) # 3'UTR
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.375856 1.636514
#sample estimates:
#  odds ratio 
#1.500416
fisher.test(data.frame(c(28,4558-28),c(41,3972-41))) # 5'UTR
#p-value = 0.03877
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3522574 0.9838893
#sample estimates:
#  odds ratio 
#0.5926667
     
# In Nucleus: Are  more likely to be in a adult over prenatal or vice versa?
fisher.test(data.frame(c(150,8044-150),c(114,4688-114))) # CDS
#p-value = 0.03326
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5919377 0.9841713
#sample estimates:
#  odds ratio 
#0.7624387
fisher.test(data.frame(c(3214,8044-3214),c(2070,4688-2070))) # Intron
#p-value = 3.748e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7818745 0.9058926
#sample estimates:
#  odds ratio 
#0.8416017
fisher.test(data.frame(c(1460,8044-1460),c(868,4688-868))) # Other
#p-value = 0.6177
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8884259 1.0724153
#sample estimates:
#  odds ratio 
#0.9758917 
fisher.test(data.frame(c(3155,8044-3155),c(1577,4688-1577))) # 3'UTR
#p-value = 3.072e-10
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.179875 1.373741
#sample estimates:
#  odds ratio 
#1.273027
fisher.test(data.frame(c(65,8044-65),c(59,4688-59))) # 5'UTR
#p-value = 0.01478
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4414743 0.9270061
#sample estimates:
#  odds ratio 
#0.639198


## Gene ontology of different groups of editing sites

entrez.editing = lapply(unique, function(x) as.character(na.omit(x$EntrezID)))

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

cytosolOnly = lapply(cytosolOnly, function(x) as.character(na.omit(x$EntrezID)))
nucleusOnly = lapply(nucleusOnly, function(x) as.character(na.omit(x$EntrezID)))
adultOnly = lapply(adultOnly, function(x) as.character(na.omit(x$EntrezID)))
prenatalOnly = lapply(prenatalOnly, function(x) as.character(na.omit(x$EntrezID)))
ANnotAC = lapply(ANnotAC, function(x) as.character(na.omit(x$EntrezID)))
ACnotAN = lapply(ACnotAN, function(x) as.character(na.omit(x$EntrezID)))
ANnotPN = lapply(ANnotPN, function(x) as.character(na.omit(x$EntrezID)))
PNnotAN = lapply(PNnotAN, function(x) as.character(na.omit(x$EntrezID))) 
ACnotPC = lapply(ACnotPC, function(x) as.character(na.omit(x$EntrezID)))
PCnotAC = lapply(PCnotAC, function(x) as.character(na.omit(x$EntrezID)))
PCnotPN = lapply(PCnotPN, function(x) as.character(na.omit(x$EntrezID)))
PNnotPC = lapply(PNnotPC, function(x) as.character(na.omit(x$EntrezID)))

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
pdf(file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/interaction.kegg.GO.DO_unsplit_by_annotation.pdf", width=18, height=10)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO,colorBy="p.adjust",  showCategory = 45, title= "Disease Ontology Enrichment")
dev.off()

## Compare the enriched terms between the split groups
# KEGG
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Kegg_split_by_annotation.pdf")
Kegg.cytosolOnly = compareCluster(cytosolOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: cytosolOnly")
Kegg.nucleusOnly = compareCluster(nucleusOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.adultOnly = compareCluster(adultOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: adultOnly")
Kegg.prenatalOnly = compareCluster(prenatalOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: prenatalOnly")
Kegg.ANnotAC = compareCluster(ANnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotAC")
Kegg.ACnotAN = compareCluster(ACnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotAN")
Kegg.ANnotPN = compareCluster(ANnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotPN")
Kegg.PNnotAN = compareCluster(PNnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotAN")
Kegg.ACnotPC = compareCluster(ACnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotPC")
Kegg.PCnotAC = compareCluster(PCnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotAC")
Kegg.PCnotPN = compareCluster(PCnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotPN")
Kegg.PNnotPC = compareCluster(PNnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotPC")
dev.off()
save(Kegg.cytosolOnly,Kegg.adultOnly,Kegg.prenatalOnly,Kegg.ANnotAC,Kegg.ACnotAN,
     Kegg.ANnotPN,Kegg.PNnotAN,Kegg.ACnotPC,Kegg.PCnotAC,Kegg.PCnotPN,Kegg.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/kegg.objects.RNAediting.SplitByAnnotation.rda")
# Biological Process
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/BP_split_by_annotation.pdf")
BP.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: cytosolOnly")
BP.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: nucleusOnly")
BP.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: prenatalOnly")
BP.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotAC")
BP.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ACnotAN")
BP.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotAC")
BP.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotPN")
BP.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
dev.off()
save(BP.cytosolOnly,BP.nucleusOnly,BP.prenatalOnly,BP.ANnotAC,BP.ACnotAN,BP.PCnotAC,BP.PCnotPN, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/BP.objects.RNAediting.SplitByAnnotation.rda")
# Molecular Function
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/MF_split_by_annotation.pdf")
MF.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: cytosolOnly")
MF.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: nucleusOnly")
MF.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: adultOnly")
MF.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotAC")
MF.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotAN")
MF.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotPN")
MF.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotPC")
MF.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PCnotPN")
MF.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PNnotPC")
dev.off()
save(MF.cytosolOnly,MF.nucleusOnly,MF.adultOnly,MF.ANnotAC,MF.ACnotAN,MF.ANnotPN,MF.ACnotPC,MF.PCnotPN,MF.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/MF.objects.RNAediting.SplitByAnnotation.rda")
# Cellular Component
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/CC_split_by_annotation.pdf")
CC.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.cytosolOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: cytosolOnly")
CC.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: nucleusOnly")
CC.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: adultOnly")
CC.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ANnotAC")
CC.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotAN")
CC.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ANnotPN")
CC.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotPC")
CC.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PCnotAC")
CC.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PCnotPN")
CC.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotPC")
dev.off()
save(CC.cytosolOnly,CC.nucleusOnly,CC.adultOnly,CC.ANnotAC,CC.ACnotAN,CC.ANnotPN,CC.ACnotPC,CC.PCnotAC,CC.PCnotPN,CC.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/CC.objects.RNAediting.SplitByAnnotation.rda")
# Disease Ontology
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/DO_split_by_annotation.pdf")
DO.cytosolOnly = compareCluster(cytosolOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.nucleusOnly = compareCluster(nucleusOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.nucleusOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: nucleusOnly")
DO.adultOnly = compareCluster(adultOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: adultOnly")
DO.prenatalOnly = compareCluster(prenatalOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: prenatalOnly")
DO.ANnotAC = compareCluster(ANnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ACnotAN = compareCluster(ACnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ACnotAN")
DO.ANnotPN = compareCluster(ANnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotPN")
DO.PNnotAN = compareCluster(PNnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PNnotAN")
DO.ACnotPC = compareCluster(ACnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ACnotPC")
DO.PCnotAC = compareCluster(PCnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotAC")
DO.PCnotPN = compareCluster(PCnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotPN")
DO.PNnotPC = compareCluster(PNnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
dev.off()
save(DO.nucleusOnly,DO.adultOnly,DO.prenatalOnly,DO.ACnotAN,DO.ANnotPN,DO.PNnotAN,DO.ACnotPC,DO.PCnotAC,DO.PCnotPN, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/DO.objects.RNAediting.SplitByAnnotation.rda")


### Characterize the overlap with retained introns

# read in results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
nonconst = list()
for (i in 1:length(comps)){nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")}
names(nonconst) = comps
elementNROWS(nonconst)
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2043              1582              1332              2118               211               338
string = lapply(nonconst, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE))
genes = lapply(string, function(x) x[grep("ENSG", x)])
comments = lapply(string, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)]))
IR.diff = lapply(nonconst, function(y) y$A.IRratio - y$B.IRratio)
Sign = lapply(IR.diff, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult"))
IR = Map(cbind, nonconst, ensID = genes, comments = comments, IR.diff = IR.diff, Sign = Sign)
IRclean = lapply(IR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" & 
                                         y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" & y$comments=="clean"),])
lapply(IRclean, head)
IRsig = lapply(IRclean, function(x) x[which(x$p.diff<=0.05),])
elementNROWS(IRsig)
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#162                93                80               179                32                37
# Find overlaps with RNA editing sites
IRranges = lapply(IRclean, function(x) makeGRangesFromDataFrame(x, start.field="Start",end.field="End",strand.field="Direction",keep.extra.columns = T))
editing_ranges = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
IR_hits = lapply(IRranges, function(x) findOverlaps(editing_ranges, x)) # No overlaps!


## Characterize the overlap with differentially expressed genes

DEG_hits = findOverlaps(geneMapGR, reduce(editing_ranges))
length(reduce(editing_ranges)) #23354
length(unique(subjectHits(DEG_hits))) # 20157

elementNROWS(sig)
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad 
#975           1010            354           3427            350           2442             16              4 
#interacting 
#1194
editing_anno = as.data.frame(editing_anno)
DEG_editing = lapply(sig, function(x) editing_anno[which(editing_anno$collapsedconversion=="A:G / T:C" & 
                                                           editing_anno$nearestID %in% as.character(x$geneID)),])
elementNROWS(DEG_editing)
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad 
#2283           2582            405           8301            332           2319            154              0 
#interacting 
#2746
elementNROWS(lapply(DEG_editing, function(x) unique(x$editingID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad 
#1141            994            234           4021            146           1024             51              0 
#interacting 
#1371
elementNROWS(lapply(DEG_editing, function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad 
#209            198             62            726             28            219              5              0 
#interacting 
#272 
# in DEG with a LFC greater than 1
DEG_editing.1 = lapply(sig.1, function(x) editing_anno[which(editing_anno$collapsedconversion=="A:G / T:C" &
                                                               editing_anno$nearestID %in% as.character(x$geneID)),])
elementNROWS(DEG_editing.1)
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#786            435             27           1464              0            206              0              0            777
elementNROWS(lapply(DEG_editing.1, function(x) unique(x$editingID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#407            202             21            652              0             99              0              0            391
elementNROWS(lapply(DEG_editing.1, function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported ret_Ad_exp_Fet ret_Fet_exp_Ad    interacting 
#92             48              3            131              0             22              0              0             94


## Are cytosolic-specific editing sites enriched for DEG Fraction?
names(unique)
cyt = unique[["cytosolOnly"]]
cytonly.deg = lapply(sig, function(x) cyt[which(cyt$collapsedconversion=="A:G / T:C" & cyt$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cytonly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.056533 2.167616
#sample estimates:
#  odds ratio 
#1.508694 
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.206004 1.761161
#sample estimates:
#  odds ratio 
#1.455631
fisher.test(data.frame(c(58+16,975+354), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.01582
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4920598 0.9349386
#sample estimates:
#  odds ratio 
#0.6796239

## Are nuclear-specific editing sites enriched for DEG Fraction?
nuc = unique[["nucleusOnly"]]
nuconly.deg = lapply(sig, function(x) nuc[which(nuc$collapsedconversion=="A:G / T:C" & nuc$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuconly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#167           112            47           546            19           136

fisher.test(data.frame(c(167,975-167), c(112,1010-112))) # both ages: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value = 0.0001323
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.271346 2.164859
#sample estimates:
#  odds ratio 
#1.656712
fisher.test(data.frame(c(167+546,975+3427-(167+546)), c(112+136,1010+2442-(112+136)))) # adult: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.140855 2.919265
#sample estimates:
#  odds ratio 
#2.496752
fisher.test(data.frame(c(167+47,975+354-(167+47)), c(112+19,1010+350-(112+19)))) # prenatal: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value = 6.311e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.419394 2.289464
#sample estimates:
#  odds ratio 
#1.800193

## Are adult-specific editing sites enriched for DEG Fraction?
ad = unique[["adultOnly"]]
adonly.deg = lapply(sig, function(x) ad[which(ad$collapsedconversion=="A:G / T:C" & ad$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(adonly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#154           127            44           509            18           157 

fisher.test(data.frame(c(154,975-154), c(127,1010-127))) # both ages: retained or exported DEG and presence or absence of adult-specific editing site
#p-value = 0.04576
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.004371 1.695185
#sample estimates:
#  odds ratio 
#1.303989
fisher.test(data.frame(c(154+509,975+3427-(154+509)), c(127+157,1010+2442-(127+157)))) # adult: retained or exported DEG and presence or absence of adult-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.704538 2.299308
#sample estimates:
#  odds ratio 
#1.977819 
fisher.test(data.frame(c(154+44,975+354-(154+44)), c(127+18,1010+350-(127+18)))) # prenatal: retained or exported DEG and presence or absence of adult-specific editing site
#p-value = 0.001182
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.159608 1.858403
#sample estimates:
#  odds ratio 
#1.466712

## Are prenatal-specific editing sites enriched for DEG Fraction?
pren = unique[["prenatalOnly"]]
prenonly.deg = lapply(sig, function(x) pren[which(pren$collapsedconversion=="A:G / T:C" & pren$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(prenonly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#104           105            24           380            15           106 

fisher.test(data.frame(c(104,975-104), c(105,1010-105))) # both ages: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 0.8838
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7643917 1.3854113
#sample estimates:
#  odds ratio 
#1.029133
fisher.test(data.frame(c(104+380,975+3427-(104+380)), c(105+106,1010+2442-(105+106)))) # prenatal: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 2.018e-14
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.599683 2.256263
#sample estimates:
#  odds ratio 
#1.897326
fisher.test(data.frame(c(104+24,975+354-(104+24)), c(105+15,1010+350-(105+15)))) # prenatal: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 0.5052
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8406777 1.4434121
#sample estimates:
#  odds ratio 
#1.101264
 














