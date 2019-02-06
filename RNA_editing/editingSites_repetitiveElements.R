library(GenomicRanges)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


### Check repetitive element overlap

rpmsk = read.table("./Dropbox/sorted_figures/github_controlled/RepeatMasker_genomewide_HG19.txt", header=T)
rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE)
editing_annogr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
hits = findOverlaps(editing_annogr, rpmskgr)
length(unique(subjectHits(hits)))

editing_anno = rbind(cbind(editing_anno[queryHits(hits),], repName = rpmsk$repName[subjectHits(hits)],
                     repClass = rpmsk$repClass[subjectHits(hits)], repFamily = rpmsk$repFamily[subjectHits(hits)]),
                     cbind(editing_anno[-unique(queryHits(hits)),], repName = "NA", repClass = "NA", repFamily = "NA"))


### Isolate the sites present in all samples in each group

unique_bySamp_all = list(cytosolOnly = unique_bySamp$cytosolOnly[-grep("no", unique_bySamp$cytosolOnly$cytosolAll),],
                         nucleusOnly = unique_bySamp$nucleusOnly[-grep("no", unique_bySamp$nucleusOnly$nucleusAll),], 
                         adultOnly = unique_bySamp$adultOnly[-grep("no", unique_bySamp$adultOnly$adultAll),], 
                         prenatalOnly = unique_bySamp$prenatalOnly[-grep("no", unique_bySamp$prenatalOnly$prenatalAll),], 
                         ANnotAC = unique_bySamp$ANnotAC[-grep("no", unique_bySamp$ANnotAC$allAN),], 
                         ACnotAN = unique_bySamp$ACnotAN[-grep("no", unique_bySamp$ACnotAN$allAC),], 
                         ANnotPN = unique_bySamp$ANnotPN[-grep("no", unique_bySamp$ANnotPN$allAN),], 
                         PNnotAN = unique_bySamp$PNnotAN[-grep("no", unique_bySamp$PNnotAN$allPN),],
                         ACnotPC = unique_bySamp$ACnotPC[-grep("no", unique_bySamp$ACnotPC$allAC),], 
                         PCnotAC = unique_bySamp$PCnotAC[-grep("no", unique_bySamp$PCnotAC$allPC),], 
                         PCnotPN = unique_bySamp$PCnotPN[-grep("no", unique_bySamp$PCnotPN$allPC),], 
                         PNnotPC = unique_bySamp$PNnotPC[-grep("no", unique_bySamp$PNnotPC$allPN),])
elementNROWS(unique_bySamp_all)

unique_all = lapply(unique_bySamp_all[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = Map(cbind, unique_all, ensID = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]))
all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]))


### What proportion of editing sites overlap a repetitive element, and what kind?
class = editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age","repClass")]
family = editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age","repFamily")]
family$alu = ifelse(family$repFamily=="Alu","Alu","non-Alu")
family[repFamily=="NA","alu"] = "non-repeat"
editing_anno[repClass != "NA" & collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age")][,list(V1),]/
  editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age")][,list(V1),]
#Fraction  Age    Proportion overlapping a repeat
#Cytosol  Adult  0.4271956
#Nucleus    Adult 0.4468188
#Nucleus Prenatal  0.4441917
#Cytosol Prenatal 0.4416968

editing_anno[repFamily=="Alu" & collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age")][,list(V1),]/
  editing_anno[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age")][,list(V1),]*100
#1: 40.00979
#2: 41.67870
#3: 41.95933
#4: 41.75090

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/numberEditingSites_overlapping_repeats.pdf",width=10,height=6)
ggplot(class, aes(x = Fraction, y = V1, fill = repClass)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Repetitive Element Classes Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(class[repClass != "NA",,], aes(x = Fraction, y = V1, fill = repClass)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Repetitive Element Classes Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(family, aes(x = Fraction, y = V1, fill = repFamily)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Repetitive Element Families Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(family[repFamily != "NA",,], aes(x = Fraction, y = V1, fill = repFamily)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Repetitive Element Families Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

family = as.data.frame(family)
family$alu = factor(family$alu, levels = c("Alu", "non-Alu","non-repeat"))
family$Fraction = gsub("Cytosol","Cytoplasm", family$Fraction)
pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/numberEditingSites_overlapping_alu.pdf",width=5.25,height=3.5)
ggplot(as.data.frame(family), aes(x = Fraction, y = V1, fill = alu)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) + scale_fill_brewer(palette = "Accent") +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Editing Site Distribution\nin Alu Sequence") +
  theme(title = element_text(size = 20), text = element_text(size = 20),
        axis.text.x = element_text(angle = 15, hjust = .5, vjust=.9))
dev.off()


# Unique sites shared in all samples in a group and repeats
class = lapply(unique_all, function(x) x[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age","repClass")])
family = lapply(unique_all, function(x) x[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("Fraction","Age","repFamily")])
for (i in 1:length(family)){
  tmp = family[[i]]
  tmp$alu = ifelse(tmp$repFamily=="Alu","Alu","non-Alu")
  tmp[repFamily=="NA","alu"] = "non-repeat"
  family[[i]] = tmp
}
class = do.call(rbind, Map(cbind, class, group = as.list(names(class))))
class$fracage = paste0(class$Age, ":",class$Fraction)
family = do.call(rbind, Map(cbind, family, group = as.list(names(family))))
family$fracage = paste0(family$Age, ":",family$Fraction)
class$group = factor(class$group, levels= c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))
family$group = factor(family$group, levels= c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))


pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/numberEditingSites_overlapping_repeats_uniqueAll3.pdf",width=18,height=8)
ggplot(class, aes(x = fracage, y = V1, fill = repClass)) + geom_bar(stat = "identity") +
  facet_grid(. ~ group, scales = "free_x") +
  labs(fill="") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Element Classes Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(class[repClass != "NA",,], aes(x = fracage, y = V1, fill = repClass)) + geom_bar(stat = "identity") +
  facet_grid(. ~ group, scales = "free_x") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Element Classes Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(family, aes(x = fracage, y = V1, fill = repFamily)) + geom_bar(stat = "identity") +
  facet_grid(. ~ group, scales = "free_x") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Element Families Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(family[repFamily != "NA",,], aes(x = fracage, y = V1, fill = repFamily)) + geom_bar(stat = "identity") +
  facet_grid(. ~ group, scales = "free_x") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Element Families Overlapping Editing Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(family, aes(x = fracage, y = V1, fill = alu)) + geom_bar(stat = "identity") +
  facet_grid(. ~ group, scales = "free_x") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Distribution of Editing Sites in Alu Sequence") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()
