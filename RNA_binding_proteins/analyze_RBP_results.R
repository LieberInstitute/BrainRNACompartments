library(GenomicRanges)
library(data.table)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


## read in results

rbpvar = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/RBP-Var.score.input.p.txt", header = F, sep = "\t")
head(rbpvar)
dim(rbpvar)
rbpvar = rbpvar[,!colnames(rbpvar) %in% c("V6","V7","V8","V10","V13","V15","V16","V17","V19","V21","V23","V25")]
head(rbpvar)
colnames(rbpvar) = c("chromosome", "start", "end", "Original_Allele", "Edited_Allele","SNV", "Conversion","Strand","SiteName",
                     "RBPScore", "CLIP","miRNA", "Motif","Energy.Change","Energy.pval")
unique(rbpvar$RBPScore)
nrow(rbpvar[rbpvar$RBPScore %in% c("04","05","06"),]) # 1578

clip = lapply(strsplit(as.character(rbpvar$CLIP), ";", fixed = T), unique)
motif = lapply(strsplit(as.character(rbpvar$Motif), ";", fixed = T), unique)
names(clip) = names(motif) = paste0(rbpvar$SiteName, ".")
clip = data.frame(clip = unlist(clip), SiteName = gsub("\\..*","", names(unlist(clip))))
motif = data.frame(motif = unlist(motif), SiteName = gsub("\\..*","", names(unlist(motif))))
rbpclip = cbind(rbpvar[match(clip$SiteName, rbpvar$SiteName),!colnames(rbpvar) %in% c("CLIP", "Motif")], clip = clip$clip)
rbpmotif = cbind(rbpvar[match(motif$SiteName, rbpvar$SiteName),!colnames(rbpvar) %in% c("CLIP", "Motif")], motif = motif$motif)


## Annotate results

editing_anno$start = editing_anno$end
editing_anno_gr = makeGRangesFromDataFrame(editing_anno)
rbpclip_gr = makeGRangesFromDataFrame(rbpclip)
rbpmotif_gr = makeGRangesFromDataFrame(rbpmotif)
ovclip = findOverlaps(editing_anno_gr, rbpclip_gr)
ovmotif = findOverlaps(editing_anno_gr, rbpmotif_gr)

rbpclip = data.table(cbind(rbpclip[subjectHits(ovclip),], editing_anno[queryHits(ovclip),]))
rbpmotif = data.table(cbind(rbpmotif[subjectHits(ovmotif),], editing_anno[queryHits(ovmotif),]))
rbpclip$annotation = factor(rbpclip$annotation, levels = c("CDS", "3'UTR", "5'UTR", "Intron", "Intergenic"))
rbpmotif$annotation = factor(rbpmotif$annotation, levels = c("CDS", "3'UTR", "5'UTR", "Intron", "Intergenic"))

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/reported_RBPVar_editingSites_Annotation.pdf", width = 8.5, height = 6)
df = rbpclip[,length(unique(editingID)), by = "annotation"]
df$percent = round(df$V1/sum(df$V1)*100,2)
ggplot(df, aes(x = annotation, y = V1)) + geom_bar(stat = "identity") +
  geom_text(aes(label=paste0(percent, "%")), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("Annotation of Reported Editing Sites by RBPVar2") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

### Isolate the sites present in all samples in each group

unique_bySamp_all = list()
comp = list(c("cytosolOnly","cytosolAll"), c("nucleusOnly","nucleusAll"), c("adultOnly","adultAll"), c("prenatalOnly","prenatalAll"), 
            c("ANnotAC","allAN"), c("ACnotAN","allAC"), c("PCnotPN","allPC"), c("PNnotPC","allPN"),
            c("ACnotPC","allAC"), c("PCnotAC","allPC"), c("ANnotPN","allAN"), c("PNnotAN","allPN"))
for (i in 1:length(comp)) {
  unique_bySamp_all[[i]] = unique_bySamp[[comp[[i]][1]]][-grep("no", unique_bySamp[[comp[[i]][1]]][,comp[[i]][2]]),]
}
names(unique_bySamp_all) = lapply(comp, function(x) x[1])
elementNROWS(unique_bySamp_all)

unique_all = lapply(unique_bySamp_all[3:12], function(x) editing_anno[which(editing_anno$editingID %in% x$EditingID),,])
unique_all = Map(cbind, unique_all, 
                 geneID = lapply(unique_all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
unique_all = Map(cbind, unique_all, ensID = lapply(unique_all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique_all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]))
elementNROWS(unique_all)
lapply(unique_all, head)

all = Map(cbind, all, geneID = lapply(all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]))

unique_clip = lapply(unique_all, function(x) rbpclip[rbpclip$editingID %in% x$editingID,])
round(elementNROWS(lapply(unique_clip, function(x) unique(x$editingID)))/
        elementNROWS(lapply(unique_all, function(x) unique(x$editingID)))*100,2)
# adultOnly prenatalOnly      ANnotAC      ACnotAN      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN 
#      3.85         3.37         5.66        10.00        25.00         3.08         4.40         4.60         4.87         5.59 
unique_motif = lapply(unique_all, function(x) rbpmotif[rbpmotif$editingID %in% x$editingID,])
round(elementNROWS(lapply(unique_motif, function(x) unique(x$editingID)))/
        elementNROWS(lapply(unique_all, function(x) unique(x$editingID)))*100,2)
# adultOnly prenatalOnly      ANnotAC      ACnotAN      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN 
#      3.85         3.37         5.66        10.00        25.00         3.08         4.40         4.60         4.87         5.59 



## look into CLIP binding by annotation

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/CLIP_RBPs_byAnnotation.pdf", width = 10, height = 6)

df = as.data.frame(rbpclip[clip!=".",length(unique(editingID)), by=c("clip", "annotation")])
for (i in 1:length(unique(df$annotation))) { 
  df[df$annotation==unique(df$annotation)[i],"total"] = sum(df[df$annotation==unique(df$annotation)[i],"V1"]) }
df$percent = round(df$V1/df$total*100,2)
df[order(df$annotation),]

ggplot(df, aes(x = annotation, y = V1, fill = clip)) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

merge(aggregate(V1 ~ annotation, df, max), df) # PTBP1 is the most frequent RBP found to bind to these sites
merge(aggregate(V1 ~ annotation, df, function(x) max( x[x!=max(x)] )), df) # TARBP2, AGO2, and eIF4AIII are second most abundant
df[order(df$percent, decreasing =T),]

ggplot(df[df$clip=="PTBP1",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation: PTBP1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$clip=="TARBP2",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation: TARBP2") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$clip=="eIF4AIII",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation: eIF4AIII") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

total = as.data.frame(editing_anno[,length(unique(editingID)), by="annotation"])
tables = list()
for (i in 1:length(unique(df$annotation))) {
  tables[[i]] = data.frame(PTBP1 = c(df[df$clip=="PTBP1" & df$annotation==unique(df$annotation)[i],"V1"],
                                     sum(df[df$clip=="PTBP1" & df$annotation!=unique(df$annotation)[i],"V1"])), 
                           noPTBP1 = c(total[total$annotation==unique(df$annotation)[i],"V1"]-
                                         df[df$clip=="PTBP1" & df$annotation==unique(df$annotation)[i],"V1"],
                                       sum(total[total$annotation!=unique(df$annotation)[i],"V1"])-
                                         sum(df[df$clip=="PTBP1" & df$annotation!=unique(df$annotation)[i],"V1"])), row.names = c("inAnno","notAnno"))
}
names(tables) = unique(df$annotation)
fisher = cbind(pval = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$p.value), recursive=F)),
               OR = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$estimate), recursive=F)))
colnames(fisher) = c("p.value", "odds.ratio")
write.csv(fisher, quote=F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/fisher_PTBP1_byAnnotation.csv")



## look into CLIP binding by unique location or age

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/CLIP_RBPs_byUniqueGroup.pdf", width = 17, height = 7)

df = do.call(rbind, Map(cbind, unique_clip, Comparison = as.list(names(unique_clip))))
df$Comparison = factor(df$Comparison, levels = c("adultOnly","prenatalOnly","ANnotAC","ACnotAN","PCnotPN","PNnotPC",
                                                 "ACnotPC","PCnotAC","ANnotPN","PNnotAN"))
df = as.data.frame(df[clip!=".",length(unique(editingID)), by=c("clip", "annotation", "Comparison")])
for (i in 1:length(unique(df$Comparison))) { 
  df[df$Comparison==unique(df$Comparison)[i],"total"] = sum(df[df$Comparison==unique(df$Comparison)[i],"V1"]) }
df$percent = round(df$V1/df$total*100,2)
df[order(df$Comparison),]

ggplot(df, aes(x = annotation, y = V1, fill = clip)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation in Unique Groups") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

merge(aggregate(V1 ~ Comparison, df, max), df) # PTBP1 is the most frequent RBP found to bind to these sites
df[order(df$percent, decreasing =T),]

ggplot(df[df$clip=="PTBP1",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation in Unique Groups: PTBP1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$clip=="TARBP2",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation in Unique Groups: TARBP2") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$clip=="eIF4AIII",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,50) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation in Unique Groups: eIF4AIII") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

df = do.call(rbind, Map(cbind, unique_clip, Comparison = as.list(names(unique_clip))))
df$Comparison = factor(df$Comparison, levels = c("adultOnly","prenatalOnly","ANnotAC","ACnotAN","PCnotPN","PNnotPC",
                                                 "ACnotPC","PCnotAC","ANnotPN","PNnotAN"))
df = as.data.frame(df[clip!=".",length(unique(editingID)), by=c("clip", "Comparison")])
totalunique = do.call(rbind, Map(cbind, unique_all, Comparison = as.list(names(unique_all))))
total = as.data.frame(totalunique[,length(unique(editingID)), by="Comparison"])
for (i in 1:length(unique(df$Comparison))) {
  if (nrow(df[df$clip=="PTBP1" & df$Comparison==unique(df$Comparison)[i],])<1) {
    df = rbind(df, data.frame(clip = "PTBP1", Comparison = unique(df$Comparison)[i], V1 = 0))
  }
  tables[[i]] = data.frame(PTBP1 = c(df[df$clip=="PTBP1" & df$Comparison==unique(df$Comparison)[i],"V1"],
                                     length(unique(as.data.frame(rbpclip)[rbpclip$clip=="PTBP1","editingID"]))-
                                       df[df$clip=="PTBP1" & df$Comparison==unique(df$Comparison)[i],"V1"]), 
                           noPTBP1 = c(total[total$Comparison==unique(df$Comparison)[i],"V1"]-
                                         df[df$clip=="PTBP1" & df$Comparison==unique(df$Comparison)[i],"V1"],
                                       length(unique(editing_anno$editingID))-total[total$Comparison==unique(df$Comparison)[i],"V1"]-
                                         (length(unique(as.data.frame(rbpclip)[rbpclip$clip=="PTBP1","editingID"]))-
                                            df[df$clip=="PTBP1" & df$Comparison==unique(df$Comparison)[i],"V1"])), row.names = c("inGroup","notGroup"))
}
names(tables) = unique(df$Comparison)
fisher = cbind(pval = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$p.value), recursive=F)),
               OR = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$estimate), recursive=F)))
colnames(fisher) = c("p.value", "odds.ratio")
write.csv(fisher, quote=F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/fisher_PTBP1_byUniqueGroup.csv")



## look into RBP binding motifs by annotation

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/motif_RBPs_byAnnotation.pdf", width = 17, height = 6)

df = as.data.frame(rbpmotif[motif!=".",length(unique(editingID)), by=c("motif", "annotation")])
for (i in 1:length(unique(df$annotation))) { 
  df[df$annotation==unique(df$annotation)[i],"total"] = sum(df[df$annotation==unique(df$annotation)[i],"V1"]) }
df$percent = round(df$V1/df$total*100,2)
df[order(df$annotation),]

ggplot(df, aes(x = annotation, y = V1, fill = motif)) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

merge(aggregate(V1 ~ annotation, df, max), df) # SNRPA is the most frequent RBP found to bind to intronic, intergenic, and 5'UTR sites
                                               # HNRNPH1, HNRNPF, and HNRNPH2 are all abundant in 3'UTR  
merge(aggregate(V1 ~ annotation, df, function(x) max( x[x!=max(x)] )), df) # ELAVL1, ELAVL2, and ELAVL3, are second most abundant in introns, 5'UTR and CDS
df[order(df$percent, decreasing =T),]

ggplot(df[df$motif=="SNRPA",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,20) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation: SNRPA") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$motif=="HNRNPH1",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,20) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation: HNRNPH1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$motif=="ELAVL1",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,20) + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation: ELAVL1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

total = as.data.frame(editing_anno[,length(unique(editingID)), by="annotation"])
tables = list()
for (i in 1:length(unique(df$annotation))) {
  if (nrow(df[df$motif=="SNRPA" & df$annotation==unique(df$annotation)[i],])<1) {
    df = rbind(df, data.frame(motif = "SNRPA", annotation = unique(df$annotation)[i], V1 = 0, total =0, percent = 0))
  }
  tables[[i]] = data.frame(SNRPA = c(df[df$motif=="SNRPA" & df$annotation==unique(df$annotation)[i],"V1"],
                                     sum(df[df$motif=="SNRPA" & df$annotation!=unique(df$annotation)[i],"V1"])), 
                           noSNRPA = c(total[total$annotation==unique(df$annotation)[i],"V1"]-
                                         df[df$motif=="SNRPA" & df$annotation==unique(df$annotation)[i],"V1"],
                                       sum(total[total$annotation!=unique(df$annotation)[i],"V1"])-
                                         sum(df[df$motif=="SNRPA" & df$annotation!=unique(df$annotation)[i],"V1"])), row.names = c("inAnno","notAnno"))
}
names(tables) = unique(df$annotation)
fisher = cbind(pval = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$p.value), recursive=F)),
               OR = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$estimate), recursive=F)))
colnames(fisher) = c("p.value", "odds.ratio")
write.csv(fisher, quote=F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/fisher_SNRPA_byAnnotation.csv")
fisher[which(fisher$p.value<=(0.05/5)),]
## SNRPA is enriched in introns and depleted in CDS sequence


## look into motif binding by unique location or age

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/motif_RBPs_byUniqueGroup.pdf", width = 17, height = 7)

df = do.call(rbind, Map(cbind, unique_motif, Comparison = as.list(names(unique_motif))))
df$Comparison = factor(df$Comparison, levels = c("adultOnly","prenatalOnly","ANnotAC","ACnotAN","PCnotPN","PNnotPC",
                                                 "ACnotPC","PCnotAC","ANnotPN","PNnotAN"))
df = as.data.frame(df[motif!=".",length(unique(editingID)), by=c("motif", "annotation", "Comparison")])
for (i in 1:length(unique(df$Comparison))) { 
  df[df$Comparison==unique(df$Comparison)[i],"total"] = sum(df[df$Comparison==unique(df$Comparison)[i],"V1"]) }
df$percent = round(df$V1/df$total*100,2)
df[order(df$Comparison),]

ggplot(df, aes(x = annotation, y = V1, fill = motif)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation in Unique Groups") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

merge(aggregate(V1 ~ Comparison, df, max), df) # SNRPA is the most frequent RBP found to bind to these sites
df[order(df$percent, decreasing =T),]

ggplot(df[df$motif=="SNRPA",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation in Unique Groups: SNRPA") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$motif=="SNRPB2",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation in Unique Groups: SNRPB2") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df[df$motif=="HNRNPH1",], aes(x = annotation, y = percent)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Comparison, scales = "free") +
  geom_text(aes(label=V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") +  
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("motif RBPs By Annotation in Unique Groups: HNRNPH1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

df = do.call(rbind, Map(cbind, unique_motif, Comparison = as.list(names(unique_motif))))
df$Comparison = factor(df$Comparison, levels = c("adultOnly","prenatalOnly","ANnotAC","ACnotAN","PCnotPN","PNnotPC",
                                                 "ACnotPC","PCnotAC","ANnotPN","PNnotAN"))
df = as.data.frame(df[motif!=".",length(unique(editingID)), by=c("motif", "Comparison")])
totalunique = do.call(rbind, Map(cbind, unique_all, Comparison = as.list(names(unique_all))))
total = as.data.frame(totalunique[,length(unique(editingID)), by="Comparison"])
for (i in 1:length(unique(df$Comparison))) {
  if (nrow(df[df$motif=="SNRPA" & df$Comparison==unique(df$Comparison)[i],])<1) {
    df = rbind(df, data.frame(motif = "SNRPA", Comparison = unique(df$Comparison)[i], V1 = 0))
  }
  tables[[i]] = data.frame(SNRPA = c(df[df$motif=="SNRPA" & df$Comparison==unique(df$Comparison)[i],"V1"],
                                     length(unique(as.data.frame(rbpmotif)[rbpmotif$motif=="SNRPA","editingID"]))-
                                       df[df$motif=="SNRPA" & df$Comparison==unique(df$Comparison)[i],"V1"]), 
                           noSNRPA = c(total[total$Comparison==unique(df$Comparison)[i],"V1"]-
                                         df[df$motif=="SNRPA" & df$Comparison==unique(df$Comparison)[i],"V1"],
                                       length(unique(editing_anno$editingID))-total[total$Comparison==unique(df$Comparison)[i],"V1"]-
                                         (length(unique(as.data.frame(rbpmotif)[rbpmotif$motif=="SNRPA","editingID"]))-
                                            df[df$motif=="SNRPA" & df$Comparison==unique(df$Comparison)[i],"V1"])), row.names = c("inGroup","notGroup"))
}
names(tables) = unique(df$Comparison)
fisher = cbind(pval = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$p.value), recursive=F)),
               OR = data.frame(unlist(lapply(lapply(tables,fisher.test), function(x) x$estimate), recursive=F)))
colnames(fisher) = c("p.value", "odds.ratio")
write.csv(fisher, quote=F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/fisher_SNRPA_byUniqueGroup.csv")



### Read in results table from the RBPVar-2 website

x = read.table("/Users/amanda/Downloads/RNAedit.RBP-Var.score.p.txt", header=T)
gr = makeGRangesFromDataFrame(x, seqnames.field = "chromosome", start.field = "position", end.field = "position", keep.extra.columns = T)
ov = findOverlaps(editing_anno_gr, gr)
RBPVar_anno = cbind(editing_anno[queryHits(ov),], x[subjectHits(ov),])


### Identify the unique to a group editing sites with RBPVar info

RBPVar_unique = lapply(unique_all, function(x) RBPVar_anno[RBPVar_anno$nearestID %in% x$nearestID,])
elementNROWS(RBPVar_unique)

table(RBPVar_anno$RBP.Var_score)
hist(RBPVar_anno$P_value)


