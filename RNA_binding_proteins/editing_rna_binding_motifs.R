library(GenomicRanges)
library(data.table)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_results_editing_sites.rda")


### Are all editing sites reflected here?

colnames(rbpmap) = c("Seq_Pos", "Coordinate", "Motif","Kmer","Zscore","Pval","editingRegion","proteinID")
length(unique(rbpmap$editingRegion)) # 18907 of 18907 coordinates


### FDR correct pvalues for motif enrichment and choose the motif with the lowest pvalue per editing site and protein ID 

rbpmap$Pval = as.numeric(as.character(rbpmap$Pval))
rbpmap$Zscore = as.numeric(as.character(rbpmap$Zscore))
rbpmap$Motif = gsub(" ", "", as.character(rbpmap$Motif))
rbpmap$Coordinate = gsub(" ", "", as.character(rbpmap$Coordinate))
rbpmap$editingRegion = as.character(rbpmap$editingRegion)
rbpmap$motifID = paste0(rbpmap$Coordinate, "-", as.numeric(unlist(lapply(strsplit(rbpmap$Coordinate, ":", fixed=T), function(x) x[2])))+
                               nchar(rbpmap$Motif), ":", unlist(lapply(strsplit(rbpmap$editingRegion, ":", fixed=T), function(x) x[3])))
rbpmap$padj = p.adjust(rbpmap$Pval, method = "fdr", n = 2155398) # here n = 18907 * 114


rbpmap = data.table(rbpmap)
sigrbp = rbpmap[rbpmap[padj<=0.05,.I[padj == min(padj)], by=c("editingRegion","proteinID")]$V1]
length(unique(sigrbp$proteinID)) # 94 RBPs are represented out of 114 motifs tested


### Map the RBP motif results to the corresponding editing sites

sites = data.frame(GRanges(as.character(sigrbp$editingRegion)))
sites$start = sites$start+25
sites$end = sites$end-25
sites = makeGRangesFromDataFrame(sites)
motifsites = GRanges(as.character(sigrbp$motifID))
ov = findOverlaps(motifsites, sites)
sigrbp = sigrbp[queryHits(ov)[which(queryHits(ov)==subjectHits(ov))],] # limit to motifs that directly overlap the editing site
sites = sites[queryHits(ov)[which(queryHits(ov)==subjectHits(ov))],]

editing_anno$start = editing_anno$end
editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns =T)
ov = findOverlaps(sites, editing_anno_gr)
editing_anno = cbind(editing_anno[subjectHits(ov),,], sigrbp[queryHits(ov),,])
length(unique(editing_anno$editingID)) # 6027 editing sites significantly overlap an RBP binding motif



### Isolate the editing sites present in all samples in each group

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



### Plot the distribution of RBPs by group and in total

total = editing_anno[,length(unique(editingID)), by = "proteinID"]
total = total[order(total$V1, decreasing = T),]

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/RBPMap_protein_counts.pdf", width = 10, height = 6)

df = as.data.frame(rbpclip[clip!=".",length(unique(editingID)), by=c("clip", "annotation")])
for (i in 1:length(unique(df$annotation))) { 
  df[df$annotation==unique(df$annotation)[i],"total"] = sum(df[df$annotation==unique(df$annotation)[i],"V1"]) }
df$percent = round(df$V1/df$total*100,2)
df[order(df$annotation),]

ggplot(total, aes(x = annotation, y = V1, fill = clip)) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("CLIP RBPs By Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Map RBP table to differentially retained intron list
lapply(proteins, head)
head(rbp)
head(introns[[4]])
rbpIntrons = lapply(introns, function(x) rbp[which(rbp$intronID %in% x$intronID),])
for (i in 1:length(rbpIntrons)){
  rbp = rbpIntrons[[i]]
  int = introns[[i]]
  rbpIntrons[[i]] = cbind(rbp, int[match(rbp$intronID, int$intronID),])
}

### Do differentially retained introns by age have more RBP sites?


### What is the distribution of RBPs by different groups of introns?
length(unique(rbp$proteinID)) # 95 represented in entire list
proteinTypes_byIR = lapply(rbpIntrons, function(x) count(x, vars = "proteinID"))


### Any RBPs stick out?







# more information can be found at http://rbpmap.technion.ac.il/