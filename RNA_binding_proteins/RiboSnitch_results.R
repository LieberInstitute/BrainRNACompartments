library(GenomicRanges)
library(data.table)
library(ggplot2)


load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")

## Scan in master list of editing sites

master = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/rna_editing/data/BEDfiles/all_editingSites_includingstrand.bed")
colnames(master) = c("chromosome","start","end","siteName","space","strand")
master_gr = makeGRangesFromDataFrame(master, keep.extra.columns = T)


## Scan in RiboSnitch results 

rbsnitch = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/Transcript_RBP_SNP.filtered.strand.energy.p.txt",sep = "\t")
colnames(rbsnitch) = c("Transcript_chrom","Transcript_type","Transcript_start","Transcript_end","Transcript_strand","Transcript_information",
                       "Transcript_wildtype_sequence","Transcript_mutant_sequence","Transcript_wildtype_energy","Transcript_mutant_energy",
                       "Transcript_change_energy","EditingSite_mutant_type","EditingSite_chrom","EditingSite_start","EditingSite_end","EditingSite_ID","pval")

dim(rbsnitch) # 14624 transcripts affected
length(unique(rbsnitch$EditingSite_ID)) # 8537 editing sites associated with changes in energy

rbsnitch = cbind(rbsnitch, padj = p.adjust(rbsnitch$pval, method = "fdr")) ## how big should N be? wait for Fengbiao's response
nrow(rbsnitch[which(rbsnitch$pval<=0.05),]) # 546
nrow(rbsnitch[which(rbsnitch$padj<=0.05),]) # 0
hist(rbsnitch$pval)


## Match the results to editing site annotation

editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
rbsnitch_gr = makeGRangesFromDataFrame(rbsnitch,seqnames.field="EditingSite_chrom", start.field="EditingSite_start", end.field="EditingSite_end",
                                       strand.field="Transcript_strand", keep.extra.columns = T)
ov = findOverlaps(editing_anno_gr, rbsnitch_gr)
rbsnitch = cbind(rbsnitch[subjectHits(ov),-grep("sequence", colnames(rbsnitch))],editing_anno[queryHits(ov),])
head(rbsnitch)
hist(rbsnitch$Transcript_change_energy)


## Identify unique editing sites present in all samples in a group

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

unique_rbsnitch = lapply(unique_all, function(x) rbsnitch[rbsnitch$editingID %in% x$editingID,])
round(elementNROWS(lapply(unique_rbsnitch, function(x) unique(x$editingID)))/
        elementNROWS(lapply(unique_all, function(x) unique(x$editingID)))*100,2)
#   adultOnly prenatalOnly      ANnotAC      ACnotAN      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN 
#       65.38        60.67        54.72        63.33        58.33        33.85        61.17        50.57        60.11        48.95 
