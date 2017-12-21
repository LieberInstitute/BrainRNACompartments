library(GenomicRanges)
library(data.table)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


### For RNA Binding protein analysis, write coordinates in a way recognized by RBPMap

atog = editing_anno[collapsedconversion=="A:G / T:C",,]
rnacoord = paste0(atog$seqnames, ":", atog$end-25,"-", atog$end+25, ":", atog$strand)
rnacoord = rnacoord[!duplicated(rnacoord)]
length(rnacoord) # 18907 editing sites
int = list(c(1:5000), c(5001:10000),c(10001:15000), c(15001:length(rnacoord))) 
for (i in 1:length(int)) {
  write.table(rnacoord[int[[i]]], quote = F, sep = "\n", row.names = F, col.names = F,
              file=paste0("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/editingSiteCoordinates.RBPmap.format",i,".txt"))
}



### load in RBPMap results

x = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/RBPMap_editingSites_Batch1_results.txt", 
               header = FALSE, sep = "\t", col.names = paste0("V",seq_len(6)), fill = TRUE)
ast = grep("***", x$V1, fixed = T)
sep = list()
sep[[1]] = x[c(1:ast[1]-1),]
for (i in 1:length(ast)) { sep[[i+1]] = x[c((ast[i]+1):(ast[i+1]-1)),] }
sep = Map(cbind, sep[2:length(sep)], editingID = lapply(sep[2:length(sep)], function(x) x$V1[1]))
sep = lapply(sep, function(x) x[3:nrow(x),])
sep = do.call(rbind, sep)
proteinIDs = grep("Protein", sep$V1, fixed = T)
batch1 = list()
for (i in 1:length(proteinIDs)) { batch1[[i]] = sep[c((proteinIDs[i]):(proteinIDs[i+1]-1)),] }
batch1 = Map(cbind, batch1, proteinID = lapply(batch1, function(x) x$V1[1]))
batch1 = do.call(rbind, batch1)
batch1 = batch1[-c(grep("Protein", batch1$V1),grep("Sequence", batch1$V1)),]
batch1$proteinID = gsub("Protein: ","", batch1$proteinID)
batch1$proteinID = gsub("(Hs/Mm)","", batch1$proteinID, fixed = T)
colnames(batch1) = c("Sequence Position","Genomic Coordinate","Motif","K-mer","Z-score", "P-value", "editingID", "proteinID")

x = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/RBPMap_editingSites_Batch2_results.txt", 
               header = FALSE, sep = "\t", col.names = paste0("V",seq_len(6)), fill = TRUE)
ast = grep("***", x$V1, fixed = T)
sep = list()
sep[[1]] = x[c(1:ast[1]-1),]
for (i in 1:length(ast)) { sep[[i+1]] = x[c((ast[i]+1):(ast[i+1]-1)),] }
sep = Map(cbind, sep[2:length(sep)], editingID = lapply(sep[2:length(sep)], function(x) x$V1[1]))
sep = lapply(sep, function(x) x[3:nrow(x),])
sep = do.call(rbind, sep)
proteinIDs = grep("Protein", sep$V1, fixed = T)
batch2 = list()
for (i in 1:length(proteinIDs)) { batch2[[i]] = sep[c((proteinIDs[i]):(proteinIDs[i+1]-1)),] }
batch2 = Map(cbind, batch2, proteinID = lapply(batch2, function(x) x$V1[1]))
batch2 = do.call(rbind, batch2)
batch2 = batch2[-c(grep("Protein", batch2$V1),grep("Sequence", batch2$V1)),]
batch2$proteinID = gsub("Protein: ","", batch2$proteinID)
batch2$proteinID = gsub("(Hs/Mm)","", batch2$proteinID, fixed = T)
colnames(batch2) = c("Sequence Position","Genomic Coordinate","Motif","K-mer","Z-score", "P-value", "editingID", "proteinID")
head(batch2)

x = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/RBPMap_editingSites_Batch3_results.txt", 
               header = FALSE, sep = "\t", col.names = paste0("V",seq_len(6)), fill = TRUE)
ast = grep("***", x$V1, fixed = T)
sep = list()
sep[[1]] = x[c(1:ast[1]-1),]
for (i in 1:length(ast)) { sep[[i+1]] = x[c((ast[i]+1):(ast[i+1]-1)),] }
sep = Map(cbind, sep[2:length(sep)], editingID = lapply(sep[2:length(sep)], function(x) x$V1[1]))
sep = lapply(sep, function(x) x[3:nrow(x),])
sep = do.call(rbind, sep)
proteinIDs = grep("Protein", sep$V1, fixed = T)
batch3 = list()
for (i in 1:length(proteinIDs)) { batch3[[i]] = sep[c((proteinIDs[i]):(proteinIDs[i+1]-1)),] }
batch3 = Map(cbind, batch3, proteinID = lapply(batch3, function(x) x$V1[1]))
batch3 = do.call(rbind, batch3)
batch3 = batch3[-c(grep("Protein", batch3$V1),grep("Sequence", batch3$V1)),]
batch3$proteinID = gsub("Protein: ","", batch3$proteinID)
batch3$proteinID = gsub("(Hs/Mm)","", batch3$proteinID, fixed = T)
colnames(batch3) = c("Sequence Position","Genomic Coordinate","Motif","K-mer","Z-score", "P-value", "editingID", "proteinID")
head(batch3)

x = read.table("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/RBPMap_editingSites_Batch4_results.txt", 
               header = FALSE, sep = "\t", col.names = paste0("V",seq_len(6)), fill = TRUE)
ast = grep("***", x$V1, fixed = T)
sep = list()
sep[[1]] = x[c(1:ast[1]-1),]
for (i in 1:length(ast)) { sep[[i+1]] = x[c((ast[i]+1):(ast[i+1]-1)),] }
sep = Map(cbind, sep[2:length(sep)], editingID = lapply(sep[2:length(sep)], function(x) x$V1[1]))
sep = lapply(sep, function(x) x[3:nrow(x),])
sep = do.call(rbind, sep)
proteinIDs = grep("Protein", sep$V1, fixed = T)
batch4 = list()
for (i in 1:length(proteinIDs)) { batch4[[i]] = sep[c((proteinIDs[i]):(proteinIDs[i+1]-1)),] }
batch4 = Map(cbind, batch4, proteinID = lapply(batch4, function(x) x$V1[1]))
batch4 = do.call(rbind, batch4)
batch4 = batch4[-c(grep("Protein", batch4$V1),grep("Sequence", batch4$V1)),]
batch4$proteinID = gsub("Protein: ","", batch4$proteinID)
batch4$proteinID = gsub("(Hs/Mm)","", batch4$proteinID, fixed = T)
colnames(batch4) = c("Sequence Position","Genomic Coordinate","Motif","K-mer","Z-score", "P-value", "editingID", "proteinID")
head(batch4)

rbpmap = rbind(batch1, batch2, batch3, batch4)
rbpmap = rbind(rbpmap[1:445058,], batch4)
dim(rbpmap)
save(rbpmap, file = "./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_results_editing_sites.rda")