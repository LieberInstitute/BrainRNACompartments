library(GenomicRanges)
library(data.table)
library(ggplot2)


load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_results_editing_sites_highStringency.rda")


### Are all editing sites reflected here?

colnames(rbpmap) = c("Seq_Pos", "Coordinate", "Motif","Kmer","Zscore","Pval","editingRegion","proteinID")
length(unique(rbpmap$editingRegion)) # 18906 of 18907 coordinates


### FDR correct pvalues for motif enrichment and choose the motif with the lowest pvalue per editing site and protein ID 

rbpmap$Pval = as.numeric(as.character(rbpmap$Pval))
rbpmap$Zscore = as.numeric(as.character(rbpmap$Zscore))
rbpmap$Motif = gsub(" ", "", as.character(rbpmap$Motif))
rbpmap$Coordinate = gsub(" ", "", as.character(rbpmap$Coordinate))
rbpmap$editingRegion = as.character(rbpmap$editingRegion)
rbpmap$motifID = paste0(rbpmap$Coordinate, "-", as.numeric(unlist(lapply(strsplit(rbpmap$Coordinate, ":", fixed=T), function(x) x[2])))+
                               nchar(rbpmap$Motif)-1, ":", unlist(lapply(strsplit(rbpmap$editingRegion, ":", fixed=T), function(x) x[3])))
rbpmap$padj = p.adjust(rbpmap$Pval, method = "BH", n = 2155398) # here n = 18907 * 114

rbpmap = data.table(rbpmap)
hist(-log(rbpmap$Pval))
sigrbp = rbpmap[rbpmap[Pval<=0.05,.I[Pval == min(Pval)], by=c("editingRegion","proteinID")]$V1]
length(unique(sigrbp$proteinID)) # 94 RBPs are represented out of 114 motifs tested


### Map the RBP motif results to the corresponding editing sites

sites = data.frame(GRanges(as.character(sigrbp$editingRegion)))
sites$start = sites$end-10
sites$end = sites$end-10
sites = makeGRangesFromDataFrame(sites)
motifsites = GRanges(as.character(sigrbp$motifID))
ov = findOverlaps(motifsites, sites)
sigrbp = sigrbp[queryHits(ov)[which(queryHits(ov)==subjectHits(ov))],] # limit to motifs that directly overlap the editing site
sites = sites[queryHits(ov)[which(queryHits(ov)==subjectHits(ov))],]

editing_anno$start = editing_anno$end
editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns =T)
ov = findOverlaps(sites, editing_anno_gr)
editing_anno = rbind(cbind(editing_anno[subjectHits(ov),,], sigrbp[queryHits(ov),,]),
                     data.frame(editing_anno[-unique(subjectHits(ov)),,], Seq_Pos = "NoMatch", Coordinate = "NoMatch", 
                                Motif = "NoMatch", Kmer = "NoMatch",Zscore = "NoMatch", Pval = "NoMatch", 
                                editingRegion = "NoMatch", proteinID = "NoMatch", motifID = "NoMatch", padj = "NoMatch"))
length(unique(editing_anno[proteinID!="NoMatch",,]$editingID)) # 10336 editing sites significantly overlap an RBP binding motif (before multiple testing correction)



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
unique_all = do.call(rbind, Map(cbind, unique_all, Comparison = as.list(names(unique_all))))

all = lapply(all, function(x) editing_anno[which(editing_anno$editingID %in% x$editingID),,])
all = Map(cbind, all, 
          geneID = lapply(all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
all = do.call(rbind, Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
                         Type = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]), Comparison = as.list(names(all))))



### Plot the counts of RBPs by group and in total

total = data.frame(editing_anno[collapsedconversion == "A:G / T:C",length(unique(editingID)), by = "proteinID"])
total = total[order(total$V1, decreasing = T),]
total$proteinID = factor(total$proteinID, 
                         levels = c("SRSF5","SRSF3","MBNL1","TRA2B","CUG-BP","PTBP1","SRSF2","SRSF7","TARDBP","G3BP2","PABPC4",   
                                    "SNRNP70","SART3","PABPC1","HNRNPC","HNRNPCL1","SRSF1","PABPC3","KHDRBS1","RALY", "PABPC5","TIA1", "HNRNPL",   
                                    "HNRPLL","LIN28A","HNRNPF","HNRNPA1","IGF2BP2","IGF2BP3","YBX1","SRSF9","SRSF10","HNRNPA2B1","NOVA1","RBFOX1",   
                                    "CNOT4","SNRPA","PABPN1","PUM2","RBMS3","QKI","KHDRBS2","DAZAP1","CPEB4","HNRNPH2","RBMS1","MSI1",     
                                    "ZC3H14","RBM24","MATR3","U2AF2","HNRNPU","CPEB2","HuR","ZCRB1","RBM38","FMR1","KHDRBS3","HNRNPM",  
                                    "RBM28","ENOX1","RBM6","ZNF638","A1CF","RBM46","SAMD4A","RBM42","BRUNOL4","HNRNPK","HNRNPA1L2","BRUNOL5",  
                                    "RBM5","TRA2A","HNRNPH1","ANKHD1","YBX2","FXR1","TUT1","PCBP1","SFPQ","ESRP2","RBM41","RBM8A",    
                                    "FXR2","BRUNOL6","PCBP2","PCBP3","ZC3H10","FUS","RBM45","RBM4","RBM3","PPRC1","SRSF6", "NoMatch"))
unique = data.frame(unique_all[,length(unique(editingID)), by = c("Comparison","proteinID")])
unique$proteinID = factor(unique$proteinID, 
                          levels = c("SRSF5","SRSF3","MBNL1","TRA2B","CUG-BP","PTBP1","SRSF2","SRSF7","TARDBP","G3BP2","PABPC4",   
                                     "SNRNP70","SART3","PABPC1","HNRNPC","HNRNPCL1","SRSF1","PABPC3","KHDRBS1","RALY", "PABPC5","TIA1", "HNRNPL",   
                                     "HNRPLL","LIN28A","HNRNPF","HNRNPA1","IGF2BP2","IGF2BP3","YBX1","SRSF9","SRSF10","HNRNPA2B1","NOVA1","RBFOX1",   
                                     "CNOT4","SNRPA","PABPN1","PUM2","RBMS3","QKI","KHDRBS2","DAZAP1","CPEB4","HNRNPH2","RBMS1","MSI1",     
                                     "ZC3H14","RBM24","MATR3","U2AF2","HNRNPU","CPEB2","HuR","ZCRB1","RBM38","FMR1","KHDRBS3","HNRNPM",  
                                     "RBM28","ENOX1","RBM6","ZNF638","A1CF","RBM46","SAMD4A","RBM42","BRUNOL4","HNRNPK","HNRNPA1L2","BRUNOL5",  
                                     "RBM5","TRA2A","HNRNPH1","ANKHD1","YBX2","FXR1","TUT1","PCBP1","SFPQ","ESRP2","RBM41","RBM8A",    
                                     "FXR2","BRUNOL6","PCBP2","PCBP3","ZC3H10","FUS","RBM45","RBM4","RBM3","PPRC1","SRSF6", "NoMatch"))
unique$Comparison = factor(unique$Comparison, levels = c("adultOnly","prenatalOnly","ANnotAC","ACnotAN","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN"))
al = data.frame(all[,length(unique(editingID)), by = c("Comparison","proteinID")])
al$proteinID = factor(al$proteinID, 
                      levels = c("SRSF5","SRSF3","MBNL1","TRA2B","CUG-BP","PTBP1","SRSF2","SRSF7","TARDBP","G3BP2","PABPC4",   
                                 "SNRNP70","SART3","PABPC1","HNRNPC","HNRNPCL1","SRSF1","PABPC3","KHDRBS1","RALY", "PABPC5","TIA1", "HNRNPL",   
                                 "HNRPLL","LIN28A","HNRNPF","HNRNPA1","IGF2BP2","IGF2BP3","YBX1","SRSF9","SRSF10","HNRNPA2B1","NOVA1","RBFOX1",   
                                 "CNOT4","SNRPA","PABPN1","PUM2","RBMS3","QKI","KHDRBS2","DAZAP1","CPEB4","HNRNPH2","RBMS1","MSI1",     
                                 "ZC3H14","RBM24","MATR3","U2AF2","HNRNPU","CPEB2","HuR","ZCRB1","RBM38","FMR1","KHDRBS3","HNRNPM",  
                                 "RBM28","ENOX1","RBM6","ZNF638","A1CF","RBM46","SAMD4A","RBM42","BRUNOL4","HNRNPK","HNRNPA1L2","BRUNOL5",  
                                 "RBM5","TRA2A","HNRNPH1","ANKHD1","YBX2","FXR1","TUT1","PCBP1","SFPQ","ESRP2","RBM41","RBM8A",    
                                 "FXR2","BRUNOL6","PCBP2","PCBP3","ZC3H10","FUS","RBM45","RBM4","RBM3","PPRC1","SRSF6", "NoMatch"))
al$Comparison = factor(al$Comparison, levels = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN"))


pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/RBPMap_protein_counts_total.pdf", width = 26, height = 8)
ggplot(total[total$proteinID!="NoMatch",,], aes(x = proteinID, y = V1)) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("RBP Motif Totals A:G Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/figures/RBPMap_protein_counts_byGroup.pdf", width = 26, height = 16)
ggplot(unique[unique$proteinID!="NoMatch",,], aes(x = proteinID, y = V1)) + geom_bar(stat = "identity") +
  facet_grid(Comparison ~ ., scales = "free_y") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("RBP Motif Totals in A:G Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(al[al$proteinID!="NoMatch",,], aes(x = proteinID, y = V1)) + geom_bar(stat = "identity") +
  facet_grid(Comparison ~ ., scales = "free_y") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("RBP Motif Totals in A:G Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()



### Are certain RBPs enriched for unique editing site groups?

rbpIDs = as.character(total$proteinID)
comp = list(c("adultOnly","adultAll"), c("prenatalOnly","prenatalAll"), c("ANnotAC","allAN"), c("ACnotAN","allAC"),c("PCnotPN","allPC"), 
            c("PNnotPC","allPN"),c("ACnotPC","allAC"), c("PCnotAC","allPC"), c("ANnotPN","allAN"), c("PNnotAN","allPN"))
al = rbind(al, data.frame(Comparison = "allPC", proteinID = "SRSF6", V1 = 0))
tables = rep.int(list(vector("list", length = length(comp))), length(rbpIDs))

for (i in 1:length(rbpIDs)) {
  for (j in 1:length(comp)) {
    if (rbpIDs[i] %in% unique[which(unique$Comparison==comp[[j]][1]),"proteinID"]) {
    tables[[i]][[j]] = data.frame(yesRBP = c(unique[which(unique$Comparison==comp[[j]][1] & unique$proteinID==rbpIDs[i]),"V1"],
                                             al[which(al$Comparison==comp[[j]][2] & al$proteinID==rbpIDs[i]),"V1"]-
                                               unique[which(unique$Comparison==comp[[j]][1] & unique$proteinID==rbpIDs[i]),"V1"]),
                                  noRBP = c(sum(unique[which(unique$Comparison==comp[[j]][1] & unique$proteinID!=rbpIDs[i]),"V1"]),
                                            sum(al[which(al$Comparison==comp[[j]][2] & al$proteinID!=rbpIDs[i]),"V1"])-
                                              sum(unique[which(unique$Comparison==comp[[j]][1] & unique$proteinID!=rbpIDs[i]),"V1"])),
                                  row.names = c("UniqueSite","NotUnique"))
    } else {
      tables[[i]][[j]] = data.frame(yesRBP = c(0, al[which(al$Comparison==comp[[j]][2] & al$proteinID==rbpIDs[i]),"V1"]),
                                    noRBP = c(sum(unique[which(unique$Comparison==comp[[j]][1] & unique$proteinID!=rbpIDs[i]),"V1"]),
                                              sum(al[which(al$Comparison==comp[[j]][2] & al$proteinID!=rbpIDs[i]),"V1"])-
                                                sum(unique[which(unique$Comparison==comp[[j]][1] & unique$proteinID!=rbpIDs[i]),"V1"])),
                                    row.names = c("UniqueSite","NotUnique"))
  } }
   names(tables[[i]]) = lapply(comp, function(x) x[1]) 
}
names(tables) = rbpIDs
fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, RBP = as.list(names(fisher)), lapply(fisher, function(x) 
     do.call(rbind, Map(cbind, Group = as.list(names(x)),lapply(x, function(y) 
       data.frame(pval = y$p.value, OR = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote=F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_motif_enrichment_byGroup.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_motif_enrichment_byGroup.csv")
df[df$FDR<=0.05,]
# none


## are certain annotations of editing sites more likely to overlap an RBP?

x = data.frame(editing_anno[,length(unique(editingID)), by = c("proteinID","annotation")])
x$proteinID = factor(x$proteinID, 
                    levels = c("SRSF5","SRSF3","MBNL1","TRA2B","CUG-BP","PTBP1","SRSF2","SRSF7","TARDBP","G3BP2","PABPC4",   
                               "SNRNP70","SART3","PABPC1","HNRNPC","HNRNPCL1","SRSF1","PABPC3","KHDRBS1","RALY", "PABPC5","TIA1", "HNRNPL",   
                               "HNRPLL","LIN28A","HNRNPF","HNRNPA1","IGF2BP2","IGF2BP3","YBX1","SRSF9","SRSF10","HNRNPA2B1","NOVA1","RBFOX1",   
                               "CNOT4","SNRPA","PABPN1","PUM2","RBMS3","QKI","KHDRBS2","DAZAP1","CPEB4","HNRNPH2","RBMS1","MSI1",     
                               "ZC3H14","RBM24","MATR3","U2AF2","HNRNPU","CPEB2","HuR","ZCRB1","RBM38","FMR1","KHDRBS3","HNRNPM",  
                               "RBM28","ENOX1","RBM6","ZNF638","A1CF","RBM46","SAMD4A","RBM42","BRUNOL4","HNRNPK","HNRNPA1L2","BRUNOL5",  
                               "RBM5","TRA2A","HNRNPH1","ANKHD1","YBX2","FXR1","TUT1","PCBP1","SFPQ","ESRP2","RBM41","RBM8A",    
                               "FXR2","BRUNOL6","PCBP2","PCBP3","ZC3H10","FUS","RBM45","RBM4","RBM3","PPRC1","SRSF6", "NoMatch"))
tables = rep.int(list(vector("list", length = length(anno))), length(rbpIDs))
anno = as.character(unique(x$annotation))

for (i in 1:length(rbpIDs)) {
  for (j in 1:length(anno)) {
    if (rbpIDs[i] %in% x[which(x$annotation==anno[j]),"proteinID"]) {
      tables[[i]][[j]] = data.frame(yesRBP = c(x[which(x$annotation==anno[j] & x$proteinID==rbpIDs[i]),"V1"],
                                               sum(x[which(x$annotation!=anno[j] & x$proteinID==rbpIDs[i]),"V1"])),
                                    noRBP = c(sum(x[which(x$annotation==anno[j] & x$proteinID!=rbpIDs[i]),"V1"]),
                                              sum(x[which(x$annotation!=anno[j] & x$proteinID!=rbpIDs[i]),"V1"])),
                                    row.names = c("anno","Not"))
    } else {
      tables[[i]][[j]] = data.frame(yesRBP = c(0, sum(x[which(x$annotation!=anno[j] & x$proteinID==rbpIDs[i]),"V1"])),
                                    noRBP = c(sum(x[which(x$annotation==anno[j] & x$proteinID!=rbpIDs[i]),"V1"]),
                                              sum(x[which(x$annotation!=anno[j] & x$proteinID!=rbpIDs[i]),"V1"])),
                                    row.names = c("anno","Not"))
    } }
  names(tables[[i]]) = anno 
}
names(tables) = rbpIDs
fisher = lapply(tables, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, RBP = as.list(names(fisher)), lapply(fisher, function(x) 
  do.call(rbind, Map(cbind, Group = as.list(names(x)),lapply(x, function(y) 
    data.frame(pval = y$p.value, OR = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
df
write.csv(df, quote=F,file="./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_motif_enrichment_byAnnotation.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/RNA_binding_proteins/data/rbpmap_motif_enrichment_byAnnotation.csv")
x = df[df$FDR<=0.05 & df$OR>1 & df$Group!="Intergenic",]
x[order(x$Group),]
#      X     RBP  Group          pval       OR           FDR
#14   14   SRSF3  3'UTR  4.073403e-03 1.153741  4.499690e-02
#144 144 IGF2BP2  3'UTR  1.116565e-04 1.832377  3.119813e-03
#149 149 IGF2BP3  3'UTR  2.813087e-05 1.923309  1.214742e-03
#289 289    FMR1  3'UTR  2.227528e-03 2.152279  2.859664e-02
#5     5 NoMatch  5'UTR  2.340941e-43 3.447016  5.559736e-41
#2     2 NoMatch    CDS 4.864645e-115 4.355730 2.310706e-112
#8     8   SRSF5 Intron  2.632366e-06 1.225028  2.083956e-04
#53   53   G3BP2 Intron  5.111727e-04 1.532568  8.671679e-03
#128 128  LIN28A Intron  2.460095e-03 1.583692  3.075119e-02




# more information can be found at http://rbpmap.technion.ac.il/