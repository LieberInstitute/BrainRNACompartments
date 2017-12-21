library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


### Isolate the sites unique to each group

unique = lapply(unique_bySamp, function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique = Map(cbind, unique, geneID = lapply(unique, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
unique = Map(cbind, unique, ensID = lapply(unique, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]),
                 baseMean = lapply(unique, function(x) Ipres.down[match(x$geneID, rownames(Ipres.down)),"baseMean"]), 
                 P.LFC = lapply(unique, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "log2FoldChange"]), 
                 P.SE = lapply(unique, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "lfcSE"]), 
                 P.padj = lapply(unique, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "padj"]),
                 A.LFC = lapply(unique, function(x) Apres[match(x$geneID, rownames(Apres)), "log2FoldChange"]), 
                 A.SE = lapply(unique, function(x) Apres[match(x$geneID, rownames(Apres)), "lfcSE"]),
                 A.padj = lapply(unique, function(x) Apres[match(x$geneID, rownames(Apres)), "padj"]),
                 N.LFC = lapply(unique, function(x) Npres[match(x$geneID, rownames(Npres)), "log2FoldChange"]), 
                 N.SE = lapply(unique, function(x) Npres[match(x$geneID, rownames(Npres)), "lfcSE"]), 
                 N.padj = lapply(unique, function(x) Npres[match(x$geneID, rownames(Npres)), "padj"]),
                 C.LFC = lapply(unique, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "log2FoldChange"]), 
                 C.SE = lapply(unique, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "lfcSE"]),
                 C.padj = lapply(unique, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "padj"]))
DEGnames = c(lapply(sig[1:8], function(x) as.character(x$geneID)), lapply(age.sig, function(x) as.character(x$geneID))) 
unique = Map(cbind, unique, both_retained = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$both_retained, "both_retained", "no")),
                 both_exported = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$both_exported, "both_exported", "no")),
                 Fet_retained = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
                 Ad_retained = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
                 Fet_exported = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
                 Ad_exported = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
                 ret_Ad_exp_Fet = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
                 ret_Fet_exp_Ad = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
                 interacting = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$interacting, "interacting", "no")),
                 both_decreasing = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
                 both_increasing = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$both_increasing, "both_increasing", "no")),
                 Cyt_decreasing = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
                 Nuc_decreasing = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
                 Cyt_increasing = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
                 Nuc_increasing = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
                 decr_Nuc_incr_Cyt = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
                 decr_Cyt_incr_Nuc = lapply(unique, function(x) ifelse(x$geneID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))
lapply(unique, head)

all = Map(cbind, all, geneID = lapply(all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]),
          baseMean = lapply(all, function(x) Ipres.down[match(x$geneID, rownames(Ipres.down)),"baseMean"]), 
          P.LFC = lapply(all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "log2FoldChange"]), 
          P.SE = lapply(all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "lfcSE"]), 
          P.padj = lapply(all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "padj"]),
          A.LFC = lapply(all, function(x) Apres[match(x$geneID, rownames(Apres)), "log2FoldChange"]), 
          A.SE = lapply(all, function(x) Apres[match(x$geneID, rownames(Apres)), "lfcSE"]),
          A.padj = lapply(all, function(x) Apres[match(x$geneID, rownames(Apres)), "padj"]),
          N.LFC = lapply(all, function(x) Npres[match(x$geneID, rownames(Npres)), "log2FoldChange"]), 
          N.SE = lapply(all, function(x) Npres[match(x$geneID, rownames(Npres)), "lfcSE"]), 
          N.padj = lapply(all, function(x) Npres[match(x$geneID, rownames(Npres)), "padj"]),
          C.LFC = lapply(all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "log2FoldChange"]), 
          C.SE = lapply(all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "lfcSE"]),
          C.padj = lapply(all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "padj"]),
          both_retained = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$both_retained, "both_retained", "no")),
          both_exported = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$both_exported, "both_exported", "no")),
          Fet_retained = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
          Ad_retained = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
          Fet_exported = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
          Ad_exported = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
          ret_Ad_exp_Fet = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
          ret_Fet_exp_Ad = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
          interacting = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$interacting, "interacting", "no")),
          both_decreasing = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
          both_increasing = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$both_increasing, "both_increasing", "no")),
          Cyt_decreasing = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
          Nuc_decreasing = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
          Cyt_increasing = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
          Nuc_increasing = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
          decr_Nuc_incr_Cyt = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
          decr_Cyt_incr_Nuc = lapply(all, function(x) ifelse(x$geneID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))



### How many sites are there per gene, and is this different than the distribution of sites that aren't specific to a group and in all samples in a group?

## Number of sites per gene
numsites_bygene = lapply(unique, function(x) x[,length(unique(editingID)),by="geneID"])
numsites_bygene_all = lapply(all, function(x) x[,length(unique(editingID)),by="geneID"])
numsites_bygene = lapply(numsites_bygene,as.data.frame)
numsites_bygene_all = lapply(numsites_bygene_all, as.data.frame)
elementNROWS(numsites_bygene)
elementNROWS(numsites_bygene_all)

## Compare to the number of sites per gene in sites present in each compartment but not unique, and the total number of sites in the same genes as those of the unique sites
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
t.num.bygene = t.editedgenes.numsites = list()
for (i in 1:length(numsites_bygene)){
  t.num.bygene[[i]] = t.test(x = numsites_bygene[[group[i]]][,"V1"], y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"geneID"] %in% numsites_bygene[[group[i]]][,"geneID"]),"V1"])
  t.editedgenes.numsites[[i]] = t.test(x = numsites_bygene_all[[allgroup[i]]][which(numsites_bygene_all[[allgroup[i]]][,"geneID"] %in% numsites_bygene[[group[i]]][,"geneID"]),"V1"], 
                                       y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"geneID"] %in% numsites_bygene[[group[i]]][,"geneID"]),"V1"])
}
names(t.num.bygene) = names(t.editedgenes.numsites) = group
write.csv(rbind(Tstat = data.frame(lapply(t.num.bygene, function(x) round(x$statistic,3))), pval = data.frame(lapply(t.num.bygene, function(x) x$p.value)),
               confInt = data.frame(lapply(t.num.bygene, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(t.num.bygene, function(x) round(x$estimate,3))),
               max = data.frame(lapply(numsites_bygene, function(x) max(x$V1))), median = data.frame(lapply(numsites_bygene, function(x) median(x$V1))), 
               sd = data.frame(lapply(numsites_bygene, function(x) round(sd(x$V1),3)))), quote=F,
         file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_numEditingSites_byGene_uniqueVSall.csv")

write.csv(rbind(Tstat = data.frame(lapply(t.editedgenes.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedgenes.numsites, function(x) x$p.value)),
                confInt = data.frame(lapply(t.editedgenes.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedgenes.numsites, function(x) x$estimate)),
                max = mapply(function(x,y) max(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"geneID"] %in% numsites_bygene[[x]][,"geneID"]),"V1"]), group, allgroup),
                median = mapply(function(x,y) median(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"geneID"] %in% numsites_bygene[[x]][,"geneID"]),"V1"]), group, allgroup),
                sd = mapply(function(x,y) sd(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"geneID"] %in% numsites_bygene[[x]][,"geneID"]),"V1"]), group, allgroup)),
          quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_numEditingSites_byGene_uniquesitecontaininggenesVSnonuniquesitecontaininggenes_allSites.csv")


## What about at the exon level?

exonMap$exonID = rownames(exonMap)
exonMap_gr = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
ov = findOverlaps(editing_anno_gr, exonMap_gr)
tog = cbind(editing_anno[queryHits(ov),], exonMap[subjectHits(ov),])
unique_exons = lapply(unique, function(x) tog[(editingID %in% x$editingID),,])
all_exons = lapply(all, function(x) tog[(editingID %in% x$editingID),,])

numsites_byexon = lapply(unique_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon_all = lapply(all_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon = lapply(numsites_byexon,as.data.frame)
numsites_byexon_all = lapply(numsites_byexon_all, as.data.frame)
elementNROWS(numsites_byexon)
elementNROWS(numsites_byexon_all)

t.num.byexon = t.editedexons.numsites = list()
for (i in 1:length(numsites_byexon)){
  t.num.byexon[[i]] = t.test(x = numsites_byexon[[group[i]]][,"V1"], y = numsites_byexon_all[[allgroup[i]]][-which(numsites_byexon_all[[allgroup[i]]][,"exonID"] %in% numsites_byexon[[group[i]]][,"exonID"]),"V1"])
  t.editedexons.numsites[[i]] = t.test(x = numsites_byexon_all[[allgroup[i]]][which(numsites_byexon_all[[allgroup[i]]][,"exonID"] %in% numsites_byexon[[group[i]]][,"exonID"]),"V1"], 
                                       y = numsites_byexon_all[[allgroup[i]]][-which(numsites_byexon_all[[allgroup[i]]][,"exonID"] %in% numsites_byexon[[group[i]]][,"exonID"]),"V1"])
}
names(t.num.byexon) = names(t.editedexons.numsites) = group
write.csv(rbind(Tstat = data.frame(lapply(t.num.byexon, function(x) round(x$statistic,3))), pval = data.frame(lapply(t.num.byexon, function(x) x$p.value)),
                confInt = data.frame(lapply(t.num.byexon, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(t.num.byexon, function(x) round(x$estimate,3))),
                max = data.frame(lapply(numsites_byexon, function(x) max(x$V1))), median = data.frame(lapply(numsites_byexon, function(x) median(x$V1))), 
                sd = data.frame(lapply(numsites_byexon, function(x) round(sd(x$V1),3)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_numEditingSites_byExon_uniqueVSall.csv")

write.csv(rbind(Tstat = data.frame(lapply(t.editedexons.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedexons.numsites, function(x) x$p.value)),
                confInt = data.frame(lapply(t.editedexons.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedexons.numsites, function(x) x$estimate)),
                max = mapply(function(x,y) max(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                median = mapply(function(x,y) median(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                sd = mapply(function(x,y) sd(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup)),
          quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_numEditingSites_byExon_uniquesitecontainingexonsVSnonuniquesitecontainingexons_allSites.csv")



### Compare the proportion of editing site in each annotation in sites unique to a specific group and those that aren't
# Using whether the site overlaps an annotation group at all rather than the CDS -> UTR -> Intron -> Intergenic heirarchy

anno.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
for (i in 1:length(unique)){
  anno.site[[i]] = list(
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("CDS", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("CDS", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("CDS", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("CDS", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("Intron", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("Intron", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("Intron", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("Intron", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("5'UTR", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("5'UTR", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("5'UTR", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("5'UTR", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("5'UTR", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("5'UTR", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("3'UTR", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("3'UTR", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("3'UTR", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("3'UTR", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("3'UTR", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("3'UTR", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][(annotation=="Intergenic"),list(unique(editingID)),]), nrow(unique[[group[i]]][(annotation!="Intergenic"),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][(annotation=="Intergenic"),list(unique(editingID)),])-nrow(unique[[group[i]]][(annotation=="Intergenic"),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][(annotation!="Intergenic"),list(unique(editingID)),])-nrow(unique[[group[i]]][(annotation!="Intergenic"),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")))
  names(anno.site[[i]]) = c("CDS","Intron","5'UTR","3'UTR","Intergenic")
}
names(anno.site) = group
fisher.anno.site = lapply(anno.site, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numEditingSites_uniqueornot_byAnnotation.csv")
pval = unlist(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
names(pval[pval<=0.0008333333]) # three pound signs indicate p value <=0.0008333333 (Bonferoni corrected threshold) 



### Characterize the overlap with differentially expressed genes

## Are editing sites specific to a group more or less likely to fall in an exported or retained gene by fraction than chance?

# number of editing sites falling within significantly retained or imported genes that are group-specific, or those found in the same group that aren't specific
frac.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
for (i in 1:length(unique)){
  frac.site[[i]] = list(
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(frac.site[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG")
}
names(frac.site) = group
fisher.frac.site = lapply(frac.site, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numEditingSites_byUnique_byFraction_exportedVSretained.csv")



# number of significantly retained or exported genes that either contain an editing site specific to a group, or don't
frac.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  frac.gene[[i]] = list(
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(geneID)),]),
                            length(unique(sig$both_retained$geneID))-nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(geneID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(geneID)),]),
                            length(unique(sig$both_exported$geneID))-nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(geneID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(geneID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Ad_retained$geneID)))-
                              nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(geneID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(geneID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Ad_exported$geneID)))-
                              nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(geneID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(geneID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Fet_retained$geneID)))-
                              nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(geneID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(geneID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Fet_exported$geneID)))-
                              nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(geneID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(frac.gene[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG")
}
names(frac.gene) = group
fisher.frac.gene = lapply(frac.gene, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numGenes_byUnique_byFraction_exportedVSretained.csv")


## Are editing sites specific to a group more or less likely to fall in an increasing or decreasing genes by age than chance?

# number of editing sites falling within significantly increasing or decreasing genes that are specific and in all samples, or those found in the same group that aren't specific
age.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  age.site[[i]] = list(
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(age.site[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(age.site) = group
fisher.age.site = lapply(age.site, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numEditingSites_byUnique_byAge_increasingVSdecreasing.csv")


# number of significantly increasing or decreasing genesby age that either contain an editing site specific to a group and in all samples, or don't
age.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  age.gene[[i]] = list(
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(geneID)),]),
                              length(unique(age.sig$both_increasing$geneID))-nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(geneID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(geneID)),]),
                              length(unique(age.sig$both_decreasing$geneID))-nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(geneID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(geneID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Cyt_increasing$geneID)))-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(geneID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(geneID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Cyt_decreasing$geneID)))-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(geneID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(geneID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Nuc_increasing$geneID)))-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(geneID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(geneID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Nuc_decreasing$geneID)))-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(geneID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(age.gene[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(age.gene) = group
fisher.age.gene = lapply(age.gene, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numGenes_byUnique_byAge_increasingVSdecreasing.csv")



### Does editing rate correlate with gene expression in the group the editing site appears?
## correlate LFC with editing rate in editing rates shared between fraction

corr.shared = list(list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
unique_df = lapply(unique, as.data.frame)
all_df = lapply(all, as.data.frame)
for (i in 1:length(unique)){
  corr.shared[[i]] = list(Adult = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                            y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"A.LFC"], use = "complete.obs"),3),
                          Prenatal = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                               y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"P.LFC"], use = "complete.obs"),3),
                          Cytosol = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                              y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"C.LFC"], use = "complete.obs"),3),
                          Nucleus = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                              y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"N.LFC"], use = "complete.obs"),3))
}
names(corr.shared) = c("sharedbyFrac:inCyt","sharedbyFrac:inNuc","sharedbyAge:inAd","sharedbyAge:inPren","sharedbyFrac:inAC","sharedbyFrac:inAN","sharedbyFrac:inPC","sharedbyFrac:inPN",
                       "sharedbyAge:inAC","sharedbyAge:inPC","sharedbyAge:inAN","sharedbyAge:inPN")
write.csv(data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F))), quote=F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/correlation_degLFC_and_editingRate_inSharedSites.csv")

max(abs(data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F))))) # 0.113
min(abs(data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F))))) # 0.002

## correlate LFC with editing rate in editing sites unique to a group

corr.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  corr.site[[i]] = list(Adult = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"A.LFC"], use = "complete.obs"),3),
                        Prenatal = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"P.LFC"], use = "complete.obs"),3),
                        Cytosol = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"C.LFC"], use = "complete.obs"),3),
                        Nucleus = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"N.LFC"], use = "complete.obs"),3))
}
names(corr.site) = group
write.csv(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F))), quote=F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/correlation_degLFC_and_editingRate_inUniqueSites.csv")

max(abs(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F))))) # 0.148
min(abs(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F))))) # 0



### Is gene expression greater in the compartment/age exhibiting the editing site than in the compared group?
## t test of LFC between unique sites in both groups

t.site = list(list(),list(),list(),list(),list(),list())
comps = list(byFraction = c("cytosolOnly","nucleusOnly"), byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.site[[i]] = list(adult = t.test(x = unique_df[[comps[[i]][1]]][,"A.LFC"], y = unique_df[[comps[[i]][2]]][,"A.LFC"]),
                     prenatal = t.test(x = unique_df[[comps[[i]][1]]][,"P.LFC"], y = unique_df[[comps[[i]][2]]][,"P.LFC"]),
                     cytosol = t.test(x = unique_df[[comps[[i]][1]]][,"C.LFC"], y = unique_df[[comps[[i]][2]]][,"C.LFC"]),
                     nucleus = t.test(x = unique_df[[comps[[i]][1]]][,"N.LFC"], y = unique_df[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.site) = names(comps)
write.csv(rbind(Tstat = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                pval = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                confInt = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                estMeans = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_LFC_ofUniqueEditingSiteGenes_byGroup.csv")



### In intronic editing sites, is the IR ratio greater in the compartment exhibiting the site than the compared group?

names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
  IRres[[i]][,"Chr"] = paste0("chr", IRres[[i]][,"Chr"])
}
names(IRres) = shortenedNames
IRres_gr = lapply(IRres, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr", start.field = "Start", end.field = "End", strand.field = "Direction", keep.extra.columns = T))
ov = lapply(IRres_gr, function(x) findOverlaps(editing_anno_gr, x))
togintron = list()
for (i in 1:length(IRres)){
  tmp = IRres[[i]]
  togintron[[i]] = cbind(editing_anno[queryHits(ov[[i]]),], IRratio = tmp$IRratio[subjectHits(ov[[i]])],
                   sampleIDintron = names(IRres)[i])
}
togintron = do.call(rbind, togintron)
IRratio_editing = lapply(unique_df, function(x) togintron[which(togintron$editingID %in% x$editingID),])

t.IRratio = list()
for (i in 1:length(comps)){
  t.IRratio[[i]] = t.test(x = IRratio_editing[[comps[[i]][1]]][,"IRratio"], y = IRratio_editing[[comps[[i]][2]]][,"IRratio"])
}
names(t.IRratio) = names(comps)
t.IRratio.editing = rbind(Tstat = data.frame(lapply(t.IRratio, function(x) x$statistic)), pval = data.frame(lapply(t.IRratio, function(x) x$p.value)),
                          confInt = data.frame(lapply(t.IRratio, function(x) x$conf.int)), estMeans = data.frame(lapply(t.IRratio, function(x) x$estimate)))
write.csv(t.IRratio.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_IRratio_ofUniqueEditingSiteIntrons_byGroup.csv", quote=F)



### in 3'UTR editing sites, is the exon differentially expressed by fraction or group?

tog.3UTR = tog[which(tog$UTR3=="3'UTR"),]
counts3UTR = exonCounts.down[which(rownames(exonCounts.down) %in% tog.3UTR$exonID), grep("polyA", colnames(exonCounts.down))]
match(rownames(pd[grep("polyA", rownames(pd)),]), colnames(counts3UTR))

## DESeq2 on exons, is the exon the last exon with the greatest count
exonsA = DESeqDataSetFromMatrix(countData = counts3UTR[,-grep("53", colnames(counts3UTR))], 
                                    colData = pd[which(pd$Fetal == "Adult" & pd$Library =="polyA"),], design = ~ Zone)
exonsF = DESeqDataSetFromMatrix(countData = counts3UTR[,grep("53", colnames(counts3UTR))], 
                             colData = pd[which(pd$Fetal == "Prenatal" & pd$Library =="polyA"),], design = ~ Zone)
exonsC = DESeqDataSetFromMatrix(countData = counts3UTR[,grep("C", colnames(counts3UTR))], 
                                colData = pd[which(pd$Zone == "Cytosol" & pd$Library =="polyA"),], design = ~ Fetal)
exonsN = DESeqDataSetFromMatrix(countData = counts3UTR[,grep("N", colnames(counts3UTR))], 
                                colData = pd[which(pd$Zone == "Nucleus" & pd$Library =="polyA"),], design = ~ Fetal)
exonsA = DESeq(exonsA)
exonsF = DESeq(exonsF)
exonsC = DESeq(exonsC)
exonsN = DESeq(exonsN)
exonres = list(exonsA.res = results(exonsA), exonsF.res = results(exonsF), exonsC.res = results(exonsC), exonsN.res = results(exonsN))
elementNROWS(exonres)
elementNROWS(lapply(exonres, function(x) x[which(x$padj<=0.05),]))

tog.3UTR = cbind(tog.3UTR, A.LFC = exonres[["exonsA.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsA.res"]])),"log2FoldChange"],
                 A.padj = exonres[["exonsA.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsA.res"]])),"padj"],
                 P.LFC = exonres[["exonsF.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsF.res"]])),"log2FoldChange"],
                 P.padj = exonres[["exonsF.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsF.res"]])),"padj"],
                 C.LFC = exonres[["exonsC.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsC.res"]])),"log2FoldChange"],
                 C.padj = exonres[["exonsC.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsC.res"]])),"padj"],
                 N.LFC = exonres[["exonsN.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsN.res"]])),"log2FoldChange"],
                 N.padj = exonres[["exonsN.res"]][match(tog.3UTR$exonID, rownames(exonres[["exonsN.res"]])),"padj"])


## Compare LFC and significance of groups of editing sites in 3'UTR

unique_3UTR = lapply(unique, function(x) tog.3UTR[which(tog.3UTR$editingID %in% x$editingID),])
t.3UTR.site = list(list())
comps = list(byFraction = c("cytosolOnly","nucleusOnly"), byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.3UTR.site[[i]] = list(adult = t.test(x = unique_3UTR[[comps[[i]][1]]][,"A.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = unique_3UTR[[comps[[i]][1]]][,"P.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = unique_3UTR[[comps[[i]][1]]][,"C.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = unique_3UTR[[comps[[i]][1]]][,"N.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.3UTR.site) = names(comps)
write.csv(rbind(Tstat = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                pval = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                confInt = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                estMeans = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_3UTR.LFC_ofUniqueEditingSite3UTR_byGroup.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.3UTR.site = list(list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
for (i in 1:length(group)){
  corr.3UTR.site[[i]] = list(fraction = cor(x = unique_3UTR[[group[i]]][,"A.LFC"], y = unique_3UTR[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = unique_3UTR[[group[i]]][,"C.LFC"], y = unique_3UTR[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.3UTR.site) = group
data.frame(lapply(corr.3UTR.site, function(x) unlist(x, recursive=F)))
#         cytosolOnly nucleusOnly adultOnly prenatalOnly   ACnotAN   ANnotAC   PCnotPN   PNnotPC   ACnotPC   PCnotAC   ANnotPN   PNnotAN
#fraction   0.4270921   0.5724515 0.5993952    0.4833157 0.4410391 0.6397155 0.6105422 0.6139386 0.6168924 0.6260270 0.6321284 0.5066662
#age        0.9383614   0.9405513 0.9308989    0.9610020 0.9172280 0.9232515 0.9333838 0.9464980 0.9348204 0.9512715 0.9321368 0.9612866


# Correlate editing rate and 3'UTR LFC by group in shared editing sites

all_3UTR = lapply(all, function(x) tog.3UTR[which(tog.3UTR$editingID %in% x$editingID),])
corr.shared = list(list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
unique_df = lapply(unique_3UTR, as.data.frame)
all_df = lapply(all_3UTR, as.data.frame)
for (i in 1:length(unique)){
  corr.shared[[i]] = list(Adult = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                            y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"A.LFC"], use = "complete.obs"),3),
                          Prenatal = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                               y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"P.LFC"], use = "complete.obs"),3),
                          Cytosol = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                              y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"C.LFC"], use = "complete.obs"),3),
                          Nucleus = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                              y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"N.LFC"], use = "complete.obs"),3))
}
names(corr.shared) = c("sharedbyFrac:inCyt","sharedbyFrac:inNuc","sharedbyAge:inAd","sharedbyAge:inPren","sharedbyFrac:inAC","sharedbyFrac:inAN","sharedbyFrac:inPC","sharedbyFrac:inPN",
                       "sharedbyAge:inAC","sharedbyAge:inPC","sharedbyAge:inAN","sharedbyAge:inPN")
write.csv(data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F))), quote=F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/correlation_3UTR.LFC_and_editingRate_inSharedSites.csv")


# Correlate editing rate and 3'UTR LFC by group in unique editing sites

corr.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  corr.site[[i]] = list(Adult = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"A.LFC"], use = "complete.obs"),3),
                        Prenatal = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"P.LFC"], use = "complete.obs"),3),
                        Cytosol = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"C.LFC"], use = "complete.obs"),3),
                        Nucleus = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"N.LFC"], use = "complete.obs"),3))
}
names(corr.site) = group
write.csv(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F))), quote=F,
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/correlation_3UTR.LFC_and_editingRate_inUniqueSites.csv")


## Compare the number of edited 3'UTRs that are significantly and non-significantly up- and down-expressed per group

fisher.3UTR = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  fisher.3UTR[[i]] = list(
    data.frame(exported = c(nrow(unique_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(unique_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(unique_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(unique_3UTR[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_3UTR[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(unique_3UTR[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_3UTR[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(unique_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(unique_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_3UTR[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_3UTR[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.3UTR[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.3UTR) = group
fisher.3UTR.editing = lapply(fisher.3UTR, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.3UTR.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.3UTR.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))), quote=F,
          file= "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numEditingSites_byUnique_byFraction_byAge_3UTR.csv")



### in CDS editing sites, is the exon differentially expressed by fraction or group?

tog.CDS = tog[which(tog$cds=="CDS"),]
countsCDS = exonCounts.down[which(rownames(exonCounts.down) %in% tog.CDS$exonID), grep("polyA", colnames(exonCounts.down))]
match(rownames(pd[grep("polyA", rownames(pd)),]), colnames(countsCDS))

## DESeq2 on exons, is the exon the last exon with the greatest count
exonsA = DESeqDataSetFromMatrix(countData = countsCDS[,-grep("53", colnames(countsCDS))], 
                                colData = pd[which(pd$Fetal == "Adult" & pd$Library =="polyA"),], design = ~ Zone)
exonsF = DESeqDataSetFromMatrix(countData = countsCDS[,grep("53", colnames(countsCDS))], 
                                colData = pd[which(pd$Fetal == "Prenatal" & pd$Library =="polyA"),], design = ~ Zone)
exonsC = DESeqDataSetFromMatrix(countData = countsCDS[,grep("C", colnames(countsCDS))], 
                                colData = pd[which(pd$Zone == "Cytosol" & pd$Library =="polyA"),], design = ~ Fetal)
exonsN = DESeqDataSetFromMatrix(countData = countsCDS[,grep("N", colnames(countsCDS))], 
                                colData = pd[which(pd$Zone == "Nucleus" & pd$Library =="polyA"),], design = ~ Fetal)
exonsA = DESeq(exonsA)
exonsF = DESeq(exonsF)
exonsC = DESeq(exonsC)
exonsN = DESeq(exonsN)
exonres = list(exonsA.res = results(exonsA), exonsF.res = results(exonsF), exonsC.res = results(exonsC), exonsN.res = results(exonsN))
elementNROWS(exonres)
elementNROWS(lapply(exonres, function(x) x[which(x$padj<=0.05),]))

tog.CDS = cbind(tog.CDS, A.LFC = exonres[["exonsA.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsA.res"]])),"log2FoldChange"],
                A.padj = exonres[["exonsA.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsA.res"]])),"padj"],
                P.LFC = exonres[["exonsF.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsF.res"]])),"log2FoldChange"],
                P.padj = exonres[["exonsF.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsF.res"]])),"padj"],
                C.LFC = exonres[["exonsC.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsC.res"]])),"log2FoldChange"],
                C.padj = exonres[["exonsC.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsC.res"]])),"padj"],
                N.LFC = exonres[["exonsN.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsN.res"]])),"log2FoldChange"],
                N.padj = exonres[["exonsN.res"]][match(tog.CDS$exonID, rownames(exonres[["exonsN.res"]])),"padj"])


## Compare LFC and significance of groups of editing sites in CDS

unique_CDS = lapply(unique, function(x) tog.CDS[which(tog.CDS$editingID %in% x$editingID),])
t.CDS.site = list(list())
comps = list(byFraction = c("cytosolOnly","nucleusOnly"), byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.CDS.site[[i]] = list(adult = t.test(x = unique_CDS[[comps[[i]][1]]][,"A.LFC"], y = unique_CDS[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = unique_CDS[[comps[[i]][1]]][,"P.LFC"], y = unique_CDS[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = unique_CDS[[comps[[i]][1]]][,"C.LFC"], y = unique_CDS[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = unique_CDS[[comps[[i]][1]]][,"N.LFC"], y = unique_CDS[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.CDS.site) = names(comps)
write.csv(rbind(Tstat = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                pval = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                confInt = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                estMeans = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/ttest_CDS.LFC_ofUniqueEditingSiteCDSs_byGroup.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.CDS.site = list(list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
for (i in 1:length(group)){
  corr.CDS.site[[i]] = list(fraction = cor(x = unique_CDS[[group[i]]][,"A.LFC"], y = unique_CDS[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = unique_CDS[[group[i]]][,"C.LFC"], y = unique_CDS[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.CDS.site) = group
data.frame(lapply(corr.CDS.site, function(x) unlist(x, recursive=F)))
#         cytosolOnly nucleusOnly adultOnly prenatalOnly    ACnotAN   ANnotAC   PCnotPN   PNnotPC   ACnotPC   PCnotAC   ANnotPN   PNnotAN
#fraction   0.2564349   0.5747068 0.4567068    0.4021567 0.06840717 0.6816928 0.5565109 0.7186106 0.3936159 0.7456216 0.4921897 0.2411119
#age        0.9465056   0.7582785 0.8292040    0.8769978 0.88817216 0.8204965 0.9346631 0.6444958 0.8414013 0.8521293 0.8313999 0.8227004


## Compare the number of edited CDSs that are significantly and non-significantly up- and down-expressed per group

fisher.CDS = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  fisher.CDS[[i]] = list(
    data.frame(exported = c(nrow(unique_CDS[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(unique_CDS[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_CDS[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(unique_CDS[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(unique_CDS[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_CDS[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(unique_CDS[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_CDS[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_CDS[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(unique_CDS[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_CDS[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(unique_CDS[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_CDS[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_CDS[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_CDS[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_CDS[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(unique_CDS[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.CDS[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.CDS) = group
fisher.CDS.editing = lapply(fisher.CDS, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.CDS.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.CDS.editing, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_numEditingSites_byUnique_byFraction_byAge_CDS.csv")

# Are group-specific edited sequences present in an unedited state in other groups, or is expression sequestered to the group showing the site?