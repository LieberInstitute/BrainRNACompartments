library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")


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
elementNROWS(lapply(unique_bySamp_all, function(x) x[which(x$AllNos==9),]))

unique_all = lapply(unique_bySamp_all[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = Map(cbind, unique_all, ensID = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]),
                 baseMean = lapply(unique_all, function(x) Ipres.down[match(x$nearestID, rownames(Ipres.down)),"baseMean"]), 
                 Prenatal.LFC = lapply(unique_all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "log2FoldChange"]), 
                 Prenatal.SE = lapply(unique_all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "lfcSE"]), 
                 Prenatal.padj = lapply(unique_all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "padj"]),
                 Adult.LFC = lapply(unique_all, function(x) Apres[match(x$nearestID, rownames(Apres)), "log2FoldChange"]), 
                 Adult.SE = lapply(unique_all, function(x) Apres[match(x$nearestID, rownames(Apres)), "lfcSE"]),
                 Adult.padj = lapply(unique_all, function(x) Apres[match(x$nearestID, rownames(Apres)), "padj"]),
                 Nucleus.LFC = lapply(unique_all, function(x) Npres[match(x$nearestID, rownames(Npres)), "log2FoldChange"]), 
                 Nucleus.SE = lapply(unique_all, function(x) Npres[match(x$nearestID, rownames(Npres)), "lfcSE"]), 
                 Nucleus.padj = lapply(unique_all, function(x) Npres[match(x$nearestID, rownames(Npres)), "padj"]),
                 Cytosol.LFC = lapply(unique_all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "log2FoldChange"]), 
                 Cytosol.SE = lapply(unique_all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "lfcSE"]),
                 Cytosol.padj = lapply(unique_all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "padj"]))
DEGnames = c(lapply(sig[1:8], function(x) as.character(x$geneID)), lapply(age.sig, function(x) as.character(x$geneID))) 
unique_all = Map(cbind, unique_all, both_retained = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_retained, "both_retained", "no")),
                 both_exported = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_exported, "both_exported", "no")),
                 Fet_retained = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
                 Ad_retained = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
                 Fet_exported = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
                 Ad_exported = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
                 ret_Ad_exp_Fet = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
                 ret_Fet_exp_Ad = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
                 interacting = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$interacting, "interacting", "no")),
                 both_decreasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
                 both_increasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_increasing, "both_increasing", "no")),
                 Cyt_decreasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
                 Nuc_decreasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
                 Cyt_increasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
                 Nuc_increasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
                 decr_Nuc_incr_Cyt = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
                 decr_Cyt_incr_Nuc = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))
elementNROWS(unique_all)
lapply(unique_all, head)

all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]),
          baseMean = lapply(all, function(x) Ipres.down[match(x$nearestID, rownames(Ipres.down)),"baseMean"]), 
          Prenatal.LFC = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "log2FoldChange"]), 
          Prenatal.SE = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "lfcSE"]), 
          Prenatal.padj = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "padj"]),
          Adult.LFC = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "log2FoldChange"]), 
          Adult.SE = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "lfcSE"]),
          Adult.padj = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "padj"]),
          Nucleus.LFC = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "log2FoldChange"]), 
          Nucleus.SE = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "lfcSE"]), 
          Nucleus.padj = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "padj"]),
          Cytosol.LFC = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "log2FoldChange"]), 
          Cytosol.SE = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "lfcSE"]),
          Cytosol.padj = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "padj"]),
          both_retained = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_retained, "both_retained", "no")),
          both_exported = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_exported, "both_exported", "no")),
          Fet_retained = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
          Ad_retained = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
          Fet_exported = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
          Ad_exported = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
          ret_Ad_exp_Fet = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
          ret_Fet_exp_Ad = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
          interacting = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$interacting, "interacting", "no")),
          both_decreasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
          both_increasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_increasing, "both_increasing", "no")),
          Cyt_decreasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
          Nuc_decreasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
          Cyt_increasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
          Nuc_increasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
          decr_Nuc_incr_Cyt = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
          decr_Cyt_incr_Nuc = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))

elementNROWS(lapply(unique_all, function(x) unique(x$editingID)))
#adultOnly prenatalOnly      ANnotAC      ACnotAN      ANnotPN      PNnotAN      ACnotPC      PCnotAC      PCnotPN      PNnotPC 
#130           89          159           30          554          286          273          261           12           65

### How many sites are there per gene, and is this different than the distribution of sites that aren't specific to a group and in all samples in a group?

## Number of sites per gene
numsites_bygene = lapply(unique_all, function(x) x[,length(unique(editingID)),by="nearestID"])
numsites_bygene_all = lapply(all, function(x) x[,length(unique(editingID)),by="nearestID"])
numsites_bygene = lapply(numsites_bygene,as.data.frame)
numsites_bygene_all = lapply(numsites_bygene_all, as.data.frame)


## Compare to the number of sites per gene in sites present in each compartment but not unique, and the total number of sites in the same genes as those of the unique sites
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
t.num.bygene = t.editedgenes.numsites = list()
for (i in 1:length(numsites_bygene)){
  t.num.bygene[[i]] = t.test(x = numsites_bygene[[group[i]]][,"V1"], y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"nearestID"] %in% numsites_bygene[[group[i]]][,"nearestID"]),"V1"])
  t.editedgenes.numsites[[i]] = t.test(x = numsites_bygene_all[[allgroup[i]]][which(numsites_bygene_all[[allgroup[i]]][,"nearestID"] %in% numsites_bygene[[group[i]]][,"nearestID"]),"V1"], 
                                       y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"nearestID"] %in% numsites_bygene[[group[i]]][,"nearestID"]),"V1"])
}
names(t.num.bygene) = names(t.editedgenes.numsites) = group
write.csv(rbind(Tstat = data.frame(lapply(t.num.bygene, function(x) round(x$statistic,3))), pval = data.frame(lapply(t.num.bygene, function(x) x$p.value)),
                confInt = data.frame(lapply(t.num.bygene, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(t.num.bygene, function(x) round(x$estimate,3))),
                max = data.frame(lapply(numsites_bygene, function(x) max(x$V1))), median = data.frame(lapply(numsites_bygene, function(x) median(x$V1))), 
                sd = data.frame(lapply(numsites_bygene, function(x) round(sd(x$V1),3)))), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perGene_inUniqueSiteGenes_vsNonUniqueSiteGenes.csv")
write.csv(rbind(Tstat = data.frame(lapply(t.editedgenes.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedgenes.numsites, function(x) x$p.value)),
                confInt = data.frame(lapply(t.editedgenes.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedgenes.numsites, function(x) x$estimate)),
                max = mapply(function(x,y) max(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup),
                median = mapply(function(x,y) median(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup),
                sd = mapply(function(x,y) sd(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup)),quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perGene_inUniqueSiteGenes_vsNonUniqueSiteGenes_notJustUniqueSites.csv")


## What about at the exon level?

exonMap$exonID = rownames(exonMap)
exonMap_gr = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
ov = findOverlaps(editing_anno_gr, exonMap_gr)
tog = cbind(editing_anno[queryHits(ov),], exonMap[subjectHits(ov),])
unique_exons = lapply(unique_all, function(x) tog[(editingID %in% x$editingID),,])
all_exons = lapply(all, function(x) tog[(editingID %in% x$editingID),,])

numsites_byexon = lapply(unique_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon_all = lapply(all_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon = lapply(numsites_byexon,as.data.frame)
numsites_byexon_all = lapply(numsites_byexon_all, as.data.frame)

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
                             sd = data.frame(lapply(numsites_byexon, function(x) round(sd(x$V1),3)))),quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perExon_inUniqueSiteExons_vsNonUniqueSiteExons.csv")
write.csv(rbind(Tstat = data.frame(lapply(t.editedexons.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedexons.numsites, function(x) x$p.value)),
                confInt = data.frame(lapply(t.editedexons.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedexons.numsites, function(x) x$estimate)),
                max = mapply(function(x,y) max(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                median = mapply(function(x,y) median(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                sd = mapply(function(x,y) sd(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup)),quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perExon_inUniqueSiteExons_vsNonUniqueSiteExons_notJustUniqueSites.csv")



### Compare the proportion of editing sites in each annotation in sites unique to a group found in all samples of that group and those that aren’t
# Using whether the site overlaps an annotation group at all rather than the CDS -> UTR -> Intron -> Other heirarchy

anno.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
unique = lapply(unique_bySamp[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
for (i in 1:length(unique_all)){
  anno.site[[i]] = list(
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("CDS", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("CDS", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("CDS", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("Intron", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("Intron", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("Intron", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("UTR5", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("UTR5", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("UTR5", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("UTR3", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("UTR3", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("UTR3", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][(annotation=="Other"),list(unique(editingID)),]), nrow(unique_all[[group[i]]][(annotation!="Other"),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][(annotation=="Other"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(annotation=="Other"),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][(annotation!="Other"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(annotation!="Other"),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")))
  names(anno.site[[i]]) = c("CDS","Intron","UTR5","UTR3","Other")
}
names(anno.site) = group
fisher.anno.site = lapply(anno.site, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
      data.frame(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) round(y$estimate,3)), recursive=F)))),quote = F,
file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_annotation_enrichment_UniqueSitesinAll_vsUniqueSitesNotInAll.csv")
pval = unlist(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
names(pval[pval<=0.05]) # three pound signs indicate p value <=0.001 (Bonferoni corrected threshold) 
# adultOnly.CDS is enriched compared to unique sites not in all samples in a group
# adultOnly.Intron is depleted compared to unique sites not in all samples in a group
# adultOnly.UTR3 is enriched compared to unique sites not in all samples in a group
### prenatalOnly.UTR3 is enriched compared to unique sites not in all samples in a group
# prenatalOnly.Other is depleted compared to unique sites not in all samples in a group
# ANnotAC.CDS is enriched compared to unique sites not in all samples in a group
# ANnotAC.Intron is enriched compared to unique sites not in all samples in a group
# PCnotAC.UTR3 is enriched compared to unique sites not in all samples in a group
### ANnotPN.UTR3 is enriched compared to unique sites not in all samples in a group
# PNnotAN.CDS is depleted compared to unique sites not in all samples in a group
### PNnotAN.UTR3 is enriched compared to unique sites not in all samples in a group
# PNnotAN.Other is depleted compared to unique sites not in all samples in a group
lapply(anno.site, function(y) data.frame(lapply(y, function(x) c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), 
                                                                 col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2])))))



### Characterize the overlap with differentially expressed genes

## Are editing sites specific to a group and shared by all samples in a group more or less likely to fall in an exported or retained gene by fraction than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
frac.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
for (i in 1:length(unique_all)){
  frac.site[[i]] = list(
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),])),
               notretained = c(nrow(unique_all[[group[i]]][(both_retained=="no"),list(unique(editingID)),]),
                               nrow(all[[allgroup[i]]][(both_retained=="no"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_retained=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])),
               notretained = c(nrow(unique_all[[group[i]]][(both_retained=="no" & Ad_retained=="no"),list(unique(editingID)),]),
                               nrow(all[[allgroup[i]]][(both_retained=="no" & Ad_retained=="no"),list(unique(editingID)),])-
                                 nrow(unique_all[[group[i]]][(both_retained=="no" & Ad_retained=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])),
               notretained = c(nrow(unique_all[[group[i]]][(both_retained=="no" & Fet_retained=="no"),list(unique(editingID)),]),
                               nrow(all[[allgroup[i]]][(both_retained=="no" & Fet_retained=="no"),list(unique(editingID)),])-
                                 nrow(unique_all[[group[i]]][(both_retained=="no" & Fet_retained=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),])),
               notexported = c(nrow(unique_all[[group[i]]][(both_exported=="no"),list(unique(editingID)),]),
                               nrow(all[[allgroup[i]]][(both_exported=="no"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_exported=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])),
               notexported = c(nrow(unique_all[[group[i]]][(both_exported=="no" & Ad_exported=="no"),list(unique(editingID)),]),
                               nrow(all[[allgroup[i]]][(both_exported=="no" & Ad_exported=="no"),list(unique(editingID)),])-
                                 nrow(unique_all[[group[i]]][(both_exported=="no" & Ad_exported=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])),
               notexported = c(nrow(unique_all[[group[i]]][(both_exported=="no" & Fet_exported=="no"),list(unique(editingID)),]),
                               nrow(all[[allgroup[i]]][(both_exported=="no" & Fet_exported=="no"),list(unique(editingID)),])-
                                 nrow(unique_all[[group[i]]][(both_exported=="no" & Fet_exported=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(frac.site[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG","bothagesRetained","adultRetained","prenatalRetained","bothagesExported","adultExported","prenatalExported")
}
names(frac.site) = group
fisher.frac.site = lapply(frac.site, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) round(y$estimate,3)), recursive=F)))),quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editingsite_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_FracDEGs.csv")
pval = unlist(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
names(pval[pval<=0.0005555556]) # three pound signs indicate p value <=0.0005555556 (Bonferoni corrected threshold) 
# ACnotAN.adultDEG sites are depleted in retained vs exported genes
# ANnotAC.adultDEG sites are enriched in retained vs exported genes
# ANnotAC.bothagesRetained sites are enriched in retained genes vs genes that are not significantly retained
### ANnotAC.adultRetained sites are enriched in retained genes vs genes that are not significantly retained
### PNnotPC.bothagesRetained sites are enriched in retained genes vs genes that are not significantly retained
### PNnotPC.adultRetained sites are enriched in retained genes vs genes that are not significantly retained
### PNnotPC.prenatalRetained sites are enriched in retained genes vs genes that are not significantly retained
# ACnotPC.adultExported sites are enriched in exported genes vs genes that are not significantly exported
lapply(frac.site, function(y) data.frame(lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2])))))


# number of significantly retained or exported genes that either contain an editing site specific to a group and in all samples, or don't
frac.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  frac.gene[[i]] = list(
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),]),
                            length(unique(sig$both_retained$geneID))-nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),]),
                            length(unique(sig$both_exported$geneID))-nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Ad_retained$geneID)))-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Ad_exported$geneID)))-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Fet_retained$geneID)))-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Fet_exported$geneID)))-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),]),
                            length(unique(sig$both_retained$geneID))-nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),])),
               notretained = c(nrow(unique_all[[group[i]]][(both_retained=="no"),list(unique(nearestID)),]),
                               nrow(geneCounts)-length(unique(sig$both_retained$geneID))-nrow(unique_all[[group[i]]][(both_retained=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Ad_retained$geneID)))-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),])),
               notretained = c(nrow(unique_all[[group[i]]][(both_retained=="no" & Ad_retained=="no"),list(unique(nearestID)),]),
                               nrow(geneCounts)-sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Ad_retained$geneID)))-
                                 nrow(unique_all[[group[i]]][(both_retained=="no" & Ad_retained=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Fet_retained=="Fet_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Fet_retained$geneID)))-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),])),
               notretained = c(nrow(unique_all[[group[i]]][(both_retained=="no" & Fet_retained=="no"),list(unique(nearestID)),]),
                               nrow(geneCounts)-sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Fet_retained$geneID)))-
                                 nrow(unique_all[[group[i]]][(both_retained=="no" & Fet_retained=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),]),
                            length(unique(sig$both_exported$geneID))-nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),])),
               notexported = c(nrow(unique_all[[group[i]]][(both_exported=="no"),list(unique(nearestID)),]),
                               nrow(geneCounts.down)-length(unique(sig$both_exported$geneID))-nrow(unique_all[[group[i]]][(both_exported=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Ad_exported$geneID)))-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),])),
               notexported = c(nrow(unique_all[[group[i]]][(both_exported=="no" & Ad_exported=="no"),list(unique(nearestID)),]),
                               nrow(geneCounts.down)-sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Ad_exported$geneID)))-
                                 nrow(unique_all[[group[i]]][(both_exported=="no" & Ad_exported=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Fet_exported=="Fet_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Fet_exported$geneID)))-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),])),
               notexported = c(nrow(unique_all[[group[i]]][(both_exported=="no" & Fet_exported=="no"),list(unique(nearestID)),]),
                               nrow(geneCounts.down)-sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Fet_exported$geneID)))-
                                 nrow(unique_all[[group[i]]][(both_exported=="no" & Fet_exported=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
    names(frac.gene[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG","bothagesRetained","adultRetained","prenatalRetained","bothagesExported","adultExported","prenatalExported")
}
names(frac.gene) = group
fisher.frac.gene = lapply(frac.gene, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) round(y$estimate,3)), recursive=F)))),quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedgene_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_FracDEGs.csv")
pval = unlist(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
names(pval[pval<=0.05]) # three pound signs indicate p value <=0.0005555556 (Bonferoni corrected threshold) 
# adultOnly.adultRetained are enriched in retained genes vs not retained
# adultOnly.adultExported are enriched in exported genes vs not exported
# prenatalOnly.adultRetained are enriched in retained genes vs not retained
# ACnotAN.adultDEG are depleted in retained genes vs exported
### ACnotAN.adultExported are enriched in exported genes vs not exported
### ANnotAC.adultDEG are enriched in retained genes vs exported
### ANnotAC.bothagesRetained are enriched in retained genes vs not retained
### ANnotAC.adultRetained are enriched in retained genes vs not retained
# ANnotAC.prenatalRetained are enriched in retained genes vs not retained
# ANnotAC.adultExported are enriched in exported genes vs not exported
# PNnotPC.adultDEG are enriched in retained genes vs exported
### PNnotPC.bothagesRetained are enriched in retained genes vs not retained
### PNnotPC.adultRetained are enriched in retained genes vs not retained
### PNnotPC.prenatalRetained are enriched in retained genes vs not retained
# ACnotPC.adultRetained are enriched in retained genes vs not retained
# ACnotPC.adultExported are enriched in exported genes vs not exported
# PCnotAC.adultRetained are enriched in retained genes vs not retained
# ANnotPN.adultDEG are enriched in retained genes vs exported
### ANnotPN.bothagesRetained are enriched in retained genes vs not retained
### ANnotPN.adultRetained are enriched in retained genes vs not retained
### ANnotPN.prenatalRetained are enriched in retained genes vs not retained
### ANnotPN.adultExported are enriched in exported genes vs not exported
# PNnotAN.adultDEG are enriched in retained genes vs exported
### PNnotAN.bothagesRetained are enriched in retained genes vs not retained
### PNnotAN.adultRetained are enriched in retained genes vs not retained
### PNnotAN.prenatalRetained are enriched in retained genes vs not retained
lapply(frac.gene, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



## Are editing sites specific to a group and shared by all samples in a group more or less likely to fall in an increasing or decreasing genes by age than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
age.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  age.site[[i]] = list(
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])),
               notincreasing = c(nrow(unique_all[[group[i]]][(both_increasing=="no"),list(unique(editingID)),]),
                                 nrow(all[[allgroup[i]]][(both_increasing=="no"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_increasing=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])),
               notincreasing = c(nrow(unique_all[[group[i]]][(both_increasing=="no" & Cyt_increasing=="no"),list(unique(editingID)),]),
                                 nrow(all[[allgroup[i]]][(both_increasing=="no" & Cyt_increasing=="no"),list(unique(editingID)),])-
                                   nrow(unique_all[[group[i]]][(both_increasing=="no" & Cyt_increasing=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])),
               notincreasing = c(nrow(unique_all[[group[i]]][(both_increasing=="no" & Nuc_increasing=="no"),list(unique(editingID)),]),
                                 nrow(all[[allgroup[i]]][(both_increasing=="no" & Nuc_increasing=="no"),list(unique(editingID)),])-
                                   nrow(unique_all[[group[i]]][(both_increasing=="no" & Nuc_increasing=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])),
               notdecreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="no"),list(unique(editingID)),]),
                                 nrow(all[[allgroup[i]]][(both_decreasing=="no"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(both_decreasing=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])),
               notdecreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="no" & Cyt_decreasing=="no"),list(unique(editingID)),]),
                                 nrow(all[[allgroup[i]]][(both_decreasing=="no" & Cyt_decreasing=="no"),list(unique(editingID)),])-
                                   nrow(unique_all[[group[i]]][(both_decreasing=="no" & Cyt_decreasing=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])),
               notdecreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="no" & Nuc_decreasing=="no"),list(unique(editingID)),]),
                                 nrow(all[[allgroup[i]]][(both_decreasing=="no" & Nuc_decreasing=="no"),list(unique(editingID)),])-
                                   nrow(unique_all[[group[i]]][(both_decreasing=="no" & Nuc_decreasing=="no"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(age.site[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG","bothfractionsIncreasing","CytosolIncreasing","NucleusIncreasing","bothfractionsDecreasing","CytosolDecreasing","NucleusDecreasing")
}
names(age.site) = group
fisher.age.site = lapply(age.site, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) round(y$estimate,3)), recursive=F)))),quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editingsite_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_AgeDEGs.csv")
pval = unlist(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
names(pval[pval<=0.0005555556])[c(grep("bothfractions",names(pval[pval<=0.0005555556])))] # three pound signs indicate p value <=0.0005555556 (Bonferoni corrected threshold) 
### adultOnly.bothfractionsDEG : unique sites enriched in increasing genes
### adultOnly.bothfractionsIncreasinG : unique sites enriched in increasing genes
# adultOnly.bothfractionsDecreasinG : unique sites depleted in decreasing genes
### prenatalOnly.bothfractionsDEG : unique sites depleted in increasing genes
# prenatalOnly.bothfractionsIncreasinG : unique sites depleted in increasing genes
### prenatalOnly.bothfractionsDecreasinG : unique sites enriched in decreasing genes
# ANnotAC.NucleusIncreasinG : unique sites enriched increasing genes
# PNnotPC.CytosolDEG : unique sites depleted in increasing genes  
# PNnotPC.CytosolIncreasinG : unique sites depleted in increasing genes
### ACnotPC.bothfractionsDEG : unique sites enriched in increasing genes
### ACnotPC.bothfractionsIncreasinG : unique sites enriched in increasing genes
### ACnotPC.bothfractionsDecreasinG : unique sites depleted in decreasing genes
### PCnotAC.bothfractionsDEG : unique sites depleted in increasing genes
# PCnotAC.bothfractionsIncreasinG : unique sites depleted in increasing genes
### PCnotAC.bothfractionsDecreasinG : unique sites enriched in decreasing genes
### ANnotPN.bothfractionsDEG : unique sites enriched in increasing genes
### ANnotPN.bothfractionsIncreasinG : unique sites enriched in increasing genes
### ANnotPN.bothfractionsDecreasinG : unique sites depleted in decreasing genes
### PNnotAN.bothfractionsDEG : unique sites depleted in increasing genes
# PNnotAN.bothfractionsIncreasinG : unique sites depleted in increasing genes
### PNnotAN.bothfractionsDecreasinG : unique sites enriched in decreasing genes
lapply(age.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))


# number of significantly increasing or decreasing genes by age that either contain an editing site specific to a group and in all samples, or don't
age.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  age.gene[[i]] = list(
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),]),
                              length(unique(age.sig$both_decreasing$geneID))-nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),])),
               increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),]),
                              length(unique(age.sig$both_increasing$geneID))-nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Cyt_decreasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),])),
               increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Cyt_increasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Nuc_decreasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),])),
               increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Nuc_increasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),]),
                              length(unique(age.sig$both_decreasing$geneID))-nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),])),
               notdecreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="no"),list(unique(nearestID)),]),
                                 nrow(geneCounts)-length(unique(age.sig$both_decreasing$geneID))-nrow(unique_all[[group[i]]][(both_decreasing=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Cyt_decreasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),])),
               notdecreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="no" & Cyt_decreasing=="no"),list(unique(nearestID)),]),
                                 nrow(geneCounts.down)-sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Cyt_decreasing$geneID)))-
                                   nrow(unique_all[[group[i]]][(both_decreasing=="no" & Cyt_decreasing=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Nuc_decreasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),])),
               notdecreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="no" & Nuc_decreasing=="no"),list(unique(nearestID)),]),
                                 nrow(geneCounts)-sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Nuc_decreasing$geneID)))-
                                   nrow(unique_all[[group[i]]][(both_decreasing=="no" & Nuc_decreasing=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),]),
                              length(unique(age.sig$both_increasing$geneID))-nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),])),
               notincreasing = c(nrow(unique_all[[group[i]]][(both_increasing=="no"),list(unique(nearestID)),]),
                                 nrow(geneCounts)-length(unique(age.sig$both_increasing$geneID))-nrow(unique_all[[group[i]]][(both_increasing=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Cyt_increasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),])),
               notincreasing = c(nrow(unique_all[[group[i]]][(both_increasing=="no" & Cyt_increasing=="no"),list(unique(nearestID)),]),
                                 nrow(geneCounts.down)-sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Cyt_increasing$geneID)))-
                                   nrow(unique_all[[group[i]]][(both_increasing=="no" & Cyt_increasing=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Nuc_increasing$geneID)))-
                                nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),])),
               notincreasing = c(nrow(unique_all[[group[i]]][(both_increasing=="no" & Nuc_increasing=="no"),list(unique(nearestID)),]),
                                 nrow(geneCounts)-sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Nuc_increasing$geneID)))-
                                   nrow(unique_all[[group[i]]][(both_increasing=="no" & Nuc_increasing=="no"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(age.gene[[i]]) = c("bothagesDEG","CytosolDEG","NucleusDEG","bothagesdecreasing","Cytosoldecreasing","Nucleusdecreasing","bothagesincreasing","Cytosolincreasing","Nucleusincreasing")
}
names(age.gene) = group
fisher.age.gene = lapply(age.gene, function(x) lapply(x, fisher.test))
write.csv(rbind(data.frame(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                data.frame(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) round(y$estimate,3)), recursive=F)))),quote = F,
          file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedgene_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_AgeDEGs.csv")
pval = unlist(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
names(pval[pval<=0.0005555556]) # three pound signs indicate p value <=0.0005555556 (Bonferoni corrected threshold) 
### adultOnly.bothagesDEG depleted in decreasing
### adultOnly.bothagesincreasing enriched in increasing
### prenatalOnly.bothagesDEG enriched in decreasing
### prenatalOnly.bothagesdecreasing enriched in decreasing
### ANnotAC.NucleusDEG depleted in decreasing
### ANnotAC.bothagesincreasing enriched in increasing
### PNnotPC.bothagesDEG enriched in decreasing
### PNnotPC.bothagesdecreasing enriched in decreasing
### ACnotPC.bothagesDEG depleted in decreasing
### ACnotPC.bothagesincreasing enriched in increasing
### PCnotAC.bothagesDEG enriched in decreasing
### PCnotAC.bothagesdecreasing enriched in decreasing
### ANnotPN.bothagesDEG depleted in decreasing
### ANnotPN.bothagesincreasing enriched in increasing
### PNnotAN.bothagesDEG enriched in decreasing
### PNnotAN.bothagesdecreasing enriched in decreasing
lapply(age.gene, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### Does editing rate correlate with gene expression in the group the editing site appears?
## correlate LFC with editing rate in editing sites unique to a group

corr.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
unique_all_df = lapply(unique_all, as.data.frame)
for (i in 1:length(unique_all)){
  corr.site[[i]] = list(Adult = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Adult.LFC"], use = "complete.obs"),3),
                             Prenatal = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Prenatal.LFC"], use = "complete.obs"),3),
                             Cytosol = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Cytosol.LFC"], use = "complete.obs"),3),
                             Nucleus = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Nucleus.LFC"], use = "complete.obs"),3))
}
names(corr.site) = group
data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))
#         adultOnly prenatalOnly ACnotAN ANnotAC PCnotPN PNnotPC ACnotPC PCnotAC ANnotPN PNnotAN
#Adult        0.141       -0.225  -0.192   0.034   0.138   0.071   0.168  -0.114   0.080  -0.204
#Prenatal     0.221       -0.022  -0.019   0.060   0.044   0.079   0.257   0.009   0.123  -0.022
#Cytosol     -0.416       -0.376  -0.038  -0.014   0.289  -0.032  -0.268  -0.239  -0.190  -0.237
#Nucleus     -0.417       -0.343   0.044  -0.007   0.346  -0.027  -0.276  -0.222  -0.196  -0.167
max(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))) # 0.346
min(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))) # -0.417
min(abs(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F))))) # 0.007


### Is gene expression greater in the compartment/age exhibiting the editing site than in the compared group?
## t test of LFC between unique sites in both groups

t.frac.site = t.age.site = list(list())
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.frac.site[[i]] = list(adult = t.test(x = unique_all_df[[comps[[i]][1]]][,"Adult.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Adult.LFC"]),
                          prenatal = t.test(x = unique_all_df[[comps[[i]][1]]][,"Prenatal.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Prenatal.LFC"]))
  t.age.site[[i]] = list(cytosol = t.test(x = unique_all_df[[comps[[i]][1]]][,"Cytosol.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Cytosol.LFC"]),
                         nucleus = t.test(x = unique_all_df[[comps[[i]][1]]][,"Nucleus.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Nucleus.LFC"]))
}
names(t.frac.site) = names(t.age.site) = names(comps)
t.LFC.editing = rbind(frac.Tstat = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      age.Tstat = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      frac.pval = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      age.pval = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      frac.confInt = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                      age.confInt = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))),
                      frac.estMeans = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))),
                      age.estMeans = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


## Plots?


### In intronic editing sites, is the IR ratio greater in the compartment exhibiting the site than the compared group?
## t test of IR ratio of introns with editing sites by group

names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
  IRres[[i]][,"Chr"] = paste0("chr", IRres[[i]][,"Chr"])
}
names(IRres) = shortenedNames
IRres_gr = lapply(IRres, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr",start.field = "Start",end.field = "End",strand.field = "Direction",keep.extra.columns = T))
ov = lapply(IRres_gr, function(x,y) findOverlaps(x, editing_anno_gr))
togintron = list()
for (i in 1:length(IRres)){
  tmp = IRres[[i]]
  togintron[[i]] = cbind(editing_anno[subjectHits(ov[[i]]),], IRratio = tmp$IRratio[queryHits(ov[[i]])],
                   sampleIDintron = names(IRres)[i], intronID = paste0(tmp$Chr[queryHits(ov[[i]])],":",tmp$Start[queryHits(ov[[i]])],"-",tmp$End[queryHits(ov[[i]])],":",tmp$Direction[queryHits(ov[[i]])]))
}
togintron = do.call(rbind, togintron)
IRratio_editing = lapply(unique_all_df, function(x) togintron[which(togintron$editingID %in% x$editingID),])

t.IRratio = list()
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.IRratio[[i]] = t.test(x = IRratio_editing[[comps[[i]][1]]][,"IRratio"], y = IRratio_editing[[comps[[i]][2]]][,"IRratio"])
}
names(t.IRratio) = names(comps)
t.IRratio.editing = rbind(Tstat = data.frame(lapply(t.IRratio, function(x) x$statistic)), pval = data.frame(lapply(t.IRratio, function(x) x$p.value)),
                          confInt = data.frame(lapply(t.IRratio, function(x) x$conf.int)), estMeans = data.frame(lapply(t.IRratio, function(x) x$estimate)))
write.csv(t.IRratio.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.IRratio.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


## Are the introns that are edited expressed in the other group where they are not edited?
IRratio_editing_dt = lapply(IRratio_editing, data.table)
IRratio_editing_dt = lapply(IRratio_editing_dt, function(x) as.data.frame(x[,list(unique(IRratio)),by=c("editingID","sampleIDintron","intronID")]))
IRratio_editing_dt = Map(cbind, IRratio_editing_dt, 
                         Fraction = lapply(IRratio_editing_dt, function(x) ifelse(x$sampleIDintron %in% c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1"),"Cytosol","Nucleus")),
                         Age = lapply(IRratio_editing_dt, function(x) ifelse(x$sampleIDintron %in% c("Br1113C1","Br2046C","Br2074C","Br1113N1","Br2046N","Br2074N"),"Adult","Prenatal")))
IRratio_editing_dt = Map(cbind, IRratio_editing_dt, Group = lapply(IRratio_editing_dt, function(x) as.character(paste(x$Age, x$Fraction, sep=":"))))
IRratio_editing_dt = lapply(IRratio_editing_dt, data.table)
lapply(IRratio_editing_dt, function(x) x[,mean(V1),by="Group"])

# Count introns with 0% retention in a group
x = lapply(IRratio_editing_dt, function(x) x[, sum(V1), by=c("editingID","intronID","Group")])
x = lapply(x, function(x) x[V1==0, length(unique(intronID)), by="Group"])
y = lapply(IRratio_editing_dt, function(x) x[, length(unique(intronID)), by="Group"])
y = do.call(rbind,Map(cbind,y,SiteGroup=as.list(names(y))))
y$Zero = NA
x = do.call(rbind,Map(cbind,x[3:8],SiteGroup=as.list(names(x)[3:8])))
class(y) = class(x) = "data.frame"
x$combo = paste0(x$Group,":",x$SiteGroup) 
y$combo = paste0(y$Group,":",y$SiteGroup)
for (i in 1:nrow(y)){
  if (y$combo[i] %in% x$combo) {
  y[i,"Zero"] = x[which(x$Group==y$Group[i] & x$SiteGroup==y$SiteGroup[i]),"V1"]
} else
  y[i,"Zero"] = 0
}
df = data.frame(percent = c(round(y$Zero/y$V1*100,2),round((y$V1-y$Zero)/y$V1*100,2)), Group = rep.int(y$Group,2), SiteGroup = rep.int(y$SiteGroup,2), 
                IR = c(rep.int("Not Expressed",nrow(y)),rep.int("Expressed",nrow(y))), TotalIntrons = c(rep.int(NA, nrow(y)),y$V1))
df$SiteGroup = factor(df$SiteGroup, levels = c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))


# Plot intron retention ratios
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/introns_overlappingUniqueInAll3_editingSites_PercentNotExpresseded_inOtherGroups.pdf",width=18,height=9)
ggplot(df, aes(x = Group, y = percent, fill=IR)) + geom_bar(stat = "identity") +
  facet_grid(. ~ SiteGroup) +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("Percent Expressed Introns Containing a Unique Editing Site\nPresent In All Group Samples") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()







### in 3'UTR editing sites, is the exon differentially expressed by group?

tog.3UTR = tog[which(tog$UTR3=="UTR3"),]
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
tog.3UTR = cbind(tog.3UTR, A.LFC = exonsA.res[match(tog.3UTR$exonID, rownames(exonsA.res)),"log2FoldChange"],
                 A.padj = exonsA.res[match(tog.3UTR$exonID, rownames(exonsA.res)),"padj"],
                 P.LFC = exonsF.res[match(tog.3UTR$exonID, rownames(exonsF.res)),"log2FoldChange"],
                 P.padj = exonsF.res[match(tog.3UTR$exonID, rownames(exonsF.res)),"padj"],
                 C.LFC = exonsC.res[match(tog.3UTR$exonID, rownames(exonsC.res)),"log2FoldChange"],
                 C.padj = exonsC.res[match(tog.3UTR$exonID, rownames(exonsC.res)),"padj"],
                 N.LFC = exonsN.res[match(tog.3UTR$exonID, rownames(exonsN.res)),"log2FoldChange"],
                 N.padj = exonsN.res[match(tog.3UTR$exonID, rownames(exonsN.res)),"padj"])


## Compare LFC and significance of groups of editing sites in 3'UTR

uniqueAll_3UTR = lapply(unique_all, function(x) tog.3UTR[which(tog.3UTR$editingID %in% x$editingID),])
t.3UTR.site = list(list())
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.3UTR.site[[i]] = list(adult = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"A.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"P.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"C.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"N.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.3UTR.site) = names(comps)
t.LFC.3UTR.editing = rbind(Tstat = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      pval = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      confInt = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                      estMeans = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.3UTR.editing, 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.3UTR.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.3UTR.site = list(list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
for (i in 1:length(group)){
  corr.3UTR.site[[i]] = list(fraction = cor(x = uniqueAll_3UTR[[group[i]]][,"A.LFC"], y = uniqueAll_3UTR[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = uniqueAll_3UTR[[group[i]]][,"C.LFC"], y = uniqueAll_3UTR[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.3UTR.site) = group
data.frame(lapply(corr.3UTR.site, function(x) unlist(x, recursive=F)))
#         adultOnly prenatalOnly   ACnotAN   ANnotAC   PCnotPN   PNnotPC   ACnotPC   PCnotAC   ANnotPN   PNnotAN
#fraction 0.8730890    0.4996052 0.3140912 0.7291933 0.2640446 0.4316696 0.8124232 0.6049925 0.8140242 0.5064677
#age      0.9043124    0.9664685 0.9115609 0.9598261 0.9853123 0.8788988 0.8788835 0.9607036 0.9293007 0.9573088


## Compare the number of edited 3'UTRs that are significantly and non-significantly up- and down-expressed per group

fisher.3UTR = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  fisher.3UTR[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_3UTR[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.3UTR[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.3UTR) = group
fisher.3UTR.editing = lapply(fisher.3UTR, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.3UTR.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                   adultOnly prenatalOnly   ACnotAN    ANnotAC PCnotPN   PNnotPC    ACnotPC    PCnotAC    ANnotPN      PNnotAN
#bothagesDEG      0.08259109    1.0000000 1.0000000 1.00000000       1 1.0000000 0.05520353 1.00000000 0.02656958 1.000000e+00
#AdultDEG         1.00000000    0.3106864 0.2854701 0.30098834       1 0.5764706 0.70051452 0.58521437 0.57841961 7.833541e-01
#PrenatalDEG      0.23972603    1.0000000 1.0000000 1.00000000       1 1.0000000 0.09005497 1.00000000 0.02875768 1.000000e+00
#bothfractionsDEG 0.43214195    1.0000000 1.0000000 1.00000000       1 0.2727273 0.57219722 1.00000000 0.68865682 1.156235e-02
#CytosolDEG       0.12935880    1.0000000 0.7063310 0.06733649       1 0.2450980 0.24832078 1.00000000 0.41151685 6.845838e-02
#NucleusDEG       0.06206983    1.0000000 0.6945843 0.14689203       1 0.2745098 0.02171784 0.02936896 0.02182513 4.126946e-05
fisher.3UTR.props = lapply(fisher.3UTR, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### in 3'UTR editing sites, is the exon differentially expressed by group?

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
tog.CDS = cbind(tog.CDS, A.LFC = exonsA.res[match(tog.CDS$exonID, rownames(exonsA.res)),"log2FoldChange"],
                A.padj = exonsA.res[match(tog.CDS$exonID, rownames(exonsA.res)),"padj"],
                P.LFC = exonsF.res[match(tog.CDS$exonID, rownames(exonsF.res)),"log2FoldChange"],
                P.padj = exonsF.res[match(tog.CDS$exonID, rownames(exonsF.res)),"padj"],
                C.LFC = exonsC.res[match(tog.CDS$exonID, rownames(exonsC.res)),"log2FoldChange"],
                C.padj = exonsC.res[match(tog.CDS$exonID, rownames(exonsC.res)),"padj"],
                N.LFC = exonsN.res[match(tog.CDS$exonID, rownames(exonsN.res)),"log2FoldChange"],
                N.padj = exonsN.res[match(tog.CDS$exonID, rownames(exonsN.res)),"padj"],
                A.basemean = exonsA.res[match(tog.CDS$exonID, rownames(exonsA.res)),"baseMean"],
                P.basemean = exonsF.res[match(tog.CDS$exonID, rownames(exonsF.res)),"baseMean"],
                C.basemean = exonsC.res[match(tog.CDS$exonID, rownames(exonsC.res)),"baseMean"],
                N.basemean = exonsN.res[match(tog.CDS$exonID, rownames(exonsN.res)),"baseMean"])


## Compare LFC and significance of groups of editing sites in 3'UTR

uniqueAll_CDS = lapply(unique_all, function(x) tog.CDS[which(tog.CDS$editingID %in% x$editingID),])
t.CDS.site = list(list())
comps = list(byFractionInAdult = c("ACnotAN","ANnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN")) # not enough observations for byAge = c("adultOnly","prenatalOnly"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), byAgeInCytosol = c("ACnotPC","PCnotAC")
for (i in 1:length(comps)){
  t.CDS.site[[i]] = list(adult = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"A.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"P.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"C.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"N.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.CDS.site) = names(comps)
t.LFC.CDS.editing = rbind(Tstat = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                           pval = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                           confInt = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                           estMeans = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.CDS.editing, 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.CDS.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.CDS.site = list(list())
for (i in 1:length(group)){
  corr.CDS.site[[i]] = list(fraction = cor(x = uniqueAll_CDS[[group[i]]][,"A.LFC"], y = uniqueAll_CDS[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = uniqueAll_CDS[[group[i]]][,"C.LFC"], y = uniqueAll_CDS[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.CDS.site) = group
data.frame(lapply(corr.CDS.site, function(x) unlist(x, recursive=F)))
#         adultOnly ACnotAN   ANnotAC PNnotPC   ACnotPC   ANnotPN PNnotAN
#fraction 0.9463874      NA 0.1542199      NA 0.9444664 0.8277404      NA
#age      0.8829940      NA 0.8139241      NA 0.8841078 0.6899047      NA


## Compare the number of edited CDSs that are significantly and non-significantly up- and down-expressed per group

fisher.CDS = list(list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","ACnotAN","ANnotAC","PNnotPC","ACnotPC","ANnotPN","PNnotAN") #not enough obs for "prenatalOnly","PCnotAC","PCnotPN"
for (i in 1:length(group)){
  fisher.CDS[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_CDS[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_CDS[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_CDS[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_CDS[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.CDS[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.CDS) = group
fisher.CDS.editing = lapply(fisher.CDS, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.CDS.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                  adultOnly ACnotAN    ANnotAC PNnotPC    ACnotPC    ANnotPN PNnotAN
#bothagesDEG      0.10000000       1 1.00000000       1 0.10000000 0.06666667       1
#AdultDEG         0.04274402       1 1.00000000       1 0.01805986 0.06636156       1
#PrenatalDEG      1.00000000       1 1.00000000       1 1.00000000 0.49090909       1
#bothfractionsDEG 0.30769231       1 0.33333333       1 0.40000000 0.35294118       1
#CytosolDEG       0.47058824       1 0.04761905       1 1.00000000 0.10521739       1
#NucleusDEG       0.23529412       1 1.00000000       1 0.31578947 0.12000000       1
fisher.CDS.props = lapply(fisher.CDS, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



#showing:
# that genes w editing site in a group affect expression of that gene in that group
# # of editing sites per gene
# annotation of distribution of sites per group
# gene set enrichment of edited genes per group
# mechanisms behind results
#questions:
# is nuclear editing sequestering transcripts
#are nuclear and cytoplasmic editing serving different purposes, or reflecting different strategies
# is this affected by age?

## Are cytosolic/nuclear/adult/prenatal-specific editing sites enriched for DEG Fraction/Age overall, and broken down by group?
# retained or exported DEG and presence or absence of cytosolic-specific editing site
# sites within retained or exported DEG and cytosolic-specific or non-specific status

### Check gene enrichment by annotation,
# retained or exported DEG and presence or absence of cytosolic-specific editing site
# Is a nuclear-specific editing site higher expressed in nucleus, for both intronic and other annotations?
# Is there a relationship between having an editing site and expression by fraction and age?
# is this affected by where in the gene the editing site falls (annotation)?
