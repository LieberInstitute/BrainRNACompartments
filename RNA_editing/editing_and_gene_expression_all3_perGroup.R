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
unique_all = Map(cbind, unique_all, geneID = lapply(unique_all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
unique_all = Map(cbind, unique_all, ensID = lapply(unique_all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique_all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]),
                 baseMean = lapply(unique_all, function(x) Ipres.down[match(x$geneID, rownames(Ipres.down)),"baseMean"]), 
                 Prenatal.LFC = lapply(unique_all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "log2FoldChange"]), 
                 Prenatal.SE = lapply(unique_all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "lfcSE"]), 
                 Prenatal.padj = lapply(unique_all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "padj"]),
                 Adult.LFC = lapply(unique_all, function(x) Apres[match(x$geneID, rownames(Apres)), "log2FoldChange"]), 
                 Adult.SE = lapply(unique_all, function(x) Apres[match(x$geneID, rownames(Apres)), "lfcSE"]),
                 Adult.padj = lapply(unique_all, function(x) Apres[match(x$geneID, rownames(Apres)), "padj"]),
                 Nucleus.LFC = lapply(unique_all, function(x) Npres[match(x$geneID, rownames(Npres)), "log2FoldChange"]), 
                 Nucleus.SE = lapply(unique_all, function(x) Npres[match(x$geneID, rownames(Npres)), "lfcSE"]), 
                 Nucleus.padj = lapply(unique_all, function(x) Npres[match(x$geneID, rownames(Npres)), "padj"]),
                 Cytosol.LFC = lapply(unique_all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "log2FoldChange"]), 
                 Cytosol.SE = lapply(unique_all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "lfcSE"]),
                 Cytosol.padj = lapply(unique_all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "padj"]))
DEGnames = c(lapply(sig[1:8], function(x) as.character(x$geneID)), lapply(age.sig, function(x) as.character(x$geneID))) 
unique_all = Map(cbind, unique_all, both_retained = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$both_retained, "both_retained", "no")),
                 both_exported = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$both_exported, "both_exported", "no")),
                 Fet_retained = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
                 Ad_retained = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
                 Fet_exported = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
                 Ad_exported = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
                 ret_Ad_exp_Fet = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
                 ret_Fet_exp_Ad = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
                 interacting = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$interacting, "interacting", "no")),
                 both_decreasing = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
                 both_increasing = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$both_increasing, "both_increasing", "no")),
                 Cyt_decreasing = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
                 Nuc_decreasing = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
                 Cyt_increasing = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
                 Nuc_increasing = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
                 decr_Nuc_incr_Cyt = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
                 decr_Cyt_incr_Nuc = lapply(unique_all, function(x) ifelse(x$geneID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))
for (i in 1:length(unique_all)) { 
  unique_all[[i]]$adult_ret = ifelse((unique_all[[i]]$both_retained=="both_retained" | unique_all[[i]]$Ad_retained=="Ad_retained"), "adult_ret","no")
  unique_all[[i]]$pren_ret = ifelse((unique_all[[i]]$both_retained=="both_retained" | unique_all[[i]]$Fet_retained=="Fet_retained"), "pren_ret","no")
  unique_all[[i]]$adult_exp = ifelse((unique_all[[i]]$both_exported=="both_exported" | unique_all[[i]]$Ad_exported=="Ad_exported"), "adult_exp","no")
  unique_all[[i]]$pren_exp = ifelse((unique_all[[i]]$both_exported=="both_exported" | unique_all[[i]]$Fet_exported=="Fet_exported"), "pren_exp","no")
  unique_all[[i]]$cyt_decr = ifelse((unique_all[[i]]$both_decreasing=="both_decreasing" | unique_all[[i]]$Cyt_decreasing=="Cyt_decreasing"), "cyt_decr","no")
  unique_all[[i]]$nuc_decr = ifelse((unique_all[[i]]$both_decreasing=="both_decreasing" | unique_all[[i]]$Nuc_decreasing=="Nuc_decreasing"), "nuc_decr","no")
  unique_all[[i]]$cyt_incr = ifelse((unique_all[[i]]$both_increasing=="both_increasing" | unique_all[[i]]$Cyt_increasing=="Cyt_increasing"), "cyt_incr","no")
  unique_all[[i]]$nuc_incr = ifelse((unique_all[[i]]$both_increasing=="both_increasing" | unique_all[[i]]$Nuc_increasing=="Nuc_increasing"), "nuc_incr","no")
  

}
elementNROWS(unique_all)
lapply(unique_all, head)

all = Map(cbind, all, geneID = lapply(all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$geneID,geneMap$gencodeID),"gene_type"]),
          baseMean = lapply(all, function(x) Ipres.down[match(x$geneID, rownames(Ipres.down)),"baseMean"]), 
          Prenatal.LFC = lapply(all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "log2FoldChange"]), 
          Prenatal.SE = lapply(all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "lfcSE"]), 
          Prenatal.padj = lapply(all, function(x) Fpres.down[match(x$geneID, rownames(Fpres.down)), "padj"]),
          Adult.LFC = lapply(all, function(x) Apres[match(x$geneID, rownames(Apres)), "log2FoldChange"]), 
          Adult.SE = lapply(all, function(x) Apres[match(x$geneID, rownames(Apres)), "lfcSE"]),
          Adult.padj = lapply(all, function(x) Apres[match(x$geneID, rownames(Apres)), "padj"]),
          Nucleus.LFC = lapply(all, function(x) Npres[match(x$geneID, rownames(Npres)), "log2FoldChange"]), 
          Nucleus.SE = lapply(all, function(x) Npres[match(x$geneID, rownames(Npres)), "lfcSE"]), 
          Nucleus.padj = lapply(all, function(x) Npres[match(x$geneID, rownames(Npres)), "padj"]),
          Cytosol.LFC = lapply(all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "log2FoldChange"]), 
          Cytosol.SE = lapply(all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "lfcSE"]),
          Cytosol.padj = lapply(all, function(x) Cpres.down[match(x$geneID, rownames(Cpres.down)), "padj"]),
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
for (i in 1:length(all)) {
    all[[i]]$adult_ret = ifelse((all[[i]]$both_retained=="both_retained" | all[[i]]$Ad_retained=="Ad_retained"), "adult_ret","no")
    all[[i]]$pren_ret = ifelse((all[[i]]$both_retained=="both_retained" | all[[i]]$Fet_retained=="Fet_retained"), "pren_ret","no")
    all[[i]]$adult_exp = ifelse((all[[i]]$both_exported=="both_exported" | all[[i]]$Ad_exported=="Ad_exported"), "adult_exp","no")
    all[[i]]$pren_exp = ifelse((all[[i]]$both_exported=="both_exported" | all[[i]]$Fet_exported=="Fet_exported"), "pren_exp","no")
    all[[i]]$cyt_decr = ifelse((all[[i]]$both_decreasing=="both_decreasing" | all[[i]]$Cyt_decreasing=="Cyt_decreasing"), "cyt_decr","no")
    all[[i]]$nuc_decr = ifelse((all[[i]]$both_decreasing=="both_decreasing" | all[[i]]$Nuc_decreasing=="Nuc_decreasing"), "nuc_decr","no")
    all[[i]]$cyt_incr = ifelse((all[[i]]$both_increasing=="both_increasing" | all[[i]]$Cyt_increasing=="Cyt_increasing"), "cyt_incr","no")
    all[[i]]$nuc_incr = ifelse((all[[i]]$both_increasing=="both_increasing" | all[[i]]$Nuc_increasing=="Nuc_increasing"), "nuc_incr","no")
}


elementNROWS(lapply(unique_all, function(x) unique(x$editingID)))
#adultOnly prenatalOnly      ANnotAC      ACnotAN      ANnotPN      PNnotAN      ACnotPC      PCnotAC      PCnotPN      PNnotPC 
#130           89          159           30          554          286          273          261           12           65

### How many sites are there per gene, and is this different than the distribution of sites that aren't specific to a group and in all samples in a group?

## Number of sites per gene
numsites_bygene = lapply(unique_all, function(x) x[,length(unique(editingID)),by="geneID"])
numsites_bygene_all = lapply(all, function(x) x[,length(unique(editingID)),by="geneID"])
numsites_bygene = lapply(numsites_bygene,as.data.frame)
numsites_bygene_all = lapply(numsites_bygene_all, as.data.frame)


## Compare to the number of sites per gene in sites present in each compartment but not unique, and the total number of sites in the same genes as those of the unique sites
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
t.num.bygene = t.editedgenes.numsites = list()
for (i in 1:length(numsites_bygene)){
  t.num.bygene[[i]] = t.test(x = numsites_bygene[[group[i]]][,"V1"], y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"geneID"] %in% numsites_bygene[[group[i]]][,"geneID"]),"V1"])
  t.editedgenes.numsites[[i]] = t.test(x = numsites_bygene_all[[allgroup[i]]][which(numsites_bygene_all[[allgroup[i]]][,"geneID"] %in% numsites_bygene[[group[i]]][,"geneID"]),"V1"], 
                                       y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"geneID"] %in% numsites_bygene[[group[i]]][,"geneID"]),"V1"])
}
names(t.num.bygene) = names(t.editedgenes.numsites) = group
numsites_bygene = numsites_bygene[match(group, names(numsites_bygene))]
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(t.num.bygene, function(x) data.frame(Tstat = x$statistic, pval = x$p.value, confInt1 = x$conf.int[1], confInt2 = x$conf.int[2], 
                                                                                                   estMeansX = x$estimate[1], estMeansY = x$estimate[2], row.names = NULL)),
                        max = lapply(numsites_bygene, function(x) max(x$V1)), median = lapply(numsites_bygene, function(x) median(x$V1)), sd = lapply(numsites_bygene, function(x) sd(x$V1))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perGene_inUniqueSiteGenes_vsNonUniqueSiteGenes.csv")

df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(t.editedgenes.numsites, function(x) data.frame(Tstat = x$statistic, pval = x$p.value, confInt1 = x$conf.int[1], confInt2 = x$conf.int[2], 
                                                                                                             estMeansX = x$estimate[1], estMeansY = x$estimate[2], row.names = NULL)),
                        max = mapply(function(x,y) max(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"geneID"] %in% numsites_bygene[[x]][,"geneID"]),"V1"]), group, allgroup),
                        median = mapply(function(x,y) median(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"geneID"] %in% numsites_bygene[[x]][,"geneID"]),"V1"]), group, allgroup),
                        sd = mapply(function(x,y) sd(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"geneID"] %in% numsites_bygene[[x]][,"geneID"]),"V1"]), group, allgroup)))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perGene_inUniqueSiteGenes_vsNonUniqueSiteGenes_notJustUniqueSites.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perGene_inUniqueSiteGenes_vsNonUniqueSiteGenes_notJustUniqueSites.csv")
df[df$FDR<=0.05,]


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
numsites_byexon = numsites_byexon[match(group, names(numsites_byexon))]
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(t.num.byexon, function(x) data.frame(Tstat = x$statistic, pval = x$p.value, confInt1 = x$conf.int[1], confInt2 = x$conf.int[2], 
                                                                                                   estMeansX = x$estimate[1], estMeansY = x$estimate[2], row.names = NULL)),
                        max = lapply(numsites_byexon, function(x) max(x$V1)), median = lapply(numsites_byexon, function(x) median(x$V1)), sd = lapply(numsites_byexon, function(x) sd(x$V1))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perExon_inUniqueSiteExons_vsNonUniqueSiteExons.csv")

df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(t.editedexons.numsites, function(x) data.frame(Tstat = x$statistic, pval = x$p.value, confInt1 = x$conf.int[1], confInt2 = x$conf.int[2], 
                                                                                                             estMeansX = x$estimate[1], estMeansY = x$estimate[2], row.names = NULL)),
                        max = mapply(function(x,y) max(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                        median = mapply(function(x,y) median(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                        sd = mapply(function(x,y) sd(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup)))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perExon_inUniqueSiteExons_vsNonUniqueSiteExons_notJustUniqueSites.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Ttest_results_numEditingSites_perExon_inUniqueSiteExons_vsNonUniqueSiteExons_notJustUniqueSites.csv")
df[df$FDR<=0.05,]


### Compare the proportion of editing sites in each annotation in sites unique to a group found in all samples of that group and those that aren’t
# Using whether the site overlaps an annotation group at all rather than the CDS -> UTR -> Intron -> Intergenic heirarchy

anno.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
unique = lapply(unique_bySamp[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
an = c("CDS","Intron","5'UTR","3'UTR","NA:NA:NA:NA")
for (i in 1:length(unique_all)){
  for (j in 1:length(an)) {
    anno.site[[i]][[j]] = data.frame(inAll = c(nrow(unique_all[[group[i]]][grep(an[j], anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep(an[j], anno),list(unique(editingID)),])),
                                     notInAll = c(nrow(unique[[group[i]]][grep(an[j], anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep(an[j], anno),list(unique(editingID)),]),
                                                  nrow(unique[[group[i]]][-grep(an[j], anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep(an[j], anno),list(unique(editingID)),])), 
                                     row.names = c("inAnno","notInAnno")) }
  names(anno.site[[i]]) = c("CDS","Intron","5'UTR","3'UTR","Intergenic")
}
names(anno.site) = group
fisher.anno.site = lapply(anno.site, lapply, fisher.test)
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.anno.site, function(x) 
     do.call(rbind, Map(cbind, Annotation = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_annotation_enrichment_UniqueSitesinAll_vsUniqueSitesNotInAll.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_annotation_enrichment_UniqueSitesinAll_vsUniqueSitesNotInAll.csv")
df[df$FDR<=0.05,] 



### Characterize the overlap with differentially expressed genes

## Are editing sites specific to a group and shared by all samples in a group more or less likely to fall in an exported or retained gene by fraction than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
frac.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
comp = list(c("both_retained", "both_exported"),c("adult_ret","adult_exp"), c("pren_ret","pren_exp"),"both_retained","adult_ret", "pren_ret","both_exported","adult_exp","pren_exp")
for (i in 1:length(unique_all)){
  ua = as.data.frame(unique_all[[group[i]]])
  al = as.data.frame(all[[allgroup[i]]])
  for (j in 1:3) {
    frac.site[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"editingID"])),length(unique(al[al[,colnames(al)==comp[[j]][1]]==comp[[j]][1],"editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"editingID"]))),
                                     c(length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"editingID"])),length(unique(al[al[,colnames(al)==comp[[j]][2]]==comp[[j]][2],"editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"editingID"]))), row.names = c("unique","notUnique")) }
  for (j in 4:9) {
    frac.site[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]]]==comp[[j]],"editingID"])),length(unique(al[al[,colnames(al)==comp[[j]]]==comp[[j]],"editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]]]==comp[[j]],"editingID"]))),
                                     c(length(unique(ua[ua[,colnames(ua)==comp[[j]]]=="no","editingID"])),length(unique(al[al[,colnames(al)==comp[[j]]]=="no","editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]]]=="no","editingID"]))), row.names = c("unique","notUnique")) }
  names(frac.site[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG","bothagesRetained","adultRetained","prenatalRetained","bothagesExported","adultExported","prenatalExported")
}
names(frac.site) = group
fisher.frac.site = lapply(frac.site, lapply, fisher.test)
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.frac.site, function(x) do.call(rbind, Map(cbind, Comparison = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editingsite_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_FracDEGs.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editingsite_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_FracDEGs.csv")
df[df$FDR<=0.05,]


# number of significantly retained or exported genes that either contain an editing site specific to a group and in all samples, or don't

frac.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
comp = list(c("both_retained", "both_exported","both_retained", "both_exported","both_retained", "both_exported"),c("adult_ret","adult_exp","both_retained", "both_exported","Ad_retained", "Ad_exported"), 
            c("pren_ret","pren_exp","both_retained", "both_exported","Fet_retained", "Fet_exported"),
            c("both_retained", "both_retained", "both_retained"),c("adult_ret","both_retained","Ad_retained"), c("pren_ret","both_retained","Fet_retained"),
            c("both_exported","both_exported","both_exported"), c("adult_exp","both_exported","Ad_exported"),c("pren_exp","both_exported","Fet_exported"))
for (i in 1:length(unique_all)){
  ua = as.data.frame(unique_all[[group[i]]])
  al = as.data.frame(all[[allgroup[i]]])
  for (j in 1:3) {
    frac.gene[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"])),length(unique(c(sig[[comp[[j]][3]]]$geneID,sig[[comp[[j]][5]]]$geneID)))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"]))),
                                     c(length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"geneID"])),length(unique(c(sig[[comp[[j]][4]]]$geneID,sig[[comp[[j]][6]]]$geneID)))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"geneID"])))) }
  for (j in 4:9) {
    frac.gene[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"])),length(unique(c(sig[[comp[[j]][2]]]$geneID,sig[[comp[[j]][3]]]$geneID)))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"]))),
                                     c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]=="no","geneID"])),
                                       nrow(geneCounts)-length(unique(c(sig[[comp[[j]][2]]]$geneID,sig[[comp[[j]][3]]]$geneID)))-length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]=="no","geneID"])))) }
  names(frac.gene[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG","bothagesRetained","adultRetained","prenatalRetained","bothagesExported","adultExported","prenatalExported")
}
names(frac.gene) = group
fisher.frac.gene = lapply(frac.gene, lapply, fisher.test)
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.frac.gene, function(x) do.call(rbind, Map(cbind, Comparison = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedgene_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_FracDEGs.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedgene_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_FracDEGs.csv")
df[df$FDR<=0.05,]



## Are editing sites specific to a group and shared by all samples in a group more or less likely to fall in an increasing or decreasing genes by age than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
age.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
comp = list(c("both_decreasing", "both_increasing"),c("cyt_decr","cyt_incr"), c("nuc_decr","nuc_incr"),"both_decreasing","cyt_decr", "nuc_decr","both_increasing","cyt_incr","nuc_incr")
for (i in 1:length(unique_all)){
  ua = as.data.frame(unique_all[[group[i]]])
  al = as.data.frame(all[[allgroup[i]]])
  for (j in 1:3) {
    age.site[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"editingID"])),length(unique(al[al[,colnames(al)==comp[[j]][1]]==comp[[j]][1],"editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"editingID"]))),
                                     c(length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"editingID"])),length(unique(al[al[,colnames(al)==comp[[j]][2]]==comp[[j]][2],"editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"editingID"]))), row.names = c("unique","notUnique")) }
  for (j in 4:9) {
    age.site[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]]]==comp[[j]],"editingID"])),length(unique(al[al[,colnames(al)==comp[[j]]]==comp[[j]],"editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]]]==comp[[j]],"editingID"]))),
                                     c(length(unique(ua[ua[,colnames(ua)==comp[[j]]]=="no","editingID"])),length(unique(al[al[,colnames(al)==comp[[j]]]=="no","editingID"]))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]]]=="no","editingID"]))), row.names = c("unique","notUnique")) }
  names(age.site[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG","bothfractionsDecreasing","CytosolDecreasing","NucleusDecreasing","bothfractionsIncreasing","CytosolIncreasing","NucleusIncreasing")
}
names(age.site) = group
fisher.age.site = lapply(age.site, lapply, fisher.test)
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.age.site, function(x) do.call(rbind, Map(cbind, Comparison = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editingsite_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_AgeDEGs.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editingsite_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_AgeDEGs.csv")
df[df$FDR<=0.05,]


# number of significantly increasing or decreasing genes by age that either contain an editing site specific to a group and in all samples, or don't
age.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
comp = list(c("both_increasing", "both_decreasing","both_increasing", "both_decreasing","both_increasing", "both_decreasing"),
            c("cyt_incr","cyt_decr","both_increasing", "both_decreasing","Cyt_increasing", "Cyt_decreasing"), 
            c("nuc_incr","nuc_decr","both_increasing", "both_decreasing","Nuc_increasing", "Nuc_decreasing"),
            c("both_increasing", "both_increasing", "both_increasing"),c("cyt_incr","both_increasing","Cyt_increasing"), c("nuc_incr","both_increasing","Nuc_increasing"),
            c("both_decreasing","both_decreasing","both_decreasing"), c("cyt_decr","both_decreasing","Cyt_decreasing"),c("nuc_decr","both_decreasing","Nuc_decreasing"))
for (i in 1:length(unique_all)){
  ua = as.data.frame(unique_all[[group[i]]])
  al = as.data.frame(all[[allgroup[i]]])
  for (j in 1:3) {
    age.gene[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"])),length(unique(c(age.sig[[comp[[j]][3]]]$geneID,age.sig[[comp[[j]][5]]]$geneID)))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"]))),
                                    c(length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"geneID"])),length(unique(c(age.sig[[comp[[j]][4]]]$geneID,age.sig[[comp[[j]][6]]]$geneID)))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][2]]==comp[[j]][2],"geneID"])))) }
  for (j in 4:9) {
    age.gene[[i]][[j]] = data.frame(c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"])),length(unique(c(age.sig[[comp[[j]][2]]]$geneID,age.sig[[comp[[j]][3]]]$geneID)))-
                                         length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]==comp[[j]][1],"geneID"]))),
                                    c(length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]=="no","geneID"])),
                                       nrow(geneCounts)-length(unique(c(age.sig[[comp[[j]][2]]]$geneID,age.sig[[comp[[j]][3]]]$geneID)))-length(unique(ua[ua[,colnames(ua)==comp[[j]][1]]=="no","geneID"])))) }
  names(age.gene[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG","bothfractionsIncreasing","CytosolIncreasing","NucleusIncreasing","bothfractionsDecreasing","CytosolDecreasing","NucleusDecreasing")
}
names(age.gene) = group
fisher.age.gene = lapply(age.gene, lapply, fisher.test)
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.age.gene, function(x) do.call(rbind, Map(cbind, Comparison = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedgene_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_AgeDEGs.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedgene_enrichment_UniqueSitesinAll_vsAllOthers_inOrout_AgeDEGs.csv")
df[df$FDR<=0.05,]



### Does editing rate correlate with gene expression in the group the editing site appears?
## correlate LFC with editing rate in editing sites unique to a group

corr.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
unique_all_df = lapply(unique_all, as.data.frame)
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/correlation_degLFC_and_editingRate_inUniqueSites_inAllSamps.pdf")
for (i in 1:length(unique_all)){
  corr.site[[i]] = list(Adult = cor.test(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Adult.LFC"], use = "complete.obs"),
                        Prenatal = cor.test(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Prenatal.LFC"], use = "complete.obs"),
                        Cytosol = cor.test(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Cytosol.LFC"], use = "complete.obs"),
                        Nucleus = cor.test(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Nucleus.LFC"], use = "complete.obs"))
  print(ggplot(unique_all_df[[group[i]]], aes(x=rate, y = Adult.LFC)) + geom_point() + geom_smooth(method='lm') + 
          ylab("log2 fold change") + xlab("editing rate") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20)) +
          ggtitle(paste0("Editing Rate (", group[i],") vs. Expression\nBy Fraction in Adult")))
  print(ggplot(unique_all_df[[group[i]]], aes(x=rate, y = Prenatal.LFC)) + geom_point() + geom_smooth(method='lm') + 
          ylab("log2 fold change") + xlab("editing rate") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20)) +
          ggtitle(paste0("Editing Rate (", group[i],") vs. Expression\nBy Fraction in Prenatal")))
  print(ggplot(unique_all_df[[group[i]]], aes(x=rate, y = Cytosol.LFC)) + geom_point() + geom_smooth(method='lm') + 
          ylab("log2 fold change") + xlab("editing rate") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20)) +
          ggtitle(paste0("Editing Rate (", group[i],") vs. Expression\nBy Age in Cytoplasm")))
  print(ggplot(unique_all_df[[group[i]]], aes(x=rate, y = Nucleus.LFC)) + geom_point() + geom_smooth(method='lm') + 
          ylab("log2 fold change") + xlab("editing rate") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20)) +
          ggtitle(paste0("Editing Rate (", group[i],") vs. Expression\nBy Age in Nucleus")))
}
dev.off()
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(corr.site, function(x) 
     cbind(LFC = names(x), do.call(rbind, lapply(x, function(y) 
       data.frame(Tstat = y$statistic, pval = y$p.value, cor = y$estimate)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
row.names(df) = NULL
write.csv(df, quote=F,file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/correlation_degLFC_and_editingRate_inUniqueSites_inAllSamps.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/correlation_degLFC_and_editingRate_inUniqueSites_inAllSamps.csv")
df[df$FDR<=0.05 & abs(df$cor)>0.3,]



### Is gene expression greater in the compartment/age exhibiting the editing site than in the compared group?
## t test of LFC between unique sites in both groups

t.site = list(list())
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.site[[i]] = list(adult = t.test(x = unique_all_df[[comps[[i]][1]]][,"Adult.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Adult.LFC"]),
                          prenatal = t.test(x = unique_all_df[[comps[[i]][1]]][,"Prenatal.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Prenatal.LFC"]),
                          cytosol = t.test(x = unique_all_df[[comps[[i]][1]]][,"Cytosol.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Cytosol.LFC"]),
                          nucleus = t.test(x = unique_all_df[[comps[[i]][1]]][,"Nucleus.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Nucleus.LFC"]))
}  
names(t.site) = names(t.age.site) = names(comps)
df = do.call(rbind, Map(cbind, comparison = as.list(names(comps)), group1 = lapply(comps, function(x) x[1]), group2 = lapply(comps, function(x) x[2]),
                        lapply(t.site, function(x) do.call(rbind, Map(cbind, LFC = as.list(names(x)), lapply(x, function(y) 
                    data.frame(Tstat = y$statistic, pval = y$p.value, estMean1 = y$estimate[1], estMean2 = y$estimate[2])))))))
df$FDR = p.adjust(df$pval, method = "fdr")
row.names(df) = NULL
write.csv(df, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)
df[df$FDR<=0.05,]

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/t.test.of.LFC.between.unique.editing.sites.shared.by.all.in.group.pdf", width = 10,height = 5)
for (i in 1:length(comps)){
  tmp = rbind(cbind(comp = comps[[i]][1], melt(unique_all_df[[comps[[i]][1]]][,colnames(unique_all_df[[comps[[i]][1]]]) %in% c("Adult.LFC","Prenatal.LFC","Cytosol.LFC","Nucleus.LFC")]),
                    melt(unique_all_df[[comps[[i]][1]]][,colnames(unique_all_df[[comps[[i]][1]]]) %in% c("Adult.padj","Prenatal.padj","Cytosol.padj","Nucleus.padj")])),
              cbind(comp = comps[[i]][2], melt(unique_all_df[[comps[[i]][2]]][,colnames(unique_all_df[[comps[[i]][2]]]) %in% c("Adult.LFC","Prenatal.LFC","Cytosol.LFC","Nucleus.LFC")]),
                    melt(unique_all_df[[comps[[i]][2]]][,colnames(unique_all_df[[comps[[i]][2]]]) %in% c("Adult.padj","Prenatal.padj","Cytosol.padj","Nucleus.padj")])))
  tmp = tmp[,-4]
  colnames(tmp) = c("comp","Group","LFC","padj")
  tmp$sig = ifelse(tmp$padj<=0.05,"FDR < 0.05", "FDR > 0.05")
  tmp[is.na(tmp$padj),"sig"] = "FDR > 0.05"
  tmp$Group = gsub("Adult.LFC", "LFC by Fraction\nin Adult", tmp$Group)
  tmp$Group = gsub("Prenatal.LFC", "LFC by Fraction\nin Prenatal", tmp$Group)
  tmp$Group = gsub("Cytosol.LFC", "LFC by Age\n in Cytoplasm", tmp$Group)
  tmp$Group = gsub("Nucleus.LFC", "LFC by Age\n in Nucleus", tmp$Group)
  tmp$comp = gsub("adultOnly", "Adult\nOnly", tmp$comp)
  tmp$comp = gsub("prenatalOnly", "Prenatal\nOnly", tmp$comp)
  tmp$comp = gsub("ANnotPN", "AN\nnot PN", tmp$comp)
  tmp$comp = gsub("PNnotAN", "PN\nnot AN", tmp$comp)
  tmp$comp = gsub("ACnotPC", "AC\nnot PC", tmp$comp)
  tmp$comp = gsub("PCnotPN", "PC\nnot PN", tmp$comp)
  tmp$comp = gsub("PNnotPC", "PN\nnot PC", tmp$comp)
  tmp$comp = gsub("ACnotAN", "AC\nnot AN", tmp$comp)
  tmp$comp = gsub("ANnotAC", "AN\nnot AC", tmp$comp)
  tmp$comp = gsub("PCnotAC", "PC\nnot AC", tmp$comp)
  print(ggplot(tmp, aes(x=comp, y = LFC, fill = sig)) + geom_boxplot() + facet_grid(. ~ Group) +
          geom_hline(yintercept=0, linetype="dotted") + ylab("log2 fold change") + xlab("") + 
          theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank()) +
          ggtitle("Gene Expression by Editing Group") + scale_fill_manual(values=c("red3","gray47")))
  print(ggplot(tmp, aes(x=comp, y = LFC)) + geom_boxplot() + geom_jitter() + facet_grid(. ~ Group) +
          geom_hline(yintercept=0, linetype="dotted") +
          ylab("log2 fold change") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20)) +
          ggtitle("Gene Expression by Editing Group"))
  
}
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/fig4D_E_F.LFC.between.unique.editing.sites.shared.by.all.in.group.pdf", width = 7,height = 4)
tmp = rbind(cbind(comp = "adultOnly", melt(unique_all_df$adultOnly[,colnames(unique_all_df$adultOnly) %in% c("Cytosol.LFC","Nucleus.LFC")]),
                  melt(unique_all_df$adultOnly[,colnames(unique_all_df$adultOnly) %in% c("Cytosol.padj","Nucleus.padj")])),
            cbind(comp = "prenatalOnly", melt(unique_all_df$prenatalOnly[,colnames(unique_all_df$prenatalOnly) %in% c("Cytosol.LFC","Nucleus.LFC")]),
                  melt(unique_all_df$prenatalOnly[,colnames(unique_all_df$prenatalOnly) %in% c("Cytosol.padj","Nucleus.padj")])))
tmp = tmp[,-4]
colnames(tmp) = c("comp","Group","LFC","padj")
tmp$sig = ifelse(tmp$padj<=0.05,"FDR < 0.05", "FDR > 0.05")
tmp[is.na(tmp$padj),"sig"] = "FDR > 0.05"
tmp$Group = gsub("Cytosol.LFC", "LFC by Age\n in Cytoplasm", tmp$Group)
tmp$Group = gsub("Nucleus.LFC", "LFC by Age\n in Nucleus", tmp$Group)
tmp$comp = gsub("adultOnly", "Adult\nOnly", tmp$comp)
tmp$comp = gsub("prenatalOnly", "Prenatal\nOnly", tmp$comp)
ggplot(tmp, aes(x=comp, y = LFC, fill = sig)) + geom_boxplot() + facet_grid(. ~ Group) +
       geom_hline(yintercept=0, linetype="dotted") + ylab("log2 fold change") + xlab("") + 
       theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank()) +
       scale_fill_manual(values=c("red3","gray47")) + ylim(-8,8)
tmp = rbind(cbind(comp = "ACnotAN", melt(unique_all_df$ACnotAN[,colnames(unique_all_df$ACnotAN) %in% c("Adult.LFC","Prenatal.LFC")]),
                  melt(unique_all_df$ACnotAN[,colnames(unique_all_df$ACnotAN) %in% c("Adult.padj","Prenatal.padj")])),
            cbind(comp = "ANnotAC", melt(unique_all_df$ANnotAC[,colnames(unique_all_df$ANnotAC) %in% c("Adult.LFC","Prenatal.LFC")]),
                  melt(unique_all_df$ANnotAC[,colnames(unique_all_df$ANnotAC) %in% c("Adult.padj","Prenatal.padj")])))
tmp = tmp[,-4]
colnames(tmp) = c("comp","Group","LFC","padj")
tmp$sig = ifelse(tmp$padj<=0.05,"FDR < 0.05", "FDR > 0.05")
tmp[is.na(tmp$padj),"sig"] = "FDR > 0.05"
tmp$Group = gsub("Adult.LFC", "LFC by Fraction\nin Adult", tmp$Group)
tmp$Group = gsub("Prenatal.LFC", "LFC by Fraction\nin Prenatal", tmp$Group)
tmp$comp = gsub("ACnotAN", "AC\nnot AN", tmp$comp)
tmp$comp = gsub("ANnotAC", "AN\nnot AC", tmp$comp)
ggplot(tmp, aes(x=comp, y = LFC, fill = sig)) + geom_boxplot() + facet_grid(. ~ Group) +
  geom_hline(yintercept=0, linetype="dotted") + ylab("log2 fold change") + xlab("") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank()) +
  scale_fill_manual(values=c("red3","gray47")) + ylim(-2.5,2.5)
tmp = rbind(cbind(comp = "PCnotPN", melt(unique_all_df$PCnotPN[,colnames(unique_all_df$PCnotPN) %in% c("Adult.LFC","Prenatal.LFC")]),
                  melt(unique_all_df$PCnotPN[,colnames(unique_all_df$PCnotPN) %in% c("Adult.padj","Prenatal.padj")])),
            cbind(comp = "PNnotPC", melt(unique_all_df$PNnotPC[,colnames(unique_all_df$PNnotPC) %in% c("Adult.LFC","Prenatal.LFC")]),
                  melt(unique_all_df$PNnotPC[,colnames(unique_all_df$PNnotPC) %in% c("Adult.padj","Prenatal.padj")])))
tmp = tmp[,-4]
colnames(tmp) = c("comp","Group","LFC","padj")
tmp$sig = ifelse(tmp$padj<=0.05,"FDR < 0.05", "FDR > 0.05")
tmp[is.na(tmp$padj),"sig"] = "FDR > 0.05"
tmp$Group = gsub("Adult.LFC", "LFC by Fraction\nin Adult", tmp$Group)
tmp$Group = gsub("Prenatal.LFC", "LFC by Fraction\nin Prenatal", tmp$Group)
tmp$comp = gsub("PCnotPN", "PC\nnot PN", tmp$comp)
tmp$comp = gsub("PNnotPC", "PN\nnot PC", tmp$comp)
ggplot(tmp, aes(x=comp, y = LFC, fill = sig)) + geom_boxplot() + facet_grid(. ~ Group) +
  geom_hline(yintercept=0, linetype="dotted") + ylab("log2 fold change") + xlab("") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank()) +
  scale_fill_manual(values=c("red3","gray47")) + ylim(-2.5,2.5)
dev.off()


### In intronic editing sites, is the IR ratio greater in the compartment exhibiting the site than the compared group?
## t test of IR ratio of introns with editing sites by group (annotation of site == "Intron", not just overlapping an intron)

names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
  IRres[[i]][,"Chr"] = paste0("chr", IRres[[i]][,"Chr"])
}
names(IRres) = shortenedNames
IRres_gr = lapply(IRres, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr",start.field = "Start",end.field = "End",strand.field = "Direction",keep.extra.columns = T))
editing_anno_gr = makeGRangesFromDataFrame(cbind(editing_anno, geneID = ifelse(as.character(editing_anno$overlappingGene)!="NA", 
                                                                               as.character(editing_anno$overlappingGene), as.character(editing_anno$nearestID))), keep.extra.columns=T)
ov = lapply(IRres_gr, function(x) findOverlaps(x, editing_anno_gr))
togintron = do.call(rbind,mapply(function(oo, IR, n) cbind(editing_anno[subjectHits(oo),], IRratio = IR$IRratio[queryHits(oo)], sampleIDintron = n, 
                                             intronID = paste0(IR$Chr[queryHits(oo)],":",IR$Start[queryHits(oo)],"-",
                                                               IR$End[queryHits(oo)],":",IR$Direction[queryHits(oo)])), ov, IRres, names(IRres), SIMPLIFY = F))
IRratio_editing = lapply(unique_all_df, function(x) togintron[which(togintron$editingID %in% x$editingID & togintron$annotation=="Intron"),])
rbind(number.introns = unlist(lapply(IRratio_editing, function(x) length(unique(x$intronID)))),
      number.sites = unlist(lapply(IRratio_editing, function(x) length(unique(x$editingID)))),
      num.total.sites = unlist(lapply(unique_all_df, function(x) length(unique(x[x$annotation=="Intron","editingID"]))))) 
#                adultOnly prenatalOnly ANnotAC ACnotAN ANnotPN PNnotAN ACnotPC PCnotAC PCnotPN PNnotPC
#number.introns          5           17      43       3      64      77      15      37       1      34
#number.sites           11           15      39       2      81      61      23      39       1      25
#num.total.sites        30           23      76       4     184      97      57      73       2      35

IRratio_editing = lapply(IRratio_editing, data.table)
IRratio_editing = lapply(IRratio_editing, function(x) x[,unique(IRratio), by = c("intronID","sampleIDintron","editingID")])
for (i in 1:length(IRratio_editing)) {
  IRratio_editing[[i]]$Fraction = ifelse(IRratio_editing[[i]]$sampleIDintron %in% c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1"), "Cytoplasm","Nucleus") 
  IRratio_editing[[i]]$Age = ifelse(IRratio_editing[[i]]$sampleIDintron %in% c("Br1113C1","Br1113N1","Br2046C","Br2046N","Br2074C","Br2074N"),"Adult","Prenatal")
}

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/IRratio.between.unique.editing.sites.shared.by.all.in.group.pdf", width = 6,height = 4)
for (i in 1:length(comps)){
  if (nrow(IRratio_editing[[comps[[i]][1]]])>0 & nrow(IRratio_editing[[comps[[i]][2]]])>0) {
  tmp = rbind(cbind(comp = comps[[i]][1], IRratio_editing[[comps[[i]][1]]]),cbind(comp = comps[[i]][2], IRratio_editing[[comps[[i]][2]]]))
  } else if (nrow(IRratio_editing[[comps[[i]][1]]])>0 & nrow(IRratio_editing[[comps[[i]][2]]])<1) {
    tmp = cbind(comp = comps[[i]][1], IRratio_editing[[comps[[i]][1]]])
  } else if (nrow(IRratio_editing[[comps[[i]][1]]])<1 & nrow(IRratio_editing[[comps[[i]][2]]])>0) {
    tmp = cbind(comp = comps[[i]][2], IRratio_editing[[comps[[i]][2]]])
  } 
  tmp$comp = gsub("adultOnly", "Adult\nOnly", tmp$comp)
  tmp$comp = gsub("prenatalOnly", "Prenatal\nOnly", tmp$comp)
  tmp$comp = gsub("ANnotPN", "AN\nnot PN", tmp$comp)
  tmp$comp = gsub("PNnotAN", "PN\nnot AN", tmp$comp)
  tmp$comp = gsub("ACnotPC", "AC\nnot PC", tmp$comp)
  tmp$comp = gsub("PCnotPN", "PC\nnot PN", tmp$comp)
  tmp$comp = gsub("PNnotPC", "PN\nnot PC", tmp$comp)
  tmp$comp = gsub("ACnotAN", "AC\nnot AN", tmp$comp)
  tmp$comp = gsub("ANnotAC", "AN\nnot AC", tmp$comp)
  tmp$comp = gsub("PCnotAC", "PC\nnot AC", tmp$comp)
  print(ggplot(tmp, aes(x=comp, y = IRratio, fill = Fraction)) + geom_boxplot() + facet_grid(. ~ Age) +
          ylab("IR Ratio") + xlab("") + 
          theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank()) +
          ggtitle("IR Ratio by Editing Group") + scale_fill_brewer(palette="Dark2"))
}
dev.off()

t.IRratio = list()
for (i in 1:length(comps)){
  if (nrow(IRratio_editing[[comps[[i]][1]]])>0 & nrow(IRratio_editing[[comps[[i]][2]]])>0) {
    t.IRratio[[i]] = t.test(x = IRratio_editing[[comps[[i]][1]]][,"IRratio"], y = IRratio_editing[[comps[[i]][2]]][,"IRratio"])  } else if (nrow(IRratio_editing[[comps[[i]][1]]])>0 & nrow(IRratio_editing[[comps[[i]][2]]])<1) {
    } }
df = do.call(rbind, Map(cbind, comparison = as.list(names(comps)[elementNROWS(t.IRratio)>0]), group1 = lapply(comps[elementNROWS(t.IRratio)>0], function(x) x[1]), group2 = lapply(comps[elementNROWS(t.IRratio)>0], function(x) x[2]),
                        lapply(t.IRratio[elementNROWS(t.IRratio)>0], function(y) data.frame(Tstat = y$statistic, pval = y$p.value, estMean1 = y$estimate[1], estMean2 = y$estimate[2]))))
df$FDR = p.adjust(df$pval, method = "fdr")
row.names(df) = NULL
write.csv(df,quote=F, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.IRratio.between.unique.editing.sites.shared.by.all.in.group.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.IRratio.between.unique.editing.sites.shared.by.all.in.group.csv")
df[df$FDR<=0.05,]


## Are the introns that are edited expressed in the other group where they are not edited?
IRratio_editing = lapply(unique_all_df, function(x) togintron[which(togintron$editingID %in% x$editingID & togintron$annotation=="Intron"),])
IRratio_editing_dt = lapply(IRratio_editing, data.table)
IRratio_editing_dt = lapply(IRratio_editing_dt, function(x) as.data.frame(x[,list(unique(IRratio)),by=c("editingID","sampleIDintron","intronID")]))
IRratio_editing_dt = Map(cbind, IRratio_editing_dt, 
                         Fraction = lapply(IRratio_editing_dt, function(x) ifelse(x$sampleIDintron %in% c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1"),"Cytoplasm","Nucleus")),
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
for (i in 1:length(x)){ if (length(x[[i]])!=0) { x[[i]][,"SiteGroup"] = names(x)[i] }}
x = do.call(rbind,x)
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
write.csv(df, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/introns.present.inEachGroup.csv")

# Count introns with >=10% retention in a group
x = lapply(IRratio_editing_dt, function(x) x[V1>=0.1,length(unique(intronID)),by="Group"])
y = lapply(IRratio_editing_dt, function(x) x[, length(unique(intronID)), by="Group"])
y = do.call(rbind,Map(cbind,y,SiteGroup=as.list(names(y))))
y$IR = NA
for (i in 1:length(x)){ if (length(x[[i]])!=0) { x[[i]][,"SiteGroup"] = names(x)[i] }}
x = do.call(rbind,x)
class(y) = class(x) = "data.frame"
x$combo = paste0(x$Group,":",x$SiteGroup) 
y$combo = paste0(y$Group,":",y$SiteGroup)
for (i in 1:nrow(y)){
  if (y$combo[i] %in% x$combo) {
    y[i,"IR"] = x[which(x$Group==y$Group[i] & x$SiteGroup==y$SiteGroup[i]),"V1"]
  } else
    y[i,"IR"] = 0
}
df1 = data.frame(percent = c(round(y$IR/y$V1*100,2),round((y$V1-y$IR)/y$V1*100,2)), Group = rep.int(y$Group,2), SiteGroup = rep.int(y$SiteGroup,2), 
                IR = c(rep.int("Retained Intron",nrow(y)),rep.int("Not Retained Intron",nrow(y))), TotalIntrons = c(rep.int(NA, nrow(y)),y$V1))
df1$SiteGroup = factor(df1$SiteGroup, levels = c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))
df1$IR = factor(df1$IR, levels = c("Retained Intron","Not Retained Intron"))
write.csv(df1, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/introns.present.10percent.inEachGroup.csv")


# Plot intron retention ratios
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/introns_overlappingUniqueInAll3_editingSites_PercentNotExpresseded_inOtherGroups.pdf",width=18,height=9)
ggplot(df, aes(x = Group, y = percent, fill=IR)) + geom_bar(stat = "identity") +
  facet_grid(. ~ SiteGroup) +
  labs(fill="") + ylab("Percent") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("Percent Expressed Introns Containing a Unique Editing Site\nPresent In All Group Samples") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(df1, aes(x = Group, y = percent, fill=IR)) + geom_bar(stat = "identity") +
  facet_grid(. ~ SiteGroup) +
  labs(fill="") + ylab("Percent") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("Percent Introns (IR Ratio > 0.1)\nContaining a Unique Editing Site Present In All Group Samples") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()



### Edited Exon differential expression: Is the exon differentially expressed by group?

editedexons = exonCounts.down[which(rownames(exonCounts.down) %in% tog$exonID), grep("polyA", colnames(exonCounts.down))]
dim(editedexons) # 8941   12
match(rownames(pd[grep("polyA", rownames(pd)),]), colnames(editedexons))

## DESeq2 on exons, is the exon the last exon with the greatest count
exonsA = DESeqDataSetFromMatrix(countData = editedexons[,-grep("53", colnames(editedexons))], 
                                colData = pd[which(pd$Fetal == "Adult" & pd$Library =="polyA"),], design = ~ Zone)
exonsF = DESeqDataSetFromMatrix(countData = editedexons[,grep("53", colnames(editedexons))], 
                                colData = pd[which(pd$Fetal == "Prenatal" & pd$Library =="polyA"),], design = ~ Zone)
exonsC = DESeqDataSetFromMatrix(countData = editedexons[,grep("C", colnames(editedexons))], 
                                colData = pd[which(pd$Zone == "Cytosol" & pd$Library =="polyA"),], design = ~ Fetal)
exonsN = DESeqDataSetFromMatrix(countData = editedexons[,grep("N", colnames(editedexons))], 
                                colData = pd[which(pd$Zone == "Nucleus" & pd$Library =="polyA"),], design = ~ Fetal)
exonsA = DESeq(exonsA)
exonsF = DESeq(exonsF)
exonsC = DESeq(exonsC)
exonsN = DESeq(exonsN)
exonres = list(exonsA.res = results(exonsA), exonsF.res = results(exonsF), exonsC.res = results(exonsC), exonsN.res = results(exonsN))
elementNROWS(exonres)
elementNROWS(lapply(exonres, function(x) x[which(x$padj<=0.05),]))
tog = cbind(tog, A.LFC = exonres$exonsA.res[match(tog$exonID, rownames(exonres$exonsA.res)),"log2FoldChange"],
                 A.padj = exonres$exonsA.res[match(tog$exonID, rownames(exonres$exonsA.res)),"padj"],
                 P.LFC = exonres$exonsF.res[match(tog$exonID, rownames(exonres$exonsF.res)),"log2FoldChange"],
                 P.padj = exonres$exonsF.res[match(tog$exonID, rownames(exonres$exonsF.res)),"padj"],
                 C.LFC = exonres$exonsC.res[match(tog$exonID, rownames(exonres$exonsC.res)),"log2FoldChange"],
                 C.padj = exonres$exonsC.res[match(tog$exonID, rownames(exonres$exonsC.res)),"padj"],
                 N.LFC = exonres$exonsN.res[match(tog$exonID, rownames(exonres$exonsN.res)),"log2FoldChange"],
                 N.padj = exonres$exonsN.res[match(tog$exonID, rownames(exonres$exonsN.res)),"padj"])
tog = cbind(tog, A.retained = ifelse((tog$A.LFC>0 & tog$A.padj<=0.05), "retained","no"),A.exported = ifelse((tog$A.LFC<0 & tog$A.padj<=0.05), "exported","no"), 
            P.retained = ifelse((tog$P.LFC>0 & tog$P.padj<=0.05), "retained","no"),P.exported = ifelse((tog$P.LFC<0 & tog$P.padj<=0.05), "exported","no"),
            C.increasing = ifelse((tog$C.LFC<0 & tog$C.padj<=0.05), "increasing","no"),C.decreasing = ifelse((tog$C.LFC>0 & tog$C.padj<=0.05), "decreasing","no"),
            N.increasing = ifelse((tog$N.LFC<0 & tog$N.padj<=0.05), "increasing","no"),N.decreasing = ifelse((tog$N.LFC>0 & tog$N.padj<=0.05), "decreasing","no"))
save(tog, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/exons.present.inEachGroup.editing.rda")

for (i in 1:length(unique_all)) { unique_all[[i]]$start = unique_all[[i]]$end }
tog$start = tog$end
oo = lapply(unique_all, function(x) findOverlaps(makeGRangesFromDataFrame(tog[,1:5]), makeGRangesFromDataFrame(x)))
x = Map(cbind, mapply(function(ov, un) un[subjectHits(ov)], oo, unique_all, SIMPLIFY = F), 
                 lapply(oo, function(x) as.data.frame(tog)[queryHits(x), colnames(tog) %in% c("A.LFC","A.padj","P.LFC","P.padj","C.LFC","C.padj","N.LFC","N.padj")]))
unique_all_df = lapply(x, as.data.frame)

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/exons.LFC.between.unique.editing.sites.shared.by.all.in.group.pdf", width = 10,height = 5)
for (i in 1:length(comps)){
  tmp = rbind(cbind(comp = comps[[i]][1], melt(unique_all_df[[comps[[i]][1]]][,colnames(unique_all_df[[comps[[i]][1]]]) %in% c("A.LFC","P.LFC","C.LFC","N.LFC")]),
                    melt(unique_all_df[[comps[[i]][1]]][,colnames(unique_all_df[[comps[[i]][1]]]) %in% c("A.padj","P.padj","C.padj","N.padj")])),
              cbind(comp = comps[[i]][2], melt(unique_all_df[[comps[[i]][2]]][,colnames(unique_all_df[[comps[[i]][2]]]) %in% c("A.LFC","P.LFC","C.LFC","N.LFC")]),
                    melt(unique_all_df[[comps[[i]][2]]][,colnames(unique_all_df[[comps[[i]][2]]]) %in% c("A.padj","P.padj","C.padj","N.padj")])))
  tmp = tmp[,-4]
  colnames(tmp) = c("comp","Group","LFC","padj")
  tmp$sig = ifelse(tmp$padj<=0.05,"FDR < 0.05", "FDR > 0.05")
  tmp[is.na(tmp$padj),"sig"] = "FDR > 0.05"
  tmp$Group = gsub("A.LFC", "LFC by Fraction\nin Adult", tmp$Group)
  tmp$Group = gsub("P.LFC", "LFC by Fraction\nin Prenatal", tmp$Group)
  tmp$Group = gsub("C.LFC", "LFC by Age\n in Cytoplasm", tmp$Group)
  tmp$Group = gsub("N.LFC", "LFC by Age\n in Nucleus", tmp$Group)
  tmp$comp = gsub("adultOnly", "Adult\nOnly", tmp$comp)
  tmp$comp = gsub("prenatalOnly", "Prenatal\nOnly", tmp$comp)
  tmp$comp = gsub("ANnotPN", "AN\nnot PN", tmp$comp)
  tmp$comp = gsub("PNnotAN", "PN\nnot AN", tmp$comp)
  tmp$comp = gsub("ACnotPC", "AC\nnot PC", tmp$comp)
  tmp$comp = gsub("PCnotPN", "PC\nnot PN", tmp$comp)
  tmp$comp = gsub("PNnotPC", "PN\nnot PC", tmp$comp)
  tmp$comp = gsub("ACnotAN", "AC\nnot AN", tmp$comp)
  tmp$comp = gsub("ANnotAC", "AN\nnot AC", tmp$comp)
  tmp$comp = gsub("PCnotAC", "PC\nnot AC", tmp$comp)
  print(ggplot(tmp, aes(x=comp, y = LFC, fill = sig)) + geom_boxplot() + facet_grid(. ~ Group) +
          geom_hline(yintercept=0, linetype="dotted") + ylab("log2 fold change") + xlab("") + 
          theme(title = element_text(size = 20), text = element_text(size = 20), legend.title = element_blank()) +
          ggtitle("Gene Expression by Editing Group") + scale_fill_manual(values=c("red3","gray47")))
  print(ggplot(tmp, aes(x=comp, y = LFC)) + geom_boxplot() + geom_jitter() + facet_grid(. ~ Group) +
          geom_hline(yintercept=0, linetype="dotted") +
          ylab("log2 fold change") + xlab("") + theme(title = element_text(size = 20)) + theme(text = element_text(size = 20)) +
          ggtitle("Gene Expression by Editing Group"))
  
}
dev.off()



### in all exon editing sites, is the exon differentially expressed by group?

## Compare LFC and significance of groups of editing sites in all edited exons

uniqueAll_allexons = lapply(unique_all, function(x) tog[editingID %in% x$editingID,])
t.allexons.site = list(list())
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.allexons.site[[i]] = list(adult = t.test(x = uniqueAll_allexons[[comps[[i]][1]]][,"A.LFC"], y = uniqueAll_allexons[[comps[[i]][2]]][,"A.LFC"]),
                              prenatal = t.test(x = uniqueAll_allexons[[comps[[i]][1]]][,"P.LFC"], y = uniqueAll_allexons[[comps[[i]][2]]][,"P.LFC"]),
                              cytosol = t.test(x = uniqueAll_allexons[[comps[[i]][1]]][,"C.LFC"], y = uniqueAll_allexons[[comps[[i]][2]]][,"C.LFC"]),
                              nucleus = t.test(x = uniqueAll_allexons[[comps[[i]][1]]][,"N.LFC"], y = uniqueAll_allexons[[comps[[i]][2]]][,"N.LFC"]))
}
df = do.call(rbind, Map(cbind, comparison = as.list(names(comps)), group1 = lapply(comps, function(x) x[1]), group2 = lapply(comps, function(x) x[2]),
                        lapply(t.allexons.site, function(x) do.call(rbind, Map(cbind, LFC = as.list(names(x)), lapply(x, function(y) 
                          data.frame(Tstat = y$statistic, pval = y$p.value, estMean1 = y$estimate[1], estMean2 = y$estimate[2])))))))
df$FDR = p.adjust(df$pval, method = "fdr")
row.names(df) = NULL
write.csv(df, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.allexons.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.allexons.LFC.between.unique.editing.sites.shared.by.all.in.group.csv")
df[df$FDR<=0.05,]


## Correlate the LFC between age and fraction comparisons in different groups

corr.allexons.site = list(list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
for (i in 1:length(group)){
  corr.allexons.site[[i]] = list(fraction = cor.test(x = uniqueAll_allexons[[group[i]]][,"A.LFC"], y = uniqueAll_allexons[[group[i]]][,"P.LFC"], use = "complete.obs"),
                                 age = cor.test(x = uniqueAll_allexons[[group[i]]][,"C.LFC"], y = uniqueAll_allexons[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.allexons.site) = group
data.frame(lapply(corr.allexons.site, function(x) unlist(x, recursive=F)))
#         adultOnly prenatalOnly   ACnotAN   ANnotAC   PCnotPN   PNnotPC   ACnotPC   PCnotAC   ANnotPN   PNnotAN
#fraction 0.8462662    0.5131192 0.2232211 0.4769812 0.3448261 0.5970088 0.7235548 0.5832618 0.7756590 0.4907453
#age      0.9536717    0.9385681 0.9155127 0.8678858 0.9809376 0.9042264 0.9296138 0.9432402 0.9283744 0.9276533


## Compare the number of edited exons that are significantly and non-significantly up- and down-expressed per group

fisher.allexons = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  fisher.allexons[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_allexons[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_allexons[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_allexons[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_allexons[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_allexons[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_allexons[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_allexons[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_allexons[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_allexons[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_allexons[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_allexons[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_allexons[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_allexons[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.allexons[[i]]) = c("bothagesDEE","AdultDEE","PrenatalDEE","bothfractionsDEE","CytosolDEE","NucleusDEE")
}
fisher.allexons.editing = lapply(fisher.allexons, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.allexons.editing, function(x) 
  do.call(rbind, Map(cbind, exon.Cat = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedallexons_sigvsnonsig_export-retain_increase-decrease.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedallexons_sigvsnonsig_export-retain_increase-decrease.csv")
df[df$FDR<=0.05,]


## Compare being significantly DE with being unique or not

group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
all_allExons = lapply(all, function(x) tog[editingID %in% x$editingID,,])
fisher = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  fisher[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_allexons[[group[i]]][(A.exported=="exported" & P.exported=="exported"),list(unique(exonID)),]),
                            nrow(all_allExons[[allgroup[i]]][(A.exported=="exported" & P.exported=="exported"),list(unique(exonID)),])-
                              nrow(uniqueAll_allexons[[group[i]]][(A.exported=="exported" & P.exported=="exported"),list(unique(exonID)),])),
               notexported = c(nrow(uniqueAll_allexons[[group[i]]][(A.exported=="no" | P.exported=="no"),list(unique(exonID)),]),
                               nrow(all_allExons[[allgroup[i]]][(A.exported=="no" | P.exported=="no"),list(unique(exonID)),])-
                                 nrow(uniqueAll_allexons[[group[i]]][(A.exported=="no" | P.exported=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(exported = c(nrow(uniqueAll_allexons[[group[i]]][(A.exported=="exported"),list(unique(exonID)),]),
                            nrow(all_allExons[[allgroup[i]]][(A.exported=="exported"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(A.exported=="exported"),list(unique(exonID)),])),
               notexported = c(nrow(uniqueAll_allexons[[group[i]]][(A.exported=="no"),list(unique(exonID)),]),
                               nrow(all_allExons[[allgroup[i]]][(A.exported=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(A.exported=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(exported = c(nrow(uniqueAll_allexons[[group[i]]][(P.exported=="exported"),list(unique(exonID)),]),
                            nrow(all_allExons[[allgroup[i]]][(P.exported=="exported"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(P.exported=="exported"),list(unique(exonID)),])),
               notexported = c(nrow(uniqueAll_allexons[[group[i]]][(P.exported=="no"),list(unique(exonID)),]),
                               nrow(all_allExons[[allgroup[i]]][(P.exported=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(P.exported=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(retained = c(nrow(uniqueAll_allexons[[group[i]]][(A.retained=="retained" & P.retained=="retained"),list(unique(exonID)),]),
                            nrow(all_allExons[[allgroup[i]]][(A.retained=="retained" & P.retained=="retained"),list(unique(exonID)),])-
                              nrow(uniqueAll_allexons[[group[i]]][(A.retained=="retained" & P.retained=="retained"),list(unique(exonID)),])),
               notretained = c(nrow(uniqueAll_allexons[[group[i]]][(A.retained=="no" | P.retained=="no"),list(unique(exonID)),]),
                               nrow(all_allExons[[allgroup[i]]][(A.retained=="no" | P.retained=="no"),list(unique(exonID)),])-
                                 nrow(uniqueAll_allexons[[group[i]]][(A.retained=="no" | P.retained=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(retained = c(nrow(uniqueAll_allexons[[group[i]]][(A.retained=="retained"),list(unique(exonID)),]),
                            nrow(all_allExons[[allgroup[i]]][(A.retained=="retained"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(A.retained=="retained"),list(unique(exonID)),])),
               notretained = c(nrow(uniqueAll_allexons[[group[i]]][(A.retained=="no"),list(unique(exonID)),]),
                               nrow(all_allExons[[allgroup[i]]][(A.retained=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(A.retained=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(retained = c(nrow(uniqueAll_allexons[[group[i]]][(P.retained=="retained"),list(unique(exonID)),]),
                            nrow(all_allExons[[allgroup[i]]][(P.retained=="retained"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(P.retained=="retained"),list(unique(exonID)),])),
               notretained = c(nrow(uniqueAll_allexons[[group[i]]][(P.retained=="no"),list(unique(exonID)),]),
                               nrow(all_allExons[[allgroup[i]]][(P.retained=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(P.retained=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(increasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="increasing" & N.increasing=="increasing"),list(unique(exonID)),]),
                              nrow(all_allExons[[allgroup[i]]][(C.increasing=="increasing" & N.increasing=="increasing"),list(unique(exonID)),])-
                                nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="increasing" & N.increasing=="increasing"),list(unique(exonID)),])),
               notincreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="no" | N.increasing=="no"),list(unique(exonID)),]),
                                 nrow(all_allExons[[allgroup[i]]][(C.increasing=="no" | N.increasing=="no"),list(unique(exonID)),])-
                                   nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="no" | N.increasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(increasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="increasing"),list(unique(exonID)),]),
                              nrow(all_allExons[[allgroup[i]]][(C.increasing=="increasing"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="increasing"),list(unique(exonID)),])),
               notincreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="no"),list(unique(exonID)),]),
                                 nrow(all_allExons[[allgroup[i]]][(C.increasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(C.increasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(increasing = c(nrow(uniqueAll_allexons[[group[i]]][(N.increasing=="increasing"),list(unique(exonID)),]),
                              nrow(all_allExons[[allgroup[i]]][(N.increasing=="increasing"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(N.increasing=="increasing"),list(unique(exonID)),])),
               notincreasing = c(nrow(uniqueAll_allexons[[group[i]]][(N.increasing=="no"),list(unique(exonID)),]),
                                 nrow(all_allExons[[allgroup[i]]][(N.increasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(N.increasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(decreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="decreasing" & N.decreasing=="decreasing"),list(unique(exonID)),]),
                              nrow(all_allExons[[allgroup[i]]][(C.decreasing=="decreasing" & N.decreasing=="decreasing"),list(unique(exonID)),])-
                                nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="decreasing" & N.decreasing=="decreasing"),list(unique(exonID)),])),
               notdecreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="no" | N.decreasing=="no"),list(unique(exonID)),]),
                                 nrow(all_allExons[[allgroup[i]]][(C.decreasing=="no" | N.decreasing=="no"),list(unique(exonID)),])-
                                   nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="no" | N.decreasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(decreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="decreasing"),list(unique(exonID)),]),
                              nrow(all_allExons[[allgroup[i]]][(C.decreasing=="decreasing"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="decreasing"),list(unique(exonID)),])),
               notdecreasing = c(nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="no"),list(unique(exonID)),]),
                                 nrow(all_allExons[[allgroup[i]]][(C.decreasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(C.decreasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(decreasing = c(nrow(uniqueAll_allexons[[group[i]]][(N.decreasing=="decreasing"),list(unique(exonID)),]),
                              nrow(all_allExons[[allgroup[i]]][(N.decreasing=="decreasing"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(N.decreasing=="decreasing"),list(unique(exonID)),])),
               notdecreasing = c(nrow(uniqueAll_allexons[[group[i]]][(N.decreasing=="no"),list(unique(exonID)),]),
                                 nrow(all_allExons[[allgroup[i]]][(N.decreasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_allexons[[group[i]]][(N.decreasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")))
  names(fisher[[i]]) = c("bothagesDEE:Exported","AdultDEE:Exported","PrenatalDEE:Exported","bothagesDEE:Retained","AdultDEE:Retained","PrenatalDEE:Retained",
                              "bothfractionsDEE:Increasing","CytosolDEE:Increasing","NucleusDEE:Increasing","bothfractionsDEE:Decreasing","CytosolDEE:Decreasing","NucleusDEE:Decreasing")
}
fisher.editing = lapply(fisher, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.editing, function(x) 
  do.call(rbind, Map(cbind, exon.Cat = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_editedExons_sigvsnonsig_export-not_retain-not_increase-not_decrease-not.csv")
df[df$FDR<=0.05,]



### in 3'UTR editing sites, is the exon differentially expressed by group?

## Compare LFC and significance of groups of editing sites in 3'UTR

tog.3UTR = tog[annotation=="3'UTR",,]
uniqueAll_3UTR = lapply(unique_all, function(x) tog.3UTR[editingID %in% x$editingID,])
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
df = do.call(rbind, Map(cbind, comparison = as.list(names(comps)), group1 = lapply(comps, function(x) x[1]), group2 = lapply(comps, function(x) x[2]),
                        lapply(t.3UTR.site, function(x) do.call(rbind, Map(cbind, LFC = as.list(names(x)), lapply(x, function(y) 
                          data.frame(Tstat = y$statistic, pval = y$p.value, estMean1 = y$estimate[1], estMean2 = y$estimate[2])))))))
df$FDR = p.adjust(df$pval, method = "fdr")
row.names(df) = NULL
write.csv(df,file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.3UTR.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


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
#fraction 0.8382209    0.5058991 0.2414527 0.6485191 0.3552849 0.5183085 0.5853646 0.6042506 0.8046249 0.5366931
#age      0.9373044    0.9666219 0.9008494 0.9554178 0.9799071 0.8797949 0.8940680 0.9619225 0.9379955 0.9563944


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
  names(fisher.3UTR[[i]]) = c("bothagesDEE","AdultDEE","PrenatalDEE","bothfractionsDEE","CytosolDEE","NucleusDEE")
}
fisher.3UTR.editing = lapply(fisher.3UTR, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.3UTR.editing, function(x) 
  do.call(rbind, Map(cbind, exon.Cat = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_edited3UTR_sigvsnonsig_export-retain_increase-decrease.csv")
df[df$FDR<=0.05,]


## Compare being significantly DE with being unique or not

group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
all_3UTR = lapply(all, function(x) tog.3UTR[editingID %in% x$editingID,,])
fisher.3UTR = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  fisher.3UTR[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="exported" & P.exported=="exported"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(A.exported=="exported" & P.exported=="exported"),list(unique(exonID)),])-
                            nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="exported" & P.exported=="exported"),list(unique(exonID)),])),
               notexported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="no" | P.exported=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(A.exported=="no" | P.exported=="no"),list(unique(exonID)),])-
                               nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="no" | P.exported=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="exported"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(A.exported=="exported"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="exported"),list(unique(exonID)),])),
               notexported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(A.exported=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(A.exported=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(P.exported=="exported"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(P.exported=="exported"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(P.exported=="exported"),list(unique(exonID)),])),
               notexported = c(nrow(uniqueAll_3UTR[[group[i]]][(P.exported=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(P.exported=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(P.exported=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(retained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="retained" & P.retained=="retained"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(A.retained=="retained" & P.retained=="retained"),list(unique(exonID)),])-
                            nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="retained" & P.retained=="retained"),list(unique(exonID)),])),
               notretained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="no" | P.retained=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(A.retained=="no" | P.retained=="no"),list(unique(exonID)),])-
                               nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="no" | P.retained=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(retained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="retained"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(A.retained=="retained"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="retained"),list(unique(exonID)),])),
               notretained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(A.retained=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(A.retained=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(retained = c(nrow(uniqueAll_3UTR[[group[i]]][(P.retained=="retained"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(P.retained=="retained"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(P.retained=="retained"),list(unique(exonID)),])),
               notretained = c(nrow(uniqueAll_3UTR[[group[i]]][(P.retained=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(P.retained=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(P.retained=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="increasing" & N.increasing=="increasing"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(C.increasing=="increasing" & N.increasing=="increasing"),list(unique(exonID)),])-
                              nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="increasing" & N.increasing=="increasing"),list(unique(exonID)),])),
               notincreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="no" | N.increasing=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(C.increasing=="no" | N.increasing=="no"),list(unique(exonID)),])-
                                 nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="no" | N.increasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="increasing"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(C.increasing=="increasing"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="increasing"),list(unique(exonID)),])),
               notincreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(C.increasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(C.increasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.increasing=="increasing"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(N.increasing=="increasing"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(N.increasing=="increasing"),list(unique(exonID)),])),
               notincreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.increasing=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(N.increasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(N.increasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="decreasing" & N.decreasing=="decreasing"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(C.decreasing=="decreasing" & N.decreasing=="decreasing"),list(unique(exonID)),])-
                              nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="decreasing" & N.decreasing=="decreasing"),list(unique(exonID)),])),
               notdecreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="no" | N.decreasing=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(C.decreasing=="no" | N.decreasing=="no"),list(unique(exonID)),])-
                                 nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="no" | N.decreasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="decreasing"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(C.decreasing=="decreasing"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="decreasing"),list(unique(exonID)),])),
               notdecreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(C.decreasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(C.decreasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")),
    data.frame(decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.decreasing=="decreasing"),list(unique(exonID)),]),
                            nrow(all_3UTR[[allgroup[i]]][(N.decreasing=="decreasing"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(N.decreasing=="decreasing"),list(unique(exonID)),])),
               notdecreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.decreasing=="no"),list(unique(exonID)),]),
                               nrow(all_3UTR[[allgroup[i]]][(N.decreasing=="no"),list(unique(exonID)),])-nrow(uniqueAll_3UTR[[group[i]]][(N.decreasing=="no"),list(unique(exonID)),])), row.names = c("unique","notunique")))
    names(fisher.3UTR[[i]]) = c("bothagesDEE:Exported","AdultDEE:Exported","PrenatalDEE:Exported","bothagesDEE:Retained","AdultDEE:Retained","PrenatalDEE:Retained",
                                "bothfractionsDEE:Increasing","CytosolDEE:Increasing","NucleusDEE:Increasing","bothfractionsDEE:Decreasing","CytosolDEE:Decreasing","NucleusDEE:Decreasing")
}
fisher.3UTR.editing = lapply(fisher.3UTR, function(x) lapply(x, fisher.test))
df = do.call(rbind, Map(cbind, Group = as.list(group), lapply(fisher.3UTR.editing, function(x) 
  do.call(rbind, Map(cbind, exon.Cat = as.list(names(x)), lapply(x, function(y) data.frame(pval = y$p.value, OddsRatio = y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df,quote = F,file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/fisher_edited3UTR_sigvsnonsig_export-not_retain-not_increase-not_decrease-not.csv")
df[df$FDR<=0.05,]


## Are the edited 3'UTRs from the major isoform?

editedGenes = exonMap[exonMap$exonID %in% tog.3UTR$exonID,]
editedGenes = data.table(cbind(editedGenes, as.data.frame(exonCounts.down[match(editedGenes$exonID, rownames(exonCounts.down)),grep("polyA",colnames(exonCounts.down))])))
editedGenes$mean = rowMeans(as.data.frame(editedGenes)[,grep("polyA",colnames(editedGenes))])
editedGenesMax = editedGenes[editedGenes[, .I[mean == max(mean)], by=gencodeID]$V1]
tog.3UTR$dominant3UTR = ifelse(tog.3UTR$exonID %in% editedGenesMax$exonID, "Major", "Minor")
unique.3UTR = Map(cbind, lapply(unique_all, function(x) tog.3UTR[editingID %in% x$editingID,,]), name = as.list(names(unique_all)))
unique.3UTR = do.call(rbind, unique.3UTR)
unique.3UTR$name = factor(unique.3UTR$name, levels= c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/editing3UTRsites_major-minor_breakdown.pdf",width=8,height=6)
x = tog.3UTR[,length(unique(exonID)),by = c("Age","Fraction","Group","dominant3UTR")]
for (i in 1:length(unique(x$Group))){ 
  x$perc[grep(unique(x$Group)[i],x$Group)] = round(x$V1[grep(unique(x$Group)[i],x$Group)]/sum(x$V1[grep(unique(x$Group)[i],x$Group)])*100,1)}
ggplot(x, aes(x = Fraction, y = perc, fill= dominant3UTR)) + geom_bar(stat = "identity",position = position_dodge(width = 1)) +
  facet_grid(. ~ Age) +
  geom_text(aes(label = V1), vjust = -.5, position = position_dodge(width = 1)) +
  labs(fill="") +
  ylab("Percent") + ylim(0,75) + 
  xlab("") +
  ggtitle("Edited 3'UTRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
x = unique.3UTR[,length(unique(exonID)),by = c("name","dominant3UTR")]
for (i in 1:length(unique(x$name))){ 
  x$perc[grep(unique(x$name)[i],x$name)] = round(x$V1[grep(unique(x$name)[i],x$name)]/sum(x$V1[grep(unique(x$name)[i],x$name)])*100,1)}
x = rbind(x,data.frame("name"="PCnotPN","dominant3UTR"="Minor","V1" = 0, "perc" = 0))
ggplot(x, aes(x = name, y = perc, fill= dominant3UTR)) + geom_bar(stat = "identity",position = position_dodge(width = 1)) +
  geom_text(aes(label = V1), vjust = -.5, position = position_dodge(width = 1)) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Edited 3'UTRs") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()



## How many edited exons (CDS, 3'UTR, and 5'UTR) are still expressed in the fraction/age where the editing site doesn't appear?

editedExons = exonMap[exonMap$exonID %in% tog$exonID,]
editedExons = data.table(cbind(editedExons, as.data.frame(exonCounts.down[match(editedExons$exonID, rownames(exonCounts.down)),grep("polyA",colnames(exonCounts.down))])))
editedExons$meanAexpr = rowMeans(as.data.frame(editedExons)[,c(grep("Br11",colnames(editedExons)),grep("Br20",colnames(editedExons)))])
editedExons$meanPexpr = rowMeans(as.data.frame(editedExons)[,grep("53",colnames(editedExons))])
editedExons$meanCexpr = rowMeans(as.data.frame(editedExons)[,c(grep("C1_",colnames(editedExons)),grep("C_",colnames(editedExons)))])
editedExons$meanNexpr = rowMeans(as.data.frame(editedExons)[,c(grep("N1_",colnames(editedExons)),grep("N_",colnames(editedExons)))])
editedExons$meanACexpr = rowMeans(as.data.frame(editedExons)[,c("Br1113C1_polyA","Br2046C_polyA","Br2074C_polyA")])
editedExons$meanANexpr = rowMeans(as.data.frame(editedExons)[,c("Br1113N1_polyA","Br2046N_polyA","Br2074N_polyA")])
editedExons$meanPCexpr = rowMeans(as.data.frame(editedExons)[,c("Br5339C1_polyA","Br5340C1_polyA","Br5341C1_polyA")])
editedExons$meanPNexpr = rowMeans(as.data.frame(editedExons)[,c("Br5339N1_polyA","Br5340N1_polyA","Br5341N1_polyA")])
editedExons$mean = rowMeans(as.data.frame(editedExons)[,grep("polyA",colnames(editedExons))])
editedExons = cbind(tog, editedExons[match(tog$exonID, editedExons$exonID),,])

# Count exons with 0% expression in a group
df = data.table(meanExpr = c(editedExons$meanAexpr, editedExons$meanPexpr, editedExons$meanCexpr, editedExons$meanNexpr, 
                             editedExons$meanACexpr, editedExons$meanANexpr, editedExons$meanPCexpr, editedExons$meanPNexpr, editedExons$mean),
                exonID = rep.int(editedExons$exonID, 9), editingID = rep.int(editedExons$editingID, 9),
                GroupMean = c(rep.int("Adult", nrow(editedExons)),rep.int("Prenatal", nrow(editedExons)),rep.int("Cytosol", nrow(editedExons)),
                          rep.int("Nucleus", nrow(editedExons)),rep.int("Adult:Cytosol", nrow(editedExons)),rep.int("Adult:Nucleus", nrow(editedExons)),
                          rep.int("Prenatal:Cytosol", nrow(editedExons)),rep.int("Prenatal:Nucleus", nrow(editedExons)),rep.int("All", nrow(editedExons))))
df = Map(cbind, lapply(unique_all, function(x) df[editingID %in% x$editingID,,]), name = as.list(names(unique_all)))
df = do.call(rbind, df)
df$name = factor(df$name, levels= c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))

x = df[meanExpr==0, length(unique(exonID)), by=c("name","GroupMean")]
y = df[, length(unique(exonID)), by=c("name","GroupMean")]
y$Zero = NA
class(y) = class(x) = "data.frame"
x$combo = paste0(x$Group,":",x$name) 
y$combo = paste0(y$Group,":",y$name)
for (i in 1:nrow(y)){
  if (y$combo[i] %in% x$combo) {
    y[i,"Zero"] = x[which(x$name==y$name[i] & x$GroupMean==y$GroupMean[i]),"V1"]
  } else
    y[i,"Zero"] = 0
}
df = data.frame(percent = c(round(y$Zero/y$V1*100,2),round((y$V1-y$Zero)/y$V1*100,2)), Group = rep.int(y$Group,2), SiteGroup = rep.int(y$name,2), 
                IR = c(rep.int("Not Expressed",nrow(y)),rep.int("Expressed",nrow(y))), TotalExons = c(rep.int(NA, nrow(y)),y$V1))
df$SiteGroup = factor(df$SiteGroup, levels = c("adultOnly","prenatalOnly","ACnotPC","PCnotAC","ANnotPN","PNnotAN","ACnotAN","ANnotAC","PCnotPN","PNnotPC"))
write.csv(df,"./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/exons_overlappingUniqueInAll3_editingSites_PercentNotExpresseded_inOtherGroups.csv")


# Plot exon retention ratios
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/exons_overlappingUniqueInAll3_editingSites_PercentNotExpresseded_inOtherGroups.pdf",width=18,height=9)
ggplot(df[grep(":", df$Group),], aes(x = Group, y = percent, fill=IR)) + geom_bar(stat = "identity") +
  facet_grid(. ~ SiteGroup) +
  labs(fill="") +
  ylab("Percent") + 
  xlab("") + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  ggtitle("Percent Expressed Exons Containing a Unique Editing Site\nPresent In All Group Samples") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()









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
