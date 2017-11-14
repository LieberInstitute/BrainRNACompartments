

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


### write all sites in the RBP-Var2 format
editing_annogr = editing_anno[collapsedconversion=="A:G / T:C",,]
editing_annogr$start = editing_annogr$end
editing_annogr = makeGRangesFromDataFrame(editing_annogr)
editing_annogr = reduce(editing_annogr)
editing_annogr = as.data.frame(editing_annogr)
editing_annogr = paste0(editing_annogr$seqnames, ":",editing_annogr$start,"-",editing_annogr$end)
write.table(editing_annogr, quote=F, col.names= F, row.names = F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/BEDfiles/all_editingSites.txt")

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


## print file of the groups of editing sites that are unique to a group but present in all samples in a group
unique_allgr = lapply(unique_all, makeGRangesFromDataFrame)
unique_allgr = lapply(unique_allgr, reduce)
unique_alldf = lapply(unique_allgr, as.data.frame)
unique_alldf = lapply(unique_alldf, function(x) paste0(x$seqnames, ":",x$start,"-",x$end))
for (i in 1:length(unique_alldf)){
  write.table(unique_alldf[[i]], quote=F, col.names= F, row.names = F,
              file=paste0("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/BEDfiles/",names(unique_alldf)[i],"_uniqueInAllSamples_site.txt"))
}










