load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")


unique = lapply(unique_bySamp, function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique = Map(cbind, unique, geneID = lapply(unique, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))
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
unique_all = lapply(unique_bySamp_all, function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = Map(cbind, unique_all, geneID = lapply(unique_all, function(x) ifelse(as.character(x$overlappingGene)!="NA", as.character(x$overlappingGene), as.character(x$nearestID))))

aej_sets = openxlsx::read.xlsx('/Users/amanda/Dropbox/sorted_figures/new/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

## Enrichment in editing sites unique to a group

geneuniverse = as.character(na.omit(unique(editing_anno$EntrezID))) # all edited genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$EntrezGene.ID) %in% geneuniverse), ] # drop genes that are not present in the test set
aej_sets_expressed$EntrezGene.ID = as.character(aej_sets_expressed$EntrezGene.ID)
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
inGroup = lapply(unique, function(x) as.character(na.omit(unique(x$EntrezID))))
outGroup = lapply(unique, function(x) geneuniverse[!(geneuniverse %in% as.character(x$EntrezID))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$EntrezGene.ID),sum(!(inG %in% x$EntrezGene.ID)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$EntrezGene.ID), sum(!(outG %in% x$EntrezGene.ID)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)
enrich = lapply(enrich, data.frame)
enrich = rbind(do.call(rbind, Map(cbind, lapply(enrich, function(x) x["P.Value",]),Value = "P.Value")),
               do.call(rbind, Map(cbind, lapply(enrich, function(x) x["Odds Ratio",]),Value = "Odds Ratio")))
enrich = cbind(Comparison = gsub("1","",rownames(enrich)), enrich)
enrich = melt(enrich)
enrich$group = paste(enrich$Comparison, enrich$variable, sep=":")
enrich = enrich[order(enrich$Comparison),]
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Birnbaum_geneSet_enrichment_uniqueEditingSites.csv")


enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=0.05),"group"],1:4]
#     Comparison      Value          variable      value
#29      ANnotAC    P.Value      ASD.DATABASE 0.02021688
#41      ANnotAC Odds Ratio      ASD.DATABASE 0.51303840
#125     ANnotAC    P.Value Neurodegenerative 0.02313140
#137     ANnotAC Odds Ratio Neurodegenerative 0.00000000
#197     ANnotAC    P.Value      SCZ.PGC.GWAS 0.02914231
#209     ANnotAC Odds Ratio      SCZ.PGC.GWAS 0.36534121
#26  nucleusOnly    P.Value      ASD.DATABASE 0.00406509
#38  nucleusOnly Odds Ratio      ASD.DATABASE 0.48026029
#36      PNnotPC    P.Value      ASD.DATABASE 0.01684535
#48      PNnotPC Odds Ratio      ASD.DATABASE 0.45925335
enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=(0.05/120)),"group"],1:4]
# none


## Enrichment in editing sites unique to a group found in all samples of a group

inGroup = lapply(unique_all, function(x) as.character(na.omit(unique(x$EntrezID))))
outGroup = lapply(unique_all, function(x) geneuniverse[!(geneuniverse %in% as.character(x$EntrezID))])

enrich_all = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$EntrezGene.ID),sum(!(inG %in% x$EntrezGene.ID)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$EntrezGene.ID), sum(!(outG %in% x$EntrezGene.ID)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)
enrich_all = lapply(enrich_all, data.frame)
enrich_all = rbind(do.call(rbind, Map(cbind, lapply(enrich_all, function(x) x["P.Value",]),Value = "P.Value")),
               do.call(rbind, Map(cbind, lapply(enrich_all, function(x) x["Odds Ratio",]),Value = "Odds Ratio")))
enrich_all = cbind(Comparison = gsub("1","",rownames(enrich_all)), enrich_all)
enrich_all = melt(enrich_all)
enrich_all$group = paste(enrich_all$Comparison, enrich_all$variable, sep=":")
enrich_all = enrich_all[order(enrich_all$Comparison),]
write.csv(enrich_all, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Birnbaum_geneSet_enrichment_uniqueEditingSites_inAllSamps.csv")


enrich_all[enrich_all$group %in% enrich_all[which(enrich_all$Value=="P.Value" & enrich_all$value<=0.05),"group"],1:4]
#    Comparison      Value variable     value
#107    PCnotPN    P.Value      NDD  0.020434
#119    PCnotPN Odds Ratio      NDD 58.566300
enrich_all[enrich_all$group %in% enrich_all[which(enrich_all$Value=="P.Value" & enrich_all$value<=(0.05/120)),"group"],1:4]
#none

