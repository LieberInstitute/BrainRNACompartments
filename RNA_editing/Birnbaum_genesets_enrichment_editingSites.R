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

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

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
enrich = do.call(rbind, Map(cbind, Comparison = as.list(names(enrich)), lapply(enrich, function(x) 
            data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
                   row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Birnbaum_geneSet_enrichment_uniqueEditingSites.csv")
enrich = read.csv("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Birnbaum_geneSet_enrichment_uniqueEditingSites.csv")
head(enrich)

enrich[enrich$FDR<=0.05,]
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
enrich_all = do.call(rbind, Map(cbind, Comparison = as.list(names(enrich_all)), lapply(enrich_all, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich_all$FDR = p.adjust(enrich_all$P.Value, method = "fdr")
enrich_all
write.csv(enrich_all, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Birnbaum_geneSet_enrichment_uniqueEditingSites_inAllSamps.csv")

enrich_all[enrich_all$FDR<=0.05,]
# none