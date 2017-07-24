## RBP Var database analysis of RNA editing sites


load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/rna_editing/data/rna_editing_results.rda")

#format editing results for easy input
lapply(editingres, head)
editingres = editingres[c(1:6,8,10:14)]
names(editingres) = gsub(".downsampled", "", names(editingres))
editingres_df = do.call(rbind, editingres)
editingresGR = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
editingresGR$conversion = paste0(editingresGR$ref, ":", editingresGR$alt)
editingresGR = editingresGR[which(editingresGR$conversion=="A:G" | editingresGR$conversion=="T:C")]
length(editingresGR)
editing.bed = reduce(editingresGR)
length(editing.bed) # 17601
editing.bed = as.data.frame(editing.bed)
editing.bed = editing.bed[,1:3]
write.table(editing.bed, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Unique.AtoI.sites.bed",
            quote = F, row.names = F, col.names = F, sep = "\t")

# Much too large for program, so focus on subset that overlaps DEGs
# Make a bed file of editing sites overlapping DEGs
# Mapping editing sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(editingresGR, geneMapGR)
editingresGR$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
editingresGR$nearestID = names(geneMapGR)[subjectHits(dA)]
editingresGR$distToGene = mcols(dA)$distance
editingresGR$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
editingres = as.data.frame(editingresGR)

DEG_editing.1 = lapply(sig.1, function(x) editingres[which(editingres$collapsedconversion=="A:G / T:C" &
                                                             editingres$nearestID %in% as.character(x$geneID)),])
x = do.call(rbind, DEG_editing.1)
x = makeGRangesFromDataFrame(x, keep.extra.columns = T)
x = reduce(x)
x = as.data.frame(x)
x = x[,1:3]
map.split <- split(x, (as.numeric(rownames(x)) - 1) %/% 35) 
for (i in 1:length(map.split)){
  write.table(map.split[[i]], file=paste0("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/Unique.AtoI.sites.DEG.FDR0.05.LFC1_",i,".bed"),
              quote = F, row.names = F, col.names = F, sep = "\t")
}
