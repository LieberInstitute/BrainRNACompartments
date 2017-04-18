library(clusterProfiler)
require(org.Hs.eg.db)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

# Prepare significant genes in a list
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres), Crres = data.frame(Crres),Nrres = data.frame(Nrres),
               Cpres.down = data.frame(Cpres.down))

# Define universe as all genes expressed in each of the four groups
GeneUniverse = lapply(AgeList, function(x) x[which(x$baseMean > 0),])
GeneUniverse = lapply(GeneUniverse, function(x) geneMap[match(rownames(x),geneMap$gencodeID),"EntrezID"])
GeneUniverse = lapply(GeneUniverse, na.omit)
GeneUniverse = list("Cpres.Up"= GeneUniverse[["Cpres"]], "Cpres.Down"= GeneUniverse[["Cpres"]],
                    "Npres.Up"= GeneUniverse[["Npres"]], "Npres.Down"= GeneUniverse[["Npres"]],
                    "Crres.Up"= GeneUniverse[["Crres"]], "Crres.Down"= GeneUniverse[["Crres"]],
                    "Nrres.Up"= GeneUniverse[["Nrres"]], "Nrres.Down"= GeneUniverse[["Nrres"]],
                    "Cpres.Up.down"= GeneUniverse[["Cpres.down"]], "Cpres.Down.down"= GeneUniverse[["Cpres.down"]])
elementNROWS(GeneUniverse)

# Define significantly differently expressed genes
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigAgeList)
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "DownFetal"))
Sign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = unlist(DirList, recursive = F) 
sigGeneList = lapply(DirList, function(x) geneMap[match(rownames(x),geneMap$gencodeID),"EntrezID"])
sigGeneList = lapply(sigGeneList, function(x) unique(na.omit(x)))
elementNROWS(sigGeneList)
names(sigGeneList) = names = c("Cpres.Increasing", "Cpres.Decreasing", "Npres.Increasing", "Npres.Decreasing", "Crres.Increasing", "Crres.Decreasing", 
                               "Nrres.Increasing", "Nrres.Decreasing", "Cpres.Down.Increasing", "Cpres.Down.Decreasing")

# Find enriched Pathways via KEGG
keggList = mapply(function(g, bg) {
  ht=enrichKEGG(as.character(g), organism="human", universe= as.character(bg), 
                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
keggListdf = lapply(keggList, function(x) as.data.frame(x))

# Enriched Molecular Function GOs
goList_MF = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "MF", OrgDb = org.Hs.eg.db, universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_MF = lapply(goList_MF, function(x) as.data.frame(x))

# Biological Process GO enrichment
goList_BP = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "BP", OrgDb = org.Hs.eg.db, universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_BP = lapply(goList_BP, function(x) as.data.frame(x))

# Cellular Compartment GO enrichment
goList_CC = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "CC", OrgDb = org.Hs.eg.db, universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_CC = lapply(goList_CC, function(x) as.data.frame(x))

# Disease Ontology Enrichment
goList_DO = mapply(function(g, bg) {
  ht=enrichDO(as.character(g), ont = "DO", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_DO = lapply(goList_DO, function(x) as.data.frame(x))

lists = list(DO = goListdf_DO, CC = goListdf_CC, BP = goListdf_BP, MF = goListdf_MF, KEGG = keggListdf)
for (i in 1:length(lists)){
  for (j in 1:length(names)){
    if (nrow(lists[[i]][[j]])>0){
      lists[[i]][[j]] = data.frame(lists[[i]][[j]], Comparison = names[j])
    }}}
lists = lapply(lists, function(x) do.call(rbind, x))
write.csv(lists[["KEGG"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Age.GO.KEGG.csv")
write.csv(lists[["BP"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Age.GO.BP.csv")
write.csv(lists[["MF"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Age.GO.MF.csv")
write.csv(lists[["CC"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Age.GO.CC.csv")
write.csv(lists[["DO"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Age.DO.csv")

# Compare the enriched terms between four groups
# KEGG
compareKegg.Age = compareCluster(sigGeneList, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareKegg.Age.down = compareCluster(sigGeneList, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Biological Process
compareBP.Age = compareCluster(sigGeneList, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP.Age.down = compareCluster(sigGeneList, fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Molecular Function
compareMF.Age = compareCluster(sigGeneList, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF.Age.down = compareCluster(sigGeneList, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Cellular Component
compareCC.Age = compareCluster(sigGeneList, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC.Age.down = compareCluster(sigGeneList, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Disease Ontology
compareDO.Age = compareCluster(sigGeneList, fun="enrichDO",  ont = "DO", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO.Age,colorBy="p.adjust",  showCategory = 20, title= "Disease Ontology Enrichment")

compareDO.Age.down = compareCluster(sigGeneList, fun="enrichDO",  ont = "DO", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO.Age.down,colorBy="p.adjust",  showCategory = 20, title= "Disease Ontology Enrichment")

save(compareKegg.Age, compareBP.Age, compareMF.Age, compareCC.Age, compareDO.Age,
     compareKegg.Age.down, compareBP.Age.down, compareMF.Age.down, compareCC.Age.down, compareDO.Age.down,
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Age.kegg.GO.DO.objects.rda")