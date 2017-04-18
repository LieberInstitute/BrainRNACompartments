library("clusterProfiler")
require("org.Hs.eg.db")

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

# Prepare significant genes in a list
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down))

# Define universe as all genes expressed in each of the four groups
GeneUniverse = lapply(FracList, function(x) x[which(x$baseMean > 0),])
elementNROWS(GeneUniverse)
GeneUniverse = lapply(GeneUniverse, function(x) geneMap[match(rownames(x),geneMap$gencodeID),"EntrezID"])
GeneUniverse = lapply(GeneUniverse, na.omit)
GeneUniverse = list("Apres.Up"= GeneUniverse[["Apres"]], "Apres.Down"= GeneUniverse[["Apres"]],
                    "Fpres.Up"= GeneUniverse[["Fpres"]], "Fpres.Down"= GeneUniverse[["Fpres"]],
                    "Arres.Up"= GeneUniverse[["Arres"]], "Arres.Down"= GeneUniverse[["Arres"]],
                    "Frres.Up"= GeneUniverse[["Frres"]], "Frres.Down"= GeneUniverse[["Frres"]],
                    "Fpres.Up.down"= GeneUniverse[["Fpres.down"]], "Fpres.Down.Down"= GeneUniverse[["Fpres.down"]])

# Define significantly differently expressed genes
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
SigList = lapply(Sign, function(x) split(x, x$Sign))
SigList = unlist(SigList, recursive = F) 
sigGeneList = lapply(SigList, function(x) geneMap[match(rownames(x),geneMap$gencodeID),"EntrezID"])
sigGeneList = lapply(sigGeneList, function(x) unique(na.omit(x)))
elementNROWS(sigGeneList)
names(sigGeneList) = names = c("Apres.Cytosolic", "Apres.Nuclear", "Fpres.Cytosolic", "Fpres.Nuclear", "Arres.Cytosolic", "Arres.Nuclear", 
          "Frres.Cytosolic", "Frres.Nuclear", "Fpres.Down.Cytosolic", "Fpres.Down.Nuclear")

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
write.csv(lists[["KEGG"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Fraction.GO.KEGG.csv")
write.csv(lists[["BP"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Fraction.GO.BP.csv")
write.csv(lists[["MF"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Fraction.GO.MF.csv")
write.csv(lists[["CC"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Fraction.GO.CC.csv")
write.csv(lists[["DO"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Fraction.DO.csv")

# Compare the enriched terms between four groups
elementNROWS(sigGeneList)
names(sigGeneList) = c("Adult\nPolyA\nCytosolic", "Adult\nPolyA\nNuclear", "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear", 
                       "Adult\nRibozero\nCytosolic", "Adult\nRibozero\nNuclear", "Fetal\nRibozero\nCytosolic", "Fetal\nRibozero\nNuclear",
                       "Fetal\nPolyA\nCytosolic", "Fetal\nPolyA\nNuclear")
# KEGG
compareKegg = compareCluster(sigGeneList[1:8], fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareKegg,colorBy="p.adjust",  showCategory = 30, title= "KEGG Pathway Enrichment")

compareKegg.down = compareCluster(sigGeneList[c(1:2,5:10)], fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareKegg.down,colorBy="p.adjust",  showCategory = 30, title= "KEGG Pathway Enrichment")

# Biological Process
compareBP = compareCluster(sigGeneList[1:8], fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 30, title= "Biological Process GO Enrichment")
compareBPDropped = dropGO(compareBP,  level = 3)

compareBP.down = compareCluster(sigGeneList[c(1:2,5:10)], fun="enrichGO",  ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 30, title= "Biological Process GO Enrichment")
compareBPDropped = dropGO(compareBP.down,  level = 3)

# Molecular Function
compareMF = compareCluster(sigGeneList[1:8], fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 30, title= "Molecular Function GO Enrichment")
compareMFDropped = dropGO(compareMF,  level = 3)

compareMF.down = compareCluster(sigGeneList[c(1:2,5:10)], fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 30, title= "Molecular Function GO Enrichment")
compareMFDropped = dropGO(compareMF.down,  level = 3)

# Cellular Component
compareCC = compareCluster(sigGeneList[1:8], fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 30, title= "Cellular Compartment GO Enrichment")
compareCCDropped = dropGO(compareCC,  level = 3)

compareCC.down = compareCluster(sigGeneList[c(1:2,5:10)], fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 30, title= "Cellular Compartment GO Enrichment")
compareCCDropped = dropGO(compareCC.down,  level = 3)

# Disease Ontology
compareDO = compareCluster(sigGeneList[1:8], fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

compareDO.down = compareCluster(sigGeneList[c(1:2,5:10)], fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO.down,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

save(compareKegg, compareBP, compareMF, compareCC, compareDO,
     compareKegg.down, compareBP.down, compareMF.down, compareCC.down, compareDO.down, 
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Fraction.kegg.GO.DO.objects.rda")
