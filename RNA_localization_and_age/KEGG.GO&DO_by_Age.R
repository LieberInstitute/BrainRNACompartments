library("clusterProfiler")
require("S4Vectors")
require("org.Hs.eg.db")
require("Rgraphviz")
library("DOSE")

# Prepare significant genes in a list
AgeList = list(Cpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Cpres.csv"),
               Npres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Npres.csv"),
               Crres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Crres.csv"),
               Nrres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nrres.csv"))

# Define universe as all genes expressed in each of the four groups
GeneUniverse = lapply(AgeList, function(x) x[which(x$baseMean > 0),])
elementLengths(GeneUniverse)
GeneUniverse = lapply(GeneUniverse, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
GeneUniverse = list("Cpres.Up"= GeneUniverse[["Cpres"]], "Cpres.Down"= GeneUniverse[["Cpres"]],
                    "Npres.Up"= GeneUniverse[["Npres"]], "Npres.Down"= GeneUniverse[["Npres"]],
                    "Crres.Up"= GeneUniverse[["Crres"]], "Crres.Down"= GeneUniverse[["Crres"]],
                    "Nrres.Up"= GeneUniverse[["Nrres"]], "Nrres.Down"= GeneUniverse[["Nrres"]])
elementLengths(GeneUniverse)

# Define significantly differently expressed genes
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigAgeList)
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "DownFetal"))
Sign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Cpres.Up = DirList[["Cpres"]][["UpFetal"]], Cpres.Down = DirList[["Cpres"]][["DownFetal"]],
               Npres.Up = DirList[["Npres"]][["UpFetal"]], Npres.Down = DirList[["Npres"]][["DownFetal"]],
               Crres.Up = DirList[["Crres"]][["UpFetal"]], Crres.Down = DirList[["Crres"]][["DownFetal"]],
               Nrres.Up = DirList[["Nrres"]][["UpFetal"]], Nrres.Down = DirList[["Nrres"]][["DownFetal"]]) 
names(DirList) = c("Cytosol\nPolyA\nUp", "Cytosol\nPolyA\nDown", "Nucleus\nPolyA\nUp", "Nucleus\nPolyA\nDown", 
                   "Cytosol\nRibozero\nUp", "Cytosol\nRibozero\nDown", "Nucleus\nRibozero\nUp", "Nucleus\nRibozero\nDown")
elementLengths(SigAgeList)
sigGeneList = lapply(DirList, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
elementLengths(sigGeneList)

# Find enriched Pathways via KEGG
keggList = mapply(function(g, bg) {
  ht=enrichKEGG(as.character(g), organism="human", universe= as.character(bg), 
                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
keggListdf = lapply(keggList, function(x) summary(x))
dotplot(keggList[["Cytosol\nPolyA\nUp"]], showCategory =15, title = "Cytosol PolyA Up in Fetal")
dotplot(keggList[["Cytosol\nPolyA\nDown"]], showCategory =15, title = "Cytosol PolyA Up in Adult")
dotplot(keggList[["Cytosol\nRibozero\nUp"]], showCategory =15, title = "Cytosol RiboZero Up in Fetal")
dotplot(keggList[["Cytosol\nRibozero\nDown"]], showCategory =15, title = "Cytosol RiboZero Up in Adult")
dotplot(keggList[["Nucleus\nPolyA\nUp"]], showCategory =15, title = "Nucleus PolyA Up in Fetal")
dotplot(keggList[["Nucleus\nPolyA\nDown" ]], showCategory =15, title = "Nucleus PolyA Up in Adult")
dotplot(keggList[["Nucleus\nRibozero\nUp"]], showCategory =15, title = "Nucleus RiboZero Up in Fetal")
dotplot(keggList[["Nucleus\nRibozero\nDown"]], showCategory =15, title = "Nucleus RiboZero Up in Adult")

# Enriched Molecular Function GOs
goList_MF = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "MF", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_MF = lapply(goList_MF, function(x) summary(x))
plotGOgraph(goList_MF[["Cytosol\nPolyA\nUp"]])
plotGOgraph(goList_MF[["Nucleus\nPolyA\nUp"]])
dotplot(goList_MF[["Cytosol\nPolyA\nUp"]], showCategory =15, title = "Cytosol PolyA Up in Fetal")
dotplot(goList_MF[["Cytosol\nPolyA\nDown"]], showCategory =15, title = "Cytosol PolyA Up in Adult")
dotplot(goList_MF[["Cytosol\nRibozero\nUp"]], showCategory =15, title = "Cytosol RiboZero Up in Fetal")
dotplot(goList_MF[["Cytosol\nRibozero\nDown"]], showCategory =15, title = "Cytosol RiboZero Up in Adult")
dotplot(goList_MF[["Nucleus\nPolyA\nUp"]], showCategory =15, title = "Nucleus PolyA Up in Fetal")
dotplot(goList_MF[["Nucleus\nPolyA\nDown" ]], showCategory =15, title = "Nucleus PolyA Up in Adult")
dotplot(goList_MF[["Nucleus\nRibozero\nUp"]], showCategory =15, title = "Nucleus RiboZero Up in Fetal")
dotplot(goList_MF[["Nucleus\nRibozero\nDown"]], showCategory =15, title = "Nucleus RiboZero Up in Adult")


# Biological Process GO enrichment
goList_BP = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "BP", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_BP = lapply(goList_BP, function(x) summary(x))
plotGOgraph(goList_BP[["Cytosol\nPolyA\nUp"]])
plotGOgraph(goList_BP[["Nucleus\nPolyA\nUp"]])
dotplot(goList_BP[["Cytosol\nPolyA\nUp"]], showCategory =15, title = "Cytosol PolyA Up in Fetal")
dotplot(goList_BP[["Cytosol\nPolyA\nDown"]], showCategory =15, title = "Cytosol PolyA Up in Adult")
dotplot(goList_BP[["Cytosol\nRibozero\nUp"]], showCategory =15, title = "Cytosol RiboZero Up in Fetal")
dotplot(goList_BP[["Cytosol\nRibozero\nDown"]], showCategory =15, title = "Cytosol RiboZero Up in Adult")
dotplot(goList_BP[["Nucleus\nPolyA\nUp"]], showCategory =15, title = "Nucleus PolyA Up in Fetal")
dotplot(goList_BP[["Nucleus\nPolyA\nDown" ]], showCategory =15, title = "Nucleus PolyA Up in Adult")
dotplot(goList_BP[["Nucleus\nRibozero\nUp"]], showCategory =15, title = "Nucleus RiboZero Up in Fetal")
dotplot(goList_BP[["Nucleus\nRibozero\nDown"]], showCategory =15, title = "Nucleus RiboZero Up in Adult")

# Cellular Compartment GO enrichment
goList_CC = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "CC", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_CC = lapply(goList_CC, function(x) summary(x))
plotGOgraph(goList_CC[["Cytosol\nPolyA\nUp"]])
plotGOgraph(goList_CC[["Nucleus\nPolyA\nUp"]])
dotplot(goList_CC[["Cytosol\nPolyA\nUp"]], showCategory =15, title = "Cytosol PolyA Up in Fetal")
dotplot(goList_CC[["Cytosol\nPolyA\nDown"]], showCategory =15, title = "Cytosol PolyA Up in Adult")
dotplot(goList_CC[["Cytosol\nRibozero\nUp"]], showCategory =15, title = "Cytosol RiboZero Up in Fetal")
dotplot(goList_CC[["Cytosol\nRibozero\nDown"]], showCategory =15, title = "Cytosol RiboZero Up in Adult")
dotplot(goList_CC[["Nucleus\nPolyA\nUp"]], showCategory =15, title = "Nucleus PolyA Up in Fetal")
dotplot(goList_CC[["Nucleus\nPolyA\nDown" ]], showCategory =15, title = "Nucleus PolyA Up in Adult")
dotplot(goList_CC[["Nucleus\nRibozero\nUp"]], showCategory =15, title = "Nucleus RiboZero Up in Fetal")
dotplot(goList_CC[["Nucleus\nRibozero\nDown"]], showCategory =15, title = "Nucleus RiboZero Up in Adult")

# Compare the enriched terms between four groups

# KEGG
compareKegg.Age = compareCluster(sigGeneList, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareKegg.Age)
plot(compareKegg.Age,colorBy="p.adjust",  showCategory = 20, title= "KEGG Pathway Enrichment")

# Biological Process
compareBP.Age = compareCluster(sigGeneList, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareBP.Age)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 20, title= "Biological Process GO Enrichment")
compareBPDropped = dropGO(compareBP.Age,  level = 1)
compareBPDropped = dropGO(compareBPDropped,  level = 2)
compareBPDropped = dropGO(compareBPDropped,  level = 3)

# Molecular Function
compareMF.Age = compareCluster(sigGeneList, fun="enrichGO",  ont = "MF", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareMF.Age)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 20, title= "Molecular Function GO Enrichment")
compareMFDropped = dropGO(compareMF.Age,  level = 1)
compareMFDropped = dropGO(compareMFDropped,  level = 2)
compareMFDropped = dropGO(compareMFDropped,  level = 3)

# Cellular Component
compareCC.Age = compareCluster(sigGeneList, fun="enrichGO",  ont = "CC", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareCC.Age)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 20, title= "Cellular Compartment GO Enrichment")
compareCCDropped = dropGO(compareCC.Age,  level = 1)
compareCCDropped = dropGO(compareCCDropped,  level = 2)

# Disease Ontology
compareDO.Age = compareCluster(sigGeneList, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareDO.Age)
plot(compareDODropped,colorBy="p.adjust",  showCategory = 20, title= "Disease Ontology Enrichment")

save(compareKegg.Age, compareBP.Age, compareMF.Age, compareCC.Age, compareDO.Age,
     file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/kegg.GO.DO.Age.objects.rda")