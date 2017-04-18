library("clusterProfiler")
require("S4Vectors")
require("org.Hs.eg.db")
require("Rgraphviz")
library("DOSE")

# Prepare significant genes in a list
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))

# Define universe as all genes expressed in each of the four groups
GeneUniverse = lapply(FracList, function(x) x[which(x$baseMean > 0),])
elementLengths(GeneUniverse)
GeneUniverse = lapply(GeneUniverse, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
GeneUniverse = list("Apres.Up"= GeneUniverse[["Apres"]], "Apres.Down"= GeneUniverse[["Apres"]],
                    "Fpres.Up"= GeneUniverse[["Fpres"]], "Fpres.Down"= GeneUniverse[["Fpres"]],
                    "Arres.Up"= GeneUniverse[["Arres"]], "Arres.Down"= GeneUniverse[["Arres"]],
                    "Frres.Up"= GeneUniverse[["Frres"]], "Frres.Down"= GeneUniverse[["Frres"]])
elementLengths(GeneUniverse)

# Define significantly differently expressed genes
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracBySign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(sigFracBySign, function(x) split(x, x$Sign))
SigList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]]) 
elementLengths(SigList)
sigGeneList = lapply(SigList, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
elementLengths(sigGeneList)

# Find enriched Pathways via KEGG
keggList = mapply(function(g, bg) {
  ht=enrichKEGG(as.character(g), organism="human", universe= as.character(bg), 
  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
keggListdf = lapply(keggList, function(x) summary(x))
dotplot(keggList[["Apres.Up"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(keggList[["Apres.Down"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(keggList[["Arres.Up"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(keggList[["Arres.Down"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(keggList[["Fpres.Up"]], showCategory =15, title = "Fetal PolyA Up in Nucleus")
dotplot(keggList[["Fpres.Down"]], showCategory =15, title = "Fetal PolyA Up in Nucleus")
dotplot(keggList[["Frres.Up"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")
dotplot(keggList[["Frres.Down"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")

# Enriched Molecular Function GOs
goList_MF = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "MF", organism="human", universe= as.character(bg), 
  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_MF = lapply(goList_MF, function(x) summary(x))
plotGOgraph(goList_MF[["Apres.Up"]])
plotGOgraph(goList_MF[["Apres.Down"]])
dotplot(goList_MF[["Apres.Up"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(goList_MF[["Apres.Down"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(goList_MF[["Arres.Up"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(goList_MF[["Arres.Down"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(goList_MF[["Frres.Up"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")
dotplot(goList_MF[["Frres.Down"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")

# Biological Process GO enrichment
goList_BP = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "BP", organism="human", universe= as.character(bg), 
  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_BP = lapply(goList_BP, function(x) summary(x))
plotGOgraph(goList_BP[["Apres.Up"]])
plotGOgraph(goList_BP[["Apres.Down"]])
dotplot(goList_BP[["Apres.Up"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(goList_BP[["Apres.Down"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(goList_BP[["Arres.Up"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(goList_BP[["Arres.Down"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(goList_BP[["Frres.Up"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")
dotplot(goList_BP[["Frres.Down"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")

# Cellular Compartment GO enrichment
goList_CC = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "CC", organism="human", universe= as.character(bg), 
  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_CC = lapply(goList_CC, function(x) summary(x))
plotGOgraph(goList_CC[["Apres.Up"]])
plotGOgraph(goList_CC[["Apres.Down"]])
dotplot(goList_CC[["Apres.Up"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(goList_CC[["Apres.Down"]], showCategory =15, title = "Adult PolyA Up in Nucleus")
dotplot(goList_CC[["Arres.Up"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(goList_CC[["Arres.Down"]], showCategory =15, title = "Adult RiboZero Up in Nucleus")
dotplot(goList_CC[["Frres.Down"]], showCategory =15, title = "Fetal RiboZero Up in Nucleus")

# Compare the enriched terms between four groups
elementLengths(sigGeneList)
names(sigGeneList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                       "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", "Fetal\nRibozero\nDown")
# KEGG
compareKegg = compareCluster(sigGeneList, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareKegg)
plot(compareKegg,colorBy="p.adjust",  showCategory = 20, title= "KEGG Pathway Enrichment")

# Biological Process
compareBP = compareCluster(sigGeneList, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareBP)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 20, title= "Biological Process GO Enrichment")
compareBPDropped = dropGO(compareBP,  level = 1)
compareBPDropped = dropGO(compareBPDropped,  level = 2)
compareBPDropped = dropGO(compareBPDropped,  level = 3)

# Molecular Function
compareMF = compareCluster(sigGeneList, fun="enrichGO",  ont = "MF", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareMF)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 20, title= "Molecular Function GO Enrichment")
compareMFDropped = dropGO(compareMF,  level = 1)
compareMFDropped = dropGO(compareMFDropped,  level = 2)
compareMFDropped = dropGO(compareMFDropped,  level = 3)

# Cellular Component
compareCC = compareCluster(sigGeneList, fun="enrichGO",  ont = "CC", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareCC)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 20, title= "Cellular Compartment GO Enrichment")
compareCCDropped = dropGO(compareCC,  level = 1)
compareCCDropped = dropGO(compareCCDropped,  level = 2)
compareCCDropped = dropGO(compareCCDropped,  level = 3)

# Disease Ontology
compareDO = compareCluster(sigGeneList, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareDO)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/kegg.GO.DO.objects.rda")
