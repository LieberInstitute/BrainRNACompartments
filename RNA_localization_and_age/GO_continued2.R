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
develop = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/developmentally-regulated.genes.csv")

# At LFC >= 1
devel = develop[which(develop$padj<=0.05 & abs(develop$log2FoldChange)>=1),]
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigFracList)

Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]],
               Devel = devel) 
names(DirList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", 
                   "Fetal\nRibozero\nDown", "Developmentally\nRegulated")
# At lfc >=2
devel = develop[which(develop$padj<=0.05 & abs(develop$log2FoldChange)>=2),]
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=2),])
elementLengths(SigFracList)

Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
Sign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]],
               Devel = devel) 
names(DirList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                   "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", 
                   "Fetal\nRibozero\nDown", "Developmentally\nRegulated")

# Define universe as all genes expressed in each of the four groups
GeneUniverse = lapply(FracList, function(x) x[which(x$baseMean > 0),])
elementLengths(GeneUniverse)
GeneUniverse = lapply(GeneUniverse, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
GeneUniverse = list("Adult\nPolyA\nUp"= GeneUniverse[["Apres"]], "Adult\nPolyA\nDown"= GeneUniverse[["Apres"]],
                    "Fetal\nPolyA\nUp"= GeneUniverse[["Fpres"]], "Fetal\nPolyA\nDown"= GeneUniverse[["Fpres"]],
                    "Adult\nRibozero\nUp"= GeneUniverse[["Arres"]], "Adult\nRibozero\nDown"= GeneUniverse[["Arres"]],
                    "Fetal\nRibozero\nUp"= GeneUniverse[["Frres"]], "Fetal\nRibozero\nDown"= GeneUniverse[["Frres"]])
elementLengths(GeneUniverse)

# Define significantly differently expressed genes
elementLengths(DirList)
sigGeneList = lapply(DirList[1:8], function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
elementLengths(sigGeneList)

# Find enriched Pathways via KEGG
keggList = mapply(function(g, bg) {
  ht=enrichKEGG(as.character(g), organism="human", universe= as.character(bg), 
                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
keggListdf = lapply(keggList, function(x) summary(x))

# Enriched Molecular Function GOs
goList_MF = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "MF", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_MF = lapply(goList_MF, function(x) summary(x))

# Biological Process GO enrichment
goList_BP = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "BP", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_BP = lapply(goList_BP, function(x) summary(x))

# Cellular Compartment GO enrichment
goList_CC = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "CC", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_CC = lapply(goList_CC, function(x) summary(x))

# Disease Ontology
goList_DO = mapply(function(g, bg) {
  ht=enrichDO(as.character(g), ont = "DO", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneList, GeneUniverse, SIMPLIFY=FALSE)
goListdf_DO = lapply(goList_DO, function(x) summary(x))

# write to csv
names = c("Adult.PolyA.Up","Adult.PolyA.Down","Fetal.PolyA.Up","Fetal.PolyA.Down","Adult.Ribozero.Up",
          "Adult.Ribozero.Down","Fetal.Ribozero.Up","Fetal.Ribozero.Down")
n = c(1,2,3,5,6,7,8)
for (i in n){
  write.csv(goListdf_CC[[i]], file=paste0("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/GO_terms/",names[i],".GO.CC.csv"))
  write.csv(keggListdf[[i]], file=paste0("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/GO_terms/",names[i],".GO.KEGG.csv"))
  write.csv(goListdf_BP[[i]], file=paste0("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/GO_terms/",names[i],".GO.BP.csv"))
  write.csv(goListdf_MF[[i]], file=paste0("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/GO_terms/",names[i],".GO.MF.csv"))
  write.csv(goListdf_DO[[i]], file=paste0("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/GO_terms/",names[i],".GO.DO.csv"))
}

# Compare the enriched terms between four groups
elementLengths(sigGeneList)
names(sigGeneList) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                       "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", "Fetal\nRibozero\nDown")
sigGeneList.UP = list("Adult:PolyA"=sigGeneList[["Adult\nPolyA\nUp"]], 
                     "Fetal:PolyA"=sigGeneList[["Fetal\nPolyA\nUp"]],
                     "Adult:Ribozero"=sigGeneList[["Adult\nRibozero\nUp"]], 
                     "Fetal:Ribozero"=sigGeneList[["Fetal\nRibozero\nUp"]])
sigGeneList.D = list("Adult:PolyA"=sigGeneList[["Adult\nPolyA\nDown"]], 
                      "Fetal:PolyA"=sigGeneList[["Fetal\nPolyA\nDown"]],
                      "Adult:Ribozero"=sigGeneList[["Adult\nRibozero\nDown"]], 
                      "Fetal:Ribozero"=sigGeneList[["Fetal\nRibozero\nDown"]])

# KEGG
compareKegg = compareCluster(sigGeneList, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareKegg)
plot(compareKegg,colorBy="p.adjust",  showCategory = 20, title= "KEGG Pathway Enrichment")

# Biological Process
compareBP = compareCluster(sigGeneList, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareBP.UP)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 20, title= "Biological Process GO Enrichment")
compareBPDropped = dropGO(compareBP,  level = 1)
compareBPDropped = dropGO(compareBPDropped,  level = 2)
compareBPDropped = dropGO(compareBPDropped,  level = 3)

compareBP.UP = compareCluster(sigGeneList.UP, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP.D = compareCluster(sigGeneList.D, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareBP.UP,colorBy="p.adjust",  showCategory = 20, title= "Biological Processes Enriched in Nuclear Genes")
plot(compareBP.D,colorBy="p.adjust",  showCategory = 20, title= "Biological Processes Enriched in Cytosolic Genes")

compareBP.lfc2 = compareCluster(sigGeneList, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareBP.lfc2)
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

compareMF.lfc2 = compareCluster(sigGeneList, fun="enrichGO",  ont = "MF", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareMF.lfc2)
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
compareDO.lfc2 = compareCluster(sigGeneList, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareDO.lfc2)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

compareDO.UP = compareCluster(sigGeneList.UP, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO.D = compareCluster(sigGeneList.D, fun="enrichDO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO.UP,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enriched in Nuclear Genes")


save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/kegg.GO.DO.objects.rda")
