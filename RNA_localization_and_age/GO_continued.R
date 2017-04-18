library("clusterProfiler")
require("S4Vectors")
require("org.Hs.eg.db")
require("Rgraphviz")
library("DOSE")

load('/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/kegg.GO.DO.objects.rda')
load("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Last.Exon>DEResults.rda")

# Explore GO for whole gene DE by fraction
DO = summary(compareDO)
dim(DO)
kegg = summary(compareKegg)
dim(kegg)
BP = summary(compareBP)
dim(BP)
BP[which(BP$Cluster=="Adult\nPolyA\nUp"),]
MF = summary(compareMF)
dim(MF)
colnames(MF)
MF[,1:3]

# Define universe as all genes expressed in each of the four groups
LEResults = lapply(LEResults, function(x) x[order(rownames(x)),])
LastExonMap = LastExonMap[order(rownames(LastExonMap)),]
LEResults = lapply(LEResults, function(x) data.frame(x, GeneID = LastExonMap$Geneid, EntrezID = LastExonMap$EntrezID))
GeneUniverse = lapply(LEResults, function(x) x[which(x$baseMean > 0),])
elementLengths(GeneUniverse)
GeneUniverse = lapply(GeneUniverse, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
GeneUniverse = list("AP.Up"= GeneUniverse[["AP"]], "AP.Down"= GeneUniverse[["AP"]],
                    "FP.Up"= GeneUniverse[["FP"]], "FP.Down"= GeneUniverse[["FP"]],
                    "AR.Up"= GeneUniverse[["AR"]], "AR.Down"= GeneUniverse[["AR"]],
                    "FR.Up"= GeneUniverse[["FR"]], "FR.Down"= GeneUniverse[["FR"]])
elementLengths(GeneUniverse)

# Define significantly differently expressed genes
SigLEList = lapply(LEResults, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigLEList)
Sign = lapply(SigLEList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigLEBySign = Map(cbind, SigLEList, Sign = Sign)
DirLEList = lapply(sigLEBySign, function(x) split(x, x$Sign))
SigLEList = list("AP.Up"= DirLEList[["AP"]][["UpNuc"]], "AP.Down"= DirLEList[["AP"]][["DownNuc"]],
                 "FP.Up"= DirLEList[["FP"]][["UpNuc"]], "FP.Down"= DirLEList[["FP"]][["DownNuc"]],
                 "AR.Up"= DirLEList[["AR"]][["UpNuc"]], "AR.Down"= DirLEList[["AR"]][["DownNuc"]],
                 "FR.Up"= DirLEList[["FR"]][["UpNuc"]], "FR.Down"= DirLEList[["FR"]][["DownNuc"]])
 elementLengths(SigLEList)
sigGeneLE = lapply(SigLEList, function(x) unique(x$EntrezID[!is.na(x$EntrezID)]))
elementLengths(sigGeneLE)
save(sigGeneLE, GeneUniverse, file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/lastExonGO.rda")
load("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/lastExonGO.rda")

# Find enriched Pathways via KEGG
keggList = mapply(function(g, bg) {
  ht=enrichKEGG(as.character(g), organism="human", universe= as.character(bg), 
                minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneLE, GeneUniverse, SIMPLIFY=FALSE)
keggList = lapply(keggList, function(x) summary(x))

# Enriched Molecular Function GOs
goList_MF = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "MF", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneLE, GeneUniverse, SIMPLIFY=FALSE)
goList_MF = lapply(goList_MF, function(x) summary(x))

# Biological Process GO enrichment
goList_BP = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "BP", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneLE, GeneUniverse, SIMPLIFY=FALSE)
goList_BP = lapply(goList_BP, function(x) summary(x))

# Cellular Compartment GO enrichment
goList_CC = mapply(function(g, bg) {
  ht=enrichGO(as.character(g), ont = "CC", organism="human", universe= as.character(bg), 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE)}, 
  sigGeneLE, GeneUniverse, SIMPLIFY=FALSE)
goList_CC = lapply(goList_CC, function(x) summary(x))

# Compare the enriched terms between four groups
elementLengths(sigGeneLE)
names(sigGeneLE) = c("Adult\nPolyA\nUp", "Adult\nPolyA\nDown", "Fetal\nPolyA\nUp", "Fetal\nPolyA\nDown", 
                       "Adult\nRibozero\nUp", "Adult\nRibozero\nDown", "Fetal\nRibozero\nUp", "Fetal\nRibozero\nDown")
# KEGG
compareKeggLE = compareCluster(sigGeneLE, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
KeggLE = summary(compareKeggLE)
plot(compareKeggLE,colorBy="p.adjust",  showCategory = 20, title= "KEGG Pathway Enrichment\nBy Last Exon Expression")

# Biological Process
compareBPLE = compareCluster(sigGeneLE, fun="enrichGO",  ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BPLE = summary(compareBPLE)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 20, title= "Biological Process GO Enrichment\nBy Last Exon Expression")
compareBPDropped = dropGO(compareBPLE,  level = 1)
compareBPDropped = dropGO(compareBPDropped,  level = 2)

# Molecular Function
compareMFLE = compareCluster(sigGeneLE, fun="enrichGO",  ont = "MF", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MFLE = summary(compareMFLE)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 20, title= "Molecular Function GO Enrichment\nBy Last Exon Expression")
compareMFDropped = dropGO(compareMFLE,  level = 1)
compareMFDropped = dropGO(compareMFDropped,  level = 2)

# Cellular Component
compareCCLE = compareCluster(sigGeneLE, fun="enrichGO",  ont = "CC", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MFLE = summary(compareCCLE)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 20, title= "Cellular Compartment GO Enrichment\nBy Last Exon Expression")
compareCCDropped = dropGO(compareCCLE,  level = 1)
compareCCDropped = dropGO(compareCCDropped,  level = 2)
compareCCDropped = dropGO(compareCCDropped,  level = 3)

# Disease Ontology
compareDOLE = compareCluster(sigGeneLE, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DOLE = summary(compareDOLE)
plot(compareDOLE,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment\nBy Last Exon Expression")

save(GeneUniverse, sigGeneLE, keggList, goList_MF, goList_BP, goList_CC, 
     compareKeggLE, compareBPLE, compareMFLE, compareCCLE, compareDOLE, file = "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/lastExonGO.rda")



