library("clusterProfiler")
require("org.Hs.eg.db")

library(GenomicRanges)
library(ggplot2)
library(plyr)

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames
lapply(IRres, head)
# Filter introns
allIntrons = IRres[[1]]
allIntrons = allIntrons[grep("clean", allIntrons$GeneIntronDetails, fixed=T),]
string = unlist(strsplit(as.character(allIntrons$GeneIntronDetails), "/", fixed = TRUE), recursive = FALSE)
allIntrons = data.frame(allIntrons[,c(1:4,6)], genes = string[grep("ENSG", string)],
                        intronID = paste0("chr",allIntrons$Chr,":",allIntrons$Start,"-",allIntrons$End))
# read in Differential IR results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
nonconst = list()
for (i in 1:length(comps)){
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
}
names(nonconst) = comps
elementNROWS(nonconst)
string = lapply(nonconst, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE))
IR.diff = lapply(nonconst, function(y) y$A.IRratio - y$B.IRratio)
dIR = Map(cbind, nonconst, ensID = lapply(string, function(y) y[grep("ENSG", y)]), 
          comments = lapply(string, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])), 
          IR.diff = IR.diff, Sign = lapply(IR.diff, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")), 
          intronID = lapply(nonconst, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End)))
dIRclean = lapply(dIR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" &
                                             y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                                             y$comments=="clean"),])
sigdIR = lapply(dIRclean, function(x) x[which(x$p.diff<=0.05),])
sigdIR = lapply(sigdIR, function(x) split(x, x$Sign))
sigdIR = unlist(sigdIR, recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                  "Cytosol-Increased", "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)]))
introns = c("All Introns" = list(allIntrons), "Introns (Fraction)" = list(do.call(rbind, lapply(dIRclean[1:2], function(x) x[,c(1:4,6,8,30,32,33,34)]))), 
            "Introns (Age)" = list(do.call(rbind, lapply(dIRclean[3:4], function(x) x[,c(1:4,6,8,30,32,33,34)]))), sigdIR)
introns[["Introns (Fraction)"]] = introns[["Introns (Fraction)"]][unique(introns[["Introns (Fraction)"]][,"intronID"]),]
introns[["Introns (Age)"]] = introns[["Introns (Age)"]][unique(introns[["Introns (Age)"]][,"intronID"]),]

## Gene Ontology
names(sig.1) = names(sig) = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only",
                              "Exported:\nPrenatal Only", "Exported:\nAdult Only","Retained: Adult/\nExported: Prenatal",
                              "Retained: Prenatal/\nExported: Adult", "Interaction")
entrezID = lapply(sig.1, function(x) na.omit(x$EntrezID))
# Define universe as all genes expressed in each of the four groups
GeneUniverse = as.character(unique(geneMap[match(rownames(Ires.down),geneMap$gencodeID),"EntrezID"]))
GeneUniverse = na.omit(GeneUniverse)
# Find enriched Pathways via KEGG
elementNROWS(entrezID)
keggList = lapply(entrezID, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
keggListdf = lapply(keggList, function(x) as.data.frame(x))
# Enriched Molecular Function GOs
goList_MF = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
goListdf_MF = lapply(goList_MF, function(x) as.data.frame(x))
# Biological Process GO enrichment
goList_BP = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
goListdf_BP = lapply(goList_BP, function(x) as.data.frame(x))
# Cellular Compartment GO enrichment
goList_CC = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
goListdf_CC = lapply(goList_CC, function(x) as.data.frame(x))
# Disease Ontology
goList_DO = lapply(entrezID, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))
goListdf_DO = lapply(goList_DO, function(x) as.data.frame(x))

# write to csv
names = c("Retained_Both", "Exported_Both", "Retained_Prenatal.Only", "Retained_Adult.Only",
          "Exported_Prenatal.Only", "Exported_Adult.Only",
          "Retained_Adult.Exported_Prenatal", "Retained_Prenatal.Exported_Adult", "Interaction")
lists = list(DO = goListdf_DO, CC = goListdf_CC, BP = goListdf_BP, MF = goListdf_MF, KEGG = keggListdf)
for (i in 1:length(lists)){
  for (j in 1:length(names)){
    if (nrow(lists[[i]][[j]])>0){
      lists[[i]][[j]] = data.frame(lists[[i]][[j]], Comparison = names[j])
    }}}
lists = lapply(lists, function(x) do.call(rbind, x))
write.csv(lists[["KEGG"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.KEGG.downsampled.csv")
write.csv(lists[["BP"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.BP.downsampled.csv")
write.csv(lists[["MF"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.MF.downsampled.csv")
write.csv(lists[["CC"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.CC.downsampled.csv")
write.csv(lists[["DO"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.DO.downsampled.csv")

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrezID, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
# Biological Process
compareBP = compareCluster(entrezID, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBPDropped = dropGO(compareBP,  level = 3)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 30, title= "Biological Process GO Enrichment")
# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMFDropped = dropGO(compareMF,  level = 3)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 30, title= "Molecular Function GO Enrichment")
# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCCDropped = dropGO(compareCC,  level = 3)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 30, title= "Cellular Compartment GO Enrichment")
# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
entrez.noLFC = lapply(sig, function(x) na.omit(x$EntrezID))
compareDO = compareCluster(entrez.noLFC, fun="enrichDO",
                           ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 40, title= "Disease Ontology Enrichment")

save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.downsampled.rda")
