library(clusterProfiler)
require(org.Hs.eg.db)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

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
names(sigdIR) = c("Cytosolic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytosolic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                  "Adult Retention\nIn Cytosol","Prenatal Retention\nIn Cytosol","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus",
                  "Cytosolic Retention", "Nuclear Retention", "Adult Retention", "Prenatal Retention")
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)]))
sigdIR = Map(cbind, sigdIR, EntrezID = lapply(sigdIR, function(x) geneMap[match(x$ensID, geneMap$ensemblID), "EntrezID"]))
entrezID.allGroups = lapply(sigdIR, function(x) na.omit(x$EntrezID))
entrezID.collapsed = lapply(list("Cytosolic Retention" = do.call(rbind,list(sigdIR[["Cytosolic Retention\nIn Adult"]],
                                                                            sigdIR[["Cytosolic Retention\nIn Prenatal"]])),
                                 "Nuclear Retention" = do.call(rbind,list(sigdIR[["Nuclear Retention\nIn Adult"]],
                                                                          sigdIR[["Nuclear Retention\nIn Prenatal"]], 
                                                                          sigdIR[["Nuclear Retention"]])),
                                 "Adult Retention" = do.call(rbind,list(sigdIR[["Adult Retention\nIn Cytosol"]],
                                                                        sigdIR[["Adult Retention\nIn Nucleus"]], 
                                                                        sigdIR[["Adult Retention"]])),
                                 "Prenatal Retention" = do.call(rbind,list(sigdIR[["Prenatal Retention\nIn Cytosol"]],
                                                                           sigdIR[["Prenatal Retention\nIn Nucleus"]], 
                                                                           sigdIR[["Prenatal Retention"]]))), function(x) na.omit(x$EntrezID))

## Compare the enriched terms between groups
# KEGG
compareKegg = compareCluster(entrezID.allGroups, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
compareKegg.collapsed = compareCluster(entrezID.collapsed, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareKegg.collapsed,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")

# Biological Process
compareBP = compareCluster(entrezID.allGroups, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareBP,colorBy="p.adjust",  showCategory = 47, title= "Biological Process GO Enrichment")
compareBP.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareBP.collapsed,colorBy="p.adjust",  showCategory = 47, title= "Biological Process GO Enrichment")

# Molecular Function
compareMF = compareCluster(entrezID.allGroups, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareMF,colorBy="p.adjust",  showCategory = 30, title= "Molecular Function GO Enrichment")
compareMF.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareMF.collapsed,colorBy="p.adjust",  showCategory = 30, title= "Molecular Function GO Enrichment")

# Cellular Component
compareCC = compareCluster(entrezID.allGroups, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareCC,colorBy="p.adjust",  showCategory = 30, title= "Cellular Compartment GO Enrichment")
compareCC.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareCC.collapsed,colorBy="p.adjust",  showCategory = 30, title= "Cellular Compartment GO Enrichment")

# Disease Ontology
compareDO = compareCluster(entrezID.allGroups, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.kegg.GO.DO.objects.rda")