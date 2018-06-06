library(clusterProfiler)
require(org.Hs.eg.db)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_results_object")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/introns_object.rda")


## Isolate Entrez IDs

IRcomp = Map(cbind, introns, EntrezID = lapply(introns, function(x) geneMap[match(x$ensID, geneMap$ensemblID),"EntrezID"]))
entrezID.allGroups = lapply(IRcomp, function(x) unique(na.omit(x$EntrezID)))
entrezID.collapsed = lapply(list("Cytoplasmic Retention" = do.call(rbind, IRcomp[grep(":Cytoplasm-",names(IRcomp))]),
                                 "Nuclear Retention" = do.call(rbind, IRcomp[grep(":Nucleus-",names(IRcomp))]),
                                 "Adult Retention" = do.call(rbind, IRcomp[grep(":Adult-",names(IRcomp))]),
                                 "Prenatal Retention" = do.call(rbind, IRcomp[grep(":Prenatal-",names(IRcomp))])), function(x) unique(na.omit(x$EntrezID)))
entrezID.allGroups = lapply(entrezID.allGroups[names(entrezID.allGroups)!="All Introns"], as.character)
entrezID.collapsed = lapply(entrezID.collapsed[names(entrezID.collapsed)!="All Introns"], as.character)


## Compare the enriched terms between groups
# KEGG
compareKegg = compareCluster(entrezID.allGroups, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareKegg.collapsed = compareCluster(entrezID.collapsed, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrezID.allGroups, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrezID.allGroups, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrezID.allGroups, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrezID.allGroups, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO.collapsed = compareCluster(entrezID.collapsed, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/A-CTest_IRintrons_KEGG_GO_DO_plots.pdf", width=16, height=9)
plot(compareKegg.collapsed,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP.collapsed,colorBy="p.adjust",  showCategory = 45, title= "Biological Process Pathway Enrichment")
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process Pathway Enrichment")
plot(compareMF.collapsed,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function Pathway Enrichment")
plot(compareMF,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function Pathway Enrichment")
plot(compareCC.collapsed,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment Pathway Enrichment")
plot(compareCC,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment Pathway Enrichment")
plot(compareDO.collapsed,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
dev.off()

save(compareKegg, compareBP, compareMF, compareCC, compareDO,compareKegg.collapsed, compareBP.collapsed, compareMF.collapsed, compareCC.collapsed,compareDO.collapsed,
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/A-CTest.intron.kegg.GO.DO.objects.rda")

plotExample = compareBP # clusterProfiler output
plotExample@compareClusterResult$Cluster = as.character(plotExample@compareClusterResult$Cluster)
plotExample@compareClusterResult$Cluster = gsub(":",":\n", plotExample@compareClusterResult$Cluster, fixed=T)
plotExample@compareClusterResult$Cluster = factor(plotExample@compareClusterResult$Cluster, levels=
                                                    c("Prenatal:\nNucleus-Increased","Cytoplasm:\nPrenatal-Increased","Nucleus:\nPrenatal-Increased"))
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/biologicalProcesses.pdf", width = 12,height = 6)
plot(plotExample,colorBy="p.adjust",  showCategory = 10, title= "Biological Process Pathway Enrichment")
dev.off()


## Read in GLM results

IDs = lapply(dIR, function(x) unlist(strsplit(rownames(x), "/", fixed=T), recursive = F))
dIR = Map(cbind, dIR, ensID = lapply(IDs, function(x) x[grep("ENSG", x)]), coord = lapply(IDs, function(x) x[grep(":", x, fixed=T)]),
          geneSymbol = lapply(IDs, function(x) x[-c(grep("ENSG", x),grep(":", x, fixed=T),grep("clean", x))]),
          A.Sign = lapply(dIR, function(y) ifelse((y$Ad.IRratio.ZoneCytoplasm - y$Ad.IRratio.ZoneNucleus) < 0,"MoreIRInNuc", "MoreIRInCyt")),
          P.Sign = lapply(dIR, function(y) ifelse((y$Pren.IRratio.ZoneCytoplasm - y$Pren.IRratio.ZoneNucleus) < 0,"MoreIRInNuc", "MoreIRInCyt")),
          C.Sign = lapply(dIR, function(y) ifelse((y$Cyt.IRratio.FetalAdult - y$Cyt.IRratio.FetalPrenatal) < 0,"MoreIRInPrenatal", "MoreIRInAdult")),
          N.Sign = lapply(dIR, function(y) ifelse((y$Nuc.IRratio.FetalAdult - y$Nuc.IRratio.FetalPrenatal) < 0,"MoreIRInPrenatal", "MoreIRInAdult")))
head(dIR[[1]])
dIR = Map(cbind, dIR, Sign = list(ifelse(dIR$byFractionInAdult$A.Sign=="MoreIRInNuc", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byFractionInPrenatal$P.Sign=="MoreIRInNuc", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byAgeInCytoplasm$C.Sign=="MoreIRInPrenatal", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byAgeInNucleus$N.Sign=="MoreIRInPrenatal", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult")),
          EntrezID = lapply(dIR, function(x) geneMap[match(x$ensID, geneMap$ensemblID),"EntrezID"]))
sigIR = lapply(dIR, function(x) x[which(x$padj<=0.05),])
elementNROWS(sigIR)
sigIR = unlist(lapply(sigIR, function(x) split(x, x$Sign)), recursive=F)
names(sigIR) = c("Cytoplasmic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytoplasmic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                  "Adult Retention\nIn Cytoplasm","Prenatal Retention\nIn Cytoplasm","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus")
entrezID.allGroups = lapply(sigIR, function(x) na.omit(x$EntrezID))
entrezID.collapsed = lapply(list("Cytoplasmic Retention" = do.call(rbind, sigIR[grep("Cytoplasmic Retention",names(sigIR))]),
                                 "Nuclear Retention" = do.call(rbind, sigIR[grep("Nuclear Retention",names(sigIR))]),
                                 "Adult Retention" = do.call(rbind, sigIR[grep("Adult Retention",names(sigIR))]),
                                 "Prenatal Retention" = do.call(rbind, sigIR[grep("Prenatal Retention",names(sigIR))])), function(x) na.omit(x$EntrezID))
entrezID.allGroups = lapply(entrezID.allGroups, as.character)
entrezID.collapsed = lapply(entrezID.collapsed, as.character)


## Compare the enriched terms between groups
# KEGG
compareKegg = compareCluster(entrezID.allGroups, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareKegg.collapsed = compareCluster(entrezID.collapsed, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrezID.allGroups, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBP.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrezID.allGroups, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMF.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrezID.allGroups, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCC.collapsed = compareCluster(entrezID.collapsed, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrezID.allGroups, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareDO.collapsed = compareCluster(entrezID.collapsed, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_IRintrons_KEGG_GO_DO_plots.pdf", width=16, height=9)
plot(compareKegg.collapsed,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP.collapsed,colorBy="p.adjust",  showCategory = 45, title= "Biological Process Pathway Enrichment")
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process Pathway Enrichment")
plot(compareMF.collapsed,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function Pathway Enrichment")
plot(compareMF,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function Pathway Enrichment")
plot(compareCC.collapsed,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment Pathway Enrichment")
plot(compareCC,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment Pathway Enrichment")
plot(compareCC.collapsed,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
dev.off()

save(compareKegg, compareBP, compareMF, compareCC, compareDO,compareKegg.collapsed, compareBP.collapsed, compareMF.collapsed, compareCC.collapsed,
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM.intron.kegg.GO.DO.objects.rda")