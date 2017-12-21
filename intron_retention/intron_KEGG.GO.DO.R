library(clusterProfiler)
require(org.Hs.eg.db)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_results_object")


# read in Differential IR results files

comps = c("Adult_PolyA_Zone_cleanIntrons_adultShared","Fetal_PolyA_Zone_cleanIntrons_prenatalShared",
          "Cytosol_PolyA_Age_cleanIntrons_cytosolShared","Nuclear_PolyA_Age_cleanIntrons_nucleusShared")
IRcomp = list()
for (i in 1:length(comps)){
  IRcomp[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/", comps[i], ".tab"), header = TRUE, comment.char="#")
}
names(IRcomp) = c("Adult_byFraction","Fetal_byFraction","Cytosol_byAge","Nuclear_byAge")
IRcomp = Map(cbind, IRcomp,
             intronID = lapply(IRcomp, function(x) paste0(x$Intron.GeneName.GeneID,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)),
             ensID = lapply(lapply(IRcomp, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE)), function(y) y[grep("ENSG", y)]),
             IR.diff = lapply(IRcomp, function(y) y$A.IRratio - y$B.IRratio),
             Sign = lapply(IRcomp, function(y) ifelse((y$A.IRratio - y$B.IRratio) < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")))
full = list(adult = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-adultShared.txt", header = TRUE),
            prenatal = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt", header = TRUE),
            cytosol = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt", header = TRUE),
            nucleus = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt", header = TRUE))
total = as.list(elementNROWS(full))
names(total) = names(IRcomp)
IRcomp = Map(cbind, IRcomp, padj = mapply(function(p,t) p.adjust(p, method = "fdr", n = t), lapply(IRcomp, function(x) x$p.diff), total),
             EntrezID = lapply(IRcomp, function(x) geneMap[match(x$ensID, geneMap$ensemblID),"EntrezID"]))


## Isolate Entrez IDs

sigdIR = lapply(IRcomp, function(x) x[which(x$padj<=0.05),])
sigdIR = unlist(lapply(sigdIR, function(x) split(x, x$Sign)), recursive=F)
names(sigdIR) = c("Cytosolic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytosolic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                  "Adult Retention\nIn Cytosol","Prenatal Retention\nIn Cytosol","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus")
entrezID.allGroups = lapply(sigdIR, function(x) na.omit(x$EntrezID))
entrezID.collapsed = lapply(list("Cytosolic Retention" = do.call(rbind, sigdIR[grep("Cytosolic Retention",names(sigdIR))]),
                                 "Nuclear Retention" = do.call(rbind, sigdIR[grep("Nuclear Retention",names(sigdIR))]),
                                 "Adult Retention" = do.call(rbind, sigdIR[grep("Adult Retention",names(sigdIR))]),
                                 "Prenatal Retention" = do.call(rbind, sigdIR[grep("Prenatal Retention",names(sigdIR))])), function(x) na.omit(x$EntrezID))
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

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/A-CTest_IRintrons_KEGG_GO_DO_plots.pdf", width=16, height=9)
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
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/A-CTest.intron.kegg.GO.DO.objects.rda")



## Read in GLM results

IDs = lapply(dIR, function(x) unlist(strsplit(rownames(x), "/", fixed=T), recursive = F))
dIR = Map(cbind, dIR, ensID = lapply(IDs, function(x) x[grep("ENSG", x)]), coord = lapply(IDs, function(x) x[grep(":", x, fixed=T)]),
          geneSymbol = lapply(IDs, function(x) x[-c(grep("ENSG", x),grep(":", x, fixed=T),grep("clean", x))]),
          A.Sign = lapply(dIR, function(y) ifelse((y$Ad.IRratio.ZoneCytosol - y$Ad.IRratio.ZoneNucleus) < 0,"MoreIRInNuc", "MoreIRInCyt")),
          P.Sign = lapply(dIR, function(y) ifelse((y$Pren.IRratio.ZoneCytosol - y$Pren.IRratio.ZoneNucleus) < 0,"MoreIRInNuc", "MoreIRInCyt")),
          C.Sign = lapply(dIR, function(y) ifelse((y$Cyt.IRratio.FetalAdult - y$Cyt.IRratio.FetalPrenatal) < 0,"MoreIRInPrenatal", "MoreIRInAdult")),
          N.Sign = lapply(dIR, function(y) ifelse((y$Nuc.IRratio.FetalAdult - y$Nuc.IRratio.FetalPrenatal) < 0,"MoreIRInPrenatal", "MoreIRInAdult")))
head(dIR[[1]])
dIR = Map(cbind, dIR, Sign = list(ifelse(dIR$byFractionInAdult$A.Sign=="MoreIRInNuc", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byFractionInPrenatal$P.Sign=="MoreIRInNuc", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byAgeInCytosol$C.Sign=="MoreIRInPrenatal", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byAgeInNucleus$N.Sign=="MoreIRInPrenatal", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult")),
          EntrezID = lapply(dIR, function(x) geneMap[match(x$ensID, geneMap$ensemblID),"EntrezID"]))
sigIR = lapply(dIR, function(x) x[which(x$padj<=0.05),])
elementNROWS(sigIR)
sigIR = unlist(lapply(sigIR, function(x) split(x, x$Sign)), recursive=F)
names(sigIR) = c("Cytosolic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytosolic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                  "Adult Retention\nIn Cytosol","Prenatal Retention\nIn Cytosol","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus")
entrezID.allGroups = lapply(sigIR, function(x) na.omit(x$EntrezID))
entrezID.collapsed = lapply(list("Cytosolic Retention" = do.call(rbind, sigIR[grep("Cytosolic Retention",names(sigIR))]),
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


