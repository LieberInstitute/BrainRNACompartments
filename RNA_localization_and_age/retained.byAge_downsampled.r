library(ggplot2)
library(plyr)
library(clusterProfiler)
require(org.Hs.eg.db)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Plot Fraction:Age interaction genes
Ipdds.down = DESeqDataSetFromMatrix(countData = geneCounts.down[,grep("polyA", colnames(geneCounts.down))], 
                                    colData = pd[which(pd$Library=="polyA"),], design = ~ Fetal + Zone + Fetal:Zone)
Ipdds.down = DESeq(Ipdds.down)
Ipres.down = results(Ipdds.down)
sigres.down = data.frame(Ipres.down[which(Ipres.down$padj<=0.05 & abs(Ipres.down$log2FoldChange)>=1),])
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/interaction_genes.polyAonly.pdf")
plots = list()
for (i in 1:nrow(sigres.down)){
  plots[[i]] = plotCounts(Ipdds.down, as.character(rownames(sigres.down[i,])), 
                               intgroup = c("Fetal", "Zone"), returnData =TRUE)
  tmp = plots[[i]]
  x = ggplot(tmp, aes(x=Fetal, y=count, fill=Zone)) + geom_boxplot() + 
    geom_jitter() +
    scale_y_log10(breaks=c(25,100,400)) +
    ylab("Normalized Count") + 
    xlab("") +
    ggtitle(paste0(as.character(geneMap[match(rownames(sigres.down[i,]),geneMap$gencodeID),"ensemblID"]),": ",
                   as.character(geneMap[match(rownames(sigres.down[i,]),geneMap$gencodeID),"Symbol"]))) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(x)
}
dev.off()

# Make list of genes
all = list(Ipres.down = data.frame(Ipres.down[order(rownames(Ipres.down)),]), 
           Apres = data.frame(Apres[order(rownames(Apres)),]), 
           Fpres.down = data.frame(Fpres.down[order(rownames(Fpres.down)),]))
all = Map(cbind, all, Sign = lapply(all, function(x) ifelse(x$log2FoldChange>0, "Pos", "Neg")))
Ipres.down = all[["Ipres.down"]]
Apres = all[["Apres"]]
Fpres.down = all[["Fpres.down"]]

AdPos = as.character(rownames(Apres[which(Apres$Sign=="Pos"),]))
AdNeg = as.character(rownames(Apres[which(Apres$Sign=="Neg"),]))
AdSig = as.character(rownames(Apres[which(Apres$padj<=0.05),]))
AdLFC = as.character(rownames(Apres[which(abs(Apres$log2FoldChange)>=1),]))
FetPos = as.character(rownames(Fpres.down[which(Fpres.down$Sign=="Pos"),]))
FetNeg = as.character(rownames(Fpres.down[which(Fpres.down$Sign=="Neg"),]))
FetSig = as.character(rownames(Fpres.down[which(Fpres.down$padj<=0.05),]))
FetLFC = as.character(rownames(Fpres.down[which(abs(Fpres.down$log2FoldChange)>=1),]))

sig = list(both_retained = Ipres.down[which(rownames(Ipres.down)%in%AdPos & rownames(Ipres.down)%in%FetPos 
                                     & rownames(Ipres.down)%in%AdSig & rownames(Ipres.down)%in%FetSig),],
           both_exported = Ipres.down[which(rownames(Ipres.down)%in%AdNeg & rownames(Ipres.down)%in%FetNeg 
                                     & rownames(Ipres.down)%in%AdSig & rownames(Ipres.down)%in%FetSig),],
           Fet_retained = Ipres.down[which(rownames(Ipres.down)%in%FetPos & !(rownames(Ipres.down)%in%AdSig) & 
                                       rownames(Ipres.down)%in%FetSig),],
           Ad_retained = Ipres.down[which(rownames(Ipres.down)%in%AdPos & rownames(Ipres.down)%in%AdSig & !(rownames(Ipres.down)%in%FetSig)),],
           Fet_exported = Ires[which(rownames(Ipres.down)%in%FetNeg & !(rownames(Ipres.down)%in%AdSig) & rownames(Ipres.down)%in%FetSig),],
           Ad_exported = Ires[which(rownames(Ipres.down)%in%AdNeg & rownames(Ipres.down)%in%AdSig & !(rownames(Ipres.down)%in%FetSig)),],
           ret_Ad_exp_Fet = Ipres.down[which(rownames(Ipres.down)%in%AdPos & rownames(Ipres.down)%in%FetNeg 
                                      & rownames(Ipres.down)%in%AdSig & rownames(Ipres.down)%in%FetSig),],
           ret_Fet_exp_Ad = Ipres.down[which(rownames(Ipres.down)%in%FetPos & rownames(Ipres.down)%in%AdNeg 
                                      & rownames(Ipres.down)%in%AdSig & rownames(Ipres.down)%in%FetSig),],
           interacting = Ipres.down[which(Ipres.down$padj<=0.05),])

data = data.frame(geneID = rownames(Ipres.down), baseMean = Ipres.down$baseMean, 
                  Prenatal.LFC = Fpres.down[match(rownames(Ipres.down), rownames(Fpres.down)), "log2FoldChange"], 
                  Prenatal.SE = Fpres.down[match(rownames(Ipres.down), rownames(Fpres.down)), "lfcSE"], 
                  Prenatal.padj = Fpres.down[match(rownames(Ipres.down), rownames(Fpres.down)), "padj"],
                  Adult.LFC = Apres[match(rownames(Ipres.down), rownames(Apres)), "log2FoldChange"], 
                  Adult.SE = Apres[match(rownames(Ipres.down), rownames(Apres)), "lfcSE"],
                  Adult.padj = Apres[match(rownames(Ipres.down), rownames(Apres)), "padj"],
                  ensID = geneMap[match(rownames(Ipres.down),as.character(rownames(geneMap))),"ensemblID"],
                  Symbol = geneMap[match(rownames(Ipres.down),as.character(rownames(geneMap))),"Symbol"],
                  EntrezID = geneMap[match(as.character(rownames(Ipres.down)),geneMap$gencodeID),"EntrezID"],
                  Type = geneMap[match(as.character(rownames(Ipres.down)),geneMap$gencodeID),"gene_type"])
sig = lapply(sig, function(x) data[which(data$geneID %in% rownames(x)),])
sig.1 = lapply(sig, function(x) x[which(abs(x$Prenatal.LFC)>=1 | abs(x$Adult.LFC)>=1),])
save(Ipres.down,Fpres.down,Apres,sig,sig.1,geneMap, 
     file = "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Annotate genes
freq = lapply(sig, function(x) count(x$Type))
TypeFreq = do.call(rbind, freq)
TypeFreq$Group = gsub("\\..*","", rownames(TypeFreq))
colnames(TypeFreq) = c("RNA_Type", "Count", "Group")

# Graph the Frequencies
ggplot(TypeFreq, aes(x = Group, y = Count, fill = RNA_Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("DEG By Library in the Nucleus") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Condense to 4 groups
group = as.character(names(freq))
RNA.Type = c("Protein-coding", "Pseudogene", "Long Non-coding", "Short Non-coding")
type = data.frame(RNA.Type = as.character(c(rep.int(RNA.Type[1], 9), rep.int(RNA.Type[2], 9),
                               rep.int(RNA.Type[3], 9),rep.int(RNA.Type[4], 9))),
                  Count = NA, Group = factor(x=rep.int(group, 4)))

for (i in 1:length(group)){
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="protein_coding" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="TR_C_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="polymorphic_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TR_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_J_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TR_J_gene" & TypeFreq$Group==group[i]),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="transcribed_processed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="unprocessed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="transcribed_unprocessed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="processed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_C_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="transcribed_unitary_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="unitary_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_V_pseudogene" & TypeFreq$Group==group[i]),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="3prime_overlapping_ncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="antisense" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="non_coding" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="lincRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="sense_intronic" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="misc_RNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="processed_transcript" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="bidirectional_promoter_lncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TEC" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="macro_lncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="sense_overlapping" & TypeFreq$Group==group[i]),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="miRNA" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="Mt_rRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="Mt_tRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="rRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="snoRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="scRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="snRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="vaultRNA" & TypeFreq$Group==group[i]),2])
}
group = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only",
          "Exported:\nPrenatal Only", "Exported:\nAdult Only",
               "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nPrenatal Only",
                                                    "Exported:\nPrenatal Only", "Exported:\nAdult Only",
                                                    "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction"))

# Graph the Frequencies
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_fraction-age.pdf", height = 7, width = 8)
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

## Limiting to 1 LFC

freq = lapply(sig.1, function(x) count(x$Type))
TypeFreq = do.call(rbind, freq)
TypeFreq$Group = gsub("\\..*","", rownames(TypeFreq))
colnames(TypeFreq) = c("RNA_Type", "Count", "Group")
# Condense to 4 groups
group = as.character(names(freq))
RNA.Type = c("Protein-coding", "Pseudogene", "Long Non-coding", "Short Non-coding")
type = data.frame(RNA.Type = as.character(c(rep.int(RNA.Type[1], 9), rep.int(RNA.Type[2],9),
                                            rep.int(RNA.Type[3], 9),rep.int(RNA.Type[4], 9))),
                  Count = NA, Group = factor(x=rep.int(group, 4)))
for (i in 1:length(group)){
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="protein_coding" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="TR_C_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="polymorphic_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TR_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_J_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TR_J_gene" & TypeFreq$Group==group[i]),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="transcribed_processed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="unprocessed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="transcribed_unprocessed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="processed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_C_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="transcribed_unitary_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="unitary_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_V_pseudogene" & TypeFreq$Group==group[i]),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="3prime_overlapping_ncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="antisense" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="non_coding" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="lincRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="sense_intronic" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="misc_RNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="processed_transcript" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="bidirectional_promoter_lncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TEC" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="macro_lncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="sense_overlapping" & TypeFreq$Group==group[i]),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="miRNA" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="Mt_rRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="Mt_tRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="rRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="snoRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="scRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="snRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="vaultRNA" & TypeFreq$Group==group[i]),2])
}
group = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only",
          "Exported:\nPrenatal Only", "Exported:\nAdult Only",
          "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nPrenatal Only",
                                                    "Exported:\nPrenatal Only", "Exported:\nAdult Only",
                                                    "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction"))
# Graph the Frequencies
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_fraction-age_LFC1.pdf", height = 7, width = 8)
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

## Gene Ontology
names(sig.1) = names(sig) = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only",
                              "Exported:\nPrenatal Only", "Exported:\nAdult Only","Retained: Adult/\nExported: Prenatal",
                              "Retained: Prenatal/\nExported: Adult", "Interaction")
entrezID = lapply(sig.1, function(x) na.omit(x$EntrezID))
# Define universe as all genes expressed in each of the four groups
GeneUniverse = as.character(unique(geneMap[match(rownames(Ipres.down),geneMap$gencodeID),"EntrezID"]))
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
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process GO Enrichment")
# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function GO Enrichment")
# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment GO Enrichment")
# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
entrez.noLFC = lapply(sig, function(x) na.omit(x$EntrezID))
compareDO = compareCluster(entrez.noLFC, fun="enrichDO",
                           ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 40, title= "Disease Ontology Enrichment")

save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.polyAonly.downsampled.rda")

### In the different groups of fraction regulated genes, is there a relationship between direction of expression over age and significance over age?
AgebyFrac = list(Cpres.down = data.frame(Cpres.down), Npres = data.frame(Npres))
AgebyFrac = Map(cbind, AgebyFrac,lapply(AgebyFrac, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Cytosol", "Nucleus"))
AgebyFrac = do.call(rbind, AgebyFrac)
AgebyFrac$FDR = ifelse(AgebyFrac$padj<=0.05, "FDR<0.05", "FDR>0.05")
AgebyFrac = AgebyFrac[which(AgebyFrac$padj!="NA"),]

fracDevel = lapply(sig, function(x) AgebyFrac[which(AgebyFrac$gencodeID %in% x$geneID),])
fracDevel = do.call(rbind, fracDevel)
fracDevel$fracReg = gsub("\\..*","", rownames(fracDevel))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nPrenatal Only",
                                      "Exported:\nPrenatal Only", "Exported:\nAdult Only",
                                      "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction"))

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/RetainedbyAge_LFCxFDR.pdf", width=24, height=6)
ggplot(fracDevel, aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ fracReg) +
  ylab("Log2 Fold Change") + 
  xlab("") +
  ggtitle(paste0("Age Expression Changes in Gene Groups Regulated by Fraction")) + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

### Is there a relationship between direction of expression by age and significance in groups of genes differentially regulated by fraction?
# In cytosol:
fracDevel = lapply(sig, function(x) AgebyFrac[which(AgebyFrac$gencodeID %in% x$geneID),])
x = lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                                                       length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
                                                     c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                                                       length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))))))
x = unlist(lapply(x,function(x) as.character(x)[1]))
x = data.frame(x)
CytCounts = do.call(rbind, lapply(fracDevel, function(x) data.frame(Sig = c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                                                       length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
                                                Nonsig = c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                                                       length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))),
                                                Dir = c("UpPrenatal", "UpAdult"))))
CytCounts$FracGroup = gsub("\\..*","", rownames(CytCounts))
CytCounts$fisher.pval = x[match(CytCounts$FracGroup,rownames(x)),]
CytCounts$Comparison = "Cytosol"

# In Nucleus:
x = lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                                                           length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
                                                         c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                                                           length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))))))
x = unlist(lapply(x,function(x) as.character(x)[1]))
x = data.frame(x)
NucCounts = do.call(rbind, lapply(fracDevel, function(x) data.frame(Sig = c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                                                                            length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
                                                                    Nonsig = c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                                                                               length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))),
                                                                    Dir = c("UpPrenatal", "UpAdult"))))
NucCounts$FracGroup = gsub("\\..*","", rownames(NucCounts))
NucCounts$fisher.pval = x[match(NucCounts$FracGroup,rownames(x)),]
NucCounts$Comparison = "Nucleus"
write.csv(rbind(CytCounts,NucCounts), 
          file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/AgeLFCxAgepval_bysigFracGroup_fisher_pvalues.csv",
          row.names = F,quote = F)