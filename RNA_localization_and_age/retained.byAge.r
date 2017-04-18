library(ggplot2)
library(plyr)
library(clusterProfiler)
require(org.Hs.eg.db)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Plot Fraction:Age interaction genes
Idds = DESeqDataSetFromMatrix(countData = geneCounts, colData = pd, design = ~ Library + Fetal + Zone + Fetal:Zone)
Idds = DESeq(Idds)
sigres = data.frame(Ires[which(Ires$padj<=0.05 & abs(Ires$log2FoldChange)>=1),])
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/interaction_genes.pdf")
plots = list()
for (i in 1:nrow(sigres)){
  plots[[i]] = plotCounts(Idds, as.character(rownames(sigres[i,])), 
                               intgroup = c("Fetal", "Zone", "Library"), returnData =TRUE)
  tmp = plots[[i]]
  tmp$Group = paste(tmp$Fetal,tmp$Zone, sep = "\n")
  x = ggplot(tmp, aes(x=Group, y=count)) + geom_boxplot() + 
    geom_jitter() +
    scale_y_log10(breaks=c(25,100,400)) +
    ylab("Normalized Count") + 
    xlab("") +
    ggtitle(paste0(as.character(geneMap[match(rownames(sigres[i,]),geneMap$gencodeID),"ensemblID"]),":",
                   as.character(geneMap[match(rownames(sigres[i,]),geneMap$gencodeID),"Symbol"]))) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(x)
}
dev.off()

# Make list of genes
all = list(Ires = data.frame(Ires[order(rownames(Ires)),]), Ares = data.frame(Ares[order(rownames(Ares)),]), 
           Fres = data.frame(Fres[order(rownames(Fres)),]))
Sign = lapply(all, function(x) ifelse(x$log2FoldChange>0, "Pos", "Neg"))
all = Map(cbind, all, Sign = Sign)
Ires = all[["Ires"]]
Ares = all[["Ares"]]
Fres = all[["Fres"]]

AdPos = as.character(rownames(Ares[which(Ares$Sign=="Pos"),]))
AdNeg = as.character(rownames(Ares[which(Ares$Sign=="Neg"),]))
AdSig = as.character(rownames(Ares[which(Ares$padj<=0.05),]))
AdLFC = as.character(rownames(Ares[which(abs(Ares$log2FoldChange)>=1),]))
FetPos = as.character(rownames(Fres[which(Fres$Sign=="Pos"),]))
FetNeg = as.character(rownames(Fres[which(Fres$Sign=="Neg"),]))
FetSig = as.character(rownames(Fres[which(Fres$padj<=0.05),]))
FetLFC = as.character(rownames(Fres[which(abs(Fres$log2FoldChange)>=1),]))

sig = list(both_retained = Ires[which(rownames(Ires)%in%AdPos & rownames(Ires)%in%FetPos 
                                     & rownames(Ires)%in%AdSig & rownames(Ires)%in%FetSig),],
           both_exported = Ires[which(rownames(Ires)%in%AdNeg & rownames(Ires)%in%FetNeg 
                                     & rownames(Ires)%in%AdSig & rownames(Ires)%in%FetSig),],
           Fet_retained = Ires[which(rownames(Ires)%in%FetPos & !(rownames(Ires)%in%AdSig) & rownames(Ires)%in%FetSig),],
           Ad_retained = Ires[which(rownames(Ires)%in%AdPos & rownames(Ires)%in%AdSig & !(rownames(Ires)%in%FetSig)),],
           Fet_exported = Ires[which(rownames(Ires)%in%FetNeg & !(rownames(Ires)%in%AdSig) & rownames(Ires)%in%FetSig),],
           Ad_exported = Ires[which(rownames(Ires)%in%AdNeg & rownames(Ires)%in%AdSig & !(rownames(Ires)%in%FetSig)),],
           ret_Ad_exp_Fet = Ires[which(rownames(Ires)%in%AdPos & rownames(Ires)%in%FetNeg 
                                      & rownames(Ires)%in%AdSig & rownames(Ires)%in%FetSig),],
           ret_Fet_exp_Ad = Ires[which(rownames(Ires)%in%FetPos & rownames(Ires)%in%AdNeg 
                                      & rownames(Ires)%in%AdSig & rownames(Ires)%in%FetSig),],
           interacting = Ires[which(Ires$padj<=0.05),])
elementNROWS(sig)
gene = geneMap[which(rownames(geneMap)%in%rownames(Ires)),]
gene = geneMap[order(rownames(geneMap)),]
data = data.frame(geneID = rownames(Ires), baseMean = Ires$baseMean, Prenatal.LFC = Fres$log2FoldChange, 
                  Prenatal.SE = Fres$lfcSE, Prenatal.padj = Fres$padj,Adult.LFC = Ares$log2FoldChange, 
                  Adult.SE = Ares$lfcSE, Adult.padj = Ares$padj, 
                  Symbol = gene[match(rownames(Ires),as.character(rownames(gene))),"Symbol"],
                  EntrezID = gene[match(as.character(rownames(Ires)),gene$gencodeID),"EntrezID"],
                  Type = gene[match(as.character(rownames(Ires)),gene$gencodeID),"gene_type"])
sig = lapply(sig, function(x) data[which(data$geneID%in%rownames(x)),])
sig.1 = lapply(sig, function(x) x[which(abs(x$Prenatal.LFC)>=1 | abs(x$Adult.LFC)>=1),])
save(Ires,Fres,Ares,sig,sig.1, 
     file = "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.rda")

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
type = data.frame(RNA.Type = as.character(c(rep.int(RNA.Type[1], 7), rep.int(RNA.Type[2], 7),
                               rep.int(RNA.Type[3], 7),rep.int(RNA.Type[4], 7))),
                  Count = NA, Group = factor(x=rep.int(group, 4)))

for (i in 1:length(group)){
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="protein_coding" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="TR_C_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="polymorphic_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TR_V_gene" & TypeFreq$Group==group[i]),2])
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
        TypeFreq[which(TypeFreq$RNA_Type=="snRNA" & TypeFreq$Group==group[i]),2])
}
group = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only",
          "Exported:\nPrenatal Only", "Exported:\nAdult Only",
          "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nPrenatal Only",
                                                    "Exported:\nPrenatal Only", "Exported:\nAdult Only", "Retained: Adult/\nExported: Prenatal",
                                                    "Retained: Prenatal/\nExported: Adult", "Interaction"))

# Graph the Frequencies
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

## Limiting to 1 LFC

freq = lapply(sig.1, function(x) count(x$Type))
TypeFreq = do.call(rbind, freq)
TypeFreq$Group = gsub("\\..*","", rownames(TypeFreq))
colnames(TypeFreq) = c("RNA_Type", "Count", "Group")
# Condense to 4 groups
group = as.character(names(freq))
RNA.Type = c("Protein-coding", "Pseudogene", "Long Non-coding", "Short Non-coding")
type = data.frame(RNA.Type = as.character(c(rep.int(RNA.Type[1], 7), rep.int(RNA.Type[2], 7),
                                            rep.int(RNA.Type[3], 7),rep.int(RNA.Type[4], 7))),
                  Count = NA, Group = factor(x=rep.int(group, 4)))
for (i in 1:length(group)){
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$RNA_Type=="protein_coding" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$RNA_Type=="TR_C_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="polymorphic_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="IG_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$RNA_Type=="TR_V_gene" & TypeFreq$Group==group[i]),2])
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
        TypeFreq[which(TypeFreq$RNA_Type=="snRNA" & TypeFreq$Group==group[i]),2])
}
group = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only","Exported:\nPrenatal Only", "Exported:\nAdult Only",
          "Retained: Adult/\nExported: Prenatal", "Retained: Prenatal/\nExported: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nPrenatal Only",
                                                    "Exported:\nPrenatal Only", "Exported:\nAdult Only","Retained: Adult/\nExported: Prenatal",
                                                    "Retained: Prenatal/\nExported: Adult", "Interaction"))
# Graph the Frequencies
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nLog2 Fold Change >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.rda")
## Gene Ontology
names(sig.1) = names(sig) = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only",
                              "Exported:\nPrenatal Only", "Exported:\nAdult Only","Retained: Adult/\nExported: Prenatal",
                              "Retained: Prenatal/\nExported: Adult", "Interaction")
entrezID = lapply(sig.1, function(x) na.omit(x$EntrezID))
# Define universe as all genes expressed in each of the four groups
GeneUniverse = as.character(unique(geneMap[match(rownames(Ires),geneMap$gencodeID),"EntrezID"]))
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
names = c("Retained_Both", "Exported_Both", "Retained_Prenatal_Only", "Retained_Adult_Only",
          "Exported_Prenatal_Only", "Exported_Adult_Only","Retained_Adult.Exported_Prenatal",
          "Retained_Prenatal.Exported_Adult", "Interaction")
lists = list(DO = goListdf_DO, CC = goListdf_CC, BP = goListdf_BP, MF = goListdf_MF, KEGG = keggListdf)
for (i in 1:length(lists)){
  for (j in 1:length(names)){
    if (nrow(lists[[i]][[j]])>0){
      lists[[i]][[j]] = data.frame(lists[[i]][[j]], Comparison = names[j])
    }}}
lists = lapply(lists, function(x) do.call(rbind, x))
write.csv(lists[["KEGG"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.KEGG.csv")
write.csv(lists[["BP"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.BP.csv")
write.csv(lists[["MF"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.MF.csv")
write.csv(lists[["CC"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.GO.CC.csv")
write.csv(lists[["DO"]], file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Interaction.DO.csv")

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrezID, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
# Biological Process
compareBP = compareCluster(entrezID, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareBPDropped = dropGO(compareBPDropped,  level = 3)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 30, title= "Biological Process GO Enrichment")
# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareMFDropped = dropGO(compareMFDropped,  level = 3)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 30, title= "Molecular Function GO Enrichment")
# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
compareCCDropped = dropGO(compareCCDropped,  level = 3)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 30, title= "Cellular Compartment GO Enrichment")
# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
entrez.noLFC = lapply(sig, function(x) na.omit(x$EntrezID))
compareDO = compareCluster(entrez.noLFC, fun="enrichDO",
                           ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.rda")