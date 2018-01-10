library(ggplot2)
library(plyr)
library(clusterProfiler)
require(org.Hs.eg.db)
library(DESeq2)
library(GenomicRanges)
library(data.table)

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
           Fet_retained = Ipres.down[which(rownames(Ipres.down)%in%FetPos & !(rownames(Ipres.down)%in%AdSig) & rownames(Ipres.down)%in%FetSig),],
           Ad_retained = Ipres.down[which(rownames(Ipres.down)%in%AdPos & rownames(Ipres.down)%in%AdSig & !(rownames(Ipres.down)%in%FetSig)),],
           Fet_exported = Ipres.down[which(rownames(Ipres.down)%in%FetNeg & !(rownames(Ipres.down)%in%AdSig) & rownames(Ipres.down)%in%FetSig),],
           Ad_exported = Ipres.down[which(rownames(Ipres.down)%in%AdNeg & rownames(Ipres.down)%in%AdSig & !(rownames(Ipres.down)%in%FetSig)),],
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
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
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
                                                    "Exported:\nAdult Only","Exported:\nPrenatal Only",
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
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_fraction-age_LFC1.pdf", height = 6, width = 8)
ggplot(type[which(type$Group!="Retained: Adult/\nExported: Prenatal" &
                    type$Group!="Retained: Prenatal/\nExported: Adult"),], 
       aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


type = data.table(type[which(type$Group!="Retained: Adult/\nExported: Prenatal" &
                           type$Group!="Retained: Prenatal/\nExported: Adult"),])
x = data.frame(type[,sum(Count), by="Group"])
type$sum = x[match(type$Group, x$Group),"V1"]
type$perc = round(type$Count/type$sum*100,1)
type = ddply(type, .(Group), transform, pos = cumsum(perc) - (0.5 * perc))
type$pos = 100

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_fraction-age_LFC1.percent.pdf", height = 6, width = 8)
ggplot(type, aes(x = Group, y = perc, fill = RNA.Type)) + 
  geom_bar(stat = "identity") +
  geom_text(data=type[type$RNA.Type=="Long Non-coding",], 
            aes(x = Group, y = pos, label = sum), size=4, nudge_y = 5) +
  coord_flip() + labs(fill="") + ylab("Percent") + xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

r = type[grep("Retained:",type$Group),]
e = type[grep("Exported:",type$Group),]
fisher.test(data.frame(c(sum(r[r$RNA.Type=="Protein-coding","Count"]), sum(r$Count)-sum(r[r$RNA.Type=="Protein-coding","Count"])),
                       c(sum(e[e$RNA.Type=="Protein-coding","Count"]), sum(e$Count)-sum(e[r$RNA.Type=="Protein-coding","Count"]))))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1900531 0.3341147
#sample estimates:
#  odds ratio 
#0.2530526 

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
# Enriched Molecular Function GOs
goList_MF = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
# Biological Process GO enrichment
goList_BP = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
# Cellular Compartment GO enrichment
goList_CC = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",
                                                  qvalueCutoff=1))
# Disease Ontology
goList_DO = lapply(entrezID, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
              minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))


# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrezID, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrezID, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Save
save(keggList, goList_MF, goList_BP, goList_CC, goList_DO,compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.polyAonly.sig1.downsampled.rda")


## plot
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/KEGG_interaction_polyAonly.sig1.downsampled.pdf", width=12,height=12)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/MF_interaction_polyAonly.sig1.downsampled.pdf", width=14,height=16)
plot(compareMF,colorBy="p.adjust",  showCategory = 85, title= "Molecular Function GO Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/BP_interaction_polyAonly.sig1.downsampled.pdf", width=12,height=64)
plot(compareBP,colorBy="p.adjust",  showCategory = 400, title= "Biological Process GO Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/CC_interaction_polyAonly.sig1.downsampled.pdf", width=12,height=16)
plot(compareCC,colorBy="p.adjust",  showCategory = 400, title= "Cellular Compartment GO Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/DO_interaction_polyAonly.sig1.downsampled.pdf", width=6,height=4)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
dev.off()


### In the different groups of fraction regulated genes, is there a relationship between direction of expression over age and significance over age?
AgebyFrac = list(Cpres.down = data.frame(Cpres.down), Npres = data.frame(Npres))
AgebyFrac = Map(cbind, AgebyFrac,lapply(AgebyFrac, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Cytosol", "Nucleus"))
AgebyFrac = do.call(rbind, AgebyFrac)
AgebyFrac$FDR = ifelse(AgebyFrac$padj<=0.05, "FDR<0.05", "FDR>0.05")
AgebyFrac = AgebyFrac[which(AgebyFrac$padj!="NA"),]

elementNROWS(sig)
fracDevel = lapply(sig[elementNROWS(sig)>0], function(x) AgebyFrac[which(AgebyFrac$gencodeID %in% x$geneID),])
names(fracDevel) = c("Retained: Both", "Exported: Both", "Retained:\nPrenatal Only", "Retained:\nAdult Only", 
                     "Exported:\nPrenatal Only", "Exported:\nAdult Only","Retained: Adult/\nExported: Prenatal", "Interaction")
fracDevel = do.call(rbind, Map(cbind, fracDevel, fracReg = as.list(names(fracDevel))))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nPrenatal Only",
                                      "Exported:\nAdult Only","Exported:\nPrenatal Only","Retained: Adult/\nExported: Prenatal", "Interaction"))

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/RetainedbyAge_LFCxFDR.pdf", width=20, height=6)
ggplot(fracDevel[fracDevel$fracReg!="Retained: Adult/\nExported: Prenatal",], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
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

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/AdultRetainednExported_Age-LFCxFDR.pdf", width=7, height=5)
ggplot(fracDevel[which(fracDevel$fracReg=="Retained:\nAdult Only" | fracDevel$fracReg=="Exported:\nAdult Only"),], 
       aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() +
  facet_grid(. ~ fracReg) +
  ylab("Log2 Fold Change") + 
  xlab("") +
  ggtitle(paste0("Age Expression Changes\nby Fraction")) + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


### Is there a relationship between direction of expression by age and significance in groups of genes differentially regulated by fraction?

fracDevel = lapply(sig, function(x) AgebyFrac[which(AgebyFrac$gencodeID %in% x$geneID),])
elementNROWS(sig)-elementNROWS(lapply(fracDevel, function(x) x[x$Comparison=="Cytosol",]))
elementNROWS(sig)-elementNROWS(lapply(fracDevel, function(x) x[x$Comparison=="Nucleus",]))

# In cytosol:
tb = lapply(fracDevel[names(fracDevel)!="ret_Fet_exp_Ad"], function(x) 
  data.frame(Sig = c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                     length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
             Nonsig = c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                        length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))),row.names = c("Decreasing","Increasing")))
df = rbind(pvalue = unlist(lapply(lapply(tb, fisher.test), function(x) x$p.value)),OR = unlist(lapply(lapply(tb, fisher.test),function(x) x$estimate))) 
CytCounts = do.call(rbind, tb)
CytCounts$FracGroup = gsub("\\..*","", rownames(CytCounts))
CytCounts$fisher.pval = df["pvalue",match(CytCounts$FracGroup,colnames(df))]
CytCounts$fisher.OR = df["OR",match(CytCounts$FracGroup,colnames(df))]
CytCounts$Comparison = "Cytosol"
CytCounts$padj = p.adjust(CytCounts$fisher.pval, method="fdr")


# In Nucleus:
tb = lapply(fracDevel[names(fracDevel)!="ret_Fet_exp_Ad"], function(x) 
  data.frame(Sig = c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                     length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
             Nonsig = c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                        length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))),row.names = c("Decreasing","Increasing")))
df = rbind(pvalue = unlist(lapply(lapply(tb, fisher.test), function(x) x$p.value)),OR = unlist(lapply(lapply(tb, fisher.test),function(x) x$estimate))) 
NucCounts = do.call(rbind, tb)
NucCounts$FracGroup = gsub("\\..*","", rownames(NucCounts))
NucCounts$fisher.pval = df["pvalue",match(NucCounts$FracGroup,colnames(df))]
NucCounts$fisher.OR = df["OR",match(NucCounts$FracGroup,colnames(df))]
NucCounts$Comparison = "Nucleus"
NucCounts$padj = p.adjust(NucCounts$fisher.pval, method="fdr")

write.csv(rbind(CytCounts,NucCounts),quote = F, 
          file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/AgeLFCxAgepval_bysigFracGroup_fisher.csv")

CytCounts[which(CytCounts$padj<=0.05),]
#                        Sig Nonsig   FracGroup  fisher.pval fisher.OR Comparison         padj
#Ad_retained.Decreasing  596    396 Ad_retained 3.374790e-03 1.3204695    Cytosol 1.349916e-02
#Ad_retained.Increasing  449    394 Ad_retained 3.374790e-03 1.3204695    Cytosol 1.349916e-02
#Ad_exported.Decreasing  384    355 Ad_exported 1.739493e-19 0.4412216    Cytosol 1.391595e-18
#Ad_exported.Increasing 1263    515 Ad_exported 1.739493e-19 0.4412216    Cytosol 1.391595e-18
#interacting.Decreasing   84     29 interacting 8.858470e-03 0.4119931    Cytosol 2.362259e-02
#interacting.Increasing  127     18 interacting 8.858470e-03 0.4119931    Cytosol 2.362259e-02
NucCounts[which(NucCounts$padj<=0.05),]
#                         Sig Nonsig     FracGroup  fisher.pval fisher.OR Comparison         padj
#both_retained.Decreasing  37     47 both_retained 5.761127e-03 0.4615974    Nucleus 2.304451e-02
#both_retained.Increasing  89     52 both_retained 5.761127e-03 0.4615974    Nucleus 2.304451e-02
#Ad_retained.Decreasing   311    304   Ad_retained 2.511130e-11 0.5085722    Nucleus 2.008904e-10
#Ad_retained.Increasing   817    406   Ad_retained 2.511130e-11 0.5085722    Nucleus 2.008904e-10


# Make Age Sig object
age = list(Ipres.down = data.frame(Ipres.down[order(rownames(Ipres.down)),]), 
           Npres = data.frame(Npres[order(rownames(Npres)),]), 
           Cpres.down = data.frame(Cpres.down[order(rownames(Cpres.down)),]))
age = Map(cbind, age, Sign = lapply(age, function(x) ifelse(x$log2FoldChange>0, "Pos", "Neg")))
Ipres.down = age[["Ipres.down"]]
Npres = age[["Npres"]]
Cpres.down = age[["Cpres.down"]]

NucPos = as.character(rownames(Npres[which(Npres$Sign=="Pos"),]))
NucNeg = as.character(rownames(Npres[which(Npres$Sign=="Neg"),]))
NucSig = as.character(rownames(Npres[which(Npres$padj<=0.05),]))
NucLFC = as.character(rownames(Npres[which(abs(Npres$log2FoldChange)>=1),]))
CytPos = as.character(rownames(Cpres.down[which(Cpres.down$Sign=="Pos"),]))
CytNeg = as.character(rownames(Cpres.down[which(Cpres.down$Sign=="Neg"),]))
CytSig = as.character(rownames(Cpres.down[which(Cpres.down$padj<=0.05),]))
CytLFC = as.character(rownames(Cpres.down[which(abs(Cpres.down$log2FoldChange)>=1),]))

age.sig = list(both_decreasing = Ipres.down[which(rownames(Ipres.down)%in%NucPos & rownames(Ipres.down)%in%CytPos 
                                            & rownames(Ipres.down)%in%NucSig & rownames(Ipres.down)%in%CytSig),],
           both_increasing = Ipres.down[which(rownames(Ipres.down)%in%NucNeg & rownames(Ipres.down)%in%CytNeg 
                                            & rownames(Ipres.down)%in%NucSig & rownames(Ipres.down)%in%CytSig),],
           Cyt_decreasing = Ipres.down[which(rownames(Ipres.down)%in%CytPos & !(rownames(Ipres.down)%in%NucSig) & rownames(Ipres.down)%in%CytSig),],
           Nuc_decreasing = Ipres.down[which(rownames(Ipres.down)%in%NucPos & rownames(Ipres.down)%in%NucSig & !(rownames(Ipres.down)%in%CytSig)),],
           Cyt_increasing = Ipres.down[which(rownames(Ipres.down)%in%CytNeg & !(rownames(Ipres.down)%in%NucSig) & rownames(Ipres.down)%in%CytSig),],
           Nuc_increasing = Ipres.down[which(rownames(Ipres.down)%in%NucNeg & rownames(Ipres.down)%in%NucSig & !(rownames(Ipres.down)%in%CytSig)),],
           decr_Nuc_incr_Cyt = Ipres.down[which(rownames(Ipres.down)%in%NucPos & rownames(Ipres.down)%in%CytNeg 
                                             & rownames(Ipres.down)%in%NucSig & rownames(Ipres.down)%in%CytSig),],
           decr_Cyt_incr_Nuc = Ipres.down[which(rownames(Ipres.down)%in%CytPos & rownames(Ipres.down)%in%NucNeg 
                                             & rownames(Ipres.down)%in%NucSig & rownames(Ipres.down)%in%CytSig),],
           interacting = Ipres.down[which(Ipres.down$padj<=0.05),])

age.data = data.frame(geneID = rownames(Ipres.down), baseMean = Ipres.down$baseMean, 
                  Cytosol.LFC = Cpres.down[match(rownames(Ipres.down), rownames(Cpres.down)), "log2FoldChange"], 
                  Cytosol.SE = Cpres.down[match(rownames(Ipres.down), rownames(Cpres.down)), "lfcSE"], 
                  Cytosol.padj = Cpres.down[match(rownames(Ipres.down), rownames(Cpres.down)), "padj"],
                  Nucleus.LFC = Npres[match(rownames(Ipres.down), rownames(Npres)), "log2FoldChange"], 
                  Nucleus.SE = Npres[match(rownames(Ipres.down), rownames(Npres)), "lfcSE"],
                  Nucleus.padj = Npres[match(rownames(Ipres.down), rownames(Npres)), "padj"],
                  ensID = geneMap[match(rownames(Ipres.down),as.character(rownames(geneMap))),"ensemblID"],
                  Symbol = geneMap[match(rownames(Ipres.down),as.character(rownames(geneMap))),"Symbol"],
                  EntrezID = geneMap[match(as.character(rownames(Ipres.down)),geneMap$gencodeID),"EntrezID"],
                  Type = geneMap[match(as.character(rownames(Ipres.down)),geneMap$gencodeID),"gene_type"])
age.sig = lapply(age.sig, function(x) age.data[which(age.data$geneID %in% rownames(x)),])
age.sig.1 = lapply(age.sig, function(x) x[which(abs(x$Cytosol.LFC)>=1 | abs(x$Nucleus.LFC)>=1),])
save(Ipres.down,Cpres.down,Npres,age.sig,age.sig.1,sig,sig.1,geneMap, 
     file = "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


## Annotation of age genes: Limiting to 1 LFC

freq = lapply(age.sig.1, function(x) count(x$Type))
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
group = c("Decreasing: Both", "Increasing: Both", "Decreasing:\nCytosol Only", "Decreasing:\nNucleus Only",
          "Increasing:\nCytosol Only", "Increasing:\nNucleus Only",
          "Decreasing: Nucleus/\nIncreasing: Cytosol", "Decreasing: Cytosol/\nIncreasing: Nucleus", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Decreasing: Both", "Increasing: Both", "Decreasing:\nCytosol Only","Decreasing:\nNucleus Only",
                                                    "Increasing:\nCytosol Only", "Increasing:\nNucleus Only",
                                                    "Decreasing: Nucleus/\nIncreasing: Cytosol", "Decreasing: Cytosol/\nIncreasing: Nucleus", "Interaction"))
# Graph the Frequencies
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_age-fraction_LFC1.pdf", height = 6, width = 8)
ggplot(type[which(type$Group!="Decreasing: Nucleus/\nIncreasing: Cytosol" &
                        type$Group!="Decreasing: Cytosol/\nIncreasing: Nucleus"),],
       aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

type = data.table(type[which(type$Group!="Decreasing: Nucleus/\nIncreasing: Cytosol" &
                               type$Group!="Decreasing: Cytosol/\nIncreasing: Nucleus"),])
x = data.frame(type[,sum(Count), by="Group"])
type$sum = x[match(type$Group, x$Group),"V1"]
type$perc = round(type$Count/type$sum*100,1)
type = ddply(type, .(Group), transform, pos = cumsum(perc) - (0.5 * perc))
type$pos = 100
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_age-fraction_LFC1.percent.pdf", height = 6, width = 8)
ggplot(type, aes(x = Group, y = perc, fill = RNA.Type)) + 
  geom_bar(stat = "identity") +
  geom_text(data=type[type$RNA.Type=="Long Non-coding",], 
            aes(x = Group, y = pos, label = sum), size=4, nudge_y = 5) +
  coord_flip() + labs(fill="") + ylab("Percent") + xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


i = type[grep("Increasing:",type$Group),]
d = type[grep("Decreasing:",type$Group),]
fisher.test(data.frame(c(sum(d[d$RNA.Type=="Protein-coding","Count"]), sum(d$Count)-sum(d[d$RNA.Type=="Protein-coding","Count"])),
                       c(sum(i[i$RNA.Type=="Protein-coding","Count"]), sum(i$Count)-sum(i[d$RNA.Type=="Protein-coding","Count"]))))
#p-value = 9.524e-13
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.329709 1.660557
#sample estimates:
#  odds ratio 
#1.485546
i = type[grep("Cytosol",type$Group),]
d = type[grep("Nucleus",type$Group),]
fisher.test(data.frame(c(sum(d[d$RNA.Type=="Protein-coding","Count"]), sum(d$Count)-sum(d[d$RNA.Type=="Protein-coding","Count"])),
                       c(sum(i[i$RNA.Type=="Protein-coding","Count"]), sum(i$Count)-sum(i[d$RNA.Type=="Protein-coding","Count"]))))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2822201 0.4251484
#sample estimates:
#  odds ratio 
#0.3465579


## Gene Ontology
names(age.sig.1) = names(age.sig) = group
entrezID = lapply(age.sig.1, function(x) na.omit(x$EntrezID))
# Define universe as all genes expressed in each of the four groups
GeneUniverse = as.character(unique(geneMap[match(rownames(Ipres.down),geneMap$gencodeID),"EntrezID"]))
GeneUniverse = na.omit(GeneUniverse)
# Find enriched Pathways via KEGG
elementNROWS(entrezID)
keggList = lapply(entrezID, function(x) enrichKEGG(as.character(x), organism="human", universe= GeneUniverse, 
                                                   minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1))
# Enriched Molecular Function GOs
goList_MF = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "MF", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
# Biological Process GO enrichment
goList_BP = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "BP", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
# Cellular Compartment GO enrichment
goList_CC = lapply(entrezID, function(x) enrichGO(as.character(x), ont = "CC", OrgDb = org.Hs.eg.db, 
                                                  universe= GeneUniverse, minGSSize=5, pAdjustMethod="BH",qvalueCutoff=1))
# Disease Ontology
goList_DO = lapply(entrezID, function(x) enrichDO(as.character(x), ont = "DO", universe= GeneUniverse, 
                                                  minGSSize=5, pAdjustMethod="BH", qvalueCutoff=1, readable=TRUE))

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrezID, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrezID, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# Save
save(keggList, goList_MF, goList_BP, goList_CC, goList_DO,compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/interaction.kegg.GO.DO.objects.polyAonly.sig1.downsampled.Agegenes.rda")


## plot
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/KEGG_interaction_polyAonly.sig1.downsampled.Agegenes.pdf", width=12,height=18)
plot(compareKegg,colorBy="p.adjust",  showCategory = 400, title= "KEGG Pathway Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/MF_interaction_polyAonly.sig1.downsampled.Agegenes.pdf", width=16,height=30)
plot(compareMF,colorBy="p.adjust",  showCategory = 400, title= "Molecular Function GO Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/BP_interaction_polyAonly.sig1.downsampled.Agegenes.pdf", width=19,height=130)
plot(compareBP,colorBy="p.adjust",  showCategory = 500, title= "Biological Process GO Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/CC_interaction_polyAonly.sig1.downsampled.Agegenes.pdf", width=14,height=40)
plot(compareCC,colorBy="p.adjust",  showCategory = 400, title= "Cellular Compartment GO Enrichment")
dev.off()
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/DO_interaction_polyAonly.sig1.downsampled.Agegenes.pdf", width=12,height=18)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")
dev.off()


## What about age genes?

### In the different groups of fraction regulated genes, is there a relationship between direction of expression over age and significance over age?
FracbyAge = list(Fpres.down = data.frame(Fpres.down), Apres = data.frame(Apres))
FracbyAge = Map(cbind, FracbyAge,lapply(FracbyAge, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Prenatal", "Adult"))
FracbyAge = do.call(rbind, FracbyAge)
FracbyAge$FDR = ifelse(FracbyAge$padj<=0.05, "FDR<0.05", "FDR>0.05")
FracbyAge = FracbyAge[which(FracbyAge$padj!="NA"),]

elementNROWS(age.sig)
fracDevel = lapply(age.sig[elementNROWS(age.sig)>0], function(x) FracbyAge[which(FracbyAge$gencodeID %in% x$geneID),])
names(fracDevel) = c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                     "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only",
                     "Decreasing: Nucleus/\nIncreasing: Cytoplasm", 
                     "Decreasing: Cytoplasm/\nIncreasing: Nucleus","Interaction")
fracDevel = do.call(rbind, Map(cbind, fracDevel, fracReg = as.list(names(fracDevel))))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                                      "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only",
                                      "Decreasing: Nucleus/\nIncreasing: Cytoplasm", 
                                      "Decreasing: Cytoplasm/\nIncreasing: Nucleus","Interaction"))


pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/AgebyFraction_LFCxFDR.pdf", width=22, height=6)
ggplot(fracDevel[fracDevel$fracReg!="Decreasing: Cytoplasm/\nIncreasing: Nucleus",], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ fracReg) +
  ylab("Log2 Fold Change") + 
  xlab("") +
  ggtitle(paste0("Fraction Expression Changes in Gene Groups Regulated by Age")) + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


### Is there a relationship between direction of expression by age and significance in groups of genes differentially regulated by fraction?

fracDevel = lapply(age.sig, function(x) FracbyAge[which(FracbyAge$gencodeID %in% x$geneID),])
elementNROWS(age.sig)-elementNROWS(lapply(fracDevel, function(x) x[x$Comparison=="Adult",]))
elementNROWS(age.sig)-elementNROWS(lapply(fracDevel, function(x) x[x$Comparison=="Prenatal",]))

# In adult:
tb = lapply(fracDevel, function(x) 
  data.frame(Sig = c(length(unique(x[which(x$Comparison=="Adult" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                     length(unique(x[which(x$Comparison=="Adult" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
             Nonsig = c(length(unique(x[which(x$Comparison=="Adult" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                        length(unique(x[which(x$Comparison=="Adult" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))),row.names = c("Retained","Exported")))
df = rbind(pvalue = unlist(lapply(lapply(tb, fisher.test), function(x) x$p.value)),OR = unlist(lapply(lapply(tb, fisher.test),function(x) x$estimate))) 
AdCounts = do.call(rbind, tb)
AdCounts$FracGroup = gsub("\\..*","", rownames(AdCounts))
AdCounts$fisher.pval = df["pvalue",match(AdCounts$FracGroup,colnames(df))]
AdCounts$fisher.OR = df["OR",match(AdCounts$FracGroup,colnames(df))]
AdCounts$Comparison = "Adult"
AdCounts$padj = p.adjust(AdCounts$fisher.pval, method="fdr")

# In Prenatal:
tb = lapply(fracDevel, function(x) 
  data.frame(Sig = c(length(unique(x[which(x$Comparison=="Prenatal" & x$log2FoldChange>0 & x$padj<=0.05),"gencodeID"])),
                     length(unique(x[which(x$Comparison=="Prenatal" & x$log2FoldChange<0 & x$padj<=0.05),"gencodeID"]))),
             Nonsig = c(length(unique(x[which(x$Comparison=="Prenatal" & x$log2FoldChange>0 & x$padj>0.05),"gencodeID"])),
                        length(unique(x[which(x$Comparison=="Prenatal" & x$log2FoldChange<0 & x$padj>0.05),"gencodeID"]))),row.names = c("Retained","Exported")))
df = rbind(pvalue = unlist(lapply(lapply(tb, fisher.test), function(x) x$p.value)),OR = unlist(lapply(lapply(tb, fisher.test),function(x) x$estimate))) 
PrenCounts = do.call(rbind, tb)
PrenCounts$FracGroup = gsub("\\..*","", rownames(PrenCounts))
PrenCounts$fisher.pval = df["pvalue",match(PrenCounts$FracGroup,colnames(df))]
PrenCounts$fisher.OR = df["OR",match(PrenCounts$FracGroup,colnames(df))]
PrenCounts$Comparison = "Prenatal"
PrenCounts$padj = p.adjust(PrenCounts$fisher.pval, method="fdr")
write.csv(rbind(AdCounts,PrenCounts),quote = F, 
          file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/FracLFCxFracpval_bysigAgeGroup_fisher.csv")

AdCounts[which(AdCounts$fisher.pval*9<=0.05),]
#                         Sig Nonsig       FracGroup  fisher.pval   fisher.OR Comparison         padj
#both_decreasing.Retained 344   1625 both_decreasing 1.014942e-04  0.72557420      Adult 1.826895e-04
#both_decreasing.Exported 405   1388 both_decreasing 1.014942e-04  0.72557420      Adult 1.826895e-04
#Cyt_decreasing.Retained  307    496  Cyt_decreasing 4.244042e-21 37.99319795      Adult 1.909819e-20
#Cyt_decreasing.Exported    2    123  Cyt_decreasing 4.244042e-21 37.99319795      Adult 1.909819e-20
#Nuc_decreasing.Retained    4     39  Nuc_decreasing 5.816124e-11  0.07149196      Adult 1.744837e-10
#Nuc_decreasing.Exported  367    255  Nuc_decreasing 5.816124e-11  0.07149196      Adult 1.744837e-10
#Cyt_increasing.Retained    4    107  Cyt_increasing 5.222355e-25  0.03846865      Adult 4.700120e-24
#Cyt_increasing.Exported  615    632  Cyt_increasing 5.222355e-25  0.03846865      Adult 4.700120e-24
#Nuc_increasing.Retained  399    538  Nuc_increasing 2.422647e-10 19.25096733      Adult 5.450957e-10
#Nuc_increasing.Exported    2     52  Nuc_increasing 2.422647e-10 19.25096733      Adult 5.450957e-10
PrenCounts[which(PrenCounts$fisher.pval*9<=0.05),]
#                         Sig Nonsig       FracGroup  fisher.pval fisher.OR Comparison         padj
#both_decreasing.Retained  43   1386 both_decreasing 9.975632e-05  2.590442   Prenatal 4.489034e-04
#both_decreasing.Exported  27   2255 both_decreasing 9.975632e-05  2.590442   Prenatal 4.489034e-04
#both_increasing.Retained  75   1103 both_increasing 4.527069e-07  3.306150   Prenatal 4.074362e-06
#both_increasing.Exported  20    973 both_increasing 4.527069e-07  3.306150   Prenatal 4.074362e-06
#Cyt_decreasing.Retained   23    422  Cyt_decreasing 4.565649e-04  6.311237   Prenatal 1.027271e-03
#Cyt_decreasing.Exported    3    348  Cyt_decreasing 4.565649e-04  6.311237   Prenatal 1.027271e-03
#Nuc_decreasing.Retained    8    151  Nuc_decreasing 3.577172e-03  5.833436   Prenatal 6.438910e-03
#Nuc_decreasing.Exported    4    442  Nuc_decreasing 3.577172e-03  5.833436   Prenatal 6.438910e-03
#Nuc_increasing.Retained   30    295  Nuc_increasing 3.338348e-04 13.586438   Prenatal 1.001504e-03
#Nuc_increasing.Exported    1    134  Nuc_increasing 3.338348e-04 13.586438   Prenatal 1.001504e-03


## make figure for age gene fraction-association
x = read.csv("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/FracLFCxFracpval_bysigAgeGroup_fisher.csv")
ad = x[which(x$padj<=0.05 & x$Comparison=="Adult"),]
fracDevel = lapply(age.sig[elementNROWS(age.sig)>0], function(x) FracbyAge[which(FracbyAge$gencodeID %in% x$geneID),])
names(fracDevel) = c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                     "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only",
                     "Decreasing: Nucleus/\nIncreasing: Cytoplasm", 
                     "Decreasing: Cytoplasm/\nIncreasing: Nucleus","Interaction")
fracDevel = do.call(rbind, Map(cbind, fracDevel, fracReg = as.list(names(fracDevel))))
fracDevel$fracReg = factor(fracDevel$fracReg, 
                           levels = c("Decreasing: Both","Increasing: Both","Decreasing:\nCytoplasm Only",
                                      "Decreasing:\nNucleus Only","Increasing:\nCytoplasm Only","Increasing:\nNucleus Only",
                                      "Decreasing: Nucleus/\nIncreasing: Cytoplasm", 
                                      "Decreasing: Cytoplasm/\nIncreasing: Nucleus","Interaction"))

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/AdultIncrDecr_Fraction-LFCxFDR.pdf", width=16, height=4)
ggplot(fracDevel[which(fracDevel$fracReg!="Decreasing: Cytoplasm/\nIncreasing: Nucleus" &
                         fracDevel$fracReg!="Decreasing: Nucleus/\nIncreasing: Cytoplasm" & 
                         fracDevel$fracReg!="Interaction" & fracDevel$Comparison=="Adult"),], 
       aes(x=fracReg, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() +
  ylab("Log2 Fold Change") + 
  xlab("") +
  ggtitle("Expression Changes by Fraction in Adult Samples in Gene Groups Regulated by Age") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


## make table of genes for supplement

load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
names(sig) = c("Retained: Both","Exported: Both","Retained: Prenatal Only","Retained: Adult Only",
               "Exported: Prenatal Only","Exported: Adult Only",
               "Retained in Adult; Exported in Prenatal", "Retained in Prenatal; Exported in Adult","Interaction")
names(age.sig) = c("Decreasing: Both","Increasing: Both","Decreasing: Cytosol Only","Decreasing: Nucleus Only",
                   "Increasing: Cytosol Only","Increasing: Nucleus Only","Decreasing in Nucleus; Increasing in Cytosol",
                   "Decreasing in Cytosol; Increasing in Nucleus", "Interaction")
for (i in 1:length(sig)) { if (nrow(sig[[i]])>0) { sig[[i]] = cbind(sig[[i]], Group = names(sig)[i]) } }
age.sig = Map(cbind, age.sig[elementNROWS(age.sig)>0], Group = names(age.sig)[elementNROWS(age.sig)>0])

x = read.csv("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_ASDenrichedGenes_retainedInBoth.csv")
frac = do.call(rbind, sig[elementNROWS(sig)>0])
frac[which(frac$gencodeID %in% x[x$variable=="ASD.CNV","gencodeID"]),"Gene Set"] = "ASD.CNV"
frac[which(frac$gencodeID %in% x[x$variable=="ASD.DATABASE","gencodeID"]),"Gene Set"] = "ASD.DATABASE"
frac[which(frac$gencodeID %in% x[x$variable=="SCZ.CNV","gencodeID"]),"Gene Set"] = "SCZ.CNV"
frac[which(frac$gencodeID %in% x[x$variable=="BPAD.GWAS","gencodeID"]),"Gene Set"] = "BPAD.GWAS"

age = do.call(rbind, age.sig[elementNROWS(age.sig)>0])
age[which(age$gencodeID %in% x[x$variable=="ASD.CNV","gencodeID"]),"Gene Set"] = "ASD.CNV"
age[which(age$gencodeID %in% x[x$variable=="ASD.DATABASE","gencodeID"]),"Gene Set"] = "ASD.DATABASE"
age[which(age$gencodeID %in% x[x$variable=="SCZ.CNV","gencodeID"]),"Gene Set"] = "SCZ.CNV"
age[which(age$gencodeID %in% x[x$variable=="BPAD.GWAS","gencodeID"]),"Gene Set"] = "BPAD.GWAS"

write.csv(frac, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/DEGs_byFrac_0.05.csv")
write.csv(age, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/DEGs_byAge_0.05.csv")
