library(org.Hs.eg.db)
library(biomaRt)
library(DESeq2)
library(parallel)
library(plyr)
library(scales)
library("ggplot2")
library(GenomicRanges)
library("clusterProfiler")
require("S4Vectors")
require("org.Hs.eg.db")
require("Rgraphviz")
library("DOSE")

load("/Users/amandaprice/Downloads/rawCounts_nucleusVsCytosol_n24.rda")

keepPeople = which(pd$totalMapped > 1e6) # removes low sequenced adult nuclear ribozero sample
keepGenes = which(rowMeans(geneCounts) > 0)
geneCounts = geneCounts[keepGenes,keepPeople]
exonCounts = exonCounts[,keepPeople]
pd = pd[keepPeople,]
geneMap = geneMap[keepGenes,]
JCOUNTS = data.frame(jCounts)
JCOUNTS = JCOUNTS[,keepPeople]
jCounts = jCounts[,keepPeople]
pd$Label = as.factor(paste(pd$Fetal, pd$Zone,pd$Library, sep="\n"))
pd$Label = factor(pd$Label, 
                  levels = c("Adult\nCytosol\npolyA", "Fetal\nCytosol\npolyA", "Adult\nNucleus\npolyA",
                             "Fetal\nNucleus\npolyA", "Adult\nCytosol\nRiboZero", "Fetal\nCytosol\nRiboZero",
                             "Adult\nNucleus\nRiboZero", "Fetal\nNucleus\nRiboZero"))
pd$WorkingID = c("Adult1_Cytosol_polyA", "Adult1_Nucleus_polyA", "Adult2_Cytosol_polyA", "Adult2_Nucleus_polyA",
                 "Adult3_Cytosol_polyA", "Adult3_Nucleus_polyA", "Fetal1_Cytosol_polyA", "Fetal1_Nucleus_polyA",
                 "Fetal2_Cytosol_polyA", "Fetal2_Nucleus_polyA", "Fetal3_Cytosol_polyA", "Fetal3_Nucleus_polyA",
                 "Adult1_Cytosol_RiboZero", "Adult2_Cytosol_RiboZero", "Adult2_Nucleus_RiboZero",
                 "Adult3_Cytosol_RiboZero", "Adult3_Nucleus_RiboZero", "Fetal1_Cytosol_RiboZero", "Fetal1_Nucleus_RiboZero",
                 "Fetal2_Cytosol_RiboZero", "Fetal2_Nucleus_RiboZero", "Fetal3_Cytosol_RiboZero", "Fetal3_Nucleus_RiboZero")

head(geneMap)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id", "hgnc_symbol",
                           "entrezgene", "description", "gene_biotype", "transcript_count"), 
            values=rownames(geneMap), mart=ensembl)
geneMap = data.frame(geneMap, EnsID = rownames(geneMap))
geneMap$Type = as.factor(sym$gene_biotype[match(geneMap$EnsID, sym$ensembl_gene_id)])
geneMap$TranscriptCount = sym$transcript_count[match(geneMap$EnsID, sym$ensembl_gene_id)]


# DESeq2
dds = DESeqDataSetFromMatrix(countData = geneCounts, colData = pd, design = ~ Library + Fetal + Zone + Fetal:Zone)
dds = DESeq(dds)
save(dds, file = "/Users/amandaprice/Dropbox/sorted_figures/new/dds_interaction.rda")
res = results(dds)
plotMA(res, alpha=0.05, main="Fraction:Age Interaction", ylim=c(-6,6))
dat <- res[order(rownames(res)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
res <- data.frame(dat, g)
res <- res[order(res$padj),]
write.csv(res, file="/Users/amandaprice/Dropbox/sorted_figures/new/interactionRes.csv")
res = res[which(res$padj<=0.05 & abs(res$log2FoldChange)>=1),]

pdf("/Users/amandaprice/Dropbox/sorted_figures/new/interaction_genes.pdf")
plots = list()
for (i in 1:nrow(res)){
  plots[[i]] = plotCounts(dds, as.character(res$EnsID[i]), 
                               intgroup = c("Fetal", "Zone", "Library"), returnData =TRUE)
  tmp = plots[[i]]
  tmp$Group = paste(tmp$Fetal,tmp$Zone, sep = "\n")
  x = ggplot(tmp, aes(x=Group, y=count)) + geom_boxplot() + 
    geom_jitter() +
    scale_y_log10(breaks=c(25,100,400)) +
    ylab("Normalized Count") + 
    xlab("") +
    ggtitle(as.character(res$Symbol[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(x)
}
dev.off()

Adult <- pd[which(pd$Fetal=="Adult"),]
Fetal <- pd[which(pd$Fetal=="Fetal"),]
Adult.counts <- geneCounts[,which(colnames(geneCounts)%in%Adult$SampleID)]
Fetal.counts <- geneCounts[,which(colnames(geneCounts)%in%Fetal$SampleID)]

Adds = DESeqDataSetFromMatrix(countData = Adult.counts, colData = Adult, design = ~ Library + Zone)
Adds = DESeq(Adds)
Ares = results(Adds)
plotMA(Ares, alpha=0.05, main="Fraction Differences in Adults", ylim=c(-6,6))
dat <- Ares[order(rownames(Ares)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Ares <- data.frame(dat, g)
Ares <- Ares[order(Ares$padj),]
write.csv(Ares, file="/Users/amandaprice/Dropbox/sorted_figures/new/Ares.csv")
Ares.1.05 = Ares[which(Ares$padj<=0.05 & abs(Ares$log2FoldChange)>=1),]
Ares.05 = Ares[which(Ares$padj<=0.05),]

Fdds = DESeqDataSetFromMatrix(countData = Fetal.counts, colData = Fetal, design = ~ Library + Zone)
Fdds = DESeq(Fdds)
Fres = results(Fdds)
plotMA(Fres, alpha=0.05, main="Fraction Differences in Fetal", ylim=c(-6,6))
dat <- Fres[order(rownames(Fres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Fres <- data.frame(dat, g)
Fres <- Fres[order(Fres$padj),]
write.csv(Fres, file="/Users/amandaprice/Dropbox/sorted_figures/new/Fres.csv")
Fres.1.05 = Fres[which(Fres$padj<=0.05 & abs(Fres$log2FoldChange)>=1),]
Fres.05 = Fres[which(Fres$padj<=0.05),]

Cytosol <- pd[which(pd$Zone=="Cytosol"),]
Nucleus <- pd[which(pd$Zone=="Nucleus"),]
Cytosol.counts <- geneCounts[,which(colnames(geneCounts)%in%Cytosol$SampleID)]
Nucleus.counts <- geneCounts[,which(colnames(geneCounts)%in%Nucleus$SampleID)]

Cdds = DESeqDataSetFromMatrix(countData = Cytosol.counts, colData = Cytosol,
                              design = ~ Library + Fetal)
Cdds = DESeq(Cdds)
Cres = results(Cdds)
plotMA(Cres, alpha=0.05, main="Age Differences in Cytosol", ylim=c(-6,6))
dat <- Cres[order(rownames(Cres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Cres <- data.frame(dat, g)
Cres <- Cres[order(Cres$padj),]
write.csv(Cres, file="/Users/amandaprice/Dropbox/sorted_figures/new/Cres.csv")
Cres.1.05 = Cres[which(Cres$padj<=0.05 & abs(Cres$log2FoldChange)>=1),]
Cres.05 = Cres[which(Cres$padj<=0.05),]

Ndds = DESeqDataSetFromMatrix(countData = Nucleus.counts, colData = Nucleus,
                              design = ~ Library + Fetal)
Ndds = DESeq(Ndds)
Nres = results(Ndds)
plotMA(Nres, alpha=0.05, main="Age Differences in Nucleus", ylim=c(-6,6))
dat <- Nres[order(rownames(Nres)),]
g <-geneMap[which(rownames(geneMap)%in%rownames(dat)),]
g <- g[order(rownames(g)),]
Nres <- data.frame(dat, g)
Nres <- Nres[order(Nres$padj),]
write.csv(Nres, file="/Users/amandaprice/Dropbox/sorted_figures/new/Nres.csv")
Nres.1.05 = Nres[which(Nres$padj<=0.05 & abs(Nres$log2FoldChange)>=1),]
Nres.05 = Nres[which(Nres$padj<=0.05),]

# Make list of genes
all = list(res = res[order(res$EnsID),], Ares = Ares[order(Ares$EnsID),], 
           Fres = Fres[order(Fres$EnsID),])
Sign = lapply(all, function(x) ifelse(x$log2FoldChange>0, "Pos", "Neg"))
all = Map(cbind, all, Sign = Sign)
res = all[["res"]]
Ares = all[["Ares"]]
Fres = all[["Fres"]]

AdPos = Ares[which(Ares$Sign=="Pos"),]
AdNeg = Ares[which(Ares$Sign=="Neg"),]
AdSig = Ares[which(Ares$padj<=0.05),]
AdLFC = Ares[which(abs(Ares$log2FoldChange)>=1),]
FetPos = Fres[which(Fres$Sign=="Pos"),]
FetNeg = Fres[which(Fres$Sign=="Neg"),]
FetSig = Fres[which(Fres$padj<=0.05),]
FetLFC = Fres[which(abs(Fres$log2FoldChange)>=1),]
AdPos = as.character(AdPos$EnsID)
AdNeg = as.character(AdNeg$EnsID)
AdSig = as.character(AdSig$EnsID)
AdLFC = as.character(AdLFC$EnsID)
FetPos = as.character(FetPos$EnsID)
FetNeg = as.character(FetNeg$EnsID)
FetSig = as.character(FetSig$EnsID)
FetLFC = as.character(FetLFC$EnsID)

sig = list(both_retained = res[which(res$EnsID%in%AdPos & res$EnsID%in%FetPos 
                                     & res$EnsID%in%AdSig & res$EnsID%in%FetSig),],
           both_exported = res[which(res$EnsID%in%AdNeg & res$EnsID%in%FetNeg 
                                     & res$EnsID%in%AdSig & res$EnsID%in%FetSig),],
           Fet_retained = res[which(res$EnsID%in%FetPos & !(res$EnsID%in%AdSig) & res$EnsID%in%FetSig),],
           Ad_retained = res[which(res$EnsID%in%AdPos & res$EnsID%in%AdSig & !(res$EnsID%in%FetSig)),],
           ret_Ad_exp_Fet = res[which(res$EnsID%in%AdPos & res$EnsID%in%FetNeg 
                                      & res$EnsID%in%AdSig & res$EnsID%in%FetSig),],
           ret_Fet_exp_Ad = res[which(res$EnsID%in%FetPos & res$EnsID%in%AdNeg 
                                      & res$EnsID%in%AdSig & res$EnsID%in%FetSig),],
           interacting = res[which(res$padj<=0.05),])

data = data.frame(EnsID = res$EnsID, baseMean = res$baseMean, Fetal.LFC = Fres$log2FoldChange, 
                  Fetal.SE = Fres$lfcSE, Fetal.padj = Fres$padj,Adult.LFC = Ares$log2FoldChange, 
                  Adult.SE = Ares$lfcSE, Adult.padj = Ares$padj,Symbol = res$Symbol, EntrezID = res$EntrezID,
                  Type = res$Type)
sig = lapply(sig, function(x) data[which(data$EnsID%in%x$EnsID),])
sig.1 = lapply(sig, function(x) x[which(abs(x$Fetal.LFC)>=1 | abs(x$Adult.LFC)>=1),])
save(res,Fres,Ares,sig,sig.1, file = "./Dropbox/sorted_figures/new/retained.byAge.rda")


load("./Dropbox/sorted_figures/new/retained.byAge.rda")
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
for (i in 1:7){
  f = freq[[i]]
  f$x = as.character(f$x)
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="protein_coding"),2], 
        f[which(f$x=="polymorphic_pseudogene"),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="pseudogene"),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="3prime_overlapping_ncrna"),2],
        f[which(f$x=="antisense"),2], 
        f[which(f$x=="lincRNA"),2],
        f[which(f$x=="sense_intronic"),2], 
        f[which(f$x=="misc_RNA"),2],
        f[which(f$x=="processed_transcript"),2], 
        f[which(f$x=="sense_overlapping"),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="miRNA"),2], 
        f[which(f$x=="rRNA"),2],
        f[which(f$x=="snoRNA"),2], 
        f[which(f$x=="snRNA"),2])
}
group = c("Retained: Both", "Exported: Both", "Retained:\nFetal Only", "Retained:\nAdult Only",
               "Retained: Adult/\nExported: Fetal", "Retained: Fetal/\nExported: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nFetal Only",
  "Retained: Adult/\nExported: Fetal", "Retained: Fetal/\nExported: Adult", "Interaction"))

# Graph the Frequencies
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Limiting to 1 LFC
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
for (i in 1:7){
  f = freq[[i]]
  f$x = as.character(f$x)
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="protein_coding"),2], 
        f[which(f$x=="polymorphic_pseudogene"),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="pseudogene"),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="3prime_overlapping_ncrna"),2],
        f[which(f$x=="antisense"),2], 
        f[which(f$x=="lincRNA"),2],
        f[which(f$x=="sense_intronic"),2], 
        f[which(f$x=="misc_RNA"),2],
        f[which(f$x=="processed_transcript"),2], 
        f[which(f$x=="sense_overlapping"),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x=="miRNA"),2], 
        f[which(f$x=="rRNA"),2],
        f[which(f$x=="snoRNA"),2], 
        f[which(f$x=="snRNA"),2])
}
group = c("Retained: Both", "Exported: Both", "Retained:\nFetal Only", "Retained:\nAdult Only",
          "Retained: Adult/\nExported: Fetal", "Retained: Fetal/\nExported: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Retained: Both", "Exported: Both", "Retained:\nAdult Only", "Retained:\nFetal Only",
                                                    "Retained: Adult/\nExported: Fetal", "Retained: Fetal/\nExported: Adult", "Interaction"))
# Graph the Frequencies
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nLog2 Fold Change >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Gene Ontology
names(sig.1) = names(sig) = c("Retained: Both", "Exported: Both", "Retained:\nFetal Only", "Retained:\nAdult Only",
                 "Retained: Adult/\nExported: Fetal", "Retained: Fetal/\nExported: Adult", "Interaction")
entrezID = lapply(sig.1, function(x) na.omit(x$EntrezID))
# Define universe as all genes expressed in each of the four groups
GeneUniverse = as.character(unique(res$EntrezID[!is.na(res$EntrezID)]))

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
names = c("Retained_Both", "Exported_Both", "Retained_Fetal.Only", "Retained_Adult.Only",
                 "Retained_Adult.Exported_Fetal", "Retained_Fetal.Exported_Adult", "Interaction")
n = c(1,2,4,7)
for (i in n){write.csv(keggListdf[[i]], file=paste0("/Users/amandaprice/Dropbox/sorted_figures/new/GO_terms/",names[i],".GO.KEGG.csv"))}
for (i in n){write.csv(goListdf_BP[[i]], file=paste0("/Users/amandaprice/Dropbox/sorted_figures/new/GO_terms/",names[i],".GO.BP.csv"))}
for (i in n){write.csv(goListdf_MF[[i]], file=paste0("/Users/amandaprice/Dropbox/sorted_figures/new/GO_terms/",names[i],".GO.MF.csv"))}
n = c(2,4,7)
for (i in n){write.csv(goListdf_CC[[i]], file=paste0("/Users/amandaprice/Dropbox/sorted_figures/new/GO_terms/",names[i],".GO.CC.csv"))}
write.csv(goListdf_DO[[7]], file="/Users/amandaprice/Dropbox/sorted_figures/new/GO_terms/interaction_DO.csv")

# Compare the enriched terms between 7 groups
# KEGG
compareKegg = compareCluster(entrezID, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
summary(compareKegg)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
# Biological Process
compareBP = compareCluster(entrezID, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
as.data.frame(compareBP)
plot(compareBPDropped,colorBy="p.adjust",  showCategory = 20, title= "Biological Process GO Enrichment")
compareBPDropped = dropGO(compareBP,  level = 1)
compareBPDropped = dropGO(compareBPDropped,  level = 2)
compareBPDropped = dropGO(compareBPDropped,  level = 3)

# Molecular Function
compareMF = compareCluster(entrezID, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
as.data.frame(compareMF)
plot(compareMFDropped,colorBy="p.adjust",  showCategory = 20, title= "Molecular Function GO Enrichment")
compareMFDropped = dropGO(compareMF,  level = 1)
compareMFDropped = dropGO(compareMFDropped,  level = 2)
compareMFDropped = dropGO(compareMFDropped,  level = 3)

# Cellular Component
compareCC = compareCluster(entrezID, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
as.data.frame(compareCC)
plot(compareCCDropped,colorBy="p.adjust",  showCategory = 20, title= "Cellular Compartment GO Enrichment")
compareCCDropped = dropGO(compareCC,  level = 1)
compareCCDropped = dropGO(compareCCDropped,  level = 2)
compareCCDropped = dropGO(compareCCDropped,  level = 3)

# Disease Ontology
compareDO = compareCluster(entrezID, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
as.data.frame(compareDO)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")

entrez.noLFC = lapply(sig, function(x) na.omit(x$EntrezID))
compareDO = compareCluster(entrez.noLFC, fun="enrichDO",
                           ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(compareDO,colorBy="p.adjust",  showCategory = 30, title= "Disease Ontology Enrichment")


save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/kegg.GO.DO.objects.rda")

