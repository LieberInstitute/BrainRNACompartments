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
Irdds = DESeqDataSetFromMatrix(countData = geneCounts[,grep("RiboZero", colnames(geneCounts))], 
                                    colData = pd[which(pd$Library=="RiboZero"),], design = ~ Fetal + Zone + Fetal:Zone)
Irdds = DESeq(Irdds)
Irres = results(Irdds)
sigres = data.frame(Irres[which(Irres$padj<=0.05 & abs(Irres$log2FoldChange)>=1),])


# Make list of genes
all = list(Irres = data.frame(Irres[order(rownames(Irres)),]), 
           Arres = data.frame(Arres[order(rownames(Arres)),]), 
           Frres = data.frame(Frres[order(rownames(Frres)),]))
all = Map(cbind, all, Sign = lapply(all, function(x) ifelse(x$log2FoldChange>0, "Pos", "Neg")))
Irres = all[["Irres"]]
Arres = all[["Arres"]]
Frres = all[["Frres"]]

AdPos = as.character(rownames(Arres[which(Arres$Sign=="Pos"),]))
AdNeg = as.character(rownames(Arres[which(Arres$Sign=="Neg"),]))
AdSig = as.character(rownames(Arres[which(Arres$padj<=0.05),]))
AdLFC = as.character(rownames(Arres[which(abs(Arres$log2FoldChange)>=1),]))
FetPos = as.character(rownames(Frres[which(Frres$Sign=="Pos"),]))
FetNeg = as.character(rownames(Frres[which(Frres$Sign=="Neg"),]))
FetSig = as.character(rownames(Frres[which(Frres$padj<=0.05),]))
FetLFC = as.character(rownames(Frres[which(abs(Frres$log2FoldChange)>=1),]))

sig = list(both_retained = Irres[which(rownames(Irres)%in%AdPos & rownames(Irres)%in%FetPos 
                                            & rownames(Irres)%in%AdSig & rownames(Irres)%in%FetSig),],
           both_exported = Irres[which(rownames(Irres)%in%AdNeg & rownames(Irres)%in%FetNeg 
                                            & rownames(Irres)%in%AdSig & rownames(Irres)%in%FetSig),],
           Fet_retained = Irres[which(rownames(Irres)%in%FetPos & !(rownames(Irres)%in%AdSig) & rownames(Irres)%in%FetSig),],
           Ad_retained = Irres[which(rownames(Irres)%in%AdPos & rownames(Irres)%in%AdSig & !(rownames(Irres)%in%FetSig)),],
           Fet_exported = Irres[which(rownames(Irres)%in%FetNeg & !(rownames(Irres)%in%AdSig) & rownames(Irres)%in%FetSig),],
           Ad_exported = Irres[which(rownames(Irres)%in%AdNeg & rownames(Irres)%in%AdSig & !(rownames(Irres)%in%FetSig)),],
           ret_Ad_exp_Fet = Irres[which(rownames(Irres)%in%AdPos & rownames(Irres)%in%FetNeg 
                                             & rownames(Irres)%in%AdSig & rownames(Irres)%in%FetSig),],
           ret_Fet_exp_Ad = Irres[which(rownames(Irres)%in%FetPos & rownames(Irres)%in%AdNeg 
                                             & rownames(Irres)%in%AdSig & rownames(Irres)%in%FetSig),],
           interacting = Irres[which(Irres$padj<=0.05),])

data = data.frame(geneID = rownames(Irres), baseMean = Irres$baseMean, 
                  Prenatal.LFC = Frres[match(rownames(Irres), rownames(Frres)), "log2FoldChange"], 
                  Prenatal.SE = Frres[match(rownames(Irres), rownames(Frres)), "lfcSE"], 
                  Prenatal.padj = Frres[match(rownames(Irres), rownames(Frres)), "padj"],
                  Adult.LFC = Arres[match(rownames(Irres), rownames(Arres)), "log2FoldChange"], 
                  Adult.SE = Arres[match(rownames(Irres), rownames(Arres)), "lfcSE"],
                  Adult.padj = Arres[match(rownames(Irres), rownames(Arres)), "padj"],
                  ensID = geneMap[match(rownames(Irres),as.character(rownames(geneMap))),"ensemblID"],
                  Symbol = geneMap[match(rownames(Irres),as.character(rownames(geneMap))),"Symbol"],
                  EntrezID = geneMap[match(as.character(rownames(Irres)),geneMap$gencodeID),"EntrezID"],
                  Type = geneMap[match(as.character(rownames(Irres)),geneMap$gencodeID),"gene_type"])
sig = lapply(sig, function(x) data[which(data$geneID %in% rownames(x)),])
sig.1 = lapply(sig, function(x) x[which(abs(x$Prenatal.LFC)>=1 | abs(x$Adult.LFC)>=1),])

# Make Age Sig object
age = list(Irres = data.frame(Irres[order(rownames(Irres)),]), 
           Nrres = data.frame(Nrres[order(rownames(Nrres)),]), 
           Crres = data.frame(Crres[order(rownames(Crres)),]))
age = Map(cbind, age, Sign = lapply(age, function(x) ifelse(x$log2FoldChange>0, "Pos", "Neg")))
Irres = age[["Irres"]]
Nrres = age[["Nrres"]]
Crres = age[["Crres"]]

NucPos = as.character(rownames(Nrres[which(Nrres$Sign=="Pos"),]))
NucNeg = as.character(rownames(Nrres[which(Nrres$Sign=="Neg"),]))
NucSig = as.character(rownames(Nrres[which(Nrres$padj<=0.05),]))
NucLFC = as.character(rownames(Nrres[which(abs(Nrres$log2FoldChange)>=1),]))
CytPos = as.character(rownames(Crres[which(Crres$Sign=="Pos"),]))
CytNeg = as.character(rownames(Crres[which(Crres$Sign=="Neg"),]))
CytSig = as.character(rownames(Crres[which(Crres$padj<=0.05),]))
CytLFC = as.character(rownames(Crres[which(abs(Crres$log2FoldChange)>=1),]))

age.sig = list(both_decreasing = Irres[which(rownames(Irres)%in%NucPos & rownames(Irres)%in%CytPos 
                                                  & rownames(Irres)%in%NucSig & rownames(Irres)%in%CytSig),],
               both_increasing = Irres[which(rownames(Irres)%in%NucNeg & rownames(Irres)%in%CytNeg 
                                                  & rownames(Irres)%in%NucSig & rownames(Irres)%in%CytSig),],
               Cyt_decreasing = Irres[which(rownames(Irres)%in%CytPos & !(rownames(Irres)%in%NucSig) & rownames(Irres)%in%CytSig),],
               Nuc_decreasing = Irres[which(rownames(Irres)%in%NucPos & rownames(Irres)%in%NucSig & !(rownames(Irres)%in%CytSig)),],
               Cyt_increasing = Irres[which(rownames(Irres)%in%CytNeg & !(rownames(Irres)%in%NucSig) & rownames(Irres)%in%CytSig),],
               Nuc_increasing = Irres[which(rownames(Irres)%in%NucNeg & rownames(Irres)%in%NucSig & !(rownames(Irres)%in%CytSig)),],
               decr_Nuc_incr_Cyt = Irres[which(rownames(Irres)%in%NucPos & rownames(Irres)%in%CytNeg 
                                                    & rownames(Irres)%in%NucSig & rownames(Irres)%in%CytSig),],
               decr_Cyt_incr_Nuc = Irres[which(rownames(Irres)%in%CytPos & rownames(Irres)%in%NucNeg 
                                                    & rownames(Irres)%in%NucSig & rownames(Irres)%in%CytSig),],
               interacting = Irres[which(Irres$padj<=0.05),])

age.data = data.frame(geneID = rownames(Irres), baseMean = Irres$baseMean, 
                      Cytoplasm.LFC = Crres[match(rownames(Irres), rownames(Crres)), "log2FoldChange"], 
                      Cytoplasm.SE = Crres[match(rownames(Irres), rownames(Crres)), "lfcSE"], 
                      Cytoplasm.padj = Crres[match(rownames(Irres), rownames(Crres)), "padj"],
                      Nucleus.LFC = Nrres[match(rownames(Irres), rownames(Nrres)), "log2FoldChange"], 
                      Nucleus.SE = Nrres[match(rownames(Irres), rownames(Nrres)), "lfcSE"],
                      Nucleus.padj = Nrres[match(rownames(Irres), rownames(Nrres)), "padj"],
                      ensID = geneMap[match(rownames(Irres),as.character(rownames(geneMap))),"ensemblID"],
                      Symbol = geneMap[match(rownames(Irres),as.character(rownames(geneMap))),"Symbol"],
                      EntrezID = geneMap[match(as.character(rownames(Irres)),geneMap$gencodeID),"EntrezID"],
                      Type = geneMap[match(as.character(rownames(Irres)),geneMap$gencodeID),"gene_type"])
age.sig = lapply(age.sig, function(x) age.data[which(age.data$geneID %in% rownames(x)),])
age.sig.1 = lapply(age.sig, function(x) x[which(abs(x$Cytoplasm.LFC)>=1 | abs(x$Nucleus.LFC)>=1),])
save(Irres,Crres,Nrres,age.sig,age.sig.1,sig,sig.1,geneMap, 
     file = "./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.Ribozero.rda")


### Annotate by fraction
## Limiting to 1 LFC

freq = lapply(sig.1, function(x) data.frame(table(x$Type)))
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
group = c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nPrenatal Only", "Nuclear:\nAdult Only",
          "Cytoplasmic:\nPrenatal Only", "Cytoplasmic:\nAdult Only",
          "Nuclear: Adult/\nCytoplasmic: Prenatal", "Nuclear: Prenatal/\nCytoplasmic: Adult", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Nuclear: Both", "Cytoplasmic: Both", "Nuclear:\nAdult Only", "Nuclear:\nPrenatal Only",
                                                    "Cytoplasmic:\nAdult Only","Cytoplasmic:\nPrenatal Only", 
                                                    "Nuclear: Adult/\nCytoplasmic: Prenatal", "Nuclear: Prenatal/\nCytoplasmic: Adult", "Interaction"))

# Graph the Frequencies
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_fraction-age_LFC1.ribozero.pdf", height = 6, width = 8)
ggplot(type[which(type$Group!="Nuclear: Adult/\nCytoplasmic: Prenatal" &
                    type$Group!="Nuclear: Prenatal/\nCytoplasmic: Adult"),], 
       aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


type = data.table(type[which(type$Group!="Nuclear: Adult/\nCytoplasmic: Prenatal" &
                               type$Group!="Nuclear: Prenatal/\nCytoplasmic: Adult"),])
x = data.frame(type[,sum(Count), by="Group"])
type$sum = x[match(type$Group, x$Group),"V1"]
type$perc = round(type$Count/type$sum*100,1)
type = ddply(type, .(Group), transform, pos = cumsum(perc) - (0.5 * perc))
type$pos = 100

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_fraction-age_LFC1.percent.ribozero.pdf", height = 6, width = 8)
ggplot(type, aes(x = Group, y = perc, fill = RNA.Type)) + 
  geom_bar(stat = "identity") +
  geom_text(data=type[type$RNA.Type=="Long Non-coding",], 
            aes(x = Group, y = pos, label = sum), size=4, nudge_y = 5) +
  coord_flip() + labs(fill="") + ylab("Percent") + xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Annotation of age genes: Limiting to 1 LFC

freq = lapply(age.sig.1, function(x) data.frame(table(x$Type)))
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
group = c("Decreasing: Both", "Increasing: Both", "Decreasing:\nCytoplasm Only", "Decreasing:\nNucleus Only",
          "Increasing:\nCytoplasm Only", "Increasing:\nNucleus Only",
          "Decreasing: Nucleus/\nIncreasing: Cytoplasm", "Decreasing: Cytoplasm/\nIncreasing: Nucleus", "Interaction")
type$Group = factor(x=rep.int(group, 4), levels = c("Decreasing: Both", "Increasing: Both", "Decreasing:\nNucleus Only", "Decreasing:\nCytoplasm Only",
                                                    "Increasing:\nNucleus Only","Increasing:\nCytoplasm Only", 
                                                    "Decreasing: Nucleus/\nIncreasing: Cytoplasm", "Decreasing: Cytoplasm/\nIncreasing: Nucleus", "Interaction"))

# Graph the Frequencies
pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_age-fraction_LFC1.ribozero.pdf", height = 6, width = 8)
ggplot(type[which(type$Group!="Decreasing: Nucleus/\nIncreasing: Cytoplasm" &
                    type$Group!="Decreasing: Cytoplasm/\nIncreasing: Nucleus"),],
       aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

type = data.table(type[which(type$Group!="Decreasing: Nucleus/\nIncreasing: Cytoplasm" &
                               type$Group!="Decreasing: Cytoplasm/\nIncreasing: Nucleus"),])
x = data.frame(type[,sum(Count), by="Group"])
type$sum = x[match(type$Group, x$Group),"V1"]
type$perc = round(type$Count/type$sum*100,1)
type = ddply(type, .(Group), transform, pos = cumsum(perc) - (0.5 * perc))
type$pos = 100

pdf("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/annotation_DEG_interaction_age-fraction_LFC1.percent.ribozero.pdf", height = 6, width = 8)
ggplot(type, aes(x = Group, y = perc, fill = RNA.Type)) + 
  geom_bar(stat = "identity") +
  geom_text(data=type[type$RNA.Type=="Long Non-coding",], 
            aes(x = Group, y = pos, label = sum), size=4, nudge_y = 5) +
  coord_flip() + labs(fill="") + ylab("Percent") + xlab("") +
  ggtitle("Gene Annotation:\nabs(Log2 Fold Change) >1") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()





