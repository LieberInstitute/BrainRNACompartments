library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Make list of length of significant genes in a list
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres), 
                Arres = data.frame(Arres), Frres = data.frame(Frres),
                Fpres_down = data.frame(Fpres.down))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracBySign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(sigFracBySign, function(x) split(x, x$Sign))
SigList = unlist(DirList, recursive=F)
elementNROWS(SigList)
maps = lapply(SigList, function(x) geneMap[match(rownames(x),geneMap$gencodeID),])
elementNROWS(maps)

type = lapply(maps, function(x) as.data.frame(x$gene_type))
freq = lapply(type, function(x) count(x))
TypeFreq = do.call(rbind, freq)
TypeFreq[grep("Apres.DownNuc",rownames(TypeFreq)),"Group"] = "Adult\nPolyA\nCytosolic"
TypeFreq[grep("Apres.UpNuc",rownames(TypeFreq)),"Group"] = "Adult\nPolyA\nNuclear"
TypeFreq[grep("Fpres.DownNuc",rownames(TypeFreq)),"Group"] = "Prenatal\nPolyA\nCytosolic"
TypeFreq[grep("Fpres.UpNuc",rownames(TypeFreq)),"Group"] = "Prenatal\nPolyA\nNuclear"
TypeFreq[grep("Arres.DownNuc",rownames(TypeFreq)),"Group"] = "Adult\nRiboZero\nCytosolic"
TypeFreq[grep("Arres.UpNuc",rownames(TypeFreq)),"Group"] = "Adult\nRiboZero\nNuclear"
TypeFreq[grep("Frres.DownNuc",rownames(TypeFreq)),"Group"] = "Prenatal\nRiboZero\nCytosolic"
TypeFreq[grep("Frres.UpNuc",rownames(TypeFreq)),"Group"] = "Prenatal\nRiboZero\nNuclear"
TypeFreq[grep("Fpres_down.DownNuc",rownames(TypeFreq)),"Group"] = "Prenatal\nPolyA\nCytosolic"
TypeFreq[grep("Fpres_down.UpNuc",rownames(TypeFreq)),"Group"] = "Prenatal\nPolyA\nNuclear"
colnames(TypeFreq) = c("RNA_Type", "Count", "Group")
TypeFreq$Group = factor(TypeFreq$Group, levels = c("Prenatal\nRiboZero\nNuclear","Prenatal\nPolyA\nNuclear",
                           "Adult\nRiboZero\nNuclear","Adult\nPolyA\nNuclear",
                           "Prenatal\nRiboZero\nCytosolic","Prenatal\nPolyA\nCytosolic",
                           "Adult\nRiboZero\nCytosolic","Adult\nPolyA\nCytosolic"))

# Graph the Frequencies
ggplot(TypeFreq, aes(x = Group, y = Count, fill = RNA_Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation (abs(LFC) >1; FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Graph as 100% Stacked
ggplot(TypeFreq, aes(x = Group, y = Count, fill = RNA_Type)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format()) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation (abs(LFC) >1; FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Condense to 4 categories
RNA.Type = c("Protein-coding", "Pseudogene", "Long Non-coding", "Short Non-coding")
type = data.frame(RNA.Type = rep.int(c(rep.int(RNA.Type[1], 2), rep.int(RNA.Type[2], 2),rep.int(RNA.Type[3], 2),rep.int(RNA.Type[4], 2)),4),
                  Count = NA, 
                  Group = factor(x=c(rep.int(c("Prenatal\nRiboZero\nNuclear","Prenatal\nRiboZero\nCytosolic"), 4),
                                     rep.int(c("Prenatal\nPolyA\nNuclear","Prenatal\nPolyA\nCytosolic"), 4),
                                     rep.int(c("Adult\nRiboZero\nNuclear","Adult\nRiboZero\nCytosolic"), 4),
                                     rep.int(c("Adult\nPolyA\nNuclear","Adult\nPolyA\nCytosolic"), 4))))
type.down = data.frame(RNA.Type = rep.int(c(rep.int(RNA.Type[1], 2), rep.int(RNA.Type[2], 2),rep.int(RNA.Type[3], 2),rep.int(RNA.Type[4], 2)),4),
                  Count = NA, 
                  Group = factor(x=c(rep.int(c("Prenatal\nRiboZero\nNuclear","Prenatal\nRiboZero\nCytosolic"), 4),
                                     rep.int(c("Prenatal\nPolyA\nNuclear","Prenatal\nPolyA\nCytosolic"), 4),
                                     rep.int(c("Adult\nRiboZero\nNuclear","Adult\nRiboZero\nCytosolic"), 4),
                                     rep.int(c("Adult\nPolyA\nNuclear","Adult\nPolyA\nCytosolic"), 4))))
TypeFreq.down = TypeFreq[c(1:27,34:81),]
TypeFreq = TypeFreq[1:77,]
group = unique(type$Group)
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
ifelse(sum(type$Count)==sum(TypeFreq$Count), "YES","NO")

for (i in 1:length(group)){
  type.down[which(type.down$RNA.Type=="Protein-coding" & type.down$Group==group[i]),2] = 
    sum(TypeFreq.down[which(TypeFreq.down$RNA_Type=="protein_coding" & TypeFreq.down$Group==group[i]),2], 
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="TR_C_gene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="polymorphic_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="IG_V_gene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="TR_V_gene" & TypeFreq.down$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Pseudogene" & type.down$Group==group[i]),2] = 
    sum(TypeFreq.down[which(TypeFreq.down$RNA_Type=="pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="transcribed_processed_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="unprocessed_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="transcribed_unprocessed_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="processed_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="IG_C_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="transcribed_unitary_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="unitary_pseudogene" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="IG_V_pseudogene" & TypeFreq.down$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Long Non-coding" & type.down$Group==group[i]),2] = 
    sum(TypeFreq.down[which(TypeFreq.down$RNA_Type=="3prime_overlapping_ncRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="antisense" & TypeFreq.down$Group==group[i]),2], 
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="non_coding" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="lincRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="sense_intronic" & TypeFreq.down$Group==group[i]),2], 
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="misc_RNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="processed_transcript" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="bidirectional_promoter_lncRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="TEC" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="macro_lncRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="sense_overlapping" & TypeFreq.down$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Short Non-coding" & type.down$Group==group[i]),2] = 
    sum(TypeFreq.down[which(TypeFreq.down$RNA_Type=="miRNA" & TypeFreq.down$Group==group[i]),2], 
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="Mt_rRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="Mt_tRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="rRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="snoRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="scRNA" & TypeFreq.down$Group==group[i]),2],
        TypeFreq.down[which(TypeFreq.down$RNA_Type=="snRNA" & TypeFreq.down$Group==group[i]),2])
}
ifelse(sum(type.down$Count)==sum(TypeFreq.down$Count), "YES","NO")
ifelse(type$Group==type.down$Group, "YES","NO")

type$Group = type.down$Group = factor(type$Group, 
                                      levels = c("Prenatal\nRiboZero\nNuclear","Prenatal\nPolyA\nNuclear",
                                                 "Adult\nRiboZero\nNuclear","Adult\nPolyA\nNuclear",
                                                 "Prenatal\nRiboZero\nCytosolic","Prenatal\nPolyA\nCytosolic",
                                                 "Adult\nRiboZero\nCytosolic","Adult\nPolyA\nCytosolic"))

# Graph the Frequencies
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

ggplot(type.down, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))