library(org.Hs.eg.db)
library(biomaRt)
library(parallel)
library('derfinder')
library(plyr)
library(scales)
library("ggplot2")

# Make list of length of significant genes in a list
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracBySign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(sigFracBySign, function(x) split(x, x$Sign))
SigList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]]) 
SigList = lapply(SigList, function(x) as.character(x$X))
names(SigList) = c("Adult\nPolyA\nNucleus", "Adult\nPolyA\nCytosol", "Fetal\nPolyA\nNucleus", "Fetal\nPolyA\nCytosol", 
                   "Adult\nRibozero\nNucleus", "Adult\nRibozero\nCytosol", "Fetal\nRibozero\nNucleus", "Fetal\nRibozero\nCytosol")
elementLengths(SigList)
maps = lapply(SigList, function(x) geneMap[which(geneMap$EnsID %in% x),])
elementLengths(maps)
type = lapply(maps, function(x) as.data.frame(x$Type))
freq = lapply(type, function(x) count(x))
TypeFreq = do.call(rbind, freq)
TypeFreq$Group = gsub("\\..*","", rownames(TypeFreq))
colnames(TypeFreq) = c("RNA_Type", "Count", "Group")

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
group = names(SigList)
type = data.frame(RNA.Type = c(rep.int(RNA.Type[1], 8), rep.int(RNA.Type[2], 8),rep.int(RNA.Type[3], 8),rep.int(RNA.Type[4], 8)),
                  Count = NA, 
                  Group = factor(x=rep.int(group, 4), levels=c("Fetal\nRibozero\nNucleus","Fetal\nPolyA\nNucleus","Adult\nRibozero\nNucleus","Adult\nPolyA\nNucleus",
                                                               "Fetal\nRibozero\nCytosol","Fetal\nPolyA\nCytosol","Adult\nRibozero\nCytosol","Adult\nPolyA\nCytosol")))

for (i in 1:length(group)){
  f = freq[[i]]
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = sum(f[which(f$x.Type=="protein_coding"),2], f[which(f$x.Type=="TR_C_gene"),2],f[which(f$x.Type=="polymorphic_pseudogene"),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = sum(f[which(f$x.Type=="pseudogene"),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = sum(f[which(f$x.Type=="antisense"),2], f[which(f$x.Type=="lincRNA"),2],
                                                                               f[which(f$x.Type=="sense_intronic"),2], f[which(f$x.Type=="misc_RNA"),2],
                                                                               f[which(f$x.Type=="processed_transcript"),2], f[which(f$x.Type=="sense_overlapping"),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = sum(f[which(f$x.Type=="miRNA"),2], f[which(f$x.Type=="rRNA"),2],
                                                                                f[which(f$x.Type=="snoRNA"),2], f[which(f$x.Type=="snRNA"),2])
}

# Graph the Frequencies
ggplot(type, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) 