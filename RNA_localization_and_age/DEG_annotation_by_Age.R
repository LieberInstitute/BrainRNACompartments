library(parallel)
library('derfinder')
library(plyr)
library(scales)
library("ggplot2")

# Make list of length of significant genes in a list
AgeList = list(Cpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Cpres.csv"),
               Npres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Npres.csv"),
               Crres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Crres.csv"),
               Nrres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Nrres.csv"))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "DownFetal"))
Sign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(Sign, function(x) split(x, x$Sign))
DirList = list(Cpres.Up = DirList[["Cpres"]][["UpFetal"]], Cpres.Down = DirList[["Cpres"]][["DownFetal"]],
               Npres.Up = DirList[["Npres"]][["UpFetal"]], Npres.Down = DirList[["Npres"]][["DownFetal"]],
               Crres.Up = DirList[["Crres"]][["UpFetal"]], Crres.Down = DirList[["Crres"]][["DownFetal"]],
               Nrres.Up = DirList[["Nrres"]][["UpFetal"]], Nrres.Down = DirList[["Nrres"]][["DownFetal"]]) 
names(DirList) = c("Cytosol\nPolyA\nUp", "Cytosol\nPolyA\nDown", "Nucleus\nPolyA\nUp", "Nucleus\nPolyA\nDown", 
                   "Cytosol\nRibozero\nUp", "Cytosol\nRibozero\nDown", "Nucleus\nRibozero\nUp", "Nucleus\nRibozero\nDown")
elementLengths(SigAgeList)

DirList = lapply(DirList, function(x) as.character(x$X))
elementLengths(DirList)
maps = lapply(DirList, function(x) geneMap[which(geneMap$EnsID %in% x),])
elementLengths(maps)
type = lapply(maps, function(x) as.data.frame(x$Type))
freq = lapply(type, function(x) count(x))
TypeFreq = do.call(rbind, freq)
TypeFreq$Group = as.factor(gsub("\\..*","", rownames(TypeFreq)))
colnames(TypeFreq) = c("RNA_Type", "Count", "Group")
TypeFreq = TypeFreq[order(TypeFreq$RNA_Type),]
rownames(TypeFreq) = NULL

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