library(org.Hs.eg.db)
library(biomaRt)
library(parallel)
library('derfinder')
library(plyr)
library(scales)
library("ggplot2")

# Library differences in the nucleus
Lnres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Lnres.csv")
SigLn = Lnres[which(Lnres$padj<=0.05 & abs(Lnres$log2FoldChange) >=1),]
dim(SigLn)
SigLn = data.frame(SigLn, Sign= ifelse(SigLn$log2FoldChange > 0,"UpRibo", "UpPolyA"))
Ln = list(Ribozero = SigLn[which(SigLn$Sign=="UpRibo"),], PolyA = SigLn[which(SigLn$Sign=="UpPolyA"),])
Ln = lapply(Ln, function(x) as.data.frame(x$Type))
Lnfreq = lapply(Ln, function(x) count(x))
TypeFreq = do.call(rbind, Lnfreq)
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

# Condense to 4 categories
RNA.Type = c("Protein-coding", "Pseudogene", "Long Non-coding", "Short Non-coding")
type = data.frame(RNA.Type = c(rep.int(RNA.Type[1], 2), rep.int(RNA.Type[2], 2),rep.int(RNA.Type[3], 2),rep.int(RNA.Type[4], 2)),
                  Count = NA, 
                  Group = factor(x=rep.int(c("Ribozero", "PolyA"), 4)))
for (i in 1:2){
  f = Lnfreq[[i]]
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x.Type=="protein_coding"),2], 
        f[which(f$x.Type=="TR_C_gene"),2],
        f[which(f$x.Type=="polymorphic_pseudogene"),2],
        f[which(f$x.Type=="IG_V_gene"),2],
        f[which(f$x.Type=="TR_V_gene"),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = 
    sum(f[which(f$x.Type=="pseudogene"),2],
        f[which(f$x.Type=="IG_C_pseudogene"),2],
        f[which(f$x.Type=="IG_V_pseudogene"),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x.Type=="3prime_overlapping_ncrna"),2],
        f[which(f$x.Type=="antisense"),2], 
        f[which(f$x.Type=="lincRNA"),2],
        f[which(f$x.Type=="sense_intronic"),2], 
        f[which(f$x.Type=="misc_RNA"),2],
        f[which(f$x.Type=="processed_transcript"),2], 
        f[which(f$x.Type=="sense_overlapping"),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = 
    sum(f[which(f$x.Type=="miRNA"),2], 
        f[which(f$x.Type=="rRNA"),2],
        f[which(f$x.Type=="snoRNA"),2], 
        f[which(f$x.Type=="snRNA"),2])
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