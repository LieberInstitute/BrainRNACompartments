library(ggplot2)
library(plyr)

load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/characterize_transcriptome/data/DESeq2_results.rda")
load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Library differences in the nucleus
Ln = list(Lnres = data.frame(Lnres), Lnres.down = data.frame(Lnres.down))
SigLn = lapply(Ln, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
SigLn = lapply(SigLn, function(x) data.frame(x, Sign=ifelse(x$log2FoldChange > 0,"UpRibo", "UpPolyA")))
Ln = lapply(SigLn, function(x) list(RiboZero = x[which(x$Sign=="UpRibo"),], PolyA = x[which(x$Sign=="UpPolyA"),]))
type = lapply(Ln, function(x) lapply(x, function(y) geneMap[match(rownames(y),geneMap$gencodeID),]))
Lnfreq = list()
for (i in 1:length(Ln)){Lnfreq[[i]] = Map(cbind,Ln[[i]],type[[i]])}
names(Lnfreq) = names(Ln)
Lnfreq = lapply(Lnfreq, function(x) lapply(x, function(y) count(y$gene_type)))
TypeFreq = lapply(Lnfreq, function(x) do.call(rbind, x))
Group = lapply(TypeFreq, function(x) gsub("\\..*","", rownames(x)))
TypeFreq = Map(cbind, TypeFreq, Group = Group)

# Condense to 4 categories
RNA.Type = c("Protein-coding", "Pseudogene", "Long Non-coding", "Short Non-coding")
type = data.frame(RNA.Type = c(rep.int(RNA.Type[1], 2), rep.int(RNA.Type[2], 2),rep.int(RNA.Type[3], 2),rep.int(RNA.Type[4], 2)),
                  Count = NA, 
                  Group = factor(x=rep.int(c("RiboZero", "PolyA"), 4)))
type.down = data.frame(RNA.Type = c(rep.int(RNA.Type[1], 2), rep.int(RNA.Type[2], 2),rep.int(RNA.Type[3], 2),rep.int(RNA.Type[4], 2)),
                  Count = NA, 
                  Group = factor(x=rep.int(c("RiboZero", "PolyA"), 4)))
Lnres = TypeFreq[["Lnres"]]
Lnres.down = TypeFreq[["Lnres.down"]]
group = c("PolyA", "RiboZero")
for (i in 1:2){
  type[which(type$RNA.Type=="Protein-coding" & type$Group==group[i]),2] = 
    sum(Lnres[which(Lnres$x=="protein_coding" & Lnres$Group==group[i]),2], 
        Lnres[which(Lnres$x=="TR_C_gene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="polymorphic_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="IG_V_gene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="TR_V_gene" & Lnres$Group==group[i]),2])
  type[which(type$RNA.Type=="Pseudogene" & type$Group==group[i]),2] = 
    sum(Lnres[which(Lnres$x=="pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="transcribed_processed_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="unprocessed_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="transcribed_unprocessed_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="processed_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="IG_C_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="transcribed_unitary_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="unitary_pseudogene" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="IG_V_pseudogene" & Lnres$Group==group[i]),2])
  type[which(type$RNA.Type=="Long Non-coding" & type$Group==group[i]),2] = 
    sum(Lnres[which(Lnres$x=="3prime_overlapping_ncRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="antisense" & Lnres$Group==group[i]),2], 
        Lnres[which(Lnres$x=="non_coding" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="lincRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="sense_intronic" & Lnres$Group==group[i]),2], 
        Lnres[which(Lnres$x=="misc_RNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="processed_transcript" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="bidirectional_promoter_lncRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="TEC" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="macro_lncRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="sense_overlapping" & Lnres$Group==group[i]),2])
  type[which(type$RNA.Type=="Short Non-coding" & type$Group==group[i]),2] = 
    sum(Lnres[which(Lnres$x=="miRNA" & Lnres$Group==group[i]),2], 
        Lnres[which(Lnres$x=="Mt_rRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="Mt_tRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="rRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="snoRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="scRNA" & Lnres$Group==group[i]),2],
        Lnres[which(Lnres$x=="snRNA" & Lnres$Group==group[i]),2])
}
for (i in 1:2){
  type.down[which(type.down$RNA.Type=="Protein-coding" & type.down$Group==group[i]),2] = 
    sum(Lnres.down[which(Lnres.down$x=="protein_coding" & Lnres.down$Group==group[i]),2], 
        Lnres.down[which(Lnres.down$x=="TR_C_gene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="polymorphic_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="IG_V_gene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="TR_V_gene" & Lnres.down$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Pseudogene" & type.down$Group==group[i]),2] = 
    sum(Lnres.down[which(Lnres.down$x=="pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="transcribed_processed_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="unprocessed_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="transcribed_unprocessed_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="processed_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="IG_C_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="transcribed_unitary_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="unitary_pseudogene" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="IG_V_pseudogene" & Lnres.down$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Long Non-coding" & type.down$Group==group[i]),2] = 
    sum(Lnres.down[which(Lnres.down$x=="3prime_overlapping_ncRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="antisense" & Lnres.down$Group==group[i]),2], 
        Lnres.down[which(Lnres.down$x=="non_coding" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="lincRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="sense_intronic" & Lnres.down$Group==group[i]),2], 
        Lnres.down[which(Lnres.down$x=="misc_RNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="processed_transcript" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="bidirectional_promoter_lncRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="TEC" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="macro_lncRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="sense_overlapping" & Lnres.down$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Short Non-coding" & type.down$Group==group[i]),2] = 
    sum(Lnres.down[which(Lnres.down$x=="miRNA" & Lnres.down$Group==group[i]),2], 
        Lnres.down[which(Lnres.down$x=="Mt_rRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="Mt_tRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="rRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="snoRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="scRNA" & Lnres.down$Group==group[i]),2],
        Lnres.down[which(Lnres.down$x=="snRNA" & Lnres.down$Group==group[i]),2])
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

ggplot(type.down, aes(x = Group, y = Count, fill = RNA.Type)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
