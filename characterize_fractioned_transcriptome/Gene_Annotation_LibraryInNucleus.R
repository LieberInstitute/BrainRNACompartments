library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Library differences in the nucleus

Ln = data.frame(Lnres.down)
SigLn = Ln[which(Ln$padj<=0.05 & abs(Ln$log2FoldChange) >=1),]
SigLn = data.frame(SigLn, 
                   Sign = ifelse(SigLn$log2FoldChange > 0,"UpRibo", "UpPolyA"))
Ln = list(RiboZero = SigLn[which(SigLn$Sign=="UpRibo"),], 
          PolyA = SigLn[which(SigLn$Sign=="UpPolyA"),])
type = lapply(Ln, function(y) geneMap[match(rownames(y),geneMap$gencodeID),])

Lnfreq = Map(cbind, Ln, type)
Lnfreq = lapply(Lnfreq, function(y) count(y$gene_type))
TypeFreq = do.call(rbind, Lnfreq)
TypeFreq$Group = gsub("\\..*","", rownames(TypeFreq))
TypeFreq$Group = ifelse(TypeFreq$Group=="PolyA", "poly(A)", "Ribo-Zero")


# Condense to 4 categories

RNA.Type = c("Protein-coding", "Pseudogene", 
             "Long Non-coding", "Short Non-coding")
type.down = data.frame(RNA.Type = c(rep.int(RNA.Type[1], 2), 
                                    rep.int(RNA.Type[2], 2),
                                    rep.int(RNA.Type[3], 2),
                                    rep.int(RNA.Type[4], 2)),
                  Count = NA, 
                  Group = factor(x=rep.int(c("Ribo-Zero", "poly(A)"), 4)))
group = c("poly(A)", "Ribo-Zero")

for (i in 1:2){
  type.down[which(type.down$RNA.Type=="Protein-coding" & type.down$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$x=="protein_coding" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$x=="TR_C_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="polymorphic_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="IG_V_gene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="TR_V_gene" & TypeFreq$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Pseudogene" & type.down$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$x=="pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="transcribed_processed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="unprocessed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="transcribed_unprocessed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="processed_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="IG_C_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="transcribed_unitary_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="unitary_pseudogene" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="IG_V_pseudogene" & TypeFreq$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Long Non-coding" & type.down$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$x=="3prime_overlapping_ncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="antisense" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$x=="non_coding" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="lincRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="sense_intronic" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$x=="misc_RNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="processed_transcript" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="bidirectional_promoter_lncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="TEC" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="macro_lncRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="sense_overlapping" & TypeFreq$Group==group[i]),2])
  type.down[which(type.down$RNA.Type=="Short Non-coding" & type.down$Group==group[i]),2] = 
    sum(TypeFreq[which(TypeFreq$x=="miRNA" & TypeFreq$Group==group[i]),2], 
        TypeFreq[which(TypeFreq$x=="Mt_rRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="Mt_tRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="rRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="snoRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="scRNA" & TypeFreq$Group==group[i]),2],
        TypeFreq[which(TypeFreq$x=="snRNA" & TypeFreq$Group==group[i]),2])
}
type.down


# Graph the Frequencies

pdf(paste0("./Dropbox/BrainRNACompartments/paper/Genome_Research_revision/",
           "Figures/FigS2/",
           "FigS2A.nuclear_DEG_byLibrary_annotation_downsampled.pdf"),
    width = 5.5, height = 4.5)
ggplot(type.down, aes(x = Group, y = Count, fill = RNA.Type)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Accent") +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle("Gene Annotation") +
  theme(title = element_text(size = 20),
        text = element_text(size = 20),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
dev.off()
