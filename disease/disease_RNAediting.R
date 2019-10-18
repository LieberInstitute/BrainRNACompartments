library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)

path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "rna_editing/data/unique_editingSites_bySample.rda"))
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))
load(paste0(path, "updated_gene_sets.rda"))

## Get disease gene editing sites

geneuniverse <- na.omit(unique(as.character(
  geneMap[which(geneMap$gencodeID %in% 
                  rownames(Ipres.down)),"gencodeID"])))
splitSets <- lapply(updated, function(f) 
  f[which(f$gencodeID %in% geneuniverse), ])

di <- lapply(splitSets, function(x) 
  editing_anno[overlappingGene %in% x$gencodeID & 
                 collapsedconversion=="A:G / T:C",,])
names(di) = c("Intellectual\nDisability","Neuro-\ndevel.","Neuro-\ndegen.",
              "SCZ\n(GWAS)","BPAD\n(GWAS)","ASD\n(SFARI)","ASD\n(CNV)",
              "SCZ\n(CNV)","SCZ\n(SNV)")


## Percentage of genes edited

data.frame(total = elementNROWS(lapply(splitSets, function(x) unique(x$gencodeID))),
           edited = elementNROWS(lapply(di, function(x) unique(x$overlappingGene))))
#                   total edited
# ID                   88     16
# NDD                  15      6
# Neurodegenerative   173     24
# SCZ.GWAS           1107     89
# BPAD.GWAS           290     36
# ASD.SFARI           589     75
# ASD.CNV             167     26
# SCZ.CNV             103     17
# SCZ.SNV             252     31

df <- data.frame(edited = elementNROWS(lapply(di, function(x) 
                                   unique(x$overlappingGene))), group = names(di), 
                 total = elementNROWS(lapply(splitSets, function(x) 
                                   unique(x$gencodeID))))
df$percentage <- df$edited/df$total * 100
df$group <- factor(df$group, levels = c("ASD\n(CNV)","SCZ\n(CNV)",
                                        "ASD\n(SFARI)","BPAD\n(GWAS)",
                                        "SCZ\n(SNV)","Neuro-\ndegen.",
                                        "SCZ\n(GWAS)", "Neuro-\ndevel.",
                                        "Intellectual\nDisability"))
df$cat <- ifelse(df$group %in% c("ASD\n(CNV)","SCZ\n(CNV)"),
                     "Nuclear in Both", "Not Nuclear")
df[which(df$group %in% c("ASD\n(SFARI)","BPAD\n(GWAS)","SCZ\n(SNV)",
                           "Neuro-\ndegen.")),"cat"] <- "Nuclear in Adult Only"
df$cat <- factor(df$cat, levels = c("Nuclear in Both",
                                    "Nuclear in Adult Only",
                                    "Not Nuclear"))
df

pdf(paste0(path, "disease/figures/diseaseGene_percentEdited.pdf"),
    width = 11.5, height = 3.5)
ggplot(df, aes(x = group, y = percentage, fill = cat), colour = "black") + 
  geom_bar(stat = "identity", colour="black") +
  scale_fill_manual(values = c("cornsilk4", "antiquewhite3", "white")) +
  geom_text(data=df, aes(x = group, y = percentage, label = total), 
            size=4, nudge_y = 2) +
  labs(fill="") + ylab("Percent") + xlab("") + 
  ggtitle("Percentage of Genes Edited") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "bottom")
dev.off()

pdf(paste0(path, "disease/figures/diseaseGene_percentEdited2.pdf"),
    width = 4.6, height = 2.8)
ggplot(df, aes(x = cat, y = percentage, fill = cat), colour="black") + 
  geom_boxplot(colour="black") +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) + 
  geom_jitter() + coord_flip() +
  labs(fill="") + ylab("Percent") + xlab("") + 
  ggtitle("Genes Edited") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "none")
ggplot(df, aes(x = cat, y = total, fill = cat),colour="black") + 
  geom_boxplot(colour="black") +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) + 
  geom_jitter() + coord_flip() +
  labs(fill="") + ylab("Number") + xlab("") + 
  ggtitle("Genes Edited") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "none")
dev.off()


res = list("Nuclear in Both" = t.test(df[which(df$cat=="Nuclear in Both"),"percentage"], 
                                      df[which(df$cat=="Not Nuclear"),"percentage"]),
           "Nuclear in Adult Only" = t.test(df[which(df$cat=="Nuclear in Adult Only"),
                                               "percentage"], 
                                            df[which(df$cat=="Not Nuclear"),"percentage"]),
           "Nuclear in Both or Adult Only" = t.test(df[which(df$cat 
                                                             %in% c("Nuclear in Both",
                                                                    "Nuclear in Adult Only")),
                                                       "percentage"], 
                                                    df[which(df$cat=="Not Nuclear"),
                                                       "percentage"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) 
                           data.frame(Tstat = x$statistic, 
                                      mean.NucGroup = x$estimate[1], 
                                      mean.Other = x$estimate[2], 
                                      pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group       Tstat mean.NucGroup mean.Other      pval
#t                Nuclear in Both  1.00995447      16.03686    13.0416 0.4144495
#t1         Nuclear in Adult Only -0.04174808      12.91840    13.0416 0.9703922
#t2 Nuclear in Both or Adult Only  0.30415117      13.95789    13.0416 0.7870138
#FDR
#t  0.9703922
#t1 0.9703922
#t2 0.9703922
 
res = list("Nuclear in Both" = t.test(df[which(df$cat=="Nuclear in Both"),"edited"], 
                                      df[which(df$cat=="Not Nuclear"),"edited"]),
           "Nuclear in Adult Only" = t.test(df[which(df$cat=="Nuclear in Adult Only"),"edited"], 
                                            df[which(df$cat=="Not Nuclear"),"edited"]),
           "Nuclear in Both or Adult Only" = t.test(df[which(df$cat %in% c("Nuclear in Both",
                                                                           "Nuclear in Adult Only")),
                                                       "edited"], 
                                                    df[which(df$cat=="Not Nuclear"),"edited"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, 
                                                            mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], 
                                                            pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group       Tstat mean.NucGroup mean.Other      pval
#t                Nuclear in Both -0.55063252      21.50000   36.33333 0.6345606
#t1         Nuclear in Adult Only  0.32595163      46.50000   36.33333 0.7631445
#t2 Nuclear in Both or Adult Only  0.06319815      38.16667   36.33333 0.9538404
#FDR
#t  0.9538404
#t1 0.9538404
#t2 0.9538404
