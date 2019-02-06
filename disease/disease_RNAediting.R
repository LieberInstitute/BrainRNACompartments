library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)


load("./Dropbox/sorted_figures/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


## Get disease gene editing sites

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
pgc2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"Symbol"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$Gene.Symbol) %in% geneuniverse ), ] # drop genes that are not present in the test set
aej_sets_expressed$Gene.Symbol = as.character(aej_sets_expressed$Gene.Symbol)
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
splitSets$PGC2 = data.frame(Gene.Symbol = unique(pgc2$Symbol))
splitSets = splitSets[which(names(splitSets)!="SCZ PGC GWAS")]


di = lapply(AEJmap, function(x) editing_anno[overlappingGene %in% x$gencodeID & collapsedconversion=="A:G / T:C",,])
names(di) = c("ASD\n(CNV)","ASD\n(Database)","BPAD\n(GWAS)","Intellectual\nDisability","Neuro-\ndevel.","Neuro-\ndegen.",
              "SCZ\n(CNV)","SCZ\n(Meta\nanalysis)","SCZ\n(SNV)","SCZ\n(PGC2)")


## Percentage of genes edited

elementNROWS(lapply(splitSets, function(x) unique(x$Gene.Symbol)))
#ASD CNV      ASD DATABASE         BPAD GWAS                ID               NDD Neurodegenerative           SCZ CNV SCZ Meta-analysis           SCZ SNV 
#    156               235               117                88                30                46                98                35               197 
#PGC2 
# 404 

elementNROWS(lapply(di, function(x) unique(x$nearestSymbol)))
#ASD\n(CNV)          ASD\n(Database)             BPAD\n(GWAS) Intellectual\nDisability           Neuro-\ndevel.           Neuro-\ndegen. 
#        10                       26                       10                       16                        6                        2 
#SCZ\n(CNV)    SCZ\n(Meta\nanalysis)               SCZ\n(SNV)              SCZ\n(PGC2) 
#        15                        2                       22                       45 

df = data.frame(edited = elementNROWS(lapply(di, function(x) unique(x$nearestSymbol))), group = names(di), 
                total = elementNROWS(lapply(splitSets, function(x) unique(x$Gene.Symbol))))
df$percentage = df$edited/df$total * 100
df$group = factor(df$group, levels = c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)","BPAD\n(GWAS)",
                                       "SCZ\n(SNV)","SCZ\n(PGC2)","SCZ\n(Meta\nanalysis)",
                                       "Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))
df$cat = ifelse(df$group %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)"), "Nuclear in Both", "Not Nuclear")
df[which(df$group %in% c("BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")),"cat"] = "Nuclear in Adult Only"
df$cat = factor(df$cat, levels = c("Nuclear in Both","Nuclear in Adult Only","Not Nuclear"))
df

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/diseaseGene_percentEdited.pdf",width=11.5,height=3.5)
ggplot(df, aes(x = group, y = percentage, fill = cat),colour="black") + geom_bar(stat = "identity", colour="black") +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  geom_text(data=df, aes(x = group, y = percentage, label = total), size=4, nudge_y = 2) +
  labs(fill="") + ylab("Percent") + xlab("") + ggtitle("Percentage of Genes Edited") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "bottom")
# scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "gray43"))
dev.off()

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/diseaseGene_percentEdited2.pdf",width=4.6,height=2.8)
ggplot(df, aes(x = cat, y = percentage, fill = cat),colour="black") + geom_boxplot(colour="black") +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) + geom_jitter() + coord_flip() +
  labs(fill="") + ylab("Percent") + xlab("") + ggtitle("Genes Edited") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
ggplot(df, aes(x = cat, y = total, fill = cat),colour="black") + geom_boxplot(colour="black") +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) + geom_jitter() + coord_flip() +
  labs(fill="") + ylab("Number") + xlab("") + ggtitle("Genes Edited") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none")
dev.off()


res = list("Nuclear in Both" = t.test(df[which(df$cat=="Nuclear in Both"),"percentage"], df[which(df$cat!="Nuclear in Both"),"percentage"]),
           "Nuclear in Adult Only" = t.test(df[which(df$cat=="Nuclear in Adult Only"),"percentage"], df[which(df$cat!="Nuclear in Adult Only"),"percentage"]),
           "Nuclear in Both or Adult Only" = t.test(df[which(df$cat!="Not Nuclear"),"percentage"], df[which(df$cat=="Not Nuclear"),"percentage"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group      Tstat mean.NucGroup mean.Other      pval       FDR
#t                Nuclear in Both -0.1094507      10.92674   11.29958 0.9169525 0.9169525
#t1         Nuclear in Adult Only -0.5054679      10.28438   11.57488 0.6281351 0.9169525
#t2 Nuclear in Both or Adult Only -0.3413419      10.60556   12.06098 0.7521112 0.9169525

res = list("Nuclear in Both" = t.test(df[which(df$cat=="Nuclear in Both"),"edited"], df[which(df$cat!="Nuclear in Both"),"edited"]),
           "Nuclear in Adult Only" = t.test(df[which(df$cat=="Nuclear in Adult Only"),"edited"], df[which(df$cat!="Nuclear in Adult Only"),"edited"]),
           "Nuclear in Both or Adult Only" = t.test(df[which(df$cat!="Not Nuclear"),"edited"], df[which(df$cat=="Not Nuclear"),"edited"]))
res = do.call(rbind, Map(cbind, Group = as.list(names(res)), 
                         lapply(res, function(x) data.frame(Tstat = x$statistic, mean.NucGroup = x$estimate[1], 
                                                            mean.Other = x$estimate[2], pval = x$p.value))))
res$FDR = p.adjust(res$pval, method="fdr")
res
#                           Group    Tstat mean.NucGroup mean.Other       pval       FDR
#t                Nuclear in Both 0.306786      17.00000   14.71429 0.76778548 0.7677855
#t1         Nuclear in Adult Only 1.360278      25.66667   11.00000 0.28646065 0.4296910
#t2 Nuclear in Both or Adult Only 2.338707      21.33333    6.50000 0.04892579 0.1467774
