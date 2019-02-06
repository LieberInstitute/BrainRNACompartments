library(reshape2)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

## Get PGC2 genes

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
pgc2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]


## Enrichment in genes differentially expressed by fraction

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"Symbol"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$Gene.Symbol) %in% geneuniverse ), ] # drop genes that are not present in the test set
aej_sets_expressed$Gene.Symbol = as.character(aej_sets_expressed$Gene.Symbol)
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
splitSets$PGC2 = data.frame(Gene.Symbol = unique(pgc2$Symbol))
splitSets = splitSets[which(names(splitSets)!="SCZ PGC GWAS")]

inGroup = lapply(sig, function(x) as.character(na.omit(unique(x$Symbol))))
outGroup = lapply(sig, function(x) geneuniverse[!(geneuniverse %in% as.character(x$Symbol))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$Gene.Symbol),sum(!(inG %in% x$Gene.Symbol)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$Gene.Symbol), sum(!(outG %in% x$Gene.Symbol)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)
enrich = lapply(enrich, data.frame)
enrich = do.call(rbind, Map(cbind, Comparison = as.list(names(enrich)), lapply(enrich, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")


write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")
enrich = read.csv("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")

enrich[enrich$FDR<=0.05,]
#    X    Comparison           GeneSet      P.Value OddsRatio          FDR
#1   1 both_retained           ASD.CNV 7.298168e-06  7.484935 0.0006568351
#2   2 both_retained      ASD.DATABASE 1.785385e-04  4.853614 0.0040171163
#7   7 both_retained           SCZ.CNV 1.446259e-03  6.463768 0.0144625894
#14 14 both_exported                ID 1.094979e-03 16.104118 0.0123185143
#32 32   Ad_retained      ASD.DATABASE 5.400726e-04  2.145052 0.0077731442
#33 33   Ad_retained         BPAD.GWAS 3.361111e-05  3.137994 0.0015125000
#39 39   Ad_retained           SCZ.SNV 2.758216e-04  2.315756 0.0049647884
#40 40   Ad_retained              PGC2 4.767475e-03  1.677378 0.0429072748
#54 54   Ad_exported                ID 6.045779e-04  2.677269 0.0077731442
#56 56   Ad_exported Neurodegenerative 5.153724e-05  4.251756 0.0015461171


# These are the genes in the retained sets:

genes = list(both_retained.ASD.CNV = cbind(geneMap[which(geneMap$Symbol %in% inGroup$both_retained[(inGroup$both_retained %in% splitSets$"ASD CNV"$Gene.Symbol)]),], 
                                           Comparison = "both_retained", variable = "ASD.CNV"),
             both_retained.ASD.DATABASE = cbind(geneMap[which(geneMap$Symbol %in% inGroup$both_retained[(inGroup$both_retained %in% splitSets$"ASD DATABASE"$Gene.Symbol)]),],
                                                Comparison = "both_retained", variable = "ASD.DATABASE"),
             both_retained.SCZ.CNV = cbind(geneMap[which(geneMap$Symbol %in% inGroup$both_retained[(inGroup$both_retained %in% splitSets$"SCZ CNV"$Gene.Symbol)]),],
                                           Comparison = "both_retained", variable = "SCZ.CNV"),
             Ad_retained.BPAD.GWAS = cbind(geneMap[which(geneMap$Symbol %in% inGroup$Ad_retained[(inGroup$Ad_retained %in% splitSets$"BPAD GWAS"$Gene.Symbol)]),],
                                           Comparison = "Ad_retained", variable = "BPAD.GWAS"),
             Ad_retained.ASD.DATABASE = cbind(geneMap[which(geneMap$Symbol %in% inGroup$Ad_retained[(inGroup$Ad_retained %in% splitSets$"ASD DATABASE"$Gene.Symbol)]),],
                                              Comparison = "Ad_retained", variable = "ASD.DATABASE"),
             Ad_retained.SCZ.SNV = cbind(geneMap[which(geneMap$Symbol %in% inGroup$Ad_retained[(inGroup$Ad_retained %in% splitSets$"SCZ SNV"$Gene.Symbol)]),],
                                         Comparison = "Ad_retained", variable = "SCZ.SNV"),
             Ad_retained.PGC2 = cbind(geneMap[which(geneMap$Symbol %in% inGroup$Ad_retained[which(inGroup$Ad_retained %in% splitSets$"PGC2"$Gene.Symbol[which(splitSets$"PGC2"$Gene.Symbol!="")])]),],
                                      Comparison = "Ad_retained", variable = "PGC2"),
             both_exported.ID = cbind(geneMap[which(geneMap$Symbol %in% inGroup$both_exported[which(inGroup$both_exported %in% splitSets$"ID"$Gene.Symbol)]),],
                                      Comparison = "both_exported", variable = "ID"),
             Ad_exported.ID = cbind(geneMap[which(geneMap$Symbol %in% inGroup$Ad_exported[which(inGroup$Ad_exported %in% splitSets$"ID"$Gene.Symbol)]),],
                                    Comparison = "Ad_exported", variable = "ID"),
             Ad_exported.Neurodegenerative = cbind(geneMap[which(geneMap$Symbol %in% inGroup$Ad_exported[which(inGroup$Ad_exported %in% splitSets$"Neurodegenerative"$Gene.Symbol)]),],
                                                   Comparison = "Ad_exported", variable = "Neurodegenerative"))

write.csv(do.call(rbind, genes), quote=F,
          file="./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_ASDenrichedGenes_retainedInBoth.csv")


### plot RPKM of these genes

genes = read.csv("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_ASDenrichedGenes_retainedInBoth.csv")

## ASD: retained in both
counts = melt(geneRpkm.down[which(rownames(geneRpkm.down) %in% genes[which(genes$variable %in% c("ASD.CNV","ASD.DATABASE") & genes$Comparison=="both_retained"),"gencodeID"]),
                       grep("poly", colnames(geneRpkm.down))])
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts[grep("53", counts$Var2),"Age"] = "Prenatal"
counts[-grep("53", counts$Var2),"Age"] = "Adult"
counts$Age = factor(counts$Age, levels = c("Prenatal", "Adult"))
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]
length(unique(counts$sym))

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/ASD_rpkm_plots_nucBoth.pdf", width=11.5, height=3.25)
ggplot(counts, aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") +
  facet_wrap(. ~ sym, nrow =2) + ylab("Log(RPKM+1)") + xlab("") +
  ggtitle("Nuclear-enriched genes associated with Autism Spectrum Disorder (both ages)") + 
  theme(title = element_text(size = 16), text = element_text(size = 16)) + labs(fill="") + 
  theme(legend.background = element_rect(fill = "transparent"), axis.text.x = element_text(angle = 15, hjust = .5, vjust=.9),
        legend.key = element_rect(fill = "transparent", color = "transparent"),legend.position = c(0.95, 0.25))
dev.off()


## SCZ: retained in both

counts = melt(geneRpkm.down[which(rownames(geneRpkm.down) %in% genes[which(genes$variable=="SCZ.CNV" & genes$Comparison=="both_retained"),"gencodeID"]),
                       grep("poly", colnames(geneRpkm.down))])
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts[grep("53", counts$Var2),"Age"] = "Prenatal"
counts[-grep("53", counts$Var2),"Age"] = "Adult"
counts$Age = factor(counts$Age, levels = c("Prenatal", "Adult"))
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/SZ_rpkm_plots_nucBoth.pdf", width=8, height=2.5)
ggplot(counts, aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") +
  facet_grid(. ~ sym) + ylab("Log(RPKM+1)") + xlab("") +
  ggtitle("Nuclear-enriched genes associated with\nSchizophrenia (both ages)") + 
  theme(title = element_text(size = 16), text = element_text(size = 16)) + labs(fill="") + 
  theme(legend.background = element_rect(fill = "transparent"), axis.text.x = element_text(angle = 15, hjust = .5, vjust=.9),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


## SCZ: Retained in adult 

counts = melt(geneRpkm.down[which(rownames(geneRpkm.down) %in% genes[which(genes$variable %in% c("SCZ.SNV","PGC2") & 
                                                                             genes$Comparison=="Ad_retained"),"gencodeID"]),
                            grep("poly", colnames(geneRpkm.down))])
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts = counts[-grep("53", counts$Var2),]
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/Ad_retained_SCZ_rpkm.pdf", width=38, height=2.25)
ggplot(counts, aes(x=sym, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + ylab("Log(RPKM+1)") + xlab("") +
  ggtitle(paste0("Nuclear-enriched genes associated with Schizophrenia (in Adult)")) + 
  theme(title = element_text(size = 16), text = element_text(size = 16), legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        axis.text.x = element_text(angle = 15, hjust = .5, vjust=.9)) + labs(fill="")
dev.off()


## SCZ: Retained in adult 

counts = melt(geneRpkm.down[which(rownames(geneRpkm.down) %in% genes[which(genes$variable=="BPAD.GWAS"),"gencodeID"]),
                            grep("poly", colnames(geneRpkm.down))])
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts = counts[-grep("53", counts$Var2),]
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/Ad_retained_BPAD_rpkm.pdf", width=13, height=2.25)
ggplot(counts, aes(x=sym, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + ylab("Log(RPKM+1)") + xlab("") +
  ggtitle(paste0("Nuclear-enriched genes associated with Bipolar Affective Disorder (in Adult)")) + 
  theme(title = element_text(size = 16), text = element_text(size = 16), legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        axis.text.x = element_text(angle = 15, hjust = .5, vjust=.9)) + labs(fill="")
dev.off()


## Check diseases with other etiologies
unique(enrich$GeneSet)
enrich[enrich$GeneSet %in% c("ID","NDD","Neurodegenerative") & enrich$FDR<=0.05,]
#    X    Comparison           GeneSet      P.Value OddsRatio         FDR
#14 14 both_exported                ID 1.094979e-03 16.104118 0.012318514
#54 54   Ad_exported                ID 6.045779e-04  2.677269 0.007773144
#56 56   Ad_exported Neurodegenerative 5.153724e-05  4.251756 0.001546117


## update table S2 (Fraction DEGs)

genes = read.csv("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_ASDenrichedGenes_retainedInBoth.csv")

genes$variable = gsub("ASD.CNV", "Autism Spectrum Disorder (from CNVs)", genes$variable)
genes$variable = gsub("ASD.DATABASE", "Autism Spectrum Disorder (from Safari database)", genes$variable)
genes$variable = gsub("BPAD.GWAS", "Bipolar Affective Disorder (from GWAS)", genes$variable)
genes$variable = gsub("ID", "Intellecual Disability", genes$variable)
genes$variable = gsub("NDD", "Syndromal Neurodevelopmental Disorder", genes$variable)
genes$variable = gsub("Neurodegenerative", "Neurodegenerative Disorder", genes$variable)
genes$variable = gsub("SCZ.Meta.analysis", "Schizophrenia (from meta-analysis)", genes$variable)
genes$variable = gsub("SCZ.SNV", "Schizophrenia (from rare SNVs)", genes$variable)
genes$variable = gsub("PGC2", "Schizophrenia (from PGC2)", genes$variable)
genes$variable = gsub("SCZ.CNV", "Schizophrenia (from CNVs)", genes$variable)

tab = openxlsx::read.xlsx("./Dropbox/BrainRNACompartments/paper/Tables_Nuc_Cyt.xlsx", sheet = 2, startRow = 3)
head(tab)
tab$Disease = genes[match(tab$geneID, genes$gencodeID),"variable"]
write.csv(tab, file = "./Dropbox/BrainRNACompartments/paper/newTableS2.csv", quote = F)


## Make Odds Ratio enrichment plot

enrich = read.csv("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")
head(enrich)
enrich$GeneSet = gsub("ASD.CNV", "ASD\n(CNV)", enrich$GeneSet)
enrich$GeneSet = gsub("ASD.DATABASE", "ASD\n(Database)", enrich$GeneSet)
enrich$GeneSet = gsub("BPAD.GWAS", "BPAD\n(GWAS)", enrich$GeneSet)
enrich$GeneSet = gsub("ID", "Intellectual\nDisability", enrich$GeneSet)
enrich$GeneSet = gsub("NDD", "Neuro-\ndevel.", enrich$GeneSet)
enrich$GeneSet = gsub("Neurodegenerative", "Neuro-\ndegen.", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.Meta.analysis", "SCZ\n(Meta\nanalysis)", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.SNV", "SCZ\n(SNV)", enrich$GeneSet)
enrich$GeneSet = gsub("PGC2", "SCZ\n(PGC2)", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.CNV", "SCZ\n(CNV)", enrich$GeneSet)
enrich$GeneSet = factor(enrich$GeneSet, levels = 
                        c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)","BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)","SCZ\n(Meta\nanalysis)",
                          "Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))
enrich$cat = ifelse(enrich$GeneSet %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)"), "Nuclear in Both", "Not Nuclear")
enrich[which(enrich$GeneSet %in% c("BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")),"cat"] = "Nuclear in Adult Only"
enrich$cat = factor(enrich$cat, levels = c("Nuclear in Both","Nuclear in Adult Only","Not Nuclear"))

pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/birnbaum_OR_plot.pdf", width = 6, height = 2.5)
ggplot(enrich[enrich$cat!="Not Nuclear" & enrich$FDR<=0.05,], aes(GeneSet, OddsRatio, fill = cat)) +
  geom_col() + ylab("Odds Ratio") +  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "gray43")) +
  xlab("") + geom_hline(yintercept=1, linetype="dotted") + ggtitle("Nuclear-Enriched Genes") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.title=element_blank(),
        legend.position = c(0.75, 0.8), legend.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()



#### Enrichment in genes differentially expressed by age

inGroup = lapply(age.sig, function(x) as.character(na.omit(unique(x$Symbol))))
outGroup = lapply(age.sig, function(x) geneuniverse[!(geneuniverse %in% as.character(x$Symbol))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$Gene.Symbol),sum(!(inG %in% x$Gene.Symbol)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$Gene.Symbol), sum(!(outG %in% x$Gene.Symbol)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)
enrich = lapply(enrich, data.frame)
enrich = do.call(rbind, Map(cbind, Comparison = as.list(names(enrich)), lapply(enrich, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_AgeDEGs.csv")

enrich[enrich$FDR<=0.05,]
#        Comparison           GeneSet      P.Value OddsRatio          FDR
#2  both_decreasing      ASD.DATABASE 8.179849e-12  2.816400 7.361864e-10
#3  both_decreasing         BPAD.GWAS 5.402583e-03  1.872504 4.420295e-02
#4  both_decreasing                ID 6.222424e-11  4.543574 2.800091e-09
#5  both_decreasing               NDD 6.561070e-07  6.802815 1.180993e-05
#9  both_decreasing           SCZ.SNV 5.877166e-08  2.497970 1.322362e-06
#10 both_decreasing              PGC2 1.402493e-04  1.659534 1.577805e-03
#11 both_increasing           ASD.CNV 4.754962e-08  2.776990 1.322362e-06
#12 both_increasing      ASD.DATABASE 3.476458e-03  1.633693 3.476458e-02
#13 both_increasing         BPAD.GWAS 2.918922e-05  2.476658 3.752900e-04
#16 both_increasing Neurodegenerative 4.604769e-03  2.636742 4.144292e-02
#17 both_increasing           SCZ.CNV 3.357161e-06  2.931991 5.035741e-05


## Limit to 1 LFC difference
inGroup = lapply(sig.1, function(x) as.character(na.omit(unique(x$Symbol))))
outGroup = lapply(sig.1, function(x) geneuniverse[!(geneuniverse %in% as.character(x$Symbol))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$Gene.Symbol),sum(!(inG %in% x$Gene.Symbol)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$Gene.Symbol), sum(!(outG %in% x$Gene.Symbol)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate)
  names(dat) <- c("P.Value","Odds Ratio")
  return(dat)
}), inGroup, outGroup, SIMPLIFY =F)
enrich = lapply(enrich, data.frame)
enrich = do.call(rbind, Map(cbind, Comparison = as.list(names(enrich)), lapply(enrich, function(x) 
  data.frame(GeneSet = colnames(x), P.Value = as.numeric(x["P.Value",]), OddsRatio = as.numeric(x["Odds Ratio",]), 
             row.names=NULL))))
enrich$FDR = p.adjust(enrich$P.Value, method = "fdr")

enrich[enrich$FDR<=0.05,]
# none :(