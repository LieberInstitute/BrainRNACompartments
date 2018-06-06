library(reshape2)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

## Enrichment in genes differentially expressed by fraction

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"EntrezID"]))) # all edited genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$EntrezGene.ID) %in% geneuniverse), ] # drop genes that are not present in the test set
aej_sets_expressed$EntrezGene.ID = as.character(aej_sets_expressed$EntrezGene.ID)
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
inGroup = lapply(sig, function(x) as.character(na.omit(unique(x$EntrezID))))
outGroup = lapply(sig, function(x) geneuniverse[!(geneuniverse %in% as.character(x$EntrezID))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$EntrezGene.ID),sum(!(inG %in% x$EntrezGene.ID)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$EntrezGene.ID), sum(!(outG %in% x$EntrezGene.ID)))
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
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")


enrich[enrich$FDR<=0.05,]
#      Comparison           GeneSet      P.Value OddsRatio         FDR
#1  both_retained           ASD.CNV 1.741394e-05  6.701774 0.001567255
#2  both_retained      ASD.DATABASE 4.417327e-04  4.278013 0.013251982
#7  both_retained           SCZ.CNV 2.568253e-03  5.655885 0.037575888
#14 both_exported                ID 1.798478e-03 13.516748 0.034433895
#33   Ad_retained         BPAD.GWAS 1.912994e-03  2.358939 0.034433895
#54   Ad_exported                ID 2.922569e-03  2.278997 0.037575888
#56   Ad_exported Neurodegenerative 2.791624e-04  3.567415 0.012562308


# what the what??? These are the genes in the both retained set:

genes = list(cnv = cbind(geneMap[which(as.character(geneMap$EntrezID) %in% 
                                         inGroup$both_retained[(inGroup$both_retained %in% splitSets$"ASD CNV"$EntrezGene.ID)]),], 
                         Comparison = "both_retained", variable = "ASD.CNV"),
             db = cbind(geneMap[which(as.character(geneMap$EntrezID) %in% 
                                        inGroup$both_retained[(inGroup$both_retained %in% splitSets$"ASD DATABASE"$EntrezGene.ID)]),],
                        Comparison = "both_retained", variable = "ASD.DATABASE"),
             sz.cnv = cbind(geneMap[which(as.character(geneMap$EntrezID) %in% 
                                            inGroup$both_retained[(inGroup$both_retained %in% splitSets$"SCZ CNV"$EntrezGene.ID)]),],
                            Comparison = "both_retained", variable = "SCZ.CNV"),
             BPAD = cbind(geneMap[which(as.character(geneMap$EntrezID) %in% 
                                          inGroup$Ad_retained[(inGroup$both_retained %in% splitSets$"BPAD GWAS"$EntrezGene.ID)]),],
                          Comparison = "Ad_retained", variable = "BPAD.GWAS"))

write.csv(do.call(rbind, genes), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_ASDenrichedGenes_retainedInBoth.csv")

## plot RPKM of these genes
counts = geneRpkm.down[rownames(geneRpkm.down) %in% c(rownames(genes$cnv),rownames(genes$db)),grep("poly", colnames(geneRpkm.down))]
counts = melt(counts)
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts[grep("53", counts$Var2),"Age"] = "P"
counts[-grep("53", counts$Var2),"Age"] = "A"
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/ASD_rpkm_plots.pdf", width=24, height=3)
ggplot(counts, aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") +
  facet_grid(. ~ sym) +
  ylab("Log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("Expression of Genes Associated with Autism")) + 
  theme(title = element_text(size = 16)) +
  theme(text = element_text(size = 16)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

counts = geneRpkm.down[rownames(geneRpkm.down) %in% rownames(genes$sz.cnv),grep("poly", colnames(geneRpkm.down))]
counts = melt(counts)
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts[grep("53", counts$Var2),"Age"] = "P"
counts[-grep("53", counts$Var2),"Age"] = "A"
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/SZ.CNV_rpkm_plots.pdf", width=10, height=3)
ggplot(counts, aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") +
  facet_grid(. ~ sym) +
  ylab("Log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("Expression of Genes Associated with Schizophrenia")) + 
  theme(title = element_text(size = 16)) +
  theme(text = element_text(size = 16)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

counts = geneRpkm.down[rownames(geneRpkm.down) %in% rownames(genes$BPAD),grep("poly", colnames(geneRpkm.down))]
counts = melt(counts)
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts[grep("53", counts$Var2),"Age"] = "P"
counts[-grep("53", counts$Var2),"Age"] = "A"
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]
pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/BPAD_rpkm_plots.pdf", width=12, height=3)
ggplot(counts, aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() + scale_fill_brewer(palette="Dark2") +
  facet_grid(. ~ sym) +
  ylab("Log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("Expression of Genes Associated with Bipolar")) + 
  theme(title = element_text(size = 16)) +
  theme(text = element_text(size = 16)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()




## Enrichment in genes differentially expressed by age

inGroup = lapply(age.sig, function(x) as.character(na.omit(unique(x$EntrezID))))
outGroup = lapply(age.sig, function(x) geneuniverse[!(geneuniverse %in% as.character(x$EntrezID))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$EntrezGene.ID),sum(!(inG %in% x$EntrezGene.ID)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$EntrezGene.ID), sum(!(outG %in% x$EntrezGene.ID)))
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
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_AgeDEGs.csv")


enrich[enrich$FDR<=0.05,]
#        Comparison      GeneSet      P.Value OddsRatio          FDR
#2  both_decreasing ASD.DATABASE 3.848874e-07  2.185606 1.731994e-05
#4  both_decreasing           ID 1.584933e-08  3.669451 1.426440e-06
#5  both_decreasing          NDD 5.583036e-06  5.646020 1.004946e-04
#10 both_decreasing      SCZ.SNV 3.750855e-06  2.165057 1.004946e-04
#11 both_increasing      ASD.CNV 5.164925e-06  2.333495 1.004946e-04
#13 both_increasing    BPAD.GWAS 1.165816e-03  2.010292 1.498906e-02
#17 both_increasing      SCZ.CNV 9.408016e-05  2.463959 1.411202e-03


# Are the ASD genes the same as the ones that were in "both retained"?

age = rbind(cbind(geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_decreasing[(inGroup$both_decreasing %in% splitSets$"ASD DATABASE"$EntrezGene.ID)]),], 
                  Comparison = "both_decreasing", variable = "ASD.DATABASE"),
            cbind(geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_decreasing[(inGroup$both_decreasing %in% splitSets$"ID"$EntrezGene.ID)]),], 
                  Comparison = "both_decreasing", variable = "ID"),
            cbind(geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_decreasing[(inGroup$both_decreasing %in% splitSets$"NDD"$EntrezGene.ID)]),], 
                  Comparison = "both_decreasing", variable = "NDD"),
            cbind(geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_decreasing[(inGroup$both_decreasing %in% splitSets$"SCZ SNV"$EntrezGene.ID)]),], 
                  Comparison = "both_decreasing", variable = "SCZ.SNV"),
            cbind(geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_increasing[(inGroup$both_increasing %in% splitSets$"ASD CNV"$EntrezGene.ID)]),], 
                  Comparison = "both_increasing", variable = "ASD.CNV"),
            cbind(geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_increasing[(inGroup$both_increasing %in% splitSets$"SCZ CNV"$EntrezGene.ID)]),], 
                  Comparison = "both_increasing", variable = "SCZ.CNV"))
head(age)
genes = list("Age\nASD CNVs\n(Increasing)" = unique(as.character(age[age$Comparison=="both_increasing" & age$variable=="ASD.CNV","gencodeID"])),
             "Age\nASD Database\n(Decreasing)" = unique(as.character(age[age$Comparison=="both_decreasing" & age$variable=="ASD.DATABASE","gencodeID"])),
             "Fraction\nASD CNVs" = unique(as.character(cnv$gencodeID)), "Fraction\nASD Database" = unique(as.character(db$gencodeID)))

venn.diagram(genes, paste0("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/figures/Birnbaum_geneSet_enrichedGenes_overlap.jpeg"), 
             main="Age- vs. Fraction-associated ASD genes", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"), alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)

write.csv(age[which(age$gencodeID %in% cnv$gencodeID | age$gencodeID %in% db$gencodeID),colnames(age)!="gencodeTx"], quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_enrichedGenes_retainedInBoth_and_ageAssociated.csv")


## Limit to 1 LFC difference
inGroup = lapply(sig.1, function(x) as.character(na.omit(unique(x$EntrezID))))
outGroup = lapply(sig.1, function(x) geneuniverse[!(geneuniverse %in% as.character(x$EntrezID))])

enrich = mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$EntrezGene.ID),sum(!(inG %in% x$EntrezGene.ID)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$EntrezGene.ID), sum(!(outG %in% x$EntrezGene.ID)))
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