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
enrich = rbind(do.call(rbind, Map(cbind, lapply(enrich, function(x) x["P.Value",]),Value = "P.Value")),
               do.call(rbind, Map(cbind, lapply(enrich, function(x) x["Odds Ratio",]),Value = "Odds Ratio")))
enrich = cbind(Comparison = gsub("1","",rownames(enrich)), enrich)
enrich = melt(enrich)
enrich$group = paste(enrich$Comparison, enrich$variable, sep=":")
enrich = enrich[order(enrich$Comparison),]
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")


enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=0.05),"group"],1:4]
#       Comparison      Value          variable        value
#60    Ad_exported    P.Value                ID 2.922569e-03
#69    Ad_exported Odds Ratio                ID 2.278997e+00
#96    Ad_exported    P.Value Neurodegenerative 2.791624e-04
#105   Ad_exported Odds Ratio Neurodegenerative 3.567415e+00
#22    Ad_retained    P.Value      ASD.DATABASE 6.205271e-03
#31    Ad_retained Odds Ratio      ASD.DATABASE 1.818759e+00
#40    Ad_retained    P.Value         BPAD.GWAS 1.912994e-03
#49    Ad_retained Odds Ratio         BPAD.GWAS 2.358939e+00
#166   Ad_retained    P.Value           SCZ.SNV 7.555499e-03
#175   Ad_retained Odds Ratio           SCZ.SNV 1.846415e+00
#56  both_exported    P.Value                ID 1.798478e-03
#65  both_exported Odds Ratio                ID 1.351675e+01
#1   both_retained    P.Value           ASD.CNV 1.741394e-05
#10  both_retained Odds Ratio           ASD.CNV 6.701774e+00
#19  both_retained    P.Value      ASD.DATABASE 4.417327e-04
#28  both_retained Odds Ratio      ASD.DATABASE 4.278013e+00
#109 both_retained    P.Value           SCZ.CNV 2.568253e-03
#118 both_retained Odds Ratio           SCZ.CNV 5.655885e+00
#45    interacting    P.Value         BPAD.GWAS 4.395678e-02
#54    interacting Odds Ratio         BPAD.GWAS 3.129840e+00

enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=(0.05/90)),"group"],1:4]
#       Comparison      Value          variable        value
#96    Ad_exported    P.Value Neurodegenerative 2.791624e-04
#105   Ad_exported Odds Ratio Neurodegenerative 3.567415e+00
#1   both_retained    P.Value           ASD.CNV 1.741394e-05
#10  both_retained Odds Ratio           ASD.CNV 6.701774e+00
#19  both_retained    P.Value      ASD.DATABASE 4.417327e-04
#28  both_retained Odds Ratio      ASD.DATABASE 4.278013e+00

p = enrich[enrich$Value=="P.Value",]
p$padj = p.adjust(p$value, method = "fdr")
p[p$padj<=0.05,c(-2,-5)]
enrich[enrich$group %in% p[p$padj<=0.05,"group"],]
#       Comparison      Value          variable        value                         group
#60    Ad_exported    P.Value                ID 2.922569e-03                Ad_exported:ID
#69    Ad_exported Odds Ratio                ID 2.278997e+00                Ad_exported:ID
#96    Ad_exported    P.Value Neurodegenerative 2.791624e-04 Ad_exported:Neurodegenerative
#105   Ad_exported Odds Ratio Neurodegenerative 3.567415e+00 Ad_exported:Neurodegenerative
#40    Ad_retained    P.Value         BPAD.GWAS 1.912994e-03         Ad_retained:BPAD.GWAS
#49    Ad_retained Odds Ratio         BPAD.GWAS 2.358939e+00         Ad_retained:BPAD.GWAS
#56  both_exported    P.Value                ID 1.798478e-03              both_exported:ID
#65  both_exported Odds Ratio                ID 1.351675e+01              both_exported:ID
#1   both_retained    P.Value           ASD.CNV 1.741394e-05         both_retained:ASD.CNV
#10  both_retained Odds Ratio           ASD.CNV 6.701774e+00         both_retained:ASD.CNV
#19  both_retained    P.Value      ASD.DATABASE 4.417327e-04    both_retained:ASD.DATABASE
#28  both_retained Odds Ratio      ASD.DATABASE 4.278013e+00    both_retained:ASD.DATABASE
#109 both_retained    P.Value           SCZ.CNV 2.568253e-03         both_retained:SCZ.CNV
#118 both_retained Odds Ratio           SCZ.CNV 5.655885e+00         both_retained:SCZ.CNV

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
counts = geneRpkm.down[rownames(geneRpkm.down) %in% c(rownames(cnv),rownames(db)),grep("poly", colnames(geneRpkm.down))]
counts = melt(counts)
counts[grep("N", counts$Var2),"Fraction"] = "Nucleus"
counts[grep("C", counts$Var2),"Fraction"] = "Cytoplasm"
counts[grep("53", counts$Var2),"Age"] = "P"
counts[-grep("53", counts$Var2),"Age"] = "A"
counts$sym = geneMap[match(counts$Var1, geneMap$gencodeID),"Symbol"]

pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/ASD_rpkm_plots.pdf", width=24, height=3)
ggplot(counts, aes(x=Age, y=log(value+1), fill=Fraction), color=Fraction) + 
  geom_boxplot() +
  facet_grid(. ~ sym) +
  ylab("Log(RPKM+1)") + 
  xlab("") +
  ggtitle(paste0("Expression of Genes Associated with ASD")) + 
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
enrich = rbind(do.call(rbind, Map(cbind, lapply(enrich, function(x) x["P.Value",]),Value = "P.Value")),
                   do.call(rbind, Map(cbind, lapply(enrich, function(x) x["Odds Ratio",]),Value = "Odds Ratio")))
enrich = cbind(Comparison = gsub("1","",rownames(enrich)), enrich)
enrich = melt(enrich)
enrich$group = paste(enrich$Comparison, enrich$variable, sep=":")
enrich = enrich[order(enrich$Comparison),]
write.csv(enrich, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_AgeDEGs.csv")


enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=(0.05/100)),"group"],1:4]
#         Comparison      Value     variable        value
#19  both_decreasing    P.Value ASD.DATABASE 3.848874e-07
#28  both_decreasing Odds Ratio ASD.DATABASE 2.185606e+00
#55  both_decreasing    P.Value           ID 1.584933e-08
#64  both_decreasing Odds Ratio           ID 3.669451e+00
#73  both_decreasing    P.Value          NDD 5.583036e-06
#82  both_decreasing Odds Ratio          NDD 5.646020e+00
#163 both_decreasing    P.Value      SCZ.SNV 3.750855e-06
#172 both_decreasing Odds Ratio      SCZ.SNV 2.165057e+00
#2   both_increasing    P.Value      ASD.CNV 5.164925e-06
#11  both_increasing Odds Ratio      ASD.CNV 2.333495e+00
#110 both_increasing    P.Value      SCZ.CNV 9.408016e-05
#119 both_increasing Odds Ratio      SCZ.CNV 2.463959e+00


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
enrich = rbind(do.call(rbind, Map(cbind, lapply(enrich, function(x) x["P.Value",]),Value = "P.Value")),
               do.call(rbind, Map(cbind, lapply(enrich, function(x) x["Odds Ratio",]),Value = "Odds Ratio")))
enrich = cbind(Comparison = gsub("1","",rownames(enrich)), enrich)
enrich = melt(enrich)
enrich$group = paste(enrich$Comparison, enrich$variable, sep=":")
enrich = enrich[order(enrich$Comparison),]


enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=0.05),"group"],1:4]
#       Comparison      Value          variable        value
#40    Ad_retained    P.Value         BPAD.GWAS  0.005944548
#49    Ad_retained Odds Ratio         BPAD.GWAS  2.898079796
#148   Ad_retained    P.Value      SCZ.PGC.GWAS  0.026527536
#157   Ad_retained Odds Ratio      SCZ.PGC.GWAS  2.532200879
#92  both_exported    P.Value Neurodegenerative  0.034813888
#101 both_exported Odds Ratio Neurodegenerative 30.571697580
#19  both_retained    P.Value      ASD.DATABASE  0.004371475
#28  both_retained Odds Ratio      ASD.DATABASE  4.143867756
#99    interacting    P.Value Neurodegenerative  0.038050368
#108   interacting Odds Ratio Neurodegenerative  6.818302145

p = enrich[enrich$Value=="P.Value",]
p$padj = p.adjust(p$value, method = "fdr")
enrich[enrich$group %in% enrich[which(enrich$Value=="P.Value" & enrich$value<=(0.05/90)),"group"],1:4]
# none :(
enrich[enrich$group %in% p[p$padj<=0.05,"group"],]
# none :(

