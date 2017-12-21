library(reshape2)
library(VennDiagram)

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

# what the what??? These are the genes in the both retained set:

cnv = geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_retained[(inGroup$both_retained %in% splitSets$"ASD CNV"$EntrezGene.ID)]),]
#Chr    Start      End Strand Length            gencodeID       ensemblID      gene_type   Symbol EntrezID Class meanExprs NumTx
#ENSG00000205363.5_1  chr15 74027772 74045088      -   6573  ENSG00000205363.5_1 ENSG00000205363 protein_coding C15orf59   388135 InGen  6.857281     4
#ENSG00000184702.17_2 chr22 19701987 19712295      +   5914 ENSG00000184702.17_2 ENSG00000184702 protein_coding    SEPT5     5413 InGen 60.710733    16
#ENSG00000099889.13_1 chr22 19957419 20004331      -   5621 ENSG00000099889.13_1 ENSG00000099889 protein_coding    ARVCF      421 InGen 10.841653    11
#ENSG00000099899.14_2 chr22 20099389 20104915      -   3963 ENSG00000099899.14_2 ENSG00000099899 protein_coding   TRMT2A    27037 InGen  5.143360    20
#ENSG00000099904.15_1 chr22 20116979 20135530      +   5542 ENSG00000099904.15_1 ENSG00000099904 protein_coding   ZDHHC8    29801 InGen  6.627343     7
#ENSG00000099949.18_2 chr22 21333751 21353327      +   7556 ENSG00000099949.18_2 ENSG00000099949 protein_coding    LZTR1     8216 InGen  9.224593    18
#ENSG00000100241.20_2 chr22 50883429 50913453      -   9441 ENSG00000100241.20_2 ENSG00000100241 protein_coding     SBF1     6305 InGen 23.156485     8
#ENSG00000100258.17_2 chr22 50941378 50946120      -   3054 ENSG00000100258.17_2 ENSG00000100258 protein_coding     LMF2    91289 InGen 10.247724     6
#ENSG00000251322.7_2  chr22 51112843 51171726      +   7455  ENSG00000251322.7_2 ENSG00000251322 protein_coding   SHANK3    85358 InGen  5.138953     3

db = geneMap[which(as.character(geneMap$EntrezID) %in% inGroup$both_retained[(inGroup$both_retained %in% splitSets$"ASD DATABASE"$EntrezGene.ID)]),]
#Chr     Start       End Strand Length            gencodeID       ensemblID      gene_type  Symbol EntrezID Class meanExprs NumTx
#ENSG00000113805.8_1   chr3  74311719  74570291      -   5109  ENSG00000113805.8_1 ENSG00000113805 protein_coding   CNTN3     5067 InGen  6.638813     2
#ENSG00000164418.19_1  chr6 101846664 102517958      +   7375 ENSG00000164418.19_1 ENSG00000164418 protein_coding   GRIK2     2898 InGen  9.277783    11
#ENSG00000107816.17_2 chr10 102756375 102767593      +   6547 ENSG00000107816.17_2 ENSG00000107816 protein_coding   LZTS2    84445 InGen  4.744216     7
#ENSG00000135127.11_2 chr12 120427673 120532298      +   5484 ENSG00000135127.11_2 ENSG00000135127 protein_coding  CCDC64    92558 InGen  5.605429    11
#ENSG00000196557.6    chr16   1203241   1271771      +   8484    ENSG00000196557.6 ENSG00000196557 protein_coding CACNA1H     8912 InGen  9.160196     9
#ENSG00000103197.16_2 chr16   2097466   2138721      +  11121 ENSG00000103197.16_2 ENSG00000103197 protein_coding    TSC2     7249 InGen 11.111742    26
#ENSG00000071655.17_2 chr19   1573595   1592800      -   7845 ENSG00000071655.17_2 ENSG00000071655 protein_coding    MBD3    53615 InGen  7.448974    11
#ENSG00000104936.17_2 chr19  46272975  46285810      -   5807 ENSG00000104936.17_2 ENSG00000104936 protein_coding    DMPK     1760 InGen  6.655554    20
#ENSG00000251322.7_2  chr22  51112843  51171726      +   7455  ENSG00000251322.7_2 ENSG00000251322 protein_coding  SHANK3    85358 InGen  5.138953     3

write.csv(rbind(cbind(cnv, Comparison = "both_retained", variable = "ASD.CNV"), cbind(db, Comparison = "both_retained", variable = "ASD.DATABASE")), quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_ASDenrichedGenes_retainedInBoth.csv")


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