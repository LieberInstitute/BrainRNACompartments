library(reshape2)
library(ggplot2)
library(GenomicRanges)

path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "QC_section/data/rpkmCounts_combined_NucVSCyt_n23.rda"))
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))

original_sets <- openxlsx::read.xlsx(paste0(path, 
                                  'Birnbaum_2013_AJP_Supplementary_table.xlsx'))
original_sets <- split(original_sets, original_sets$Gene.Set)
data.frame(elementNROWS(original_sets))
#ASD CNV                                   176
#ASD SFARI                                 241
#BPAD GWAS                                 123
#ID                                         88
#NDD                                        31
#Neurodegenerative                          46
#SCZ CNV                                   113
#SCZ Meta-analysis                          36
#SCZ PGC GWAS                              106
#SCZ SNV                                   212


################################################################################
# Go through lists and update them
################################################################################

# https://ajp.psychiatryonline.org/doi/suppl/10.1176/appi.ajp.2014.13111452/suppl_file/758_ds001.pdf

# SCZ Meta-analysis is old, let's discard it

updated <- original_sets[names(original_sets)!="SCZ Meta-analysis"]


# SCZ PGC GWAS should be updated to the genes in the 145 new GWAS loci

SCZ.GWAS <- openxlsx::read.xlsx(paste0(path, 
                                       'disease/data/41588_2018_59_MOESM3_ESM.xlsx'),
                                sheet = 5, startRow = 8, colNames = TRUE)
head(SCZ.GWAS)
SCZ.GWAS$Chromosome <- paste0("chr", SCZ.GWAS$Chromosome)
SCZ.GWAS <- makeGRangesFromDataFrame(SCZ.GWAS, seqnames.field = "Chromosome",
                                     start.field = "Start.(BP)",
                                     end.field = "End.(BP)",
                                     keep.extra.columns = TRUE)
SCZ.GWAS

updated <- updated[names(updated)!="SCZ PGC GWAS"]
updated$SCZ.GWAS <- geneMap[subjectHits(
  findOverlaps(SCZ.GWAS, 
               makeGRangesFromDataFrame(geneMap))),]


# BPAD GWAS should be update to the newest BP GWAS results

BPAD.GWAS <- openxlsx::read.xlsx(paste0(path, 
                                 'disease/data/41588_2019_397_MOESM3_ESM.xlsx'),
                                 sheet = 4, rows = 6:36, colNames = FALSE)
head(BPAD.GWAS)
BPAD.GWAS$X6 <- paste0("chr", BPAD.GWAS$X6)
BPAD.GWAS <- makeGRangesFromDataFrame(BPAD.GWAS, seqnames.field = "X6",
                                     start.field = "X2",
                                     end.field = "X3")
BPAD.GWAS

updated <- updated[names(updated)!="BPAD GWAS"]
updated$BPAD.GWAS <- geneMap[subjectHits(
  findOverlaps(BPAD.GWAS, 
               makeGRangesFromDataFrame(geneMap))),]


# NDD was based on https://www.cell.com/abstract/S0092-8674(12)00411-4
  # after looking into it I think I'll leave this one the same

updated$NDD <- geneMap[which(geneMap$Symbol %in% 
                               c("LDLRAD4", original_sets$NDD$Gene.Symbol)),]


# ASD SFARI should be refreshed.

# http://sfari.org/resources/sfari-gene
# last updated June 20, 2019

sfari <- read.csv(paste0(path, 
            'disease/data/SFARI-Gene_genes_08-29-2019release_10-06-2019export.csv'),
            colClasses = "character")
head(sfari)
sfari <- sfari[which(sfari$syndromic==0),]
sfari <- sfari[-grep("Syndromic", sfari$genetic.category),]
sfari <- sfari[which(sfari$genetic.category!=""),]
sfari <- sfari[-which(sfari$gene.symbol %in% updated$`ASD CNV`$Gene.Symbol),]
sfari$chromosome <- paste0("chr", sfari$chromosome)
head(sfari)

updated <- updated[names(updated)!="ASD DATABASE"]
updated$ASD.SFARI <- geneMap[which(geneMap$Symbol %in% sfari$gene.symbol),]


# Neurodegenerative list should be refreshed to newest database lists.

#included Alzheimer’s Disease, Parkinson’s Disease, Fronto-temporal Dementia, 
#Amyotrophic Lateral Sclerosis, Huntington’s Disease, and Multiple Sclerosis. 
# add https://www.ncbi.nlm.nih.gov/pubmed/29566793
# add https://www-nature-com.proxy1.library.jhu.edu/articles/ng.3955
# skip https://www.nature.com/articles/s41588-018-0311-9
# add https://www.nature.com/articles/s41588-019-0358-2

gwas <- list(ALS = c("KIF5A", "TNIP1", "C9orf72", "TBK1", "UNC13A", "C21orf2"),
             PD = c("GBA","NUCKS1","SLC41A1","SIPA1L2","TMEM163","CCNT2","STK39",
                    "CHMP2B","MCCC1","TMEM175","DGKQ","FAM200B","CD38","FAM47E",
                    "SNCA","HLA-DRB6","HLA-DQA1","KLHL7","NUPL2","GPNMB","MICU3",
                    "BAG3","DLG2","MIR4697","LRRK2","OGFOD2","GCH1","TMEM229B",
                    "VPS13C","ZNF646","KAT8","ARHGAP27","CRHR1","SPPL2C","MAPT",
                    "STH","KANSL1","SYT4","LSM7","DDRGK1","ITPKB","IL1R2","SCN3A",                                                     
                    "SATB1","NCKIPSD","CDC71","ALAS1","TLR9","DNAH1","BAP1",
                    "PHF7","NISCH","STAB1","ITIH3","ITIH4","ANK2","CAMK2D","ELOVL7",
                    "ZNF184","CTSB","SORBS3","PDLIM2","C8orf58","BIN3","SH3GL2",
                    "FAM171A1","GALC","COQ7","TOX3","ATP6V0A1","PSMC3IP","TUBG2"),
             AD = c("CR1","AL691452.1","RP11-78B10.2","CRIL","INPP5D","HLA-DRA",
                    "HLA-DRB5","RNU1-61P","HLA-DRB6","HLA-DRB1","HLA-DQA1",
                    "HLA-DQB1","HLA-DQB1-AS1","LOC101929555","UNC5CL","TSPO2",
                    "APOBEC2","OARD1","NFYA","ADCY10P1","TREML1","TREM2","TREML2",
                    "TREML3P","TREML4","RNA5SP207","TREML5P","TREM1","RNU6-643P",
                    "NCR2","AL136967.1","RP1-149M18.3","AL355353.1","RP11-385F7.1",
                    "CD2AP","Y_RNA","ADGRF2","PMS2P1","STAG3L5P","PVRIG2P",
                    "PILRB","STAG3L5P","PVRIG2P","MIR6840","PILRB","PILRA",
                    "ZCWPW1","MEPCE","RP11","758P17.2","PPP1R35","RP11-758P17.3",
                    "C7orf61","RN7SL161P","TSC22D4","AC092849.1","NYAP1",
                    "RP11-44M6.1","RN7SL416P","AGFG2","SAP25","RP11-44M6.7",
                    "LRCH4","ZASP","FBXO24","PCOLCE-AS1","EPHA1","EPHA1-AS1",
                    "PTK2B","CLU","MIR6843","RP11-138I18.2","MYBPC3","SPI1",
                    "RP11-750H9.5","AC090559.2","MIR4487","SLC39A13","PSMC3",
                    "RAPSN","RNU6-1302P","snoU13","MS4A2","AP001257.1",
                    "MS4A6A","MS4A4E","MIR6503","MS4A4A","PICALM","RNU6-560P",
                    "SORL1","FERMT2","SLC24A4","ADAM10","snoU13","HSP90AB4P",
                    "RN7SKP95","U3","LOC101928725","RP11-30K9.6","FAM63B",
                    "RP11-30K9.4","C16orf62","AC002550.5","KNOP1","IQCK",
                    "CTD-2380F24.1","TANC2","RP11-269G24.3","RP11-269G24.4",
                    "CYB561","RP11-269G24.6","ABCA7","ARHGAP45","HMHA1",
                    "CSTF1","CASS4"))

# ALS comes from Nicolas (2018) "gene" column in table 1
# PD comes from Chang et al (2017) "Candidate gene" columns in tables 1 & 2
# AD comes from Kunkle et al (2019) Table S8 column: "Gene(s) and Gene Features in LD Block"

gwas <- unique(unlist(gwas))
gwas <- c(gwas,"NCR2","GPR111","SNORD3A","HMHA1", 
          original_sets$Neurodegenerative$Gene.Symbol)
updated$Neurodegenerative <- geneMap[which(geneMap$Symbol %in% gwas),]


# ID list I'm leaving alone.

updated$ID <- geneMap[which(geneMap$Symbol %in% original_sets$ID$Gene.Symbol),]


# ASD CNV: Not updating as the SFARI list will reflect any missed genes.

updated$ASD.CNV <- updated$`ASD CNV`
updated <- updated[names(updated)!="ASD CNV"]
updated$ASD.CNV <- geneMap[which(geneMap$Symbol %in% c(updated$ASD.CNV$Gene.Symbol, 
                                                       "F2R","F2RL3","PKM",
                                                       "CELF6","CDC45","FAN1",
                                                       "GID8","NPAP1","HEXA-AS1",
                                                       "TANGO2","EOGT")),]


# SCZ CNV: could add from https://www.nature.com/articles/ng.3725, but it looks 
  # like many genes are already included

updated$SCZ.CNV <- updated$`SCZ CNV`
updated <- updated[names(updated)!="SCZ CNV"]
updated$SCZ.CNV <- geneMap[which(geneMap$Symbol %in% c(original_sets$`SCZ CNV`$Gene.Symbol, 
                                                       "CDC45","FAN1","TANGO2",
                                                       "CEP19","SLC51A")),]
c(original_sets$`SCZ CNV`$Gene.Symbol, 
  "CDC45","FAN1","TANGO2",
  "CEP19","SLC51A")[-which(c(original_sets$`SCZ CNV`$Gene.Symbol, 
                             "CDC45","FAN1","TANGO2",
                             "CEP19","SLC51A") %in% geneMap$Symbol)]

# SCZ SNV
# after looking at https://www-nature-com.proxy1.library.jhu.edu/articles/nature12929
# https://www.nature.com/articles/nn.4402 and https://www.nature.com/articles/ng.3903
# only add highest evidence dURV genes from table S2 (Genovese, Nat Neuro, 2016)

genoveseTabS2= c("KL","STXBP2","PCYT1A","CPVL","KDM5B","HEPHL1","AP5Z1","TTN",
                 "CSNK1G3","KCNH7","ZDBF2","ITPR1","PRSS16","CCDC30","LAMB2",
                 "PCSK6","PRDM10","KIAA1033","NDST2","TMPO","KIAA1432","CDO1",
                 "OGDHL","SORCS1","KDM5C","PLEKHG7","JARID2","ARFGEF1","ASAH1",
                 "PMS2","TRIM2","TRAF3IP1","GOLGA4","MPP4","KLHL25","WBSCR22",
                 "QRICH1","ISYNA1","VSTM2A","DBNL","THG1L","LY75","LY75-CD302",
                 "DNM1L","CHRNE","MTAP","TBC1D2B","ETFBKMT","CELF2","KDM2B",
                 "ARHGEF38","TENM1","NPRL2","LDLRAD4","DCDC1")

updated$SCZ.SNV <- geneMap[c(geneMap$Symbol %in% c(original_sets$`SCZ SNV`$Gene.Symbol, 
                                                   genoveseTabS2,"METTL20","CELF2",
                                                   "KDM2B","DNMBP","TENM1","NPRL2",
                                                   "LDLRAD4","DCDC1","METTL20")),]
updated <- updated[names(updated)!="SCZ SNV"]
elementNROWS(updated)

save(updated, geneMap,splitSets, file = paste0(path, "updated_gene_sets.rda"))


################################################################################
# Redo Fisher's exact test
################################################################################

path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "updated_gene_sets.rda"))
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))


## Enrichment in genes differentially expressed by fraction

geneuniverse <- as.character(na.omit(unique(
                  geneMap[which(geneMap$gencodeID %in% 
                                  rownames(Ipres.down)),"Symbol"])))
splitSets <- lapply(updated, function(f) 
                             f[which(f$Symbol %in% geneuniverse), ])

inGroup = lapply(sig, function(x) as.character(na.omit(unique(x$Symbol))))
outGroup = lapply(sig, function(x) geneuniverse[!(geneuniverse %in% 
                                                    as.character(na.omit(x$Symbol)))])
inGroup = lapply(inGroup, function(x) x[which(x!="")])
outGroup = lapply(outGroup, function(x) x[which(x!="")])


enrich <- mapply(function(inG,outG) lapply(splitSets, function(x) {
  INGROUP_OVERLAP = c( sum(inG %in% x$Symbol), sum(!(inG %in% x$Symbol)))
  OUTGROUP_OVERLAP= c(sum(outG %in% x$Symbol), sum(!(outG %in% x$Symbol)))
  enrich_table = cbind(INGROUP_OVERLAP, OUTGROUP_OVERLAP)
  res = fisher.test(enrich_table)
  dat=c(res$p.value, res$estimate, enrich_table[1,], enrich_table[2,],
        res$conf.int[1], res$conf.int[2])
  names(dat) <- c("P.Value","Odds Ratio","YesSig.YesSet","NoSig.YesSet",
                  "YesSig.NoSet","NoSig.NoSet","conf.int.lower","conf.int.upper")
  return(dat)
  }), inGroup, outGroup, SIMPLIFY =F)

enrich <- lapply(enrich, data.frame)
enrich <- do.call(rbind, Map(cbind, Comparison = as.list(names(enrich)), 
                             lapply(enrich, function(x) 
                               data.frame(GeneSet = colnames(x),
                                          YesSig.YesSet = as.numeric(x["YesSig.YesSet",]),
                                          YesSig.NoSet = as.numeric(x["YesSig.NoSet",]),
                                          NoSig.YesSet = as.numeric(x["NoSig.YesSet",]),
                                          NoSig.NoSet = as.numeric(x["NoSig.NoSet",]),
                                          P.Value = as.numeric(x["P.Value",]), 
                                          OddsRatio = as.numeric(x["Odds Ratio",]),
                                          conf.int.lower = as.numeric(x["conf.int.lower",]),
                                          conf.int.upper = as.numeric(x["conf.int.upper",]),
                                          row.names = NULL))))
enrich$FDR <- p.adjust(enrich$P.Value, method = "fdr")


write.csv(enrich, quote = FALSE, 
          file = paste0(path, "RNA_localization_and_age/data/",
                        "Birnbaum_geneSet_enrichment_FractionDEGs_updated.csv"))
enrich <- read.csv(paste0(path, "RNA_localization_and_age/data/",
                          "Birnbaum_geneSet_enrichment_FractionDEGs_updated.csv"))

enrich[enrich$FDR<=0.05,c("Comparison","GeneSet","OddsRatio")]
#      Comparison           GeneSet OddsRatio
#7  both_retained           ASD.CNV  6.959218
#8  both_retained           SCZ.CNV  6.132846
#10 both_exported                ID 16.104118
#30   Ad_retained Neurodegenerative  2.696315
#32   Ad_retained         BPAD.GWAS  2.038887
#33   Ad_retained         ASD.SFARI  1.746940
#36   Ad_retained           SCZ.SNV  2.465757
#46   Ad_exported                ID  2.677269

# These are the genes in the significant sets:
paste0(enrich[enrich$FDR<=0.05,"Comparison"],".",enrich[enrich$FDR<=0.05,"GeneSet"])

genes <- list(
  both_retained.ASD.CNV = cbind(geneMap[which(geneMap$Symbol %in% 
                                                inGroup$both_retained[(inGroup$both_retained %in% 
                                                                         splitSets$ASD.CNV$Symbol)]),],
                                Comparison = "both_retained", variable = "ASD.CNV"),        
  both_retained.SCZ.CNV = 
    cbind(geneMap[which(geneMap$Symbol %in% 
                          inGroup$both_retained[(inGroup$both_retained %in% 
                                                   splitSets$"SCZ.CNV"$Symbol)]),],
          Comparison = "both_retained", variable = "SCZ.CNV"),
  both_exported.ID = 
    cbind(geneMap[which(geneMap$Symbol %in% 
                          inGroup$both_exported[(inGroup$both_exported %in% 
                                                   splitSets$ID$Symbol)]),], 
          Comparison = "both_exported", variable = "ID"), 
  Ad_retained.Neurodegenerative = 
    cbind(geneMap[which(geneMap$Symbol %in%
                          inGroup$Ad_retained[which(inGroup$Ad_retained %in%
                                                      splitSets$Neurodegenerative$Symbol)]),],
          Comparison = "Ad_retained", variable = "Neurodegenerative"),
  Ad_retained.BPAD.GWAS = 
    cbind(geneMap[which(geneMap$Symbol %in% 
                          inGroup$Ad_retained[which(inGroup$Ad_retained %in% 
                                                 splitSets$BPAD.GWAS$Symbol)]),],
          Comparison = "Ad_retained", variable = "BPAD.GWAS"),        
  Ad_retained.ASD.SFARI = 
    cbind(geneMap[which(geneMap$Symbol %in% 
                          inGroup$Ad_retained[(inGroup$Ad_retained %in% 
                                                 splitSets$ASD.SFARI$Symbol)]),],
          Comparison = "Ad_retained", variable = "ASD.SFARI"),         
  Ad_retained.SCZ.SNV = 
    cbind(geneMap[which(geneMap$Symbol %in% 
                          inGroup$Ad_retained[(inGroup$Ad_retained %in% 
                                                 splitSets$"SCZ.SNV"$Symbol)]),],
          Comparison = "Ad_retained", variable = "SCZ.SNV"),           
  Ad_exported.ID = 
    cbind(geneMap[which(geneMap$Symbol %in% 
                          inGroup$Ad_exported[which(inGroup$Ad_exported %in% 
                                                      splitSets$"ID"$Symbol)]),],
          Comparison = "Ad_exported", variable = "ID"))
elementNROWS(genes)
write.csv(do.call(rbind, genes), quote=F,
          file = paste0(path, "RNA_localization_and_age/data/",
                        "diseasegenes_enriched_updated.csv"))


## update table S2 (Fraction DEGs)

genes = read.csv(paste0(path, "RNA_localization_and_age/data/",
                        "diseasegenes_enriched_updated.csv"))

genes$variable = gsub("ASD.CNV", "Autism Spectrum Disorder (from CNVs)", genes$variable)
genes$variable = gsub("ASD.SFARI", "Autism Spectrum Disorder (from SAFARI)", genes$variable)
genes$variable = gsub("BPAD.GWAS", "Bipolar Affective Disorder (from GWAS)", genes$variable)
genes$variable = gsub("ID", "Intellectual Disability", genes$variable)
genes$variable = gsub("NDD", "Syndromal Neurodevelopmental Disorder", genes$variable)
genes$variable = gsub("Neurodegenerative", "Neurodegenerative Disorder", genes$variable)
genes$variable = gsub("SCZ.SNV", "Schizophrenia (from rare SNVs)", genes$variable)
genes$variable = gsub("SCZ.GWAS", "Schizophrenia (from GWAS)", genes$variable)
genes$variable = gsub("SCZ.CNV", "Schizophrenia (from CNVs)", genes$variable)

tab = openxlsx::read.xlsx(paste0("./Dropbox/BrainRNACompartments/paper/",
                                 "Genome_Research_revision/",
                                 "Tables_S1-S7_GENOME_2019_250217.xlsx"), 
                                 sheet = 2, startRow = 4, cols = 1:12, colNames = T)
head(tab)
tab$Disease = genes[match(tab$Gene.ID, genes$gencodeID),"variable"]

openxlsx::write.xlsx(tab, file = "./Dropbox/BrainRNACompartments/paper/newTableS2.xlsx")


enrich$GeneSet = gsub("ASD.CNV", "Autism Spectrum Disorder (from CNVs)", enrich$GeneSet)
enrich$GeneSet = gsub("ASD.SFARI", "Autism Spectrum Disorder (from SAFARI)", enrich$GeneSet)
enrich$GeneSet = gsub("BPAD.GWAS", "Bipolar Affective Disorder (from GWAS)", enrich$GeneSet)
enrich$GeneSet = gsub("ID", "Intellectual Disability", enrich$GeneSet)
enrich$GeneSet = gsub("NDD", "Syndromal Neurodevelopmental Disorder", enrich$GeneSet)
enrich$GeneSet = gsub("Neurodegenerative", "Neurodegenerative Disorder", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.SNV", "Schizophrenia (from rare SNVs)", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.GWAS", "Schizophrenia (from GWAS)", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.CNV", "Schizophrenia (from CNVs)", enrich$GeneSet)

enrich$Comparison = gsub("both_retained", "Nuclear: Both", enrich$Comparison)
enrich$Comparison = gsub("both_exported", "Cytoplasmic: Both", enrich$Comparison)
enrich$Comparison = gsub("Fet_retained", "Nuclear: Prenatal Only", enrich$Comparison)
enrich$Comparison = gsub("Ad_retained", "Nuclear: Adult Only", enrich$Comparison)
enrich$Comparison = gsub("Fet_exported", "Cytoplasmic: Prenatal Only", enrich$Comparison)
enrich$Comparison = gsub("Ad_exported", "Cytoplasmic: Adult Only", enrich$Comparison)
enrich$Comparison = gsub("ret_Ad_exp_Fet", "Nuclear in Adult; Cytoplasmic in Prenatal",
                         enrich$Comparison)
enrich$Comparison = gsub("ret_Fet_exp_Ad", "Nuclear in Prenatal; Cytoplasmic in Adult", 
                         enrich$Comparison)
enrich$Comparison = gsub("interacting", "Interaction", enrich$Comparison)

openxlsx::write.xlsx(enrich, file = paste0("./Dropbox/BrainRNACompartments/",
                                           "paper/Birnbaum_geneSet_enrichment_",
                                           "FractionDEGs_updated.xlsx"))

## Make New Odds Ratio enrichment plot

enrich = read.csv(paste0(path, "RNA_localization_and_age/data/",
                         "Birnbaum_geneSet_enrichment_FractionDEGs_updated.csv"))
enrich[enrich$FDR<=0.05,c("Comparison","GeneSet","OddsRatio")]

head(enrich)
enrich$GeneSet = gsub("ASD.CNV", "ASD\n(CNV)", enrich$GeneSet)
enrich$GeneSet = gsub("ASD.SFARI", "ASD\n(SFARI)", enrich$GeneSet)
enrich$GeneSet = gsub("BPAD.GWAS", "BPAD\n(GWAS)", enrich$GeneSet)
enrich$GeneSet = gsub("ID", "Intellectual\nDisability", enrich$GeneSet)
enrich$GeneSet = gsub("NDD", "Neuro-\ndevel.", enrich$GeneSet)
enrich$GeneSet = gsub("Neurodegenerative", "Neuro-\ndegen.", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.SNV", "SCZ\n(SNV)", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.GWAS", "SCZ\n(GWAS)", enrich$GeneSet)
enrich$GeneSet = gsub("SCZ.CNV", "SCZ\n(CNV)", enrich$GeneSet)
enrich$GeneSet = factor(enrich$GeneSet, levels = 
                          c("ASD\n(CNV)","SCZ\n(CNV)","ASD\n(SFARI)","BPAD\n(GWAS)",
                            "SCZ\n(SNV)","Neuro-\ndegen.","SCZ\n(GWAS)",
                            "Neuro-\ndevel.","Intellectual\nDisability"))

enrich$cat <- ifelse(enrich$GeneSet %in% c("ASD\n(CNV)","SCZ\n(CNV)"),
                     "Nuclear in Both", "Not Nuclear")
enrich[which(enrich$GeneSet %in% c("ASD\n(SFARI)","BPAD\n(GWAS)","SCZ\n(SNV)",
                                   "Neuro-\ndegen.")),"cat"] <- "Nuclear in Adult Only"
enrich$cat <- factor(enrich$cat, levels = c("Nuclear in Both",
                                            "Nuclear in Adult Only",
                                            "Not Nuclear"))

pdf(paste0(path,"disease/figures/birnbaum_OR_plot_updated.pdf"), 
    width = 6, height = 4)
ggplot(enrich[enrich$X %in% c(7,8,30,32,33,36),], 
       aes(GeneSet, OddsRatio, fill = cat)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=conf.int.lower, ymax=conf.int.upper), width=.2,
                position=position_dodge(.9)) +
  ylab("Odds Ratio") +  
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "gray43")) +
  xlab("") + geom_hline(yintercept=1, linetype="dotted") + 
  ggtitle("Nuclear-Enriched Genes") + 
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.title=element_blank(), 
        legend.position = c(0.7, 0.85), 
        legend.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


## Print gene sets for table

splitSets <- lapply(splitSets, function(x) as.character(na.omit(unique(x$Symbol))))
splitSets <- lapply(splitSets, function(x) x[which(x!="")])

names(splitSets) <- gsub("ASD.CNV", "Autism Spectrum Disorder (from CNVs)", names(splitSets))
names(splitSets) <- gsub("ASD.SFARI", "Autism Spectrum Disorder (from SAFARI)", names(splitSets))
names(splitSets) <- gsub("BPAD.GWAS", "Bipolar Affective Disorder (from GWAS)", names(splitSets))
names(splitSets) <- gsub("ID", "Intellectual Disability", names(splitSets))
names(splitSets) <- gsub("NDD", "Syndromal Neurodevelopmental Disorder", names(splitSets))
names(splitSets) <- gsub("Neurodegenerative", "Neurodegenerative Disorder", names(splitSets))
names(splitSets) <- gsub("SCZ.SNV", "Schizophrenia (from rare SNVs)", names(splitSets))
names(splitSets) <- gsub("SCZ.GWAS", "Schizophrenia (from GWAS)", names(splitSets))
names(splitSets) <- gsub("SCZ.CNV", "Schizophrenia (from CNVs)", names(splitSets))
splitSets

elementNROWS(sig)



path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))


ids <- do.call(rbind, Map(cbind, Group = as.list(names(sig[elementNROWS(sig)>0])), 
                          lapply(sig[elementNROWS(sig)>0], function(x) 
                            data.frame(Symbol = unique(as.character(x$Symbol))))))
write.csv(ids, quote = FALSE,
          file = paste0(path, "RNA_localization_and_age/data/sigGeneSymbols.csv"))

ids1 <- do.call(rbind, Map(cbind, Group = as.list(names(sig.1[elementNROWS(sig.1)>0])), 
                          lapply(sig.1[elementNROWS(sig.1)>0], function(x) 
                            data.frame(Symbol = unique(as.character(x$Symbol))))))
write.csv(ids, quote = FALSE,
          file = paste0(path, "RNA_localization_and_age/data/sig.1GeneSymbols.csv"))

