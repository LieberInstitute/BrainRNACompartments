library(ggplot2)

load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')
enrich = read.csv("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
pgc2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]


# Assign lengths to disease genes

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"Symbol"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$Gene.Symbol) %in% geneuniverse), ] # drop genes that are not present in the test set
aej_sets_expressed = aej_sets_expressed[which(aej_sets_expressed$Gene.Set!="SCZ PGC GWAS"),]
aej_sets_expressed = rbind(aej_sets_expressed[,colnames(aej_sets_expressed) %in% c("Gene.Set","Gene.Symbol")], 
                           data.frame(Gene.Set = "PGC2", Gene.Symbol = pgc2$Symbol))
length(unique(aej_sets$Gene.Symbol)) # 1007
length(unique(aej_sets_expressed$Gene.Symbol)) # 1240
len = cbind(aej_sets_expressed, geneMap[match(aej_sets_expressed$Gene.Symbol, geneMap$Symbol),colnames(geneMap) %in% c("Length", "gencodeID")])
len$cat = ifelse(len$Gene.Set %in% c("ASD CNV","ASD DATABASE","SCZ CNV"), "Nuclear in Both", "Not Nuclear")
len[which(len$Gene.Set %in% c("BPAD GWAS","SCZ SNV","PGC2")),"cat"] = "Nuclear in Adult Only"

head(len)

# Compare retained gene set lengths to others

enrich[enrich$FDR<=0.05 & enrich$Comparison %in% c("both_retained","Ad_retained"),"GeneSet"]
res = list(allNuc.vs.disease = t.test(x = len[len$cat!="Not Nuclear","Length"], y = len[len$cat=="Not Nuclear","Length"]),
           bothNuc.vs.disease = t.test(x = len[len$cat=="Nuclear in Both","Length"], y = len[len$cat!="Nuclear in Both","Length"]),
           adNuc.vs.disease = t.test(x = len[len$cat=="Nuclear in Adult Only","Length"], y = len[len$cat!="Nuclear in Adult Only","Length"]),
           
           allNuc.vs.disease.noTTN = t.test(x = len[len$cat!="Not Nuclear" & len$Gene.Symbol!="TTN","Length"], 
                                            y = len[len$cat=="Not Nuclear" & len$Gene.Symbol!="TTN","Length"]),
           bothNuc.vs.disease.noTTN = t.test(x = len[len$cat=="Nuclear in Both" & len$Gene.Symbol!="TTN","Length"], 
                                             y = len[len$cat!="Nuclear in Both" & len$Gene.Symbol!="TTN","Length"]),
           ad.Nuc.vs.disease.noTTN = t.test(x = len[len$cat=="Nuclear in Adult Only" & len$Gene.Symbol!="TTN","Length"], 
                                            y = len[len$cat!="Nuclear in Adult Only" & len$Gene.Symbol!="TTN","Length"]),
           
           allNuc.vs.allGenes = t.test(x = len[len$cat!="Not Nuclear","Length"], 
                                       y = geneMap[!geneMap$Symbol %in% len[len$cat!="Not Nuclear","Gene.Symbol"],"Length"]),
           bothNuc.vs.allGenes = t.test(x = len[len$cat=="Nuclear in Both","Length"], 
                                        y = geneMap[!geneMap$Symbol %in% len[len$cat=="Nuclear in Both","Gene.Symbol"],"Length"]),
           adNuc.vs.allGenes = t.test(x = len[len$cat=="Nuclear in Adult Only","Length"], 
                                      y = geneMap[!geneMap$Symbol %in% len[len$cat=="Nuclear in Adult Only","Gene.Symbol"],"Length"]),
           
           adNuc.vs.allGenes.noTTN = t.test(x = len[len$cat!="Not Nuclear" & len$Gene.Symbol!="TTN","Length"], 
                                            y = geneMap[!geneMap$Symbol %in% len[len$cat!="Not Nuclear" & len$Gene.Symbol!="TTN","Gene.Symbol"],"Length"]),
           bothNuc.vs.allGenes.noTTN = t.test(x = len[len$cat=="Nuclear in Both" & len$Gene.Symbol!="TTN","Length"], 
                                              y = geneMap[!geneMap$Symbol %in% len[len$cat=="Nuclear in Both" & len$Gene.Symbol!="TTN","Gene.Symbol"],"Length"]),
           adNuc.vs.allGenes.noTTN = t.test(x = len[len$cat=="Nuclear in Adult Only" & len$Gene.Symbol!="TTN","Length"], 
                                            y = geneMap[!geneMap$Symbol %in% len[len$cat=="Nuclear in Adult Only" & len$Gene.Symbol!="TTN","Gene.Symbol"],"Length"]))


res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            diseaseMean = x$estimate[1], restMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
write.csv(res, file = "./Dropbox/sorted_figures/github_controlled/disease/figures/geneSet_length_ttest.csv")

res[res$FDR<=0.05,]
#                  Comparison     Tstat         pval diseaseMean restMean          FDR
#1          allNuc.vs.disease -5.635511 4.744633e-08    5879.202 8368.462 5.693559e-08
#2         bothNuc.vs.disease  4.220835 2.754719e-05    7196.440 5749.614 2.754719e-05
#3           adNuc.vs.disease -8.044180 2.215194e-15    5200.456 7535.440 3.322792e-15
#4    allNuc.vs.disease.noTTN -5.682031 3.961024e-08    5897.107 8368.462 5.281365e-08
#5   bothNuc.vs.disease.noTTN  4.281310 2.045875e-05    6967.383 5869.879 2.231863e-05
#6    ad.Nuc.vs.disease.noTTN -8.434380 8.470616e-17    5330.012 7373.227 1.452106e-16
#7         allNuc.vs.allGenes 18.453901 1.168248e-68    5879.202 3201.681 4.672990e-68
#8        bothNuc.vs.allGenes 15.963149 1.730060e-46    7196.440 2215.128 5.190181e-46
#9          adNuc.vs.allGenes 13.576974 1.272268e-38    5200.456 3249.888 2.544537e-38
#10   adNuc.vs.allGenes.noTTN 21.882516 2.984813e-92    5897.107 3205.296 3.581776e-91
#11 bothNuc.vs.allGenes.noTTN 22.369224 6.275669e-77    6967.383 2217.079 3.765401e-76
#12   adNuc.vs.allGenes.noTTN 14.262444 5.445507e-42    5330.012 3249.888 1.306922e-41


max(len[len$Gene.Set=="ASD DATABASE","Length"])
len[len$Length==118976,] # titin
max(len[len$Gene.Set=="ASD DATABASE","Length"])-
  max(len[len$Gene.Set=="ASD DATABASE" & len$Length!=118976,"Length"]) # 88035
max(len$Length)-max(len[len$Length!=118976,"Length"]) # 73327

head(len)
len$Gene.Set = gsub("ASD CNV", "ASD\n(CNV)", len$Gene.Set)
len$Gene.Set = gsub("ASD DATABASE", "ASD\n(Database)", len$Gene.Set)
len$Gene.Set = gsub("BPAD GWAS", "BPAD\n(GWAS)", len$Gene.Set)
len$Gene.Set = gsub("ID", "Intellectual\nDisability", len$Gene.Set)
len$Gene.Set = gsub("NDD", "Neuro-\ndevel.", len$Gene.Set)
len$Gene.Set = gsub("Neurodegenerative", "Neuro-\ndegen.", len$Gene.Set)
len$Gene.Set = gsub("SCZ Meta-analysis", "SCZ\n(Meta\nanalysis)", len$Gene.Set)
len$Gene.Set = gsub("SCZ SNV", "SCZ\n(SNV)", len$Gene.Set)
len$Gene.Set = gsub("PGC2", "SCZ\n(PGC2)", len$Gene.Set)
len$Gene.Set = gsub("SCZ CNV", "SCZ\n(CNV)", len$Gene.Set)
len$Gene.Set = factor(len$Gene.Set, levels = 
                          c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)","BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)","SCZ\n(Meta\nanalysis)",
                            "Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))
len$cat = ifelse(len$Gene.Set %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)"), "Nuclear in Both", "Not Nuclear")
len[which(len$Gene.Set %in% c("BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")),"cat"] = "Nuclear in Adult Only"
len$cat = factor(len$cat, levels = c("Nuclear in Both","Nuclear in Adult Only","Not Nuclear"))
len[len$cat=="Nuclear in Adult Only",]


pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/disease_geneSet_length.pdf", width = 13, height = 3)
ggplot(len, aes(x=Gene.Set, y=Length/1000, fill = cat)) + geom_boxplot() + 
  ylab("Gene Length (Kb)") + xlab("") +
  ylim(0,46) + scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  ggtitle("Gene Lengths by Disease Gene Set") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
ggplot(len, aes(x=Gene.Set, y=Length/1000, fill = cat)) + geom_boxplot() + 
  ylab("Gene Length (Kb)") + xlab("") + scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "white")) +
  ggtitle("Gene Lengths by Disease Gene Set") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
dev.off()
