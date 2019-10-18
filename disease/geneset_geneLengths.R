library(ggplot2)

path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "updated_gene_sets.rda"))
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))

enrich <- read.csv(paste0(path, "RNA_localization_and_age/data/",
                          "Birnbaum_geneSet_enrichment_FractionDEGs_updated.csv"))


# Assign lengths to disease genes

geneuniverse <- as.character(na.omit(unique(
        geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"gencodeID"])))
splitSets <- lapply(updated, function(f) f[which(f$gencodeID %in% geneuniverse),])
splitSets <- do.call(rbind, Map(cbind, Gene.Set = as.list(names(splitSets)), splitSets))

len <- splitSets[,colnames(splitSets) %in% c("Gene.Set","Length","gencodeID")]

len$Gene.Set = gsub("ASD.CNV", "ASD\n(CNV)", len$Gene.Set)
len$Gene.Set = gsub("ASD.SFARI", "ASD\n(SFARI)", len$Gene.Set)
len$Gene.Set = gsub("BPAD.GWAS", "BPAD\n(GWAS)", len$Gene.Set)
len$Gene.Set = gsub("ID", "Intellectual\nDisability", len$Gene.Set)
len$Gene.Set = gsub("NDD", "Neuro-\ndevel.", len$Gene.Set)
len$Gene.Set = gsub("Neurodegenerative", "Neuro-\ndegen.", len$Gene.Set)
len$Gene.Set = gsub("SCZ.SNV", "SCZ\n(SNV)", len$Gene.Set)
len$Gene.Set = gsub("SCZ.GWAS", "SCZ\n(GWAS)", len$Gene.Set)
len$Gene.Set = gsub("SCZ.CNV", "SCZ\n(CNV)", len$Gene.Set)
len$Gene.Set = factor(len$Gene.Set, levels = 
                          c("ASD\n(CNV)","ASD\n(SFARI)","BPAD\n(GWAS)",
                            "SCZ\n(SNV)","SCZ\n(CNV)","Neuro-\ndegen.",
                            "SCZ\n(GWAS)", "Neuro-\ndevel.",
                            "Intellectual\nDisability"))

len$cat <- ifelse(len$Gene.Set %in% c("ASD\n(CNV)","SCZ\n(CNV)"),
                  "Nuclear in Both", "Not Nuclear")
len[which(len$Gene.Set %in% c("ASD\n(SFARI)","BPAD\n(GWAS)","SCZ\n(SNV)",
                                   "Neuro-\ndegen.")),"cat"] <- "Nuclear in Adult Only"
len$cat <- factor(len$cat, levels = c("Nuclear in Both",
                                      "Nuclear in Adult Only",
                                      "Not Nuclear"))
head(len)
unique(len$cat)

# Compare retained gene set lengths to others

spl <- split(len, len$Gene.Set)
lapply(spl, function(x) c(range = range(x$Length), mean = mean(x$Length)) )
res <- list(allNuc.vs.disease = t.test(x = len[len$cat %in% c("Nuclear in Both",
                                                              "Nuclear in Adult Only"),
                                               "Length"], 
                                       y = len[len$cat %in% c("Both Fractions",
                                                              "Not Nuclear"),"Length"]),
            bothNuc.vs.disease = t.test(x = len[len$cat=="Nuclear in Both","Length"], 
                                        y = len[len$cat %in% c("Both Fractions",
                                                               "Not Nuclear"),"Length"]),
            adNuc.vs.disease = t.test(x = len[len$cat=="Nuclear in Adult Only","Length"], 
                                      y = len[len$cat %in% c("Both Fractions",
                                                             "Not Nuclear"),"Length"]),
           
            allNuc.vs.disease.noTTN = t.test(x = len[len$cat %in% c("Nuclear in Both",
                                                                    "Nuclear in Adult Only") & 
                                                       len$Length!=118976,"Length"], 
                                             y = len[len$cat %in% c("Both Fractions",
                                                                    "Not Nuclear") & 
                                                       len$Length!=118976,"Length"]),
            bothNuc.vs.disease.noTTN = t.test(x = len[len$cat=="Nuclear in Both" & 
                                                        len$Length!=118976,"Length"], 
                                              y = len[len$cat %in% c("Both Fractions",
                                                                     "Not Nuclear") & 
                                                        len$Length!=118976,"Length"]),
            ad.Nuc.vs.disease.noTTN = t.test(x = len[len$cat=="Nuclear in Adult Only" & 
                                                       len$Length!=118976,"Length"], 
                                             y = len[len$cat %in% c("Both Fractions",
                                                                    "Not Nuclear") & 
                                                       len$Length!=118976,"Length"]),
           
            allNuc.vs.allGenes = t.test(x = len[len$cat %in% c("Nuclear in Both",
                                                               "Nuclear in Adult Only"),
                                                "Length"], 
                                        y = geneMap[!geneMap$gencodeID %in% 
                                                      len[len$cat %in% c("Nuclear in Both",
                                                                         "Nuclear in Adult Only"),
                                                          "gencodeID"],"Length"]),
            bothNuc.vs.allGenes = t.test(x = len[len$cat=="Nuclear in Both","Length"], 
                                         y = geneMap[!geneMap$gencodeID %in% 
                                                       len[len$cat=="Nuclear in Both",
                                                           "gencodeID"],"Length"]),
            adNuc.vs.allGenes = t.test(x = len[len$cat=="Nuclear in Adult Only","Length"], 
                                       y = geneMap[!geneMap$gencodeID %in% 
                                                     len[len$cat=="Nuclear in Adult Only",
                                                         "gencodeID"],"Length"]),
           
            adNuc.vs.allGenes.noTTN = t.test(x = len[len$cat %in% c("Nuclear in Both",
                                                                    "Nuclear in Adult Only") & 
                                                       len$Length!=118976,"Length"], 
                                             y = geneMap[!geneMap$gencodeID %in% 
                                                          len[len$cat %in% c("Nuclear in Both",
                                                                             "Nuclear in Adult Only") & 
                                                                len$Length!=118976,"gencodeID"],"Length"]),
            bothNuc.vs.allGenes.noTTN = t.test(x = len[len$cat=="Nuclear in Both" & 
                                                        len$Length!=118976,"Length"], 
                                               y = geneMap[!geneMap$gencodeID %in% 
                                                            len[len$cat=="Nuclear in Both" & 
                                                                  len$Length!=118976,"gencodeID"],
                                                          "Length"]),
            adNuc.vs.allGenes.noTTN = t.test(x = len[len$cat=="Nuclear in Adult Only" & 
                                                      len$gencodeID!=118976,"Length"], 
                                             y = geneMap[!geneMap$gencodeID %in% 
                                                          len[len$cat=="Nuclear in Adult Only" & 
                                                                len$Length!=118976,"gencodeID"],
                                                        "Length"]))


res = do.call(rbind, Map(cbind, Comparison = names(res), 
                         lapply(res, function(x) data.frame(
                           Tstat = x$statistic, pval = x$p.value,
                           diseaseMean = x$estimate[1], 
                           restMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
write.csv(res, file = paste0(path, "disease/figures/geneSet_length_ttest.csv"))

res
#                  Comparison     Tstat          pval diseaseMean restMean           FDR
#1          allNuc.vs.disease 16.084779  7.705842e-56    6611.735 3723.343  1.321002e-55
#2         bothNuc.vs.disease  7.222248  2.468684e-12    5694.970 3723.343  2.468684e-12
#3           adNuc.vs.disease 15.877231  2.790093e-54    6783.390 3723.343  4.185139e-54
#4    allNuc.vs.disease.noTTN 16.884378  5.696658e-61    6546.064 3723.343  1.139332e-60
#5   bothNuc.vs.disease.noTTN  7.222248  2.468684e-12    5694.970 3723.343  2.468684e-12
#6    ad.Nuc.vs.disease.noTTN 16.908745  5.275645e-61    6705.533 3723.343  1.139332e-60
#7         allNuc.vs.allGenes 33.829836 2.689239e-193    6611.735 2133.255 1.613543e-192
#8        bothNuc.vs.allGenes 14.131308  2.466567e-34    5694.970 2238.801  2.959881e-34
#9          adNuc.vs.allGenes 30.978035 2.035167e-162    6783.390 2144.460 8.007434e-162
#10   adNuc.vs.allGenes.noTTN 38.345373 3.668767e-234    6546.064 2135.246 4.402521e-233
#11 bothNuc.vs.allGenes.noTTN 14.131308  2.466567e-34    5694.970 2238.801  2.959881e-34
#12   adNuc.vs.allGenes.noTTN 30.962065 2.669145e-162    6783.390 2146.444 8.007434e-162

max(len[len$Gene.Set=="SCZ\n(SNV)","Length"])
len[len$Length==118976,] # titin
max(len[len$Gene.Set=="SCZ\n(SNV)","Length"])-
  max(len[len$Gene.Set=="SCZ\n(SNV)" & len$Length!=118976,"Length"]) # 75874
max(len$Length)-max(len[len$Length!=118976,"Length"]) # 75874

head(len)


pdf(paste0(path, "disease/figures/disease_geneSet_length.pdf"), 
    width = 13, height = 3)
ggplot(len, aes(x=Gene.Set, y=Length/1000, fill = cat)) + geom_boxplot() + 
  ylab("Gene Length (Kb)") + xlab("") +
  ylim(0,46) + 
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "grey", "white")) +
  ggtitle("Gene Lengths by Disease Gene Set") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
ggplot(len, aes(x=Gene.Set, y=Length/1000, fill = cat)) + geom_boxplot() + 
  ylab("Gene Length (Kb)") + xlab("") + 
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "grey", "white")) +
  ggtitle("Gene Lengths by Disease Gene Set") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
dev.off()
