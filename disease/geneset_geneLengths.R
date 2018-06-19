library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/new/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')
enrich = read.csv("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/Birnbaum_geneSet_enrichment_FractionDEGs.csv")


# Assign lengths to disease genes

geneuniverse = as.character(na.omit(unique(geneMap[which(geneMap$gencodeID %in% rownames(Ipres.down)),"Symbol"]))) # all expressed genes
aej_sets_expressed = aej_sets[which(as.character(aej_sets$Gene.Symbol) %in% geneuniverse), ] # drop genes that are not present in the test set
length(unique(aej_sets$Gene.Symbol)) # 1007
length(unique(aej_sets_expressed$Gene.Symbol)) # 947
len = cbind(aej_sets_expressed, geneMap[match(aej_sets_expressed$Gene.Symbol, geneMap$Symbol),colnames(geneMap) %in% c("Length", "gencodeID")])


# Compare retained gene set lengths to others

enrich[enrich$FDR<=0.05 & enrich$Comparison=="both_retained","GeneSet"]
res = list(Allsets.v.disease = t.test(x = len[len$Gene.Set %in% c("ASD CNV","ASD DATABASE","SCZ CNV","BPAD GWAS"),"Length"],
                         y = len[!len$Gene.Set %in% c("ASD CNV","ASD DATABASE","SCZ CNV"),"Length"]),
           BPAD.v.disease = t.test(x = len[len$Gene.Set == "BPAD GWAS","Length"],
                                   y = len[len$Gene.Set != "BPAD GWAS","Length"]),
           ASD.CNV.v.disease = t.test(x = len[len$Gene.Set == "ASD CNV","Length"],
                             y = len[len$Gene.Set != "ASD CNV","Length"]),
            ASD.DATABASE.v.disease = t.test(x = len[len$Gene.Set == "ASD DATABASE","Length"],
                                  y = len[len$Gene.Set != "ASD DATABASE","Length"]),
            SCZ.CNV.v.disease = t.test(x = len[len$Gene.Set == "SCZ CNV","Length"],
                             y = len[len$Gene.Set != "SCZ CNV","Length"]),
           Allsets.v.all = t.test(x = len[len$Gene.Set %in% c("ASD CNV","ASD DATABASE","SCZ CNV","BPAD GWAS"),"Length"],
                                      y = geneMap[!geneMap$gencodeID %in% len[len$Gene.Set %in% c("ASD CNV","ASD DATABASE","SCZ CNV"),"gencodeID"],"Length"]),
           BPAD.v.all = t.test(x = len[len$Gene.Set == "BPAD GWAS","Length"],
                                  y = geneMap[!geneMap$gencodeID %in% len[len$Gene.Set =="BPAD GWAS","gencodeID"],"Length"]),
           ASD.CNV.v.all = t.test(x = len[len$Gene.Set == "ASD CNV","Length"],
                                      y = geneMap[!geneMap$gencodeID %in% len[len$Gene.Set =="ASD CNV","gencodeID"],"Length"]),
           ASD.DATABASE.v.all = t.test(x = len[len$Gene.Set == "ASD DATABASE","Length"],
                                           y = geneMap[!geneMap$gencodeID %in% len[len$Gene.Set == "ASD DATABASE","gencodeID"],"Length"]),
           SCZ.CNV.v.all = t.test(x = len[len$Gene.Set == "SCZ CNV","Length"],
                                      y = geneMap[!geneMap$gencodeID %in% len[len$Gene.Set == "SCZ CNV","gencodeID"],"Length"]))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            diseaseMean = x$estimate[1], restMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
write.csv(res, file = "./Dropbox/sorted_figures/new/github_controlled/disease/figures/geneSet_length_ttest.csv")

res[res$FDR<=0.05,]
#2          BPAD.v.disease -2.663323 8.566176e-03    6070.684 7529.644 9.517973e-03
#3       ASD.CNV.v.disease -4.042939 6.845300e-05    5965.167 7607.110 9.779000e-05
#4  ASD.DATABASE.v.disease  2.962179 3.316670e-03    8763.711 6998.879 4.145838e-03
#5       SCZ.CNV.v.disease -5.172468 6.817025e-07    5398.173 7567.640 1.136171e-06
#6           Allsets.v.all 17.580704 3.038193e-56    6979.091 2215.128 3.038193e-55
#7              BPAD.v.all  7.489405 1.469118e-11    6070.684 2242.915 2.938236e-11
#8           ASD.CNV.v.all 10.650744 3.105944e-20    5965.167 2240.705 1.035315e-19
#9      ASD.DATABASE.v.all 11.484335 1.647416e-24    8763.711 2224.845 8.237081e-24
#10          SCZ.CNV.v.all  8.534498 1.905627e-13    5398.173 2245.220 4.764069e-13

max(len[len$Gene.Set=="ASD DATABASE","Length"])
len[len$Length==118976,] # titin
max(len[len$Gene.Set=="ASD DATABASE","Length"])-
  max(len[len$Gene.Set=="ASD DATABASE" & len$Length!=118976,"Length"]) # 88035
max(len$Length)-max(len[len$Length!=118976,"Length"]) # 73327

head(len)
len$cat = ifelse(len$Gene.Set %in% c("ASD CNV","ASD DATABASE","SCZ CNV"), "Nuclear in Both", "Not Nuclear")
len[len$Gene.Set=="BPAD GWAS","cat"] = "Nuclear in Adult Only"
len$cat = factor(len$cat, levels = c("Nuclear in Both","Nuclear in Adult Only","Not Nuclear"))
len$Gene.Set = gsub("NDD", "Neuro-\ndevel.", len$Gene.Set)     
len$Gene.Set = gsub("Neurodegenerative", "Neuro-\ndegen.", len$Gene.Set) 
len$Gene.Set = gsub("SCZ Meta-analysis" ,"SCZ\n(Meta\nanalysis)", len$Gene.Set)
len$Gene.Set = gsub("SCZ PGC GWAS","SCZ\n(GWAS)", len$Gene.Set)     
len$Gene.Set = gsub("SCZ CNV", "SCZ\n(CNV)", len$Gene.Set)         
len$Gene.Set = gsub("SCZ SNV", "SCZ\n(SNV)", len$Gene.Set)
len$Gene.Set = gsub("ASD CNV", "ASD\n(CNV)", len$Gene.Set)
len$Gene.Set = gsub("ASD DATABASE", "ASD\n(Database)", len$Gene.Set)
len$Gene.Set = gsub("BPAD GWAS","BPAD\n(GWAS)", len$Gene.Set)     
len$Gene.Set = gsub("ID", "Intellectual\nDisability", len$Gene.Set)
len$Gene.Set = factor(len$Gene.Set, levels = 
                        c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)","BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(Meta\nanalysis)",
                          "SCZ\n(GWAS)","Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))

pdf("./Dropbox/sorted_figures/new/github_controlled/disease/figures/disease_geneSet_length.pdf", width = 13.5, height = 3.5)
ggplot(len, aes(x=Gene.Set, y=Length/1000, fill = cat)) + geom_boxplot() + 
  ylab("Gene Length (Kb)") + xlab("") +
  ylim(0,46) + scale_fill_manual(values=c("cornsilk4", "antiquewhite2", "white")) +
  ggtitle("Gene Lengths by Disease Gene Set") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
ggplot(len, aes(x=Gene.Set, y=Length/1000, fill = cat)) + geom_boxplot() + 
  ylab("Gene Length (Kb)") + xlab("") + scale_fill_manual(values=c("cornsilk4", "antiquewhite2", "white")) +
  ggtitle("Gene Lengths by Disease Gene Set") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
dev.off()
