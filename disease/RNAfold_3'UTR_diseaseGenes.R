library(data.table)
library(ggplot2)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


## Identify the 3'UTRs

exonMap$exonID = rownames(exonMap)
dt = data.table(cbind(exonMap, as.data.frame(exonCounts.down[match(rownames(exonMap), rownames(exonCounts.down)),
                                                             grep("polyA",colnames(exonCounts.down))])))
dt$mean = rowMeans(as.data.frame(dt)[,grep("polyA",colnames(dt))])

txdb = loadDb("./Dropbox/sorted_figures/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
UTR3 = threeUTRsByTranscript(txdb, use.names=T)
UTR3 = unlist(UTR3)

exonGR = makeGRangesFromDataFrame(dt, keep = T)
ov = findOverlaps(UTR3, exonGR, type="equal")
UTR3 = exonGR[unique(subjectHits(ov))]
UTR3 = data.table(as.data.frame(UTR3))


## Annotate the "Major" and "Minor" isoform 3'UTR by labeling the highest expressed 3'UTR per gene

dtMax = UTR3[UTR3[, .I[mean == max(mean)], by=gencodeID]$V1]
UTR3$dominant3UTR = ifelse(UTR3$exonID %in% dtMax$exonID, "Major", "Minor")


## RNA secondary structure of disease genes

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

AEJmap = geneMap[which(geneMap$Symbol %in% as.character(aej_sets$Gene.Symbol)),]
AEJmap = cbind(AEJmap, Gene.Set = aej_sets[match(AEJmap$Symbol, aej_sets$Gene.Symbol),"Gene.Set"])

pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)

AEJmap = c(split(AEJmap, AEJmap$Gene.Set), list(PGC2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]))
elementNROWS(AEJmap)
# ASD CNV      ASD DATABASE         BPAD GWAS                ID               NDD Neurodegenerative           SCZ CNV SCZ Meta-analysis 
#      87               211                84                78                30                41                99                33 
# SCZ PGC GWAS           SCZ SNV              PGC2 
#          100               196               635 

df = data.frame(set = names(AEJmap), number = elementNROWS(AEJmap))
df$set = gsub("ASD CNV", "ASD\n(CNV)", df$set)
df$set = gsub("ASD DATABASE", "ASD\n(Database)", df$set)
df$set = gsub("BPAD GWAS", "BPAD\n(GWAS)", df$set)
df$set = gsub("ID", "Intellectual\nDisability", df$set)
df$set = gsub("NDD", "Neuro-\ndevel.", df$set)
df$set = gsub("Neurodegenerative", "Neuro-\ndegen.", df$set)
df$set = gsub("SCZ Meta-analysis", "SCZ\n(Meta\nanalysis)", df$set)
df$set = gsub("SCZ SNV", "SCZ\n(SNV)", df$set)
df$set = gsub("PGC2", "SCZ\n(PGC2)", df$set)
df$set = gsub("SCZ CNV", "SCZ\n(CNV)", df$set)
df = df[which(df$set!="SCZ PGC GWAS"),]
df$set = factor(df$set, levels = c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)",
                                   "BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)",
                                   "SCZ\n(Meta\nanalysis)","Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))
df$cat = ifelse(df$set %in% c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)"), "Nuclear in Both", "Not Nuclear")
df[which(df$set %in% c("BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)")),"cat"] = "Nuclear in Adult Only"
df$cat = factor(df$cat, levels = c("Nuclear in Both","Nuclear in Adult Only","Not Nuclear"))


pdf("./Dropbox/sorted_figures/github_controlled/disease/figures/birnbaum_numberGenes_plot.pdf", width = 10.6, height = 3.25)
ggplot(df, aes(set, number, fill = cat)) +
  geom_col() + ylab("Number") +  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "gray43")) +
  xlab("") + ggtitle("Number of Genes Per Set") + 
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.title=element_blank(),
        legend.position = c(0.75, 0.8), legend.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


diseaseUTR = do.call(rbind, Map(cbind, lapply(AEJmap[which(names(AEJmap)!="SCZ PGC GWAS")], function(x) UTR3[which(UTR3$gencodeID %in% x$gencodeID),,]), 
                                fracReg = as.list(names(AEJmap)[which(names(AEJmap)!="SCZ PGC GWAS")])))
head(diseaseUTR)


### Is there a significant difference between the free energy of the 3'UTRs that are regulated in each group?

## Make the input FASTA file for RNAfold

diseaseUTR = makeGRangesFromDataFrame(diseaseUTR, keep.extra.columns = T)
names(diseaseUTR) = diseaseUTR$exonID
Seq = getSeq(Hsapiens, diseaseUTR[diseaseUTR$dominant3UTR=="Major",])
Seq = RNAStringSet(reverseComplement(Seq))

writeXStringSet(Seq, filepath = "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/diseaseGenes_3UTRs.fa")
length(unique(names(Seq)))


## Run RNAfold

# /Users/amandaprice/Desktop/ViennaRNA/bin/RNAfold /Users/amandaprice/Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/diseaseGenes_3UTRs.fa > ./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/disease_3UTRs.out.txt 

## Plot the difference between gene sets

res = read.table("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/disease_3UTRs.out.txt", header = FALSE, fill = TRUE)
energy = data.frame(exonID = res$V1[grep(">", res$V1, fixed = T)], energy = res$V2[grep("(", res$V2, fixed=T)])
energy$exonID = gsub(">","", energy$exonID)
energy$energy = gsub("(","",energy$energy, fixed=T)
energy$energy = gsub(")","",energy$energy, fixed=T)
energy$energy = as.numeric(energy$energy)

major = data.frame(diseaseUTR)[diseaseUTR$dominant3UTR=="Major",]
major$energy = energy[match(major$exonID, energy$exonID),"energy"]

major$cat = ifelse(major$fracReg %in% c("ASD CNV","ASD DATABASE","SCZ CNV"), "Nuclear in Both", "Not Nuclear")
major[major$fracReg %in% c("BPAD GWAS","SCZ SNV","PGC2"),"cat"] = "Nuclear in Adult Only"
major$cat = factor(major$cat, levels = c("Nuclear in Both","Nuclear in Adult Only","Not Nuclear"))
major$fracReg = gsub("NDD", "Neuro-\ndevel.", major$fracReg)     
major$fracReg = gsub("Neurodegenerative", "Neuro-\ndegen.", major$fracReg) 
major$fracReg = gsub("SCZ Meta-analysis" ,"SCZ\n(Meta\nanalysis)", major$fracReg)
major$fracReg = gsub("PGC2","SCZ\n(PGC2)", major$fracReg)     
major$fracReg = gsub("SCZ CNV", "SCZ\n(CNV)", major$fracReg)         
major$fracReg = gsub("SCZ SNV", "SCZ\n(SNV)", major$fracReg)
major$fracReg = gsub("ASD CNV", "ASD\n(CNV)", major$fracReg)
major$fracReg = gsub("ASD DATABASE", "ASD\n(Database)", major$fracReg)
major$fracReg = gsub("BPAD GWAS","BPAD\n(GWAS)", major$fracReg)     
major$fracReg = gsub("ID", "Intellectual\nDisability", major$fracReg)
major$fracReg = factor(major$fracReg, levels = c("ASD\n(CNV)","ASD\n(Database)","SCZ\n(CNV)",
                                                 "BPAD\n(GWAS)","SCZ\n(SNV)","SCZ\n(PGC2)",
                                                 "SCZ\n(Meta\nanalysis)","Neuro-\ndevel.","Intellectual\nDisability","Neuro-\ndegen."))


pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/free_energy_diseaseGenes_3UTRs.pdf",width=11.5,height=3.5)
ggplot(major, aes(x = fracReg, y = energy, fill = cat)) + geom_jitter() + geom_boxplot() +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite2", "white")) +
  labs(fill="") + ylab("MFE") + xlab("") +
  ggtitle("3'UTR Predicted Minimum Free Energy (MFE)") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "bottom")
dev.off()


res = list(Allnuc.v.rest = t.test(x = major[which(major$cat!= "Not Nuclear"),"energy"],
                                      y = major[which(major$cat=="Not Nuclear"),"energy"]),
           NucBoth.v.rest = t.test(x = major[which(major$cat== "Nuclear in Both"),"energy"],
                                   y = major[which(major$cat!="Nuclear in Both"),"energy"]),
           ASD.CNV.v.disease = t.test(x = major[which(major$cat== "Nuclear in Adult Only"),"energy"],
                                      y = major[which(major$cat!="Nuclear in Adult Only"),"energy"]))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            diseaseMean = x$estimate[1], restMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#         Comparison      Tstat      pval diseaseMean  restMean       FDR
#      Allnuc.v.rest  0.8840318 0.3779661   -219.7034 -249.0865 0.5584784
#     NucBoth.v.rest -0.5855785 0.5584784   -235.0764 -218.5121 0.5584784
#  ASD.CNV.v.disease  1.4163411 0.1572730   -202.8592 -239.5380 0.4718189


## Check the 3'UTR length

res = list(Allnuc.v.rest = t.test(x = major[which(major$cat!= "Not Nuclear"),"Length"],
                                  y = major[which(major$cat=="Not Nuclear"),"Length"]),
           NucBoth.v.rest = t.test(x = major[which(major$cat== "Nuclear in Both"),"Length"],
                                   y = major[which(major$cat!="Nuclear in Both"),"Length"]),
           ASD.CNV.v.disease = t.test(x = major[which(major$cat== "Nuclear in Adult Only"),"Length"],
                                      y = major[which(major$cat!="Nuclear in Adult Only"),"Length"]))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            diseaseMean = x$estimate[1], restMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
#res
#         Comparison      Tstat       pval diseaseMean restMean       FDR
#      Allnuc.v.rest -1.8329255 0.06898581    711.8561 926.0093 0.1034787
#     NucBoth.v.rest  0.4934879 0.62194615    774.1878 731.4490 0.6219461
#  ASD.CNV.v.disease -1.9963148 0.04633828    669.1198 822.5357 0.1034787

