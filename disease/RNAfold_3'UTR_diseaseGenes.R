library(data.table)
library(ggplot2)
library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

path <- "./Dropbox/sorted_figures/github_controlled/"
load(paste0(path, "QC_section/data/rawCounts_combined_NucVSCyt_n23.rda"))
load(paste0(path, "RNA_localization_and_age/data/retained.byAge.downsampled.rda"))


## Identify the 3'UTRs

exonMap$exonID = rownames(exonMap)
dt = data.table(cbind(exonMap, as.data.frame(exonCounts.down
                                             [match(rownames(exonMap), 
                                                    rownames(exonCounts.down)),
                                               grep("polyA",
                                                    colnames(exonCounts.down))])))
dt$mean = rowMeans(as.data.frame(dt)[,grep("polyA",colnames(dt))])

txdb = loadDb(paste0(path, "intron_retention/data/SGSeq_out/",
                     "gencode.v25lift37.annotation.sqlite"))
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

load(paste0(path, "updated_gene_sets.rda"))

geneuniverse <- as.character(na.omit(unique(
  geneMap[which(geneMap$gencodeID %in% 
                  rownames(Ipres.down)),"Symbol"])))
splitSets <- lapply(updated, function(f) 
  f[which(f$Symbol %in% geneuniverse), ])

df = data.frame(set = names(splitSets), number = elementNROWS(splitSets))
df$set = gsub("ASD.CNV", "ASD\n(CNV)", df$set)
df$set = gsub("ASD.SFARI", "ASD\n(SFARI)", df$set)
df$set = gsub("BPAD.GWAS", "BPAD\n(GWAS)", df$set)
df$set = gsub("ID", "Intellectual\nDisability", df$set)
df$set = gsub("NDD", "Neuro-\ndevel.", df$set)
df$set = gsub("Neurodegenerative", "Neuro-\ndegen.", df$set)
df$set = gsub("SCZ.SNV", "SCZ\n(SNV)", df$set)
df$set = gsub("SCZ.GWAS", "SCZ\n(GWAS)", df$set)
df$set = gsub("SCZ.CNV", "SCZ\n(CNV)", df$set)
df$set <- factor(df$set, levels = c("ASD\n(CNV)","SCZ\n(CNV)",
                                    "ASD\n(SFARI)","BPAD\n(GWAS)",
                                    "SCZ\n(SNV)","Neuro-\ndegen.",
                                    "SCZ\n(GWAS)", "Neuro-\ndevel.",
                                    "Intellectual\nDisability"))
df$cat <- ifelse(df$set %in% c("ASD\n(CNV)","SCZ\n(CNV)"), 
                 "Nuclear in Both", "Not Nuclear")
df[which(df$set %in% c("ASD\n(SFARI)","BPAD\n(GWAS)","SCZ\n(SNV)",
                       "Neuro-\ndegen.")),"cat"] <- "Nuclear in Adult Only"
df$cat <- factor(df$cat, levels = c("Nuclear in Both",
                                   "Nuclear in Adult Only","Not Nuclear"))


pdf(paste0(path, "disease/figures/birnbaum_numberGenes_plot.pdf"), 
    width = 10.6, height = 3.25)
ggplot(df, aes(set, number, fill = cat)) +
  geom_col() + ylab("Number") +  
  scale_fill_manual(values=c("cornsilk4", "antiquewhite3", "gray43")) +
  xlab("") + ggtitle("Number of Genes Per Set") + 
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.title=element_blank(),
        legend.position = c(0.75, 0.8), 
        legend.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


diseaseUTR <- do.call(rbind, Map(cbind, lapply(splitSets, function(x) 
  UTR3[which(UTR3$gencodeID %in% x$gencodeID),,]), fracReg = 
    as.list(names(splitSets))))
head(diseaseUTR)


### Is there a significant difference between the free energy of the 3'UTRs 
## that are regulated in each group?

## Make the input FASTA file for RNAfold

diseaseUTR = makeGRangesFromDataFrame(diseaseUTR, keep.extra.columns = T)
names(diseaseUTR) = diseaseUTR$exonID
Seq = getSeq(Hsapiens, diseaseUTR[diseaseUTR$dominant3UTR=="Major",])
Seq = RNAStringSet(reverseComplement(Seq))

writeXStringSet(Seq, filepath = paste0(path, "RNA_localization_and_age/data/",
                                       "diseaseGenes_3UTRs.fa"))
length(unique(names(Seq)))


## Run RNAfold

# /Users/amandaprice/Desktop/ViennaRNA/bin/RNAfold /Users/amandaprice/Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/diseaseGenes_3UTRs.fa > ./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/disease_3UTRs.updated.txt 

## Plot the difference between gene sets

res <- read.table(paste0(path, "RNA_localization_and_age/data/RNAfold_out/",
                         "disease_3UTRs.updated.txt"), header = FALSE, fill = TRUE)
energy <- data.frame(exonID = res$V1[grep(">", res$V1, fixed = T)], 
                     energy = res$V2[grep("(", res$V2, fixed=T)])
energy$exonID <- gsub(">","", energy$exonID)
energy$energy <- gsub("(","",energy$energy, fixed=T)
energy$energy <- gsub(")","",energy$energy, fixed=T)
energy$energy <- as.numeric(energy$energy)

major <- data.frame(diseaseUTR)[diseaseUTR$dominant3UTR=="Major",]
major$energy <- energy[match(major$exonID, energy$exonID),"energy"]

major$fracReg = gsub("ASD.CNV", "ASD\n(CNV)", major$fracReg)
major$fracReg = gsub("ASD.SFARI", "ASD\n(SFARI)", major$fracReg)
major$fracReg = gsub("BPAD.GWAS", "BPAD\n(GWAS)", major$fracReg)
major$fracReg = gsub("ID", "Intellectual\nDisability", major$fracReg)
major$fracReg = gsub("NDD", "Neuro-\ndevel.", major$fracReg)
major$fracReg = gsub("Neurodegenerative", "Neuro-\ndegen.", major$fracReg)
major$fracReg = gsub("SCZ.SNV", "SCZ\n(SNV)", major$fracReg)
major$fracReg = gsub("SCZ.GWAS", "SCZ\n(GWAS)", major$fracReg)
major$fracReg = gsub("SCZ.CNV", "SCZ\n(CNV)", major$fracReg)
major$fracReg <- factor(major$fracReg, levels = c("ASD\n(CNV)","SCZ\n(CNV)",
                                    "ASD\n(SFARI)","BPAD\n(GWAS)",
                                    "SCZ\n(SNV)","Neuro-\ndegen.",
                                    "SCZ\n(GWAS)", "Neuro-\ndevel.",
                                    "Intellectual\nDisability"))
major$cat <- ifelse(major$fracReg %in% c("ASD\n(CNV)","SCZ\n(CNV)"), 
                 "Nuclear in Both", "Not Nuclear")
major[which(major$fracReg %in% c("ASD\n(SFARI)","BPAD\n(GWAS)","SCZ\n(SNV)",
                       "Neuro-\ndegen.")),"cat"] <- "Nuclear in Adult Only"
major$cat <- factor(major$cat, levels = c("Nuclear in Both",
                                    "Nuclear in Adult Only","Not Nuclear"))


pdf(paste0(path, "RNA_localization_and_age/figures/",
           "free_energy_diseaseGenes_3UTRs.pdf"), width = 11.5, height = 3.5)
ggplot(major, aes(x = fracReg, y = energy, fill = cat)) + 
  geom_jitter() + geom_boxplot() +
  scale_fill_manual(values=c("cornsilk4", "antiquewhite2", "white")) +
  labs(fill="") + ylab("MFE") + xlab("") +
  ggtitle("3'UTR Predicted Minimum Free Energy (MFE)") +
  theme(title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.position = "bottom")
dev.off()


res = list(Allnuc.v.rest = t.test(x = major[which(major$cat!= "Not Nuclear"),"energy"],
                                  y = major[which(major$cat=="Not Nuclear"),"energy"]),
           NucBoth.v.rest = t.test(x = major[which(major$cat=="Nuclear in Both"),"energy"],
                                   y = major[which(major$cat=="Not Nuclear"),"energy"]),
           NucAdult.v.rest = t.test(x = major[which(major$cat=="Nuclear in Adult Only"),"energy"],
                                      y = major[which(major$cat=="Not Nuclear"),"energy"]))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) 
  data.frame(Tstat = x$statistic, pval = x$p.value, 
             diseaseMean = x$estimate[1], restMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#       Comparison     Tstat      pval diseaseMean  restMean       FDR
#1   Allnuc.v.rest 0.2139335 0.8306736   -236.4536 -240.7078 0.9275768
#2  NucBoth.v.rest 0.0909608 0.9275768   -238.2492 -240.7078 0.9275768
#3 NucAdult.v.rest 0.2252771 0.8218346   -236.0376 -240.7078 0.9275768


## Check the 3'UTR length

res <- list(Allnuc.v.rest = t.test(x = major[which(major$cat!= 
                                                    "Not Nuclear"),"Length"],
                                  y = major[which(major$cat==
                                                    "Not Nuclear"),"Length"]),
           NucBoth.v.rest = t.test(x = major[which(major$cat==
                                                     "Nuclear in Both"),"Length"],
                                   y = major[which(major$cat==
                                                     "Not Nuclear"),"Length"]),
           NucAdult.v.rest = t.test(x = major[which(major$cat==
                                                      "Nuclear in Adult Only"),"Length"],
                                    y = major[which(major$cat==
                                                      "Not Nuclear"),"Length"]))

res <- do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) 
  data.frame(Tstat = x$statistic, pval = x$p.value,
             diseaseMean = x$estimate[1], restMean = x$estimate[2]))))
res$FDR <- p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#       Comparison      Tstat       pval diseaseMean restMean        FDR
#1   Allnuc.v.rest -1.2827897 0.20015616    747.7526 841.3444 0.30023424
#2  NucBoth.v.rest -2.3120119 0.02119584    642.8941 841.3444 0.06358753
#3 NucAdult.v.rest -0.8985553 0.36927363    773.4014 841.3444 0.36927363