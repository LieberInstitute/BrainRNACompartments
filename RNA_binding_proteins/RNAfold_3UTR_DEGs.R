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

sigUTR = do.call(rbind, Map(cbind, lapply(sig[-which(names(sig) %in% c("ret_Ad_exp_Fet", "ret_Fet_exp_Ad", "interacting"))], 
                                          function(x) UTR3[which(UTR3$gencodeID %in% x$geneID),,]), 
                            fracReg = as.list(names(sig[-which(names(sig) %in% c("ret_Ad_exp_Fet", "ret_Fet_exp_Ad", "interacting"))]))))
head(sigUTR)


### Is the nuclear one the longer one available in the transcript families?

major = as.data.frame(sigUTR)[sigUTR$dominant3UTR=="Major",]
major$Dir = major$group = NA
major[grep("retained", major$fracReg),"Dir"] = "Nuclear"
major[grep("exported", major$fracReg),"Dir"] = "Cytoplasmic"
major[grep("Ad_", major$fracReg),"group"] = "In Adult"
major[grep("Fet_", major$fracReg),"group"] = "In Prenatal"
major[grep("both", major$fracReg),"group"] = "In Both"
major$group = factor(major$group, levels = c("In Both", "In Prenatal", "In Adult"))
major$Dir = factor(major$Dir, levels = c("Cytoplasmic", "Nuclear"))

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/length_fracDEGs_3UTRs.pdf",width=4.5,height=3)
ggplot(major, aes(x = Dir, y = Length/1000, fill = Dir)) + geom_jitter() + geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  facet_grid(. ~ group) + labs(fill="") + ylab("Kb") + xlab("") +
  ggtitle("3'UTR Length") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 0.65))
dev.off()

res = list(Nuclear.v.Cytoplasmic = t.test(x = major[major$Dir=="Nuclear","Length"],
                                          y = major[major$Dir=="Cytoplasmic","Length"]),
           Nuclear.v.Cytoplasmic.both = t.test(x = major[major$fracReg=="both_retained","Length"],
                                               y = major[major$fracReg=="both_exported","Length"]),
           Nuclear.v.Cytoplasmic.AdultOnly = t.test(x = major[major$fracReg=="Ad_retained","Length"],
                                                    y = major[major$fracReg=="Ad_exported","Length"]),
           Nuclear.v.Cytoplasmic.PrenatalOnly = t.test(x = major[major$fracReg=="Fet_retained","Length"],
                                                       y = major[major$fracReg=="Fet_exported","Length"]))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            retainedMean = x$estimate[1], exportedMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#                          Comparison      Tstat         pval retainedMean exportedMean          FDR
#1              Nuclear.v.Cytoplasmic  6.9023845 7.003869e-12     805.4851     544.6010 1.400774e-11
#2         Nuclear.v.Cytoplasmic.both  0.3577633 7.222646e-01     603.1329     559.1071 7.222646e-01
#3    Nuclear.v.Cytoplasmic.AdultOnly  7.2507073 6.679374e-13     843.5193     541.3727 2.671749e-12
#4 Nuclear.v.Cytoplasmic.PrenatalOnly -3.4091896 1.323315e-02     549.3889    1400.2500 1.764420e-02
res[res$FDR<=0.05,]


### Is there a significant difference between the free energy of the 3'UTRs that are regulated in each group?

## Make the input FASTA file for RNAfold

sigUTR = makeGRangesFromDataFrame(sigUTR, keep.extra.columns = T)
names(sigUTR) = sigUTR$exonID
Seq = getSeq(Hsapiens, sigUTR[sigUTR$dominant3UTR=="Major",])
Seq = RNAStringSet(reverseComplement(Seq))

writeXStringSet(Seq, filepath = "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/fracDEGs_3UTRs.fa")
length(unique(names(Seq)))

## Run RNAfold

# /Users/amandaprice/Desktop/ViennaRNA/bin/RNAfold /Users/amandaprice/Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/fracDEGs_3UTRs.fa > ./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/fracDEGs_3UTRs.out.txt 


## Plot the difference between gene sets

res = read.table("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/fracDEGs_3UTRs.out.txt", header = FALSE, fill = TRUE)
energy = data.frame(exonID = res$V1[grep(">", res$V1, fixed = T)], energy = res$V2[grep("(", res$V2, fixed=T)])
energy$exonID = gsub(">","", energy$exonID)
energy$exonID = gsub("(","",energy$exonID, fixed=T)
energy$energy = gsub("(","",energy$energy, fixed=T)
energy$energy = gsub(")","",energy$energy, fixed=T)
energy$energy = as.numeric(energy$energy)
head(energy)

major = data.frame(sigUTR)[sigUTR$dominant3UTR=="Major",]
major$MFE = energy[match(major$exonID, energy$exonID),"energy"]
major$Dir = major$group = NA
major[grep("retained", major$fracReg),"Dir"] = "Nuclear"
major[grep("exported", major$fracReg),"Dir"] = "Cytoplasmic"
major[grep("Ad_", major$fracReg),"group"] = "In Adult"
major[grep("Fet_", major$fracReg),"group"] = "In Prenatal"
major[grep("both", major$fracReg),"group"] = "In Both"
major$group = factor(major$group, levels = c("In Both", "In Prenatal", "In Adult"))
major$Dir = factor(major$Dir, levels = c("Cytoplasmic","Nuclear"))
major$Length = major$Length/1000

major = reshape2::melt(major, measure.vars = c("MFE", "Length"))
major$variable = gsub("Length","Length (Kb)", major$variable)

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/free_energy_fracDEGs_3UTRs.pdf",width=5.4,height=3)
ggplot(major[which(major$variable=="MFE"),], aes(x = Dir, y = value, fill = Dir)) + geom_jitter() + geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  facet_grid(. ~ group) + labs(fill="") + ylab("MFE") + xlab("") +
  ggtitle("3'UTR Minimum Free Energy") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 0.65))
dev.off()

pdf("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/figures/free_energy_length_fracDEGs_3UTRs.pdf",width=6.8,height=4.5)
ggplot(major, aes(x = Dir, y = value, fill = Dir)) + geom_jitter() + geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  facet_grid(variable ~ group, scales = "free") + labs(fill="") + ylab("") + xlab("") +
  ggtitle("3'UTR Length, Minimum Free Energy") +
  theme(title = element_text(size = 20), text = element_text(size = 20), legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 0.65))
dev.off()


res = list(Nuclear.v.Cytoplasmic = t.test(x = major[major$Dir=="Nuclear" & major$variable=="MFE","value"],
                                      y = major[major$Dir=="Cytoplasmic" & major$variable=="MFE","value"]),
           Nuclear.v.Cytoplasmic.both = t.test(x = major[major$fracReg=="both_retained" & major$variable=="MFE","value"],
                                          y = major[major$fracReg=="both_exported" & major$variable=="MFE","value"]),
           Nuclear.v.Cytoplasmic.AdultOnly = t.test(x = major[major$fracReg=="Ad_retained" & major$variable=="MFE","value"],
                                                    y = major[major$fracReg=="Ad_exported" & major$variable=="MFE","value"]),
           Nuclear.v.Cytoplasmic.PrenatalOnly = t.test(x = major[major$fracReg=="Fet_retained" & major$variable=="MFE","value"],
                                                       y = major[major$fracReg=="Fet_exported" & major$variable=="MFE","value"]))
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(Tstat = x$statistic, pval = x$p.value,
                                                                                            retainedMean = x$estimate[1], exportedMean = x$estimate[2]))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
res
#                          Comparison      Tstat         pval retainedMean exportedMean          FDR
#1              Nuclear.v.Cytoplasmic -11.751156 1.306003e-30    -285.2657    -140.3806 5.224011e-30
#2         Nuclear.v.Cytoplasmic.both  -3.851094 2.129059e-04    -249.3937    -118.1600 2.838746e-04
#3    Nuclear.v.Cytoplasmic.AdultOnly -11.224747 6.275807e-28    -292.8529    -140.3322 1.255161e-27
#4 Nuclear.v.Cytoplasmic.PrenatalOnly   2.111832 6.355328e-02    -196.4833    -310.2500 6.355328e-02
res[res$FDR<=0.05,]

