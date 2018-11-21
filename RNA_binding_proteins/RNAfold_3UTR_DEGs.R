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
major$Dir = factor(major$Dir, levels = c("Nuclear", "Cytoplasmic"))

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/length_fracDEGs_3UTRs.pdf",width=8,height=4)
ggplot(major, aes(x = Dir, y = Length, fill = Dir)) + geom_boxplot() +
  facet_grid(. ~ group) + labs(fill="") + ylab("MFE") + xlab("") +
  ggtitle("3'UTR Length") +
  theme(title = element_text(size = 20), text = element_text(size = 20))
dev.off()


### Is there a significant difference between the free energy of the 3'UTRs that are regulated in each group?

## Make the input FASTA file for RNAfold

sigUTR = makeGRangesFromDataFrame(sigUTR, keep.extra.columns = T)
names(sigUTR) = sigUTR$exonID
Seq = getSeq(Hsapiens, sigUTR[sigUTR$dominant3UTR=="Major",])
Seq = RNAStringSet(reverseComplement(Seq))

writeXStringSet(Seq, filepath = "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/fracDEGs_3UTRs.fa")
length(unique(names(Seq)))

## Run RNAfold

/Users/amandaprice/Desktop/ViennaRNA/bin/RNAfold /Users/amandaprice/Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/fracDEGs_3UTRs.fa > ./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/fracDEGs_3UTRs.out.txt 

res = read.table("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/RNAfold_out/fracDEGs_3UTRs.out.txt", header = FALSE, 
                 sep = "\t", col.names = paste0("V",seq_len(4)), fill = TRUE)


## Plot the difference between gene sets

res = read.table("/Users/amandaprice/Desktop/ViennaRNA/out/out.txt", header = FALSE, fill = TRUE)
energy = data.frame(exonID = res$V1[grep(">", res$V1, fixed = T)], energy = res$V2[grep("(", res$V2, fixed=T)])
energy$exonID = gsub(">","", energy$exonID)
energy$energy = gsub("(","",energy$energy, fixed=T)
energy$energy = gsub(")","",energy$energy, fixed=T)
energy$energy = as.numeric(energy$energy)

major = sigUTR[sigUTR$dominant3UTR=="Major",]
major$energy = energy[match(major$exonID, energy$exonID),"energy"]
major$Dir = major$group = NA
major[grep("retained", major$fracReg),"Dir"] = "Nuclear"
major[grep("exported", major$fracReg),"Dir"] = "Cytoplasmic"
major[grep("Ad_", major$fracReg),"group"] = "In Adult"
major[grep("Fet_", major$fracReg),"group"] = "In Prenatal"
major[grep("both", major$fracReg),"group"] = "In Both"
major$group = factor(major$group, levels = c("In Both", "In Prenatal", "In Adult"))
major$Dir = factor(major$Dir, levels = c("Nuclear", "Cytoplasmic"))


pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/free_energy_fracDEGs_3UTRs.pdf",width=8,height=4)
ggplot(major, aes(x = Dir, y = energy, fill = Dir)) + geom_boxplot() +
  facet_grid(. ~ group) + labs(fill="") + ylab("MFE") + xlab("") +
  ggtitle("3'UTR Predicted Minimum Free Energy (MFE)") +
  theme(title = element_text(size = 20), text = element_text(size = 20))
dev.off()