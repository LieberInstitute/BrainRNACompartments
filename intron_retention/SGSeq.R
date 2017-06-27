library("GenomicFeatures")
library("GenomicRanges")
library("SGSeq")
library("ggplot2")
library("DEXSeq")
library(plyr)

load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# http://www.bioconductor.org/packages/release/bioc/vignettes/SGSeq/inst/doc/SGSeq.html#overview

###  Create the TranscriptDb object from gtf file
# http://chitka-kalyan.blogspot.com/2014/02/creating-gencode-transcript-database-in.html

# Download the latest gencode comprehensive gtf file from gencode website

gencode <- makeTxDbFromGFF("/Users/amanda/Downloads/gencode.v25lift37.annotation.gtf", 
                           dataSource=paste("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz"),
                           organism="Homo sapiens")

# Save the transcriptDb object as a sql database object
saveDb(gencode, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
loadDb("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")

# For SGSeq, create sample information dataframe

frag_length = read.table("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/MISO_out/insert-dist/means.txt", col.names = c("sample_name", "mean"))
file_bam = c(paste0("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/pipeline_output/polyA/HISAT2_out/", frag_length$sample_name[1:12], "_accepted_hits.sorted.bam"),
             paste0("/dcl01/lieber/ajaffe/Amanda/NucVsCyt/pipeline_output/downsampling/HISAT2_out/", frag_length$sample_name[13:14], "_accepted_hits.sorted.bam"))

si = data.frame(sample_name = frag_length$sample_name, file_bam = file_bam)
si$sample_name = as.character(si$sample_name)
si$file_bam = as.character(si$file_bam)
si = getBamInfo(si)
save(si, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/sample_info.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/sample_info.rda")
si$file_bam = paste0("/Users/amanda/Dropbox/NucVsCytosol/BAM/", si$sample_name, "_accepted_hits.sorted.bam")
si = si[which(si$sample_name != "Br5340C1_polyA" & si$sample_name != "Br5339C1_polyA"),]

# Extract transcript features
txf = convertToTxFeatures(gencode)
sgf = convertToSGFeatures(txf)

# Analyze splice variants in our samples
sgfc = analyzeFeatures(si, features = txf)
sgvc10 = analyzeVariants(sgfc, min_denominator = 10)
save(sgfc, sgvc, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/SGSeq_objects.rda")
save(sgfc, sgvc, sgvc10, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/SGSeq_objects10.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/SGSeq_objects10.rda")

# subset variant results by type

splicetype = list(SE = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("SE:S", x)) }), ],
                  S2E = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("S2E:S", x)) }), ],
                  RI = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("RI:R", x)) }), ],
                  MXE = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("MXE", x)) }), ],
                  A5SS.P = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("A5SS:P", x)) }), ],
                  A3SS.P = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("A3SS:P", x)) }), ],
                  A5SS.D = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("A5SS:D", x)) }), ],
                  A3SS.D = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("A3SS:D", x)) }), ],
                  AFE = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("AFE", x)) }), ],
                  ALE = sgvc[sapply(variantType(sgvc), function(x) { any(grepl("ALE", x)) }), ])

splicetype10 = list(SE = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("SE:S", x)) }), ],
                    S2E = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("S2E:S", x)) }), ],
                    RI = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("RI:R", x)) }), ],
                    MXE = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("MXE", x)) }), ],
                    A5SS.P = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("A5SS:P", x)) }), ],
                    A3SS.P = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("A3SS:P", x)) }), ],
                    A5SS.D = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("A5SS:D", x)) }), ],
                    A3SS.D = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("A3SS:D", x)) }), ],
                    AFE = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("AFE", x)) }), ],
                    ALE = sgvc10[sapply(variantType(sgvc10), function(x) { any(grepl("ALE", x)) }), ])

# Compare PSI rate of different types of splicing events by fraction and age
psi = lapply(splicetype, variantFreq)
psi10 = lapply(splicetype10, variantFreq)

psi_filt = lapply(psi, function(x) na.omit(x))
psi_filt10 = lapply(psi10, function(x) na.omit(x))

psi_df =  psi_df10 = list()
for (i in (1:length(psi_filt))){
  tmp = psi_filt[[i]]
  tmp10 = psi_filt10[[i]]
  psi_df[[i]] = data.frame(melt(tmp), VariantType = names(psi_filt)[i])
  psi_df10[[i]] = data.frame(melt(tmp10), VariantType = names(psi_filt10)[i])}
psi_df = do.call(rbind, psi_df)
psi_df10 = do.call(rbind, psi_df10)
colnames(psi_df) = colnames(psi_df10) = c("rowID", "SampleID", "PSI", "VariantType")
psi_df$rowNum = rownames(psi_df)
psi_df$Fraction = ifelse((psi_df$rowNum %in% grep("C", psi_df$SampleID)), "Cytosol", "Nucleus")
psi_df$Age = ifelse((psi_df$rowNum %in% grep("53", psi_df$SampleID)), "Prenatal", "Adult")
psi_df10$rowNum = rownames(psi_df10)
psi_df10$Fraction = ifelse((psi_df10$rowNum %in% grep("C", psi_df10$SampleID)), "Cytosol", "Nucleus")
psi_df10$Age = ifelse((psi_df10$rowNum %in% grep("53", psi_df10$SampleID)), "Prenatal", "Adult")

ggplot(psi_df, aes(x = Fraction, y = PSI, fill = Age)) + geom_boxplot() +
  facet_grid(. ~ VariantType) +
  labs(fill="") +
  ylab("PSI") + 
  xlab("") +
  ggtitle("PSI by Variant Type, Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

ggplot(psi_df10, aes(x = Fraction, y = PSI, fill = Age)) + geom_boxplot() +
  facet_grid(. ~ VariantType) +
  labs(fill="") +
  ylab("PSI") + 
  xlab("") +
  ggtitle("PSI by Variant Type, Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# PSI_by_SpliceVariantType_Fraction_Age.pdf

tfrac = tAge = list()
for (i in 1:length(names(psi))){
  tfrac[[i]] = t.test(psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Fraction=="Cytosol"),"PSI"], 
                 psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Fraction=="Nucleus"),"PSI"])
  tAge[[i]] = t.test(psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Age=="Adult"),"PSI"], 
                psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Age=="Prenatal"),"PSI"])
}
names(tfrac) = names(tAge) = names(psi)
## By fraction
#SE
#t = 1.0894, df = 80736, p-value = 0.276
#S2E
#t = 0.34473, df = 11433, p-value = 0.7303
#RI
#t = -14.228, df = 133610, p-value < 2.2e-16
#mean of x mean of y 
#0.4341496 0.4677697 
#MXE
#t = -0.0060693, df = 6968.5, p-value = 0.9952
#A5SS.P
#t = -3.5701, df = 47083, p-value = 0.0003572
#mean of x mean of y 
#0.5443473 0.5585771 
#A3SS.P
#t = -4.0372, df = 56768, p-value = 5.416e-05
#mean of x mean of y 
#0.5911044 0.6056086 
#A5SS.D
#t = 3.5539, df = 47083, p-value = 0.00038
#mean of x mean of y 
#0.4517250 0.4375905 
#A3SS.D
#t = 4.0527, df = 56767, p-value = 5.069e-05
#mean of x mean of y 
#0.4055084 0.3909861 
#AFE
#t = 0.30793, df = 34580, p-value = 0.7581
#mean of x mean of y 
#0.3857587 0.3842505 
#ALE
#t = 0.20726, df = 12790, p-value = 0.8358

## By Age
#SE
#t = -4.7852, df = 80732, p-value = 1.712e-06
#  mean of x mean of y 
#0.4690351 0.4842983 
#S2E
#t = -1.9716, df = 11433, p-value = 0.04868
#  mean of x mean of y 
#0.2012128 0.2141766 
#RI
#t = -2.4828, df = 133800, p-value = 0.01304
#  mean of x mean of y 
#0.4480242 0.4538951 
#MXE
#t = 0.77321, df = 6968.5, p-value = 0.4394
#A5SS.P
#t = 1.4025, df = 47101, p-value = 0.1608
#A3SS.P
#t = -0.50941, df = 56819, p-value = 0.6105
#A5SS.D
#t = -1.4026, df = 47101, p-value = 0.1607
#A3SS.D
#t = 0.46762, df = 56819, p-value = 0.6401
#AFE
#t = -0.1805, df = 34581, p-value = 0.8568
#ALE
#t = -0.065651, df = 12790, p-value = 0.9477

## Look for types of splicing event PSI changes in different groups of genes
lapply(sig, head)
elementNROWS(sig)
byGene = type_byGene = psi_byGene =  psi_byGene_collapsed = list()
for (i in 1:length(sig)){
  gene = sig[[i]]
  byGene[[i]] = sgvc10[sapply(geneName(sgvc10), function(x) { any(x %in% gene$geneID) }), ]
  type_byGene[[i]] = list(SE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("SE:S", x)) }), ],
                      S2E = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("S2E:S", x)) }), ],
                      RI = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("RI:R", x)) }), ],
                      MXE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("MXE", x)) }), ],
                      A5SS.P = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A5SS:P", x)) }), ],
                      A3SS.P = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A3SS:P", x)) }), ],
                      A5SS.D = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A5SS:D", x)) }), ],
                      A3SS.D = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A3SS:D", x)) }), ],
                      AFE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("AFE", x)) }), ],
                      ALE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("ALE", x)) }), ])
  psi_byGene[[i]] = lapply(type_byGene[[i]], variantFreq)
  psi_byGene[[i]] = lapply(psi_byGene[[i]], function(x) na.omit(x))
  }
names(psi_byGene) = names(type_byGene) = names(byGene) = names(sig)
elementNROWS(psi_byGene[[1]])
lapply(psi_byGene[[1]], head)

psi_byGene_df = list(list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in (1:length(psi_byGene))){
for (j in (1:length(psi_byGene[[i]]))){
  tmp = psi_byGene[[i]][[j]]
  if (nrow(tmp) > 0) {
  psi_byGene_df[[i]][[j]] = data.frame(melt(tmp), VariantType = names(psi_byGene[[i]][j]))
  }
  psi_byGene_collapsed[[i]] = do.call(rbind, psi_byGene_df[[i]])
  psi_byGene_collapsed[[i]][,"GeneGroup"] = names(psi_byGene[i])
  }}
names(psi_byGene_df) = names(psi_byGene_collapsed) = names(psi_byGene)
psi_byGene_df = do.call(rbind, psi_byGene_collapsed)

colnames(psi_byGene_df) = c("rowID", "SampleID", "PSI", "VariantType", "GeneGroup")
psi_byGene_df$rowNum = 1:nrow(psi_byGene_df)
psi_byGene_df$Fraction = ifelse((psi_byGene_df$rowNum %in% grep("C", psi_byGene_df$SampleID)), "Cytosol", "Nucleus")
psi_byGene_df$Age = ifelse((psi_byGene_df$rowNum %in% grep("53", psi_byGene_df$SampleID)), "Prenatal", "Adult")

psi_byGene_df$GeneGroup = gsub("both_retained", "Retained:\nBoth", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("both_exported", "Exported:\nBoth", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("Fet_retained", "Retained:\nPrenatal Only", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("Ad_retained", "Retained:\nAdult Only", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("Fet_exported", "Exported:\nPrenatal Only", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("Ad_exported", "Exported:\nAdult Only", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("ret_Ad_exp_Fet", "Retained: Adult/\nExported: Prenatal", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("ret_Fet_exp_Ad", "Retained: Prenatal/\nExported: Adult", psi_byGene_df$GeneGroup)
psi_byGene_df$GeneGroup = gsub("interacting", "Interaction", psi_byGene_df$GeneGroup)

ggplot(psi_byGene_df, aes(x = Fraction, y = PSI, fill = Age)) + geom_boxplot() +
  facet_grid(GeneGroup ~ VariantType) +
  labs(fill="") +
  ylab("PSI") + 
  xlab("") +
  ggtitle("PSI by Variant Type, Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# PSI_by_SpliceVariantType_Fraction_Age_in_DEGs.pdf


## Compare individual splice variant differences by variant type across fraction and age

# Construct the objects needed to test individual splice variant differences
sgv = rowRanges(sgvc10)
sgv = getSGVariantCounts(sgv, sample_info = si)
sgv
save(sgfc, sgvc, sgvc10, sgv, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/SGSeq_objects10.rda")
sgv.counts = counts(sgv)
vid = as.character(variantID(sgv))
eid = as.character(eventID(sgv))

ageFullModel = ~ sample + exon + Zone:exon + Fetal:exon
ageReducedModel = ~ sample + exon + Zone:exon
fracFullModel = ~ sample + exon + Fetal:exon + Zone:exon
fracReducedModel = ~ sample + exon + Fetal:exon

sampleData = pd[match(colnames(sgv)[1:10], rownames(pd)),]
sampleData = rbind(sampleData, pd[grep("Br5339C1_polyA", rownames(pd)),], pd[grep("Br5340C1_polyA", rownames(pd)),])
rownames(sampleData) = c(rownames(sampleData)[1:10], "Br5339C1_downsamp", "Br5340C1_downsamp")
sampleData$SampleID = rownames(sampleData)

agedxd = DEXSeqDataSet(countData = sgv.counts, featureID = vid, groupID = eid, sampleData = sampleData,
                    design= ageFullModel)
fracdxd = DEXSeqDataSet(countData = sgv.counts, featureID = vid, groupID = eid, sampleData = sampleData,
                        design= fracFullModel)
save(agedxd,fracdxd, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")

agedxd.cyt = DEXSeqDataSet(countData = sgv.counts[,c(1,3,5,9,11,12)], featureID = vid, groupID = eid, 
                           sampleData = sampleData[which(sampleData$Zone=="Cytosol"),],
                           design= ~ sample + exon + Fetal:exon)
agedxd.nuc = DEXSeqDataSet(countData = sgv.counts[,c(2,4,6,7,8,10)], featureID = vid, groupID = eid, 
                           sampleData = sampleData[which(sampleData$Zone=="Nucleus"),],
                           design= ~ sample + exon + Fetal:exon)
fracdxd.adult = DEXSeqDataSet(countData = sgv.counts[,c(1:6)], featureID = vid, groupID = eid, 
                              sampleData = sampleData[which(sampleData$Fetal=="Adult"),], 
                              design=  ~ sample + exon + Zone:exon)
fracdxd.prenatal = DEXSeqDataSet(countData = sgv.counts[,c(7:12)], featureID = vid, groupID = eid, 
                                 sampleData = sampleData[which(sampleData$Fetal=="Prenatal"),],
                                 design= ~ sample + exon + Zone:exon)
save(agedxd.cyt, agedxd.nuc, fracdxd.adult, fracdxd.prenatal, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")


# Calculate differential splicing
# by Fraction
fracdxd = estimateSizeFactors(fracdxd)
fracdxd = estimateDispersions(fracdxd, formula = fracFullModel)
plotDispEsts(fracdxd)
fracdxd = testForDEU(fracdxd, reducedModel = fracReducedModel, fullModel = fracFullModel)
fracdxd = estimateExonFoldChanges(fracdxd, fitExpToVar="Zone")
fracdxr = DEXSeqResults(fracdxd)
plotMA(fracdxr, ylim = c(-5,5), main = "Differential Splicing by Fraction", alpha=0.05)
# Differential_splicing_byFraction.pdf
save(fracdxd,fracdxr, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")

# by Age
agedxd = estimateSizeFactors(agedxd)
agedxd = estimateDispersions(agedxd, formula = ageFullModel)
plotDispEsts(agedxd)
agedxd = testForDEU(agedxd, reducedModel = ageReducedModel, fullModel = ageFullModel)
agedxd = estimateExonFoldChanges(agedxd, fitExpToVar="Fetal")
agedxr = DEXSeqResults(agedxd)
plotMA(agedxr, ylim = c(-5,5), main = "Differential Splicing by Age", alpha=0.05)
# Differential_splicing_byAge.pdf
save(agedxd,agedxr,fracdxd,fracdxr, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")

# by Fraction in Adult
fracdxd.adult = estimateSizeFactors(fracdxd.adult)
fracdxd.adult = estimateDispersions(fracdxd.adult, formula = ~ sample + exon + Zone:exon)
fracdxd.adult = testForDEU(fracdxd.adult, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Zone:exon)
fracdxd.adult = estimateExonFoldChanges(fracdxd.adult, fitExpToVar="Zone")
fracdxr.adult = DEXSeqResults(fracdxd.adult)
plotMA(fracdxr.adult, ylim = c(-5,5), main = "Differential Splicing by Fraction in Adult", alpha=0.05)
# Differential_splicing_byFraction_inAdult.pdf
save(fracdxd.adult,fracdxr.adult, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# by Fraction in Prenatal
fracdxd.prenatal = estimateSizeFactors(fracdxd.prenatal)
fracdxd.prenatal = estimateDispersions(fracdxd.prenatal, formula = ~ sample + exon + Zone:exon)
fracdxd.prenatal = testForDEU(fracdxd.prenatal, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Zone:exon)
fracdxd.prenatal = estimateExonFoldChanges(fracdxd.prenatal, fitExpToVar="Zone")
fracdxr.prenatal = DEXSeqResults(fracdxd.prenatal)
plotMA(fracdxr.prenatal, ylim = c(-5,5), main = "Differential Splicing by Fraction in Prenatal", alpha=0.05)
# Differential_splicing_byFraction_inPrenatal.pdf
save(fracdxd.adult,fracdxr.adult,fracdxd.prenatal,fracdxr.prenatal, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# by Age in Nucleus
agedxd.nuc = estimateSizeFactors(agedxd.nuc)
agedxd.nuc = estimateDispersions(agedxd.nuc, formula = ~ sample + exon + Fetal:exon)
agedxd.nuc = testForDEU(agedxd.nuc, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Fetal:exon)
agedxd.nuc = estimateExonFoldChanges(agedxd.nuc, fitExpToVar="Fetal")
agedxr.nuc = DEXSeqResults(agedxd.nuc)
plotMA(agedxr.nuc, ylim = c(-5,5), main = "Differential Splicing by Age in Nucleus", alpha=0.05)
# Differential_splicing_byFraction_inNucleus.pdf
save(fracdxd.adult,fracdxr.adult,agedxd.nuc,agedxr.nuc, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# by Age in Cytosol
agedxd.cyt = estimateSizeFactors(agedxd.cyt)
agedxd.cyt = estimateDispersions(agedxd.cyt, formula = ~ sample + exon + Fetal:exon)
agedxd.cyt = testForDEU(agedxd.cyt, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Fetal:exon)
agedxd.cyt = estimateExonFoldChanges(agedxd.cyt, fitExpToVar="Fetal")
agedxr.cyt = DEXSeqResults(agedxd.cyt)
plotMA(agedxr.cyt, ylim = c(-5,5), main = "Differential Splicing by Age in Cytosol", alpha=0.05)
# Differential_splicing_byFraction_inCytosol.pdf
save(fracdxd.adult,fracdxr.adult,fracdxd.prenatal,fracdxr.prenatal,agedxd.nuc,agedxr.nuc,agedxd.cyt,agedxr.cyt, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# Explore results
dexres = list(Fraction = fracdxr[order(fracdxr$padj),], Age = agedxr[order(agedxr$padj),], 
              frac.adult = fracdxr.adult[order(fracdxr.adult$padj),], 
              frac.prenatal = fracdxr.prenatal[order(fracdxr.prenatal$padj),], 
              age.cytosol = agedxr.cyt[order(agedxr.cyt$padj),], age.nucleus = agedxr.nuc[order(agedxr.nuc$padj),])
mcols = mcols(sgv)
dexres = lapply(dexres, function(x) DataFrame(x, vID = as.integer(x$featureID)))
dexres = Map(cbind, dexres, lapply(dexres, function(y) mcols[match(y$vID, mcols$variantID),]))
dexres = lapply(dexres, function(y) DataFrame(y, more.in.nuc.prenatal = lapply(dexres, function(z) ifelse((z[,10]>0), "Yes", "No"))))
dexres.sig = lapply(dexres, function(x) x[which(x$padj<=0.05),])

# How many of each splicing event type are significantly different by group?
elementNROWS(dexres.sig)
#Fraction           Age    frac.adult frac.prenatal   age.cytosol   age.nucleus 
#2158          4608          2512           131          2182          2226
type = lapply(dexres.sig, function(x) unlist(x$variantType))
type = lapply(type, function(x) count(x))
splice = c("SE:S","S2E:S","RI:R","MXE","A5SS:D","A3SS:D","AFE")
type = lapply(type, function(x) x[which(x$x %in% splice),])

for (i in 1:length(type)){  
  colnames(type[[i]]) = c("splice","freq")
  type[[i]] = cbind(type[[i]], comparison = names(type)[i])
}
type = do.call(rbind, type)
type$comparison = gsub("Fraction", "By Fraction", type$comparison)
type$comparison = gsub("Age", "By Age", type$comparison)
type$comparison = gsub("frac.adult", "By Fraction\nIn Adult", type$comparison)
type$comparison = gsub("frac.prenatal", "By Fraction\nIn Prenatal", type$comparison)
type$comparison = gsub("age.cytosol", "By Age\nIn Cytosol", type$comparison)
type$comparison = gsub("age.nucleus", "By Age\nIn Nucleus", type$comparison)
type$comparison = factor(type$comparison, 
                         levels = c("By Age\nIn Nucleus","By Age\nIn Cytosol","By Age",
                                    "By Fraction\nIn Prenatal","By Fraction\nIn Adult","By Fraction"))

ggplot(type[which(type$comparison!="By Fraction" & type$comparison!="By Age"),], 
       aes(x = comparison, y = freq, fill = splice)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Differentially Spliced Events (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# DSE_counts_byGroup.pdf

ggplot(type[grep("Fraction", type$comparison),], 
       aes(x = comparison, y = freq, fill = splice)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Differentially Spliced Events (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# DSE_counts_byFraction.pdf

ggplot(type[grep("Age", type$comparison),], 
       aes(x = comparison, y = freq, fill = splice)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Differentially Spliced Events (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
# DSE_counts_byAge.pdf


