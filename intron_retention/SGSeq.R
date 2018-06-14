library("GenomicFeatures")
library("GenomicRanges")
library("SGSeq")
library("ggplot2")
library("DEXSeq")
library(data.table)
library(reshape2)

load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
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

# How many variants of each type are identified?
varnames = lapply(splicetype, rowData)
varnames = lapply(varnames, function(x) length(unique(as.character(x$variantName))))
proportion = as.data.frame(unlist(varnames))
proportion$variant = rownames(proportion)
proportion = proportion[-grep("P",proportion$variant),]
proportion$variant = gsub(".D","", proportion$variant)
proportion$variant = factor(proportion$variant, levels = c("SE","S2E","RI","MXE","A5SS","A3SS","AFE","ALE"))
proportion$prop = proportion[,1] / sum(proportion[,1])
colnames(proportion)[1] = "total"
proportion$perc = paste0(round((proportion$prop*100),digits = 1),"%")

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/total_unique_splice_variants.pdf")
ggplot(proportion, aes(x = variant, y = total)) + geom_col() +
  geom_text(aes(label=perc), vjust=1.5, colour="white") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Unique Splice Variants by Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

varnames = lapply(splicetype10, rowData)
varnames = lapply(varnames, function(x) length(unique(as.character(x$variantName))))
proportion = as.data.frame(unlist(varnames))
proportion$variant = rownames(proportion)
proportion = proportion[-grep("P",proportion$variant),]
proportion$variant = gsub(".D","", proportion$variant)
proportion$prop = proportion[,1] / sum(proportion[,1])
colnames(proportion)[1] = "total"
proportion$perc = paste0(round((proportion$prop*100),digits = 1),"%")
write.csv(proportion, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/total_unique_splice_variants_10denom.csv")
proportion = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/total_unique_splice_variants_10denom.csv")
proportion$variant = gsub("RI","IR", proportion$variant)
proportion$variant = factor(proportion$variant, levels = c("SE","S2E","IR","MXE","A5SS","A3SS","AFE","ALE"))


pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/total_unique_splice_variants_10denom.pdf", width = 6, height = 5)
ggplot(proportion, aes(x = variant, y = total)) + geom_col() +
  geom_text(aes(label=perc), vjust=1.5, colour="white") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Unique Splice Variants by Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


# How many variants of each type are identified by fraction and age?
varnames = lapply(splicetype, rowData)
varnames = lapply(varnames, function(x) x$variantName)
psi = lapply(splicetype, variantFreq)
for (i in 1:length(psi)){
 tmp = psi[[i]]
 rownames(tmp) = varnames[[i]]
 psi[[i]] = tmp
}
lapply(psi, head)

psi = do.call(rbind, psi)
psi = as.data.frame(psi)
psi$suminAdNuc = rowSums(psi[,c("Br1113N1_polyA","Br2046N_polyA","Br2074N_polyA")])
psi$suminAdCyt = rowSums(psi[,c("Br1113C1_polyA","Br2046C_polyA","Br2074C_polyA")])
psi$suminFetNuc = rowSums(psi[,c("Br5339N1_polyA","Br5340N1_polyA","Br5341N1_polyA")])
psi$suminFetCyt = rowSums(psi[,c("Br5341C1_polyA","Br5339C1_downsamp","Br5340C1_downsamp")])
psi$variantID = rownames(psi)

psiByGroup = list("Adult:Cytoplasm" = psi[which(psi$suminAdCyt!="NA" & psi$suminAdCyt>0),], 
                  "Adult:Nucleus" = psi[which(psi$suminAdNuc!="NA" & psi$suminAdNuc>0),],
                  "Prenatal:Cytoplasm" = psi[which(psi$suminFetCyt!="NA" & psi$suminFetCyt>0),], 
                  "Prenatal:Nucleus" = psi[which(psi$suminFetNuc!="NA" & psi$suminFetNuc>0),])
VariantsByGroup = lapply(psiByGroup, function(x) x$variantID)
numVars = data.frame(SE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("SE",x)])))),
                     S2E = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("S2E",x)])))),
                     RI = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("RI",x)])))),
                     MXE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("MXE",x)])))),
                     A5SS = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("A5SS",x)])))),
                     A3SS = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("A3SS",x)])))),
                     AFE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("AFE",x)])))),
                     ALE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("ALE",x)])))))
numVars$Group = as.factor(rownames(numVars))
numVars = melt(numVars)
head(numVars)

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/total_unique_splice_variants_byGroup.pdf",width = 10,height = 8)
dodge <- position_dodge(width=0.9)
ggplot(numVars, aes(x = variable, y = value, fill = Group)) +
  stat_summary(position=position_dodge(),geom="bar") +
  ylim(0,21000) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Unique Splice Variants by Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

varnames = lapply(splicetype10, rowData)
varnames = lapply(varnames, function(x) x$variantName)
psi = lapply(splicetype10, variantFreq)
for (i in 1:length(psi)){
  tmp = psi[[i]]
  rownames(tmp) = varnames[[i]]
  psi[[i]] = tmp
}
lapply(psi, head)

psi = do.call(rbind, psi)
psi = as.data.frame(psi)
psi$suminAdNuc = rowSums(psi[,c("Br1113N1_polyA","Br2046N_polyA","Br2074N_polyA")])
psi$suminAdCyt = rowSums(psi[,c("Br1113C1_polyA","Br2046C_polyA","Br2074C_polyA")])
psi$suminFetNuc = rowSums(psi[,c("Br5339N1_polyA","Br5340N1_polyA","Br5341N1_polyA")])
psi$suminFetCyt = rowSums(psi[,c("Br5341C1_polyA","Br5339C1_downsamp","Br5340C1_downsamp")])
psi$suminNuc = rowSums(psi[,c("Br1113N1_polyA","Br2046N_polyA","Br2074N_polyA","Br5339N1_polyA","Br5340N1_polyA","Br5341N1_polyA")])
psi$suminCyt = rowSums(psi[,c("Br1113C1_polyA","Br2046C_polyA","Br2074C_polyA","Br5341C1_polyA","Br5339C1_downsamp","Br5340C1_downsamp")])
psi$suminFet = rowSums(psi[,c("Br5339N1_polyA","Br5340N1_polyA","Br5341N1_polyA","Br5341C1_polyA","Br5339C1_downsamp","Br5340C1_downsamp")])
psi$suminAd = rowSums(psi[,c("Br1113N1_polyA","Br2046N_polyA","Br2074N_polyA","Br1113C1_polyA","Br2046C_polyA","Br2074C_polyA")])
psi$variantID = rownames(psi)

psiByGroup = list("Adult:Cytoplasm" = psi[which(psi$suminAdCyt!="NA" & psi$suminAdCyt>0),], 
                  "Adult:Nucleus" = psi[which(psi$suminAdNuc!="NA" & psi$suminAdNuc>0),],
                  "Prenatal:Cytoplasm" = psi[which(psi$suminFetCyt!="NA" & psi$suminFetCyt>0),], 
                  "Prenatal:Nucleus" = psi[which(psi$suminFetNuc!="NA" & psi$suminFetNuc>0),],
                  "Nucleus" = psi[which(psi$suminNuc!="NA" & psi$suminNuc>0),],
                  "Cytoplasm" = psi[which(psi$suminCyt!="NA" & psi$suminCyt>0),],
                  "Prenatal" = psi[which(psi$suminFet!="NA" & psi$suminFet>0),],
                  "Adult" = psi[which(psi$suminAd!="NA" & psi$suminAd>0),])
VariantsByGroup = lapply(psiByGroup, function(x) x$variantID)
numVars = data.frame(SE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("SE",x)])))),
                     S2E = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("S2E",x)])))),
                     RI = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("RI",x)])))),
                     MXE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("MXE",x)])))),
                     A5SS = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("A5SS",x)])))),
                     A3SS = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("A3SS",x)])))),
                     AFE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("AFE",x)])))),
                     ALE = unlist(lapply(VariantsByGroup, function(x) length(unique(x[grep("ALE",x)])))))
numVars$Group = as.factor(rownames(numVars))
numVars = melt(numVars)
head(numVars)
write.csv(numVars, quote = F, file = "/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/total_unique_splice_variants_byGroup_10denom.csv")
numVars = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/total_unique_splice_variants_byGroup_10denom.csv")
(sum(numVars[numVars$Group=="Nucleus","value"])-sum(numVars[numVars$Group=="Cytoplasm","value"]))/sum(numVars[numVars$Group=="Cytoplasm","value"])*100 # 42.79054
(sum(numVars[numVars$Group=="Prenatal","value"])-sum(numVars[numVars$Group=="Adult","value"]))/sum(numVars[numVars$Group=="Adult","value"])*100 # 72.88765

numVars = numVars[grep(":", numVars$Group, fixed = T),]
numVars$variable = gsub("RI","IR", numVars$variable)
numVars$variable = factor(numVars$variable, levels = c("SE","S2E","IR","MXE","A5SS","A3SS","AFE","ALE"))

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/total_unique_splice_variants_byGroup_10denom.pdf",width = 6,height = 5)
dodge <- position_dodge(width=0.9)
ggplot(numVars, aes(x = variable, y = value, fill = Group)) +
  stat_summary(position=position_dodge(),geom="bar") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Unique Splice Variants by Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20), legend.position = "bottom") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()


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
psi_df$Fraction = ifelse((psi_df$rowNum %in% grep("C", psi_df$SampleID)), "Cytoplasm", "Nucleus")
psi_df$Age = ifelse((psi_df$rowNum %in% grep("53", psi_df$SampleID)), "Prenatal", "Adult")
psi_df10$rowNum = rownames(psi_df10)
psi_df10$Fraction = ifelse((psi_df10$rowNum %in% grep("C", psi_df10$SampleID)), "Cytoplasm", "Nucleus")
psi_df10$Age = ifelse((psi_df10$rowNum %in% grep("53", psi_df10$SampleID)), "Prenatal", "Adult")

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/PSI_by_SpliceVariantType_Fraction_Age.pdf",width = 24,height = 6)
ggplot(psi_df, aes(x = Fraction, y = PSI, fill = Age)) + geom_boxplot() +
  facet_grid(. ~ VariantType) +
  labs(fill="") +
  ylab("PSI") + 
  xlab("") +
  ggtitle("PSI by Variant Type, Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

pdf("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/PSI_by_SpliceVariantType_Fraction_Age_10denom.pdf",width = 24,height = 6)
ggplot(psi_df10, aes(x = Fraction, y = PSI, fill = Age)) + geom_boxplot() +
  facet_grid(. ~ VariantType) +
  labs(fill="") +
  ylab("PSI") + 
  xlab("") +
  ggtitle("PSI by Variant Type, Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()

tfrac = tAge = list()
for (i in 1:length(names(psi))){
  tfrac[[i]] = t.test(psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Fraction=="Cytoplasm"),"PSI"], 
                 psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Fraction=="Nucleus"),"PSI"])
  tAge[[i]] = t.test(psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Age=="Adult"),"PSI"], 
                psi_df[which(psi_df$VariantType==names(psi)[i] & psi_df$Age=="Prenatal"),"PSI"])
}
names(tfrac) = names(tAge) = names(psi)
tfracage = rbind(Tstat.frac = data.frame(lapply(tfrac, function(x) round(x$statistic,3))), pval.frac = data.frame(lapply(tfrac, function(x) x$p.value)),
                 confInt.frac = data.frame(lapply(tfrac, function(x) round(x$conf.int,3))), estMeans.frac = data.frame(lapply(tfrac, function(x) round(x$estimate,3))),
                 Tstat.age = data.frame(lapply(tAge, function(x) round(x$statistic,3))), pval.age = data.frame(lapply(tAge, function(x) x$p.value)),
                 confInt.age = data.frame(lapply(tAge, function(x) round(x$conf.int,3))), estMeans.age = data.frame(lapply(tAge, function(x) round(x$estimate,3))))
write.csv(tfracage,quote=F,file="/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/Ttest_PSI_byVariantType_CytVSNuc_AdultVSPrenatal_10denom.csv")
colnames(tfracage)[which(tfracage["pval.frac",]<=0.005)] # "RI"     "A5SS.P" "A3SS.P" "A5SS.D" "A3SS.D"
colnames(tfracage)[which(tfracage["pval.age",]<=0.005)] # "SE"

tfrac = tAge = list()
for (i in 1:length(names(psi10))){
  tfrac[[i]] = t.test(psi_df10[which(psi_df10$VariantType==names(psi10)[i] & psi_df10$Fraction=="Cytoplasm"),"PSI"], 
                      psi_df10[which(psi_df10$VariantType==names(psi10)[i] & psi_df10$Fraction=="Nucleus"),"PSI"])
  tAge[[i]] = t.test(psi_df10[which(psi_df10$VariantType==names(psi10)[i] & psi_df10$Age=="Adult"),"PSI"], 
                     psi_df10[which(psi_df10$VariantType==names(psi10)[i] & psi_df10$Age=="Prenatal"),"PSI"])
}
names(tfrac) = names(tAge) = names(psi10)
tfracage10 = rbind(Tstat.frac = data.frame(lapply(tfrac, function(x) round(x$statistic,3))), pval.frac = data.frame(lapply(tfrac, function(x) x$p.value)),
                 confInt.frac = data.frame(lapply(tfrac, function(x) round(x$conf.int,3))), estMeans.frac = data.frame(lapply(tfrac, function(x) round(x$estimate,3))),
                 Tstat.age = data.frame(lapply(tAge, function(x) round(x$statistic,3))), pval.age = data.frame(lapply(tAge, function(x) x$p.value)),
                 confInt.age = data.frame(lapply(tAge, function(x) round(x$conf.int,3))), estMeans.age = data.frame(lapply(tAge, function(x) round(x$estimate,3))))
write.csv(tfracage10,quote=F,file="/Users/amanda/Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/Ttest_PSI_byVariantType_CytVSNuc_AdultVSPrenatal.csv")
colnames(tfracage10)[which(tfracage10["pval.frac",]<=0.005)] # "RI"
colnames(tfracage10)[which(tfracage10["pval.age",]<=0.005)] # none


## Look for types of splicing event PSI changes in different groups of genes
lapply(sig, head)
elementNROWS(sig)
byGene = type_byGene = psi_byGene =  psi_byGene_collapsed = list()
for (i in 1:length(sig)){
  gene = sig[[i]]
  byGene[[i]] = sgvc10[sapply(geneName(sgvc10), function(x) { any(x %in% gene$geneID) }), ]
  if (length(byGene[[i]])>0) {
  type_byGene[[i]] = list(SE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("SE:S", x)) }), ],
                      S2E = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("S2E:S", x)) }), ],
                      RI = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("RI:R", x)) }), ],
                      MXE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("MXE", x)) }), ],
                      A5SS.P = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A5SS:P", x)) }), ],
                      A3SS.P = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A3SS:P", x)) }), ],
                      A5SS.D = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A5SS:D", x)) }), ],
                      A3SS.D = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("A3SS:D", x)) }), ],
                      AFE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("AFE", x)) }), ],
                      ALE = byGene[[i]][sapply(variantType(byGene[[i]]), function(x) { any(grepl("ALE", x)) }), ]) }
}
for (i in 1:length(sig)){ psi_byGene[[i]] = lapply(type_byGene[[i]], variantFreq) }
for (i in 1:length(sig)){ psi_byGene[[i]] = lapply(psi_byGene[[i]], na.omit) }
names(psi_byGene) = names(type_byGene) = names(byGene) = names(sig)
elementNROWS(psi_byGene[[1]])
lapply(psi_byGene[[1]], head)

psi_byGene_df = list(list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in (1:length(psi_byGene))) {
if (length(psi_byGene[[i]])>0) {
  for (j in (1:length(psi_byGene[[i]]))){
  tmp = psi_byGene[[i]][[j]]
  if (nrow(tmp) > 0) {
  psi_byGene_df[[i]][[j]] = data.frame(melt(tmp), VariantType = names(psi_byGene[[i]][j]))
  }}}}
for (i in (1:length(psi_byGene))) {
  if (length(psi_byGene_df[[i]])>0){
  psi_byGene_collapsed[[i]] = do.call(rbind, psi_byGene_df[[i]])
  psi_byGene_collapsed[[i]][,"GeneGroup"] = names(psi_byGene[i])
  }}
names(psi_byGene_df) = names(psi_byGene_collapsed) = names(psi_byGene)
psi_byGene_df = do.call(rbind, psi_byGene_collapsed)

colnames(psi_byGene_df) = c("rowID", "SampleID", "PSI", "VariantType", "GeneGroup")
psi_byGene_df$rowNum = 1:nrow(psi_byGene_df)
psi_byGene_df$Fraction = ifelse((psi_byGene_df$rowNum %in% grep("C", psi_byGene_df$SampleID)), "Cytoplasm", "Nucleus")
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

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/PSI_by_SpliceVariantType_Fraction_Age_in_DEGs.pdf", width=24,height=20)
ggplot(psi_byGene_df, aes(x = Fraction, y = PSI, fill = Age)) + geom_boxplot() +
  facet_grid(GeneGroup ~ VariantType) +
  labs(fill="") +
  ylab("PSI") + 
  xlab("") +
  ggtitle("PSI by Variant Type, Fraction and Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


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

agedxd.cyt = DEXSeqDataSet(countData = sgv.counts[,grep("C", colnames(sgv.counts))], featureID = vid, groupID = eid, 
                           sampleData = sampleData[which(sampleData$Zone=="Cytoplasm"),],
                           design= ~ sample + exon + Fetal:exon)
agedxd.nuc = DEXSeqDataSet(countData = sgv.counts[,grep("N", colnames(sgv.counts))], featureID = vid, groupID = eid, 
                           sampleData = sampleData[which(sampleData$Zone=="Nucleus"),],
                           design= ~ sample + exon + Fetal:exon)
fracdxd.adult = DEXSeqDataSet(countData = sgv.counts[,-grep("53", colnames(sgv.counts))], featureID = vid, groupID = eid, 
                              sampleData = sampleData[which(sampleData$Fetal=="Adult"),], 
                              design=  ~ sample + exon + Zone:exon)
fracdxd.prenatal = DEXSeqDataSet(countData = sgv.counts[,grep("53", colnames(sgv.counts))], featureID = vid, groupID = eid, 
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
save(fracdxd,fracdxr, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")

# by Age
agedxd = estimateSizeFactors(agedxd)
agedxd = estimateDispersions(agedxd, formula = ageFullModel)
plotDispEsts(agedxd)
agedxd = testForDEU(agedxd, reducedModel = ageReducedModel, fullModel = ageFullModel)
agedxd = estimateExonFoldChanges(agedxd, fitExpToVar="Fetal")
agedxr = DEXSeqResults(agedxd)
save(agedxd,agedxr,fracdxd,fracdxr, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")

# by Fraction in Adult
fracdxd.adult = estimateSizeFactors(fracdxd.adult)
fracdxd.adult = estimateDispersions(fracdxd.adult, formula = ~ sample + exon + Zone:exon)
fracdxd.adult = testForDEU(fracdxd.adult, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Zone:exon)
fracdxd.adult = estimateExonFoldChanges(fracdxd.adult, fitExpToVar="Zone")
fracdxr.adult = DEXSeqResults(fracdxd.adult)
save(fracdxd.adult,fracdxr.adult, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# by Fraction in Prenatal
fracdxd.prenatal = estimateSizeFactors(fracdxd.prenatal)
fracdxd.prenatal = estimateDispersions(fracdxd.prenatal, formula = ~ sample + exon + Zone:exon)
fracdxd.prenatal = testForDEU(fracdxd.prenatal, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Zone:exon)
fracdxd.prenatal = estimateExonFoldChanges(fracdxd.prenatal, fitExpToVar="Zone")
fracdxr.prenatal = DEXSeqResults(fracdxd.prenatal)
save(fracdxd.adult,fracdxr.adult,fracdxd.prenatal,fracdxr.prenatal, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# by Age in Nucleus
agedxd.nuc = estimateSizeFactors(agedxd.nuc)
agedxd.nuc = estimateDispersions(agedxd.nuc, formula = ~ sample + exon + Fetal:exon)
agedxd.nuc = testForDEU(agedxd.nuc, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Fetal:exon)
agedxd.nuc = estimateExonFoldChanges(agedxd.nuc, fitExpToVar="Fetal")
agedxr.nuc = DEXSeqResults(agedxd.nuc)
save(fracdxd.adult,fracdxr.adult,agedxd.nuc,agedxr.nuc, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# by Age in Cytoplasm
agedxd.cyt = estimateSizeFactors(agedxd.cyt)
agedxd.cyt = estimateDispersions(agedxd.cyt, formula = ~ sample + exon + Fetal:exon)
agedxd.cyt = testForDEU(agedxd.cyt, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + Fetal:exon)
agedxd.cyt = estimateExonFoldChanges(agedxd.cyt, fitExpToVar="Fetal")
agedxr.cyt = DEXSeqResults(agedxd.cyt)
save(fracdxd.adult,fracdxr.adult,fracdxd.prenatal,fracdxr.prenatal,agedxd.nuc,agedxr.nuc,agedxd.cyt,agedxr.cyt, 
     file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")

# Plot MA

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/MA_Plots_Differential_splicing_byFraction_byAge.pdf",height=6,width=6)
plotMA(fracdxr, ylim = c(-15,15), main = "Differential Splicing by Fraction", alpha=0.05)
plotMA(agedxr, ylim = c(-15,15), main = "Differential Splicing by Age", alpha=0.05)
plotMA(fracdxr.adult, ylim = c(-15,15), main = "Differential Splicing by Fraction in Adult", alpha=0.05)
plotMA(fracdxr.prenatal, ylim = c(-15,15), main = "Differential Splicing by Fraction in Prenatal", alpha=0.05)
plotMA(agedxr.nuc, ylim = c(-15,15), main = "Differential Splicing by Age in Nucleus", alpha=0.05)
plotMA(agedxr.cyt, ylim = c(-15,15), main = "Differential Splicing by Age in Cytoplasm", alpha=0.05)
dev.off()


load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/SGSeq_objects10.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_objects.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DEXSeq_singlevariable_objects.rda")


# Explore results
dexres = list(Fraction = fracdxr[order(fracdxr$padj),], Age = agedxr[order(agedxr$padj),], 
              frac.adult = fracdxr.adult[order(fracdxr.adult$padj),], 
              frac.prenatal = fracdxr.prenatal[order(fracdxr.prenatal$padj),], 
              age.cytosol = agedxr.cyt[order(agedxr.cyt$padj),], age.nucleus = agedxr.nuc[order(agedxr.nuc$padj),])
mcols = mcols(sgv)
dexres = lapply(dexres, function(x) DataFrame(x, vID = as.integer(x$featureID)))
dexres = Map(cbind, dexres, lapply(dexres, function(y) mcols[match(y$vID, mcols$variantID),]))
dexres = lapply(dexres, function(x) DataFrame(x, more.in.nuc.prenatal = ifelse(x[,10]>0, "Yes", "No")))
dexres.sig = lapply(dexres, function(x) x[which(x$padj<=0.05),])


# How many of each splicing event type are significantly different by group?

elementNROWS(dexres.sig)
#     Fraction           Age    frac.adult frac.prenatal   age.cytosol   age.nucleus 
#         2158          4608          2512           131          2182          2226 
type = lapply(dexres.sig, function(x) unlist(x$variantType))
type = lapply(type, function(x) data.frame(table(x)[which(names(table(x)) %in% c("SE:S","S2E:S","RI:R","MXE","A5SS:D","A3SS:D","AFE"))]))
type =  Map(cbind, type, comparison = list("By Fraction","By Age","By Fraction\nIn Adult","By Fraction\nIn Prenatal","By Age\nIn Cytoplasm","By Age\nIn Nucleus"))
type = do.call(rbind, type)
type$comparison = factor(type$comparison, 
                         levels = c("By Age\nIn Nucleus","By Age\nIn Cytoplasm","By Age","By Fraction\nIn Prenatal","By Fraction\nIn Adult","By Fraction"))
type$x = gsub("RI:R", "IR", type$x)

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/DSE_counts_byGroup_byFrac_byAge.pdf", width = 6,height = 5)
ggplot(type[which(type$comparison!="By Fraction" & type$comparison!="By Age"),], 
       aes(x = comparison, y = Freq, fill = x)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\n(FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(type[grep("Fraction", type$comparison),], 
       aes(x = comparison, y = Freq, fill = x)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\nBy Fraction (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(type[grep("Age", type$comparison),], 
       aes(x = comparison, y = Freq, fill = x)) + geom_bar(stat = "identity") +
  coord_flip() +
  labs(fill="") +
  ylab("Count") + xlab("") +
  ggtitle("Differentially Spliced Events\nBy Age (FDR < 0.05)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


## Which direction are the log fold changes by variant type?

lapply(dexres,head)
df = Map(cbind, lapply(dexres, function(x) data.frame(x[,c("groupID","featureID","padj")],LFC = x[,10], x[,c("more.in.nuc.prenatal","geneName","variantType","variantName")])),
         comparison = list("By Fraction","By Age","By Fraction\nIn Adult","By Fraction\nIn Prenatal","By Age\nIn Cytoplasm","By Age\nIn Nucleus"))
df = do.call(rbind, df)
df$geneName = as.character(df$geneName)
df$variantType = as.character(df$variantType)
df$variantName = as.character(df$variantName)
df$threshold = ifelse(df$padj<=0.05, "FDR < 0.05", "FDR > 0.05")
dt = data.table(df)
dt = dt[variantType %in% c("SE:S","S2E:S","RI:R","MXE","A5SS:P","A3SS:P","A5SS:D","A3SS:D","AFE","ALE"),,]
dt = dt[threshold!="NA",,]
dt$variantType = gsub("RI:R", "IR", dt$variantType)
dt$variantType = factor(dt$variantType, levels = c("SE:S","S2E:S","IR","MXE","A5SS:P","A3SS:P","A5SS:D","A3SS:D","AFE","ALE"))


pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/SGSeq_out/volcano_plots_byComparison_byVariantType.pdf",width=20,height=8)
ggplot(dt[comparison!="By Age" & comparison!="By Fraction",,], aes(x=LFC, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.5) +
  ylim(0,60) + xlim(-15,15) +
  scale_colour_manual(values=c("red3","gray47")) +
  geom_vline(xintercept=0, linetype="dotted") +
  facet_grid(comparison ~ variantType) +
  xlab("log2 Fold Change") + ylab("-log10(FDR)") +
  ggtitle("Differential Splicing by Fraction, Age and Variant Type") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 18), legend.position="none")
dev.off()


# Calculate the proportion of splice variants that are significant vs not significant

df$rnum = 1:nrow(df) 
df$SE = ifelse(df$rnum %in% grep("SE:S", df$variantType), "SE", "no")
df$S2E = ifelse(df$rnum %in% grep("S2E:S", df$variantType), "S2E", "no")
df$RI = ifelse(df$rnum %in% grep("RI:R", df$variantType), "RI", "no")
df$MXE = ifelse(df$rnum %in% grep("MXE", df$variantType), "MXE", "no")
df$A5SS.P = ifelse(df$rnum %in% grep("A5SS:P", df$variantType), "A5SS:P", "no")
df$A3SS.P = ifelse(df$rnum %in% grep("A3SS:P", df$variantType), "A3SS:P", "no")
df$A5SS.D = ifelse(df$rnum %in% grep("A5SS:D", df$variantType), "A5SS:D", "no")
df$A3SS.D = ifelse(df$rnum %in% grep("A3SS:D", df$variantType), "A3SS:D", "no")
df$AFE = ifelse(df$rnum %in% grep("AFE", df$variantType), "AFE", "no")
df$ALE = ifelse(df$rnum %in% grep("ALE", df$variantType), "ALE", "no")
colnames(df)[12:21]

prop = list(list(),list(),list(),list(),list(),list())
for (i in (1:length(unique(df$comparison)))) {
  for (j in 1:10) {
    prop[[i]][[j]] = data.frame(Anno = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$threshold=="FDR < 0.05"),"variantName"])),
                                         length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$threshold!="FDR < 0.05"),"variantName"]))),
                                notAnno = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]=="no" & df$threshold=="FDR < 0.05"),"variantName"])),
                                            length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]=="no" & df$threshold!="FDR < 0.05"),"variantName"]))),
                                row.names = c("Sig","N.S."))
  }
names(prop[[i]]) = colnames(df)[12:21]
}
names(prop) = c("ByFraction","ByAge","ByFraction:Adult","ByFraction:Prenatal","ByAge:Cytoplasm","ByAge:Nucleus")
fisher.prop = lapply(prop, function(x) lapply(x, fisher.test))
write.csv(rbind(pvalue = data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value)))),
                data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$estimate))))),quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType.csv")

d = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType.csv")
d = reshape2::melt(d)
x = d[grep("pvalue", d$X),]
x$X = gsub("pvalue.","", x$X)
y = d[grep("odds", d$X),]
y$X = gsub(".odds ratio", "", y$X)
x$comp = paste(x$X, x$variable, sep = ".")
y$comp = paste(y$X, y$variable, sep = ".")
d = cbind(x, OR = y[match(x$comp, y$comp),"value"])
d$FDR = p.adjust(d$value, method="fdr")
colnames(d) = c("VariantType","Comparison","Pvalue","comp","OddsRatio","FDR")
d = d[,colnames(d)!="comp"]
write.csv(d,quote=F,file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType.csv")


### is there a relationship between proportion of sig vs non-sig and direction of expression in a variant type?

prop = list(list(),list(),list(),list(),list(),list())
for (i in (1:length(unique(df$comparison)))) {
  for (j in 1:10) {
    prop[[i]][[j]] = data.frame(Pos = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC>0 & df$threshold=="FDR < 0.05"),"variantName"])),
                                         length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC>0 & df$threshold!="FDR < 0.05"),"variantName"]))),
                                Neg = c(length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC<0 & df$threshold=="FDR < 0.05"),"variantName"])),
                                            length(unique(df[which(df$comparison==unique(df$comparison)[i] & df[,colnames(df)[11+j]]!="no" & df$LFC<0 & df$threshold!="FDR < 0.05"),"variantName"]))),
                                row.names = c("Sig","N.S."))
  }
  names(prop[[i]]) = colnames(df)[12:21]
}
names(prop) = c("ByFraction","ByAge","ByFraction:Adult","ByFraction:Prenatal","ByAge:Cytoplasm","ByAge:Nucleus")
fisher.prop = lapply(prop, function(x) lapply(x, fisher.test))
write.csv(rbind(pvalue = data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$p.value)))),
                data.frame(lapply(fisher.prop, function(x) unlist(lapply(x, function(y) y$estimate))))),quote=F,
          file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType_byLFC_Dir.csv")

d = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType_byLFC_Dir.csv")
d = reshape2::melt(d)
x = d[grep("pvalue", d$X),]
x$X = gsub("pvalue.","", x$X)
y = d[grep("odds", d$X),]
y$X = gsub(".odds ratio", "", y$X)
x$comp = paste(x$X, x$variable, sep = ".")
y$comp = paste(y$X, y$variable, sep = ".")
d = cbind(x, OR = y[match(x$comp, y$comp),"value"])
d$FDR = p.adjust(d$value, method="fdr")
colnames(d) = c("VariantType","Comparison","Pvalue","comp","OddsRatio","FDR")
d = d[,colnames(d)!="comp"]
write.csv(d,quote=F,file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/DSE_byComparison_byVariantType_byLFC_Dir.csv")






