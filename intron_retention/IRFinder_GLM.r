library(GenomicRanges)
library(data.table)
library(ggplot2)
library(DESeq2)
source("./Dropbox/sorted_figures/IRfinder/PolyA/GLM/DESeq2Constructor.R")

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


### Load IRFinder Results
names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames


### Write results files including only the "clean" introns

lapply(lapply(IRres, function(x) unlist(strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE), recursive = F)), function(y) unique(y[seq(3,length(y),3)]))
# There are 6 types of introns: "clean", "anti-near", "anti-over", "known-exon+anti-near", "known-exon", and "known-exon+anti-near+anti-over"

IRclean = lapply(IRres, function(x) x[grep("clean", x$GeneIntronDetails, fixed=T),])
elementNROWS(IRclean) # 175475 of the 233535 measured introns don't have any overlapping antisense transcripts or serve as exons in other sense transcripts


## Filter out constitutively spliced introns by the four comparisons

IRclean = Map(cbind, IRclean, intronID = lapply(IRclean, function(x) paste0(x$GeneIntronDetails,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)))
IRratios = do.call(cbind, lapply(IRclean, function(x) x$IRratio))
rownames(IRratios) = IRclean$Br1113C1$intronID
adultIDs=rownames(IRratios[rowSums(IRratios[,-grep("53",colnames(IRratios))])!=0,])
prenatalIDs=rownames(IRratios[rowSums(IRratios[,grep("53",colnames(IRratios))])!=0,])
nucleusIDs=rownames(IRratios[rowSums(IRratios[,grep("N",colnames(IRratios))])!=0,])
cytosolIDs=rownames(IRratios[rowSums(IRratios[,grep("C",colnames(IRratios))])!=0,])
for (i in (1:length(IRclean))) { IRclean[[i]][,"intronID"] = as.character(IRclean[[i]][,"intronID"]) }

for (i in 1:length(IRclean)){
  if (names(IRclean)[i] %in% names(IRclean)[-grep("53",names(IRclean))]) {
  write.table(IRclean[[i]][which(IRclean[[i]][,"intronID"] %in% adultIDs),], 
              file=paste0(path,"PolyA/",names(IRclean)[i],"/IRFinder-IR-nondir-cleanIntrons-adultShared.txt"), row.names = F, quote = F, sep = "\t") }
  if (names(IRclean)[i] %in% names(IRclean)[grep("53",names(IRclean))]) {
    write.table(IRclean[[i]][which(IRclean[[i]][,"intronID"] %in% prenatalIDs),], 
                file=paste0(path,"PolyA/",names(IRclean)[i],"/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt"), row.names = F, quote = F, sep = "\t") }
  if (names(IRclean)[i] %in% names(IRclean)[grep("N",names(IRclean))]) {
    write.table(IRclean[[i]][which(IRclean[[i]][,"intronID"] %in% nucleusIDs),], 
                file=paste0(path,"PolyA/",names(IRclean)[i],"/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt"), row.names = F, quote = F, sep = "\t") }
  if (names(IRclean)[i] %in% names(IRclean)[grep("C",names(IRclean))]) {
    write.table(IRclean[[i]][which(IRclean[[i]][,"intronID"] %in% cytosolIDs),], 
                file=paste0(path,"PolyA/",names(IRclean)[i],"/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt"), row.names = F, quote = F, sep = "\t") }
}



### Create the file list and metadata file

pd = pd[pd$Library=="polyA",]
pd$SampleID = gsub( "_.*$", "", pd$SampleID)
pd = pd[,!colnames(pd) %in% c("leftReads","rightReads","bamFile")]
pd$Zone = factor(pd$Zone)
pd$Fetal = factor(pd$Fetal)

adultPaths = as.vector(paste0(path,"PolyA/",pd[pd$Fetal=="Adult","SampleID"],"/IRFinder-IR-nondir-cleanIntrons-adultShared.txt"))
fetalPaths = as.vector(paste0(path,"PolyA/",pd[pd$Fetal=="Prenatal","SampleID"],"/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt"))
cytPaths = as.vector(paste0(path,"PolyA/",pd[pd$Zone=="Cytosol","SampleID"],"/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt"))
nucPaths = as.vector(paste0(path,"PolyA/",pd[pd$Zone=="Nucleus","SampleID"],"/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt"))

adult = pd[pd$Fetal=="Adult",c("SampleID","Zone")]
rownames(adult) = NULL
fetal = pd[pd$Fetal=="Prenatal",c("SampleID","Zone")]
rownames(fetal) = NULL
nucleus = pd[pd$Zone=="Nucleus",c("SampleID","Fetal")]
rownames(nucleus) = NULL
cytosol = pd[pd$Zone=="Cytosol",c("SampleID","Fetal")]
rownames(cytosol) = NULL


## Prepare the DESeq objects

Ad.metaList=DESeqDataSetFromIRFinder(filePaths = adultPaths, designMatrix = adult, designFormula=~1)
Pren.metaList=DESeqDataSetFromIRFinder(filePaths = fetalPaths, designMatrix = fetal, designFormula=~1)
Nuc.metaList=DESeqDataSetFromIRFinder(filePaths = nucPaths, designMatrix=nucleus, designFormula=~1)
Cyt.metaList=DESeqDataSetFromIRFinder(filePaths = cytPaths, designMatrix=cytosol, designFormula=~1)

Ad.dds = Ad.metaList$DESeq2Object           # Extract DESeq2 Object with normalization factors ready
colData(Ad.dds)                             # Check design of matrix
Pren.dds = Pren.metaList$DESeq2Object
colData(Pren.dds)
Nuc.dds = Nuc.metaList$DESeq2Object
colData(Nuc.dds)
Cyt.dds = Cyt.metaList$DESeq2Object
colData(Cyt.dds)


## Run DESeq

design(Ad.dds) = ~ Zone + Zone:IRFinder     # Build a formula of GLM. Read below for more details. 
Ad.dds = DESeq(Ad.dds)                      # Estimate parameters and fit to model
resultsNames(Ad.dds)                        # Check the actual variable name assigned by DESeq2
design(Pren.dds) = ~ Zone + Zone:IRFinder
Pren.dds = DESeq(Pren.dds)
resultsNames(Pren.dds)
design(Cyt.dds) = ~ Fetal + Fetal:IRFinder
Cyt.dds = DESeq(Cyt.dds)
resultsNames(Cyt.dds)
design(Nuc.dds) = ~ Fetal + Fetal:IRFinder
Nuc.dds = DESeq(Nuc.dds)
resultsNames(Nuc.dds)


## Extract the results

Ad.testIR.ZoneCytosol = results(Ad.dds, name = "ZoneCytosol.IRFinderIR")
# This tests if the number of IR reads are significantly different from the sum of normal spliced reads and intronic reads, in the WT samples. 
# We might only be interested in the "log2FoldChange" column, instead of the significance.
# This is because "log2FoldChange" represents log2(number of intronic reads/the sum of normal spliced reads and intronic reads).
# It actually equals log2(IR ratio) in the WT samples.
# So the IR ratio for each intron in WT samples can be easily extracted by the following line
Ad.IRratio.ZoneCytosol = data.frame(intronIDs = rownames(Ad.testIR.ZoneCytosol), IRratio = 2^Ad.testIR.ZoneCytosol$log2FoldChange)
Ad.testIR.ZoneNucleus = results(Ad.dds, name = "ZoneNucleus.IRFinderIR")
Ad.IRratio.ZoneNucleus = data.frame(intronIDs = rownames(Ad.testIR.ZoneNucleus), IRratio = 2^Ad.testIR.ZoneNucleus$log2FoldChange)

Pren.testIR.ZoneCytosol = results(Pren.dds, name = "ZoneCytosol.IRFinderIR")
Pren.IRratio.ZoneCytosol = data.frame(intronIDs = rownames(Pren.testIR.ZoneCytosol), IRratio = 2^Pren.testIR.ZoneCytosol$log2FoldChange)    
Pren.testIR.ZoneNucleus = results(Pren.dds, name = "ZoneNucleus.IRFinderIR")
Pren.IRratio.ZoneNucleus = data.frame(intronIDs = rownames(Pren.testIR.ZoneNucleus), IRratio = 2^Pren.testIR.ZoneNucleus$log2FoldChange)    

Cyt.testIR.FetalAdult = results(Cyt.dds, name = "FetalAdult.IRFinderIR")
Cyt.IRratio.FetalAdult = data.frame(intronIDs = rownames(Cyt.testIR.FetalAdult), IRratio = 2^Cyt.testIR.FetalAdult$log2FoldChange)
Cyt.testIR.FetalPrenatal = results(Cyt.dds, name = "FetalPrenatal.IRFinderIR")
Cyt.IRratio.FetalPrenatal = data.frame(intronIDs = rownames(Cyt.testIR.FetalPrenatal), IRratio = 2^Cyt.testIR.FetalPrenatal$log2FoldChange)

Nuc.testIR.FetalAdult = results(Nuc.dds, name = "FetalAdult.IRFinderIR")
Nuc.IRratio.FetalAdult = data.frame(intronIDs = rownames(Nuc.testIR.FetalAdult), IRratio = 2^Nuc.testIR.FetalAdult$log2FoldChange)
Nuc.testIR.FetalPrenatal = results(Nuc.dds, name = "FetalPrenatal.IRFinderIR")
Nuc.IRratio.FetalPrenatal = data.frame(intronIDs = rownames(Nuc.testIR.FetalPrenatal), IRratio = 2^Nuc.testIR.FetalPrenatal$log2FoldChange)


Ad.res.diff = results(Ad.dds, contrast=list("ZoneNucleus.IRFinderIR","ZoneCytosol.IRFinderIR"))     
# This actually returns the Wald test result of each intron for differential IR analysis.
# This is because it is testing for if the two log2-transformed fold changes are different from each other of not.
# It equals to test if the log2(IRratio.WT) is the same as log2(IRratio.KO), which actually compares the IR ratio in two conditions.
Pren.res.diff = results(Pren.dds, contrast=list("ZoneNucleus.IRFinderIR","ZoneCytosol.IRFinderIR"))
Cyt.res.diff = results(Cyt.dds, contrast=list("FetalPrenatal.IRFinderIR","FetalAdult.IRFinderIR"))
Nuc.res.diff = results(Nuc.dds, contrast=list("FetalPrenatal.IRFinderIR","FetalAdult.IRFinderIR"))

dIR = list(byFractionInAdult = Ad.res.diff, byFractionInPrenatal = Pren.res.diff, byAgeInCytosol = Cyt.res.diff, byAgeInNucleus = Nuc.res.diff)
dIR = lapply(dIR, as.data.frame)
dIR = Map(cbind, dIR, intronID = lapply(dIR, rownames), 
          Ad.IRratio.ZoneCytosol = lapply(dIR, function(x) Ad.IRratio.ZoneCytosol[match(rownames(x),as.character(Ad.IRratio.ZoneCytosol$intronIDs)),"IRratio"]),
          Ad.IRratio.ZoneNucleus = lapply(dIR, function(x) Ad.IRratio.ZoneNucleus[match(rownames(x),as.character(Ad.IRratio.ZoneNucleus$intronIDs)),"IRratio"]),
          Pren.IRratio.ZoneCytosol = lapply(dIR, function(x) Pren.IRratio.ZoneCytosol[match(rownames(x),as.character(Pren.IRratio.ZoneCytosol$intronIDs)),"IRratio"]),
          Pren.IRratio.ZoneNucleus = lapply(dIR, function(x) Pren.IRratio.ZoneNucleus[match(rownames(x),as.character(Pren.IRratio.ZoneNucleus$intronIDs)),"IRratio"]),
          Cyt.IRratio.FetalAdult = lapply(dIR, function(x) Cyt.IRratio.FetalAdult[match(rownames(x),as.character(Cyt.IRratio.FetalAdult$intronIDs)),"IRratio"]),
          Cyt.IRratio.FetalPrenatal = lapply(dIR, function(x) Cyt.IRratio.FetalPrenatal[match(rownames(x),as.character(Cyt.IRratio.FetalPrenatal$intronIDs)),"IRratio"]),
          Nuc.IRratio.FetalAdult = lapply(dIR, function(x) Nuc.IRratio.FetalAdult[match(rownames(x),as.character(Nuc.IRratio.FetalAdult$intronIDs)),"IRratio"]),
          Nuc.IRratio.FetalPrenatal = lapply(dIR, function(x) Nuc.IRratio.FetalPrenatal[match(rownames(x),as.character(Nuc.IRratio.FetalPrenatal$intronIDs)),"IRratio"]))

save(dIR, pd, geneMap, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_results_object")



### Explore Differential retention results

load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_results_object")
elementNROWS(dIR)

IDs = lapply(dIR, function(x) unlist(strsplit(rownames(x), "/", fixed=T), recursive = F))
dIR = Map(cbind, dIR, ensID = lapply(IDs, function(x) x[grep("ENSG", x)]), coord = lapply(IDs, function(x) x[grep(":", x, fixed=T)]),
            geneSymbol = lapply(IDs, function(x) x[-c(grep("ENSG", x),grep(":", x, fixed=T),grep("clean", x))]),
          A.Sign = lapply(dIR, function(y) ifelse((y$Ad.IRratio.ZoneCytosol - y$Ad.IRratio.ZoneNucleus) < 0,"MoreIRInNuc", "MoreIRInCyt")),
          P.Sign = lapply(dIR, function(y) ifelse((y$Pren.IRratio.ZoneCytosol - y$Pren.IRratio.ZoneNucleus) < 0,"MoreIRInNuc", "MoreIRInCyt")),
          C.Sign = lapply(dIR, function(y) ifelse((y$Cyt.IRratio.FetalAdult - y$Cyt.IRratio.FetalPrenatal) < 0,"MoreIRInPrenatal", "MoreIRInAdult")),
          N.Sign = lapply(dIR, function(y) ifelse((y$Nuc.IRratio.FetalAdult - y$Nuc.IRratio.FetalPrenatal) < 0,"MoreIRInPrenatal", "MoreIRInAdult")))
head(dIR[[1]])
dIR = Map(cbind, dIR, Sign = list(ifelse(dIR$byFractionInAdult$A.Sign=="MoreIRInNuc", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byFractionInPrenatal$P.Sign=="MoreIRInNuc", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byAgeInCytosol$C.Sign=="MoreIRInPrenatal", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult"),
                                  ifelse(dIR$byAgeInNucleus$N.Sign=="MoreIRInPrenatal", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult")))
sigIR = lapply(dIR, function(x) x[which(x$padj<=0.05),])
elementNROWS(sigIR)
sigIRdir = unlist(lapply(sigIR, function(x) split(x, x$log2FoldChange>0)), recursive=F)



## Explore the results

counts = data.frame(rbind(elementNROWS(dIR),
                          elementNROWS(lapply(dIR, function(y) y[which(y$padj<=0.05),])),
                          elementNROWS(lapply(dIR, function(y) y[which(y$Sign=="MoreIRInNuc.Fetal"),])),
                          elementNROWS(lapply(dIR, function(y) y[which(y$Sign=="MoreIRINCyt.Adult"),])),
                          elementNROWS(lapply(dIR, function(y) y[which(y$padj<=0.05 & y$Sign=="MoreIRInNuc.Fetal"),])),
                          elementNROWS(lapply(dIR, function(y) y[which(y$padj<=0.05 & y$Sign=="MoreIRINCyt.Adult"),]))),
                    row.names = c("total measured", "FDR < 0.05", "MoreIRInNuc.Fetal", "MoreIRINCyt.Adult","FDR < 0.05 and MoreIRInNuc.Fetal","FDR < 0.05 and MoreIRINCyt.Adult"))
write.csv(counts, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_numIntrons_clean_nonconstitutivelySpliced.csv")



## Comparison of significantly vs nonsignificantly retained introns by zone/age

tables = list()
for (i in 1:4){
  tables[[i]] = data.frame(MoreIRInNuc.Fetal = c(counts["FDR < 0.05 and MoreIRInNuc.Fetal",i], counts["MoreIRInNuc.Fetal",i]-counts["FDR < 0.05 and MoreIRInNuc.Fetal",i]),
                           MoreIRINCyt.Adult = c(counts["FDR < 0.05 and MoreIRINCyt.Adult",i], counts["MoreIRINCyt.Adult",i]-counts["FDR < 0.05 and MoreIRINCyt.Adult",i]),
                           row.names = c("Sig","nonSig"))
}
names(tables) = colnames(counts)

write.csv(data.frame(rbind(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value)),unlist(lapply(lapply(tables, fisher.test), function(x) x$estimate))),row.names=c("p.value","OR")), 
          quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_sig_vs_nonSig_IR_byIntron_byFraction_byAge.csv")



# Get the DEG pval and LFC sign by Fraction for differentially retained introns
degs = list(Apres = data.frame(Apres, Sign = ifelse(Apres$log2FoldChange>0, "Pos","Neg")), Fpres = data.frame(Fpres.down, Sign = ifelse(Fpres.down$log2FoldChange>0, "Pos","Neg")), 
            Cpres = data.frame(Cpres.down), Npres = data.frame(Npres))
degs = Map(cbind, degs, lapply(degs, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
dIR = Map(cbind, dIR, 
             AP.sig = lapply(dIR, function(x) degs$Apres[match(x$ensID, degs$Apres$ensemblID),"padj"]),
             AP.LFC = lapply(dIR, function(x) degs$Apres[match(x$ensID, degs$Apres$ensemblID),"log2FoldChange"]),
             FP.sig = lapply(dIR, function(x) degs$Fpres[match(x$ensID, degs$Fpres$ensemblID),"padj"]),
             FP.LFC = lapply(dIR, function(x) degs$Fpres[match(x$ensID, degs$Fpres$ensemblID),"log2FoldChange"]),
             CP.sig = lapply(dIR, function(x) degs$Cpres[match(x$ensID, degs$Cpres$ensemblID),"padj"]),
             CP.LFC = lapply(dIR, function(x) degs$Cpres[match(x$ensID, degs$Cpres$ensemblID),"log2FoldChange"]),
             NP.sig = lapply(dIR, function(x) degs$Npres[match(x$ensID, degs$Npres$ensemblID),"padj"]),
             NP.LFC = lapply(dIR, function(x) degs$Npres[match(x$ensID, degs$Npres$ensemblID),"log2FoldChange"]))
lapply(dIR, head)



# Are genes with significantly differentially retained introns more likely to be significantly DEG by fraction?

P = c("AP", "FP", "CP", "NP")
tables = list(list(), list(), list(), list())
for (i in 1:length(dIR)){
  for (j in 1:length(P)) {
    tables[[i]][[j]] = data.frame(sigIR = c(nrow(dIR[[i]][which(dIR[[i]][,"padj"]<=0.05 & dIR[[i]][,paste0(P[j],".sig")]<=0.05),]),
                                            nrow(dIR[[i]][which(dIR[[i]][,"padj"]<=0.05 & dIR[[i]][,paste0(P[j],".sig")]>0.05),])),
                                  nonsigIR = c(nrow(dIR[[i]][which(dIR[[i]][,"padj"]>0.05 & dIR[[i]][,paste0(P[j],".sig")]<=0.05),]),
                                               nrow(dIR[[i]][which(dIR[[i]][,"padj"]>0.05 & dIR[[i]][,paste0(P[j],".sig")]>0.05),])), row.names = c("sigDEG","nonsigDEG"))
  }
  names(tables[[i]]) = paste0(P, "-DEG")
}
names(tables) = names(dIR)
fisher = lapply(tables, function(x) lapply(x, fisher.test))
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) lapply(x, function(y) y$p.value))),unlist(lapply(fisher, function(x) lapply(x, function(y) y$estimate)))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_sig_vs_nonSig_IR_byIntron_and_sig_vs_nonSig_DEG_byFraction_byAge.csv")



# Are genes with significantly differentially retained introns more likely to have LFC in one direction by fraction?

for (i in 1:length(dIR)){
  for (j in 1:length(P)) {
    tables[[i]][[j]] = data.frame(sigIR = c(nrow(dIR[[i]][which(dIR[[i]][,"padj"]<=0.05 & dIR[[i]][,paste0(P[j],".LFC")]>0),]),
                                            nrow(dIR[[i]][which(dIR[[i]][,"padj"]<=0.05 & dIR[[i]][,paste0(P[j],".LFC")]<0),])),
                                  nonsigIR = c(nrow(dIR[[i]][which(dIR[[i]][,"padj"]>0.05 & dIR[[i]][,paste0(P[j],".LFC")]>0),]),
                                               nrow(dIR[[i]][which(dIR[[i]][,"padj"]>0.05 & dIR[[i]][,paste0(P[j],".LFC")]<0),])), row.names = c("posLFC","negLFC"))
  }
  names(tables[[i]]) = paste0(P, "-DEG")
}
names(tables) = names(dIR)
fisher = lapply(tables, function(x) lapply(x, fisher.test))
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) lapply(x, function(y) y$p.value))),unlist(lapply(fisher, function(x) lapply(x, function(y) y$estimate)))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_sig_vs_nonSig_IR_byIntron_and_pos_vs_neg_geneLFC_byFraction_byAge.csv")



# Are significantly differentially retained introns by fraction more likely to have a higher IR ratio in nuclear RNA?

for (i in 1:length(dIR)){
  tables[[i]] = data.frame(sigIR = c(nrow(dIR[[i]][which(dIR[[i]][,"padj"]<=0.05 & dIR[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),]),
                                     nrow(dIR[[i]][which(dIR[[i]][,"padj"]<=0.05 & dIR[[i]][,"Sign"]=="MoreIRInCyt.Adult"),])),
                           nonsigIR = c(nrow(dIR[[i]][which(dIR[[i]][,"padj"]>0.05 & dIR[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),]),
                                        nrow(dIR[[i]][which(dIR[[i]][,"padj"]>0.05 & dIR[[i]][,"Sign"]=="MoreIRInCyt.Adult"),])), row.names = c("MoreIRInNuc.Fetal","MoreIRInCyt.Adult"))
}
names(tables) = names(dIR)
fisher = lapply(tables, fisher.test)
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) x$p.value)),unlist(lapply(fisher, function(x) x$estimate))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_sig_vs_nonSig_IR_byIntron_and_MoreIRInNuc.Fetal_byFraction_byAge.csv")



# IR regulated genes 

IRlist = unlist(lapply(dIR, function(x) split(x, x$Sign)), recursive = F)
names(IRlist) = names = c("Adult:Cytosol-Enriched","Adult:Nuclear-Enriched","Prenatal:Cytosol-Enriched","Prenatal:Nuclear-Enriched",
                          "Cytosol:Adult-Enriched","Cytosol:Prenatal-Enriched","Nucleus:Adult-Enriched","Nucleus:Prenatal-Enriched")
degs = Map(cbind, degs, Comparison = list("Adult","Prenatal","Cytosol","Nucleus"), Collapsed.Comparison = list("Fraction","Fraction","Age","Age"), 
           FDR = lapply(degs, function(x) ifelse(x$padj<=0.05, "FDR<0.05", "FDR>0.05")))
degs = do.call(rbind, degs)
degs = degs[which(degs$padj!="NA"),]
IRlist = lapply(IRlist, function(x) x[which(x$padj <= 0.05),])
IRlist = lapply(IRlist, function(x) degs[which(degs$ensemblID %in% x$ensID),])
elementNROWS(IRlist)


# Plot DEG by fraction LFC and FDR for genes that contain introns that are differentially retained

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_IRgenes_byFraction.pdf", width=8, height=8)
for (i in 1:length(IRlist)){
  g = ggplot(IRlist[[i]][which(IRlist[[i]][,"Collapsed.Comparison"]=="Fraction"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("RNA Localization by IR Ratio:\n",names[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
}
dev.off()

# Plot DEG by age LFC and FDR for genes that contain introns that are differentially retained

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_IRgenes_byAge.pdf", width=8, height=8)
for (i in 1:length(IRlist)){
  g = ggplot(IRlist[[i]][which(IRlist[[i]][,"Collapsed.Comparison"]=="Age"),], aes(x=Comparison, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("Age Expression Changes by IR Ratio:\n",names[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
}
dev.off()

## Are the Fraction or Age LFC values greater in genes containing a retained intron differentially?

t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"])
t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"])
t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"])
t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"])
tables = list(t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"]),
              t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"]),
              t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"]),
              t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Adult"),"log2FoldChange"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Prenatal"),"log2FoldChange"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"log2FoldChange"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"log2FoldChange"]))
names(tables) = c(paste0("byAgeInCytosol_", unique(degs$Comparison), "_LFC"), paste0("byAgeInNucleus_", unique(degs$Comparison), "_LFC"))
write.csv(data.frame(rbind(unlist(lapply(tables, function(x) x$statistic)),unlist(lapply(tables, function(x) x$p.value)),data.frame(lapply(tables, function(x) x$conf.int)),
                           data.frame(lapply(tables, function(x) x$estimate))), row.names=c("Tstat","p.value","conf.int1","conf.int2","mean of x","mean of y")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_ttest_nuclearIntrons_byAge_geneLFC_diff.csv")



## Are the Fraction or Age p.values lower in genes containing a retained intron differentially?

t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"padj"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"padj"])
t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"padj"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"padj"])
t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Cytosol"),"padj"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Cytosol"),"padj"])
t.test(x = IRlist[["Adult:Cytosol-Enriched"]][which(IRlist[["Adult:Cytosol-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
       y = IRlist[["Adult:Nuclear-Enriched"]][which(IRlist[["Adult:Nuclear-Enriched"]][,"Comparison"]=="Nucleus"),"padj"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Adult"),"padj"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Adult"),"padj"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Prenatal"),"padj"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Prenatal"),"padj"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Cytosol"),"padj"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Cytosol"),"padj"])
t.test(x = IRlist[["Prenatal:Cytosol-Enriched"]][which(IRlist[["Prenatal:Cytosol-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
       y = IRlist[["Prenatal:Nuclear-Enriched"]][which(IRlist[["Prenatal:Nuclear-Enriched"]][,"Comparison"]=="Nucleus"),"padj"])
tables = list(t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Adult"),"padj"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Adult"),"padj"]),
              t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Prenatal"),"padj"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Prenatal"),"padj"]),
              t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"padj"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"padj"]),
              t.test(x = IRlist[["Cytosol:Adult-Enriched"]][which(IRlist[["Cytosol:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
                     y = IRlist[["Cytosol:Prenatal-Enriched"]][which(IRlist[["Cytosol:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"padj"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Adult"),"padj"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Adult"),"padj"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Prenatal"),"padj"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Prenatal"),"padj"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Cytosol"),"padj"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Cytosol"),"padj"]),
              t.test(x = IRlist[["Nucleus:Adult-Enriched"]][which(IRlist[["Nucleus:Adult-Enriched"]][,"Comparison"]=="Nucleus"),"padj"],
                     y = IRlist[["Nucleus:Prenatal-Enriched"]][which(IRlist[["Nucleus:Prenatal-Enriched"]][,"Comparison"]=="Nucleus"),"padj"]))
names(tables) = c(paste0("byAgeInCytosol_", unique(degs$Comparison), "_padj"),paste0("byAgeInNucleus_", unique(degs$Comparison), "_padj"))
write.csv(data.frame(rbind(unlist(lapply(tables, function(x) x$statistic)),unlist(lapply(tables, function(x) x$p.value)),data.frame(lapply(tables, function(x) x$conf.int)),
                           data.frame(lapply(tables, function(x) x$estimate))), row.names=c("Tstat","p.value","conf.int1","conf.int2","mean of x","mean of y")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_ttest_nuclearIntrons_byAge_padj_diff.csv")



# plot IR ratio in cytosolic prenatal samples and LFC by age in cytosol
ggplot(dIR$Fetal_byFraction[which(dIR$Fetal_byFraction$padj<0.05),], aes(x=IR.diff, y=CP.LFC)) + 
  geom_point() +
  ylab("Log2 Fold Change") + ylim(-10,10) + 
  xlab("Difference in IR Ratio") +
  ggtitle("Age LFC and IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))


# Are genes with greater IR in cytosol in prenatal higher- or lower-expressed in prenatal than adult in cytosol?

tables = list()
for (i in 1:length(dIR)){
  tables[[i]] = list(inAdult = t.test(dIR[[i]][which(dIR[[i]][,"Sign"]!="MoreIRInNuc.Fetal"),"AP.LFC"],
                                      dIR[[i]][which(dIR[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"AP.LFC"]),
                     inPrenatal = t.test(dIR[[i]][which(dIR[[i]][,"Sign"]!="MoreIRInNuc.Fetal"),"FP.LFC"],
                                         dIR[[i]][which(dIR[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"FP.LFC"]),
                     inCytosol = t.test(dIR[[i]][which(dIR[[i]][,"Sign"]!="MoreIRInNuc.Fetal"),"CP.LFC"],
                                        dIR[[i]][which(dIR[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"CP.LFC"]),
                     inNucleus = t.test(dIR[[i]][which(dIR[[i]][,"Sign"]!="MoreIRInNuc.Fetal"),"NP.LFC"],
                                        dIR[[i]][which(dIR[[i]][,"Sign"]=="MoreIRInNuc.Fetal"),"NP.LFC"]))
}
names(tables) = names(dIR)
write.csv(data.frame(rbind(unlist(lapply(tables, function(x) lapply(x, function(y) y$statistic))),
                           unlist(lapply(tables, function(x) lapply(x, function(y) y$p.value))),
                           data.frame(lapply(tables, function(x) lapply(x, function(y) y$conf.int))),
                           data.frame(lapply(tables, function(x) lapply(x, function(y) y$estimate)))), 
                     row.names=c("Tstat","p.value","conf.int1","conf.int2","mean of x","mean of y")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_ttest_MoreIRInNuc.Fetal_IR_vs_geneLFC_inGroup.csv")



## Are dIR introns by age more likely to be dIR introns by fraction?

sigIR.intron = lapply(dIR, function(x) as.character(x[which(x$padj<=0.05),"intronID"])) 
nonsigIR.intron = lapply(dIR, function(x) as.character(x[which(x$padj>0.05),"intronID"]))
sigIR.intron.combined = list("Significantly IR\nBy Fraction" = c(sigIR.intron$byFractionInAdult, sigIR.intron$byFractionInPrenatal),
                             "Significantly IR\nBy Age" = c(sigIR.intron$byAgeInCytosol, sigIR.intron$byAgeInNucleus))
nonsigIR.intron.combined = list("Non-Significantly IR\nBy Fraction" = c(nonsigIR.intron$byFractionInAdult, nonsigIR.intron$byFractionInPrenatal),
                                "Non-Significantly IR\nBy Age" = c(nonsigIR.intron$byAgeInCytosol, nonsigIR.intron$byAgeInNucleus))
comps = list(c(1,3), c(1,4), c(1:2), c(2,3), c(2,4), c(3,4))
path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_dIR_"
ids = c("FractionbyAge_adult_cytosol", "FractionbyAge_adult_nucleus", "Fraction_adult_prenatal",
        "FractionbyAge_prenatal_cytosol","FractionbyAge_prenatal_nucleus","Age_cytosol_nucleus")

for (i in 1:length(comps)) {
  venn.diagram(c(sigIR.intron[comps[[i]]], nonsigIR.intron[comps[[i]]]), paste0(path, ids[i], ".jpeg"), 
               main=paste0("dIR_", ids[i]), col = "transparent",
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold",
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
}
venn.diagram(c(sigIR.intron.combined, nonsigIR.intron.combined), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_dIR_FractionbyAge_combined.jpeg", 
             main="dIR_FractionbyAge_combined", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
# all 7 saved together at ./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_overlap_Fraction_Age.pdf


tables = list()
for (i in 1:length(comps)) {
  tables[[i]] = data.frame(sig1 = c(length(sigIR.intron[[comps[[i]][1]]][sigIR.intron[[comps[[i]][1]]] %in% sigIR.intron[[comps[[i]][2]]]]),
                                    length(sigIR.intron[[comps[[i]][1]]][sigIR.intron[[comps[[i]][1]]] %in% nonsigIR.intron[[comps[[i]][2]]]])),
                           nonsig1 = c(length(nonsigIR.intron[[comps[[i]][1]]][nonsigIR.intron[[comps[[i]][1]]] %in% sigIR.intron[[comps[[i]][2]]]]),
                                       length(nonsigIR.intron[[comps[[i]][1]]][nonsigIR.intron[[comps[[i]][1]]] %in% nonsigIR.intron[[comps[[i]][2]]]])),row.names = c("sig2","nonsig2"))
}
names(tables) = ids
fisher = c(lapply(tables, fisher.test), list(combined_fraction_age = fisher.test(data.frame(
  sig1 = c(length(sigIR.intron.combined[[1]][sigIR.intron.combined[[1]] %in% sigIR.intron.combined[[2]]]),
           length(sigIR.intron.combined[[1]][sigIR.intron.combined[[1]] %in% nonsigIR.intron.combined[[2]]])),
  nonsig1 = c(length(nonsigIR.intron.combined[[1]][nonsigIR.intron.combined[[1]] %in% sigIR.intron.combined[[2]]]),
              length(nonsigIR.intron.combined[[1]][nonsigIR.intron.combined[[1]] %in% nonsigIR.intron.combined[[2]]])),
  row.names = c("sig2","nonsig2")))))
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) x$p.value)),unlist(lapply(fisher, function(x) x$estimate))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_overlap_sig_vs_nonSig_IR_byIntron_byGroup.csv")



## Are genes that have a dIR intron by age more likely to have a dIR intron by fraction?

sigIR.gene = lapply(dIR, function(x) unique(as.character(x[which(x$padj<=0.05),"ensID"])))
names(sigIR.gene) = paste0(names(sigIR.gene),"\nsig")
nonsigIR.gene = lapply(dIR, function(x) unique(as.character(x[which(x$padj>0.05),"ensID"])))
names(nonsigIR.gene) = paste0(names(nonsigIR.gene),"\nnonsig")
sigIR.gene.combined = list("Significantly IR\nBy Fraction" = c(sigIR.gene[["byFractionInAdult\nsig"]], sigIR.gene[["byFractionInPrenatal\nsig"]]),
                           "Significantly IR\nBy Age" = c(sigIR.gene[["byAgeInCytosol\nsig"]], sigIR.gene[["byAgeInNucleus\nsig"]]))
nonsigIR.gene.combined = list("Non-Significantly IR\nBy Fraction" = c(nonsigIR.gene[["byFractionInAdult\nnonsig"]], nonsigIR.gene[["byFractionInPrenatal\nnonsig"]]),
                              "Non-Significantly IR\nBy Age" = c(nonsigIR.gene[["byAgeInCytosol\nnonsig"]], nonsigIR.gene[["byAgeInNucleus\nnonsig"]]))

comps = list(c(1,3), c(1,4), c(1:2), c(2,3), c(2,4), c(3,4))
path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_dIR_gene_"
ids = c("FractionbyAge_adult_cytosol", "FractionbyAge_adult_nucleus", "Fraction_adult_prenatal",
        "FractionbyAge_prenatal_cytosol","FractionbyAge_prenatal_nucleus","Age_cytosol_nucleus")

for (i in 1:length(comps)) {
  venn.diagram(c(sigIR.intron[comps[[i]]], nonsigIR.intron[comps[[i]]]), paste0(path, ids[i], ".jpeg"), 
               main=paste0("dIR_", ids[i]), col = "transparent",
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold",
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
}
venn.diagram(c(sigIR.intron.combined, nonsigIR.intron.combined), 
             "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_FractionbyAge_combined.jpeg", 
             main="dIR_FractionbyAge_combined", col = "transparent",
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white",
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold",
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
# all 7 saved together at ./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/dIR_byGene_overlap_Fraction_Age.pdf


tables = list()
for (i in 1:length(comps)) {
  tables[[i]] = data.frame(sig1 = c(length(sigIR.gene[[comps[[i]][1]]][sigIR.gene[[comps[[i]][1]]] %in% sigIR.gene[[comps[[i]][2]]]]),
                                    length(sigIR.gene[[comps[[i]][1]]][sigIR.gene[[comps[[i]][1]]] %in% nonsigIR.gene[[comps[[i]][2]]]])),
                           nonsig1 = c(length(nonsigIR.gene[[comps[[i]][1]]][nonsigIR.gene[[comps[[i]][1]]] %in% sigIR.gene[[comps[[i]][2]]]]),
                                       length(nonsigIR.gene[[comps[[i]][1]]][nonsigIR.gene[[comps[[i]][1]]] %in% nonsigIR.gene[[comps[[i]][2]]]])),row.names = c("sig2","nonsig2"))
}
names(tables) = ids
fisher = c(lapply(tables, fisher.test), list(combined_fraction_age = fisher.test(data.frame(
  sig1 = c(length(sigIR.gene.combined[[1]][sigIR.gene.combined[[1]] %in% sigIR.gene.combined[[2]]]),
           length(sigIR.gene.combined[[1]][sigIR.gene.combined[[1]] %in% nonsigIR.gene.combined[[2]]])),
  nonsig1 = c(length(nonsigIR.gene.combined[[1]][nonsigIR.gene.combined[[1]] %in% sigIR.gene.combined[[2]]]),
              length(nonsigIR.gene.combined[[1]][nonsigIR.gene.combined[[1]] %in% nonsigIR.gene.combined[[2]]])),
  row.names = c("sig2","nonsig2")))))
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) x$p.value)),unlist(lapply(fisher, function(x) x$estimate))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_overlap_genesContaining_sig_vs_nonSig_IR_byGroup.csv")



## What about direction of retention in shared introns?

sigIR.dir = lapply(dIR, function(x) split(x, x$Sign))
sigIR.dir = unlist(lapply(sigIR.dir, function(x) lapply(x, function(y) unique(as.character(y[which(y$padj<=0.05),"intronID"])))), recursive = F)
sigIR.dir.combined = list("Significantly IR\nBy Fraction\nMore in Cytosol" = c(sigIR.dir$byFractionInAdult.MoreIRINCyt.Adult, sigIR.dir$byFractionInPrenatal.MoreIRINCyt.Adult),
                          "Significantly IR\nBy Fraction\nMore in Nucleus" = c(sigIR.dir$byFractionInAdult.MoreIRInNuc.Fetal, sigIR.dir$byFractionInPrenatal.MoreIRInNuc.Fetal),
                          "Significantly IR\nBy Age\nMore in Adult" = c(sigIR.dir$byAgeInCytosol.MoreIRInCyt.Adult, sigIR.dir$byAgeInNucleus.MoreIRINCyt.Adult),
                          "Significantly IR\nBy Age\nMore in Prenatal" = c(sigIR.dir$byAgeInCytosol.MoreIRInNuc.Fetal, sigIR.dir$byAgeInNucleus.MoreIRInNuc.Fetal))

comps = list(c(1,3), c(1,4), c(1,2), c(2,3), c(2,4), c(3,4))
nf = list(c(2,6), c(2,8), c(2,4), c(4,6), c(4,8), c(6,8))
ac = list(c(1,5), c(1,7), c(1,3), c(3,5), c(3,7), c(5,7))
path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/GLM_IRdir_"
ids = c("FractionbyAge_adult_cytosol", "FractionbyAge_adult_nucleus", "Fraction_adult_prenatal",
        "FractionbyAge_prenatal_cytosol","FractionbyAge_prenatal_nucleus","Age_cytosol_nucleus")

both = list()
for (i in 1:length(comps)) {
  both[[i]] = sigIR.intron[[comps[[i]][1]]][sigIR.intron[[comps[[i]][1]]] %in% sigIR.intron[[comps[[i]][2]]]]
}
both = c(both, list())
for (i in 1:length(nf)) {
  tables[[i]] = data.frame(MoreIRInNuc.Fetal1 = c(length(sigIR.dir[[nf[[i]][1]]][(sigIR.dir[[nf[[i]][1]]] %in% sigIR.dir[[nf[[i]][2]]]) & (sigIR.dir[[nf[[i]][1]]] %in% both[[i]])]),
                                                  length(sigIR.dir[[nf[[i]][1]]][sigIR.dir[[nf[[i]][1]]] %in% sigIR.dir[[ac[[i]][2]]] & (sigIR.dir[[nf[[i]][1]]] %in% both[[i]])])),
                           MoreIRInCyt.Adult1 = c(length(sigIR.dir[[ac[[i]][1]]][sigIR.dir[[ac[[i]][1]]] %in% sigIR.dir[[nf[[i]][2]]] & (sigIR.dir[[ac[[i]][1]]] %in% both[[i]])]),
                                                  length(sigIR.dir[[ac[[i]][1]]][sigIR.dir[[ac[[i]][1]]] %in% sigIR.dir[[ac[[i]][2]]] & (sigIR.dir[[ac[[i]][1]]] %in% both[[i]])])),
                           row.names = c("MoreIRInNuc.Fetal2","MoreIRInCyt.Adult2"))
}
names(tables) = names(both) = ids
fisher = c(lapply(tables, fisher.test), 
           list(combined_fraction_age = fisher.test(data.frame(sig1 = c(length(sigIR.dir.combined[[1]][sigIR.dir.combined[[1]] %in% sigIR.dir.combined[[3]]]),
                                                                        length(sigIR.dir.combined[[1]][sigIR.dir.combined[[1]] %in% sigIR.dir.combined[[4]]])),
                                                               nonsig1 = c(length(sigIR.dir.combined[[2]][sigIR.dir.combined[[2]] %in% sigIR.dir.combined[[3]]]),
                                                                           length(sigIR.dir.combined[[2]][sigIR.dir.combined[[2]] %in% sigIR.dir.combined[[4]]])),
                                                               row.names = c("sig2","nonsig2")))))
write.csv(data.frame(rbind(unlist(lapply(fisher, function(x) x$p.value)),unlist(lapply(fisher, function(x) x$estimate))),
                     row.names=c("p.value","OR")), quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GLM_fisher_overlap_sigIntrons_byDirection_byGroup.csv")

