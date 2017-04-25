library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames

# Filter introns
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings=="-"),])
clean = lapply(IRfiltered, function(x) grep("clean", x$GeneIntronDetails, fixed=T))
IRfiltered2 = list()
for (i in 1:length(IRfiltered)){
  cl = clean[[i]]
  irf = IRfiltered[[i]]
  IRfiltered2[[i]] = irf[cl,]}
names(IRfiltered2) = names(IRfiltered)
string = lapply(IRfiltered2, function(x) strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE))
string = lapply(string, function(x) unlist(x, recursive = FALSE))
y = lapply(string, function(x) grep("ENSG", x))
genes = list()
for (i in 1:length(string)){
  tmp = string[[i]]
  genes[[i]] = tmp[y[[i]]]}
names(genes) = names(string)
intronID = lapply(IRfiltered2, function(x) paste0(x$Chr,":",x$Start,"-",x$End))
IRfiltered2 = Map(cbind, IRfiltered2, genes = genes, intronID = intronID)
head(IRfiltered2[[1]]) # all introns passing QC per sample
irbyGene = lapply(IRfiltered2, function(x) data.table(x, key="genes"))
irbyGene = lapply(irbyGene,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="genes"])) # maximum IR for each gene with at least one intron passing QC

#Number of genes with introns passing filtering steps (not overlapping anything, enough coverage, etc.)
elementNROWS(irbyGene)
#Number of Genes in sig gene set
elementNROWS(sig)

#Get IR info for genes significantly regulated by fraction 
max = do.call(rbind, irbyGene)
max$SampleID = gsub("\\..*","", rownames(max))
max$Fraction = max$Age = "NA"
max[grep("C", max$SampleID), colnames(max)=="Fraction"] = "Cytosol"
max[grep("N", max$SampleID), colnames(max)=="Fraction"] = "Nucleus"
max[c(grep("1113", max$SampleID),grep("20", max$SampleID)), colnames(max)=="Age"] = "Adult"
max[grep("53", max$SampleID), colnames(max)=="Age"] = "Prenatal"
max$Group = paste(max$Age, max$Fraction, sep=":")
sigIR = lapply(sig, function(x) max[which(max$genes %in% x$ensID),])
lapply(sigIR, function(x) length(unique(x$genes)))

# genes >50% retained in at least one sample
perc = lapply(sigIR, function(x) list(sig.more=as.character(unique(x[which(x$IRratio>=0.50),colnames(x)=="genes"])),
                                            sig.less=as.character(unique(x[which(x$IRratio<0.50),colnames(x)=="genes"])),
                                            all.more=as.character(unique(max[which(max$IRratio>=0.50),colnames(max)=="genes"])),
                                            all.less=as.character(unique(max[which(max$IRratio<0.50),colnames(max)=="genes"]))))
names = names(perc)
for (i in 1:length(perc)){
venn.diagram(perc[[i]], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",names[i], ".50.percent.IRratio.jpeg"), 
             main=names[i], col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 50.percent.IRratio.overlap.pdf

#test whether genes with >50% intron retention in any sample are preferentially significantly DEG in the pattern listed
both_exported = data.frame(c(0,3+59), c(769,8486))
fisher.test(both_exported)
#p-value = 0.008925
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000000 0.6782648
#sample estimates:
#  odds ratio 
#0
both_retained = data.frame(c(16,3+43), c(474,8781))
fisher.test(both_retained)
#data:  both_retained
#p-value = 7.755e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.377414 11.692614
#sample estimates:
#  odds ratio 
#6.440472
Ad_exported =  data.frame(c(6,3+53), c(607,8648))
fisher.test(Ad_exported)
#data:  Ad_exported
#p-value = 0.2999
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5352187 3.5553318
#sample estimates:
#  odds ratio 
#1.526377
Ad_retained =  data.frame(c(1+27,2+32), c(1648,7607))
fisher.test(Ad_retained)
#data:  Ad_retained
#p-value = 6.878e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.213186 6.476791
#sample estimates:
#  odds ratio 
#3.800572
Fet_exported =  data.frame(c(0,3+59), c(83,9172))
fisher.test(Fet_exported)
#data:  Fet_exported
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.000000 6.922506
#sample estimates:
#  odds ratio 
#0
Fet_retained =  data.frame(c(1,3+58), c(113,9142))
fisher.test(Fet_retained)
#data:  Fet_retained
#p-value = 0.535
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.03275724 7.82218589
#sample estimates:
#  odds ratio 
#1.326225
ret_Ad_exp_Fet = data.frame(c(0,3+59), c(10,9245))
fisher.test(ret_Ad_exp_Fet)
#data:  ret_Ad_exp_Fet
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 68.0569
#sample estimates:
#  odds ratio 
#0
ret_Fet_exp_Ad = data.frame(c(0,3+59), c(1,9254))
fisher.test(ret_Fet_exp_Ad)
#data:  ret_Fet_exp_Ad
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00 5336.36
#sample estimates:
#  odds ratio 
#0
interacting = data.frame(c(9,3+50), c(858,8397))
fisher.test(interacting)
#data:  interacting
#p-value = 0.182
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7180923 3.4137314
#sample estimates:
#  odds ratio 
#1.661775

# Intron retention between retained and exported genes
t.test(sigIR[["both_retained"]][,"IRratio"], sigIR[["both_exported"]][,"IRratio"], alternative = "greater")
#data:  sigIR[["both_retained"]][, "IRratio"] and sigIR[["both_exported"]][, "IRratio"]
#t = 28.516, df = 4190, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.03505721        Inf
#sample estimates:
#  mean of x   mean of y 
#0.045755728 0.008552073
t.test(sigIR[["Ad_retained"]][,"IRratio"], sigIR[["Ad_exported"]][,"IRratio"], alternative = "greater")
#data:  sigIR[["Ad_retained"]][, "IRratio"] and sigIR[["Ad_exported"]][, "IRratio"]
#t = 17.033, df = 13952, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.01298768        Inf
#sample estimates:
#  mean of x  mean of y 
#0.03554051 0.02116449
t.test(sigIR[["Fet_retained"]][,"IRratio"], sigIR[["Fet_exported"]][,"IRratio"], alternative = "greater")
#data:  sigIR[["Fet_retained"]][, "IRratio"] and sigIR[["Fet_exported"]][, "IRratio"]
#t = 5.2315, df = 1434.5, p-value = 9.653e-08
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.008375712         Inf
#sample estimates:
#  mean of x  mean of y 
#0.03302164 0.02080116

# What does the distribution of IR Ratios by Sig group look like?
gr = c("Both Retained", "Both Exported", "Retained in Prenatal", "Retained in Adult",
       "Exported in Prenatal", "Exported in Adult", "Retained in Adult/\nExported in Prenatal",
       "Retained in Prenatal/\nExported in Adult", "Interaction")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density.pdf")
for (i in 1:length(sigIR)){
x = ggplot(sigIR[[i]], aes(x=IRratio)) +
  geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlim(0,0.15) +
  xlab("IR Ratio") +
  ggtitle(paste0("Intron Retention: ", gr[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
print(x)
}
dev.off() 

## look at these comparisons in genes DE by fraction in individual groups
# Get the IR ratio for Fraction genes
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres),
                Arres = data.frame(Arres),Frres = data.frame(Frres), 
                Fpres.down = data.frame(Fpres.down))
FracList = Map(cbind, FracList, lapply(FracList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
SigFracList = Map(cbind, SigFracList, Sign = Sign)
SigFracList = lapply(SigFracList, function(x) split(x, x$Sign))
SigList = unlist(SigFracList, recursive = F) 
lapply(SigList, head)
sigFracIRratio = lapply(SigList, function(x) max[which(max$genes %in% x$ensemblID),])
elementNROWS(sigFracIRratio)

# genes >50% retained in at least one sample
perc = lapply(sigFracIRratio, function(x) list(sig.more=as.character(unique(x[which(x$IRratio>=0.50),colnames(x)=="genes"])),
                                      sig.less=as.character(unique(x[which(x$IRratio<0.50),colnames(x)=="genes"])),
                                      all.more=as.character(unique(max[which(max$IRratio>=0.50),colnames(max)=="genes"])),
                                      all.less=as.character(unique(max[which(max$IRratio<0.50),colnames(max)=="genes"]))))
names = c("AP_Cytosolic","AP_Nuclear","FP_Cytosolic","FP_Nuclear","AR_Cytosolic","AR_Nuclear",
          "FR_Cytosolic","FR_Nuclear","FP_Downsampled_Cytosolic","FP_Downsampled_Nuclear")
for (i in 1:length(perc)){
  venn.diagram(perc[[i]], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",names[i], ".50.percent.IRratio.jpeg"), 
               main=names[i], col = "transparent", 
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold", 
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 50.percent.IRratio.overlap_DEG_byGroup_fraction.pdf

#test whether genes with >50% intron retention in any sample are preferentially significantly DEG in the pattern listed
AP_Cytosolic = data.frame(c(1,3+58), c(522,8733))
fisher.test(AP_Cytosolic)
#data:  AP_Cytosolic
#p-value = 0.2628
AP_Nuclear = data.frame(c(1+21,2+38), c(396,8859))
fisher.test(AP_Nuclear)
#data:  AP_Nuclear
#p-value = 8.333e-15
FP_Cytosolic =  data.frame(c(0,3+59), c(0,9255))
fisher.test(FP_Cytosolic)
#data:  FP_Cytosolic
#p-value = 1
FP_Nuclear =  data.frame(c(5,3+54), c(15,9240))
fisher.test(FP_Nuclear)
#data:  FP_Nuclear
#p-value = 1.59e-07
FP_Downsampled_Cytosolic =  data.frame(c(0,3+59), c(0,9255))
fisher.test(FP_Downsampled_Cytosolic)
#data:  FP_Downsampled_Cytosolic
#p-value = 1
FP_Downsampled_Nuclear =  data.frame(c(5,3+54), c(17,9238))
fisher.test(FP_Downsampled_Nuclear)
#data:  FP_Downsampled_Nuclear
#p-value = 2.673e-07
AR_Cytosolic =  data.frame(c(1,3+58), c(563,8692))
fisher.test(AR_Cytosolic)
#data:  AR_Cytosolic
#p-value = 0.1827
AR_Nuclear =  data.frame(c(6,3+53), c(183,9072))
fisher.test(AR_Nuclear)
#data:  AR_Nuclear
#p-value = 0.001541
FR_Cytosolic = data.frame(c(0,3+59), c(0,9255))
fisher.test(FR_Cytosolic)
#data:  FR_Cytosolic
#p-value = 1
FR_Nuclear = data.frame(c(0,3+59), c(7,9248))
fisher.test(FR_Nuclear)
#data:  FR_Nuclear
#p-value = 1

# Intron retention between retained and exported genes
t.test(sigFracIRratio[["Apres.UpNuc"]][,"IRratio"], sigFracIRratio[["Apres.DownNuc"]][,"IRratio"], alternative = "greater")
#data:  sigFracIRratio[["Apres.UpNuc"]][, "IRratio"] and sigFracIRratio[["Apres.DownNuc"]][, "IRratio"]
#t = 28.712, df = 3634, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.04987853        Inf
#sample estimates:
#  mean of x   mean of y 
#0.061425353 0.008514881
t.test(sigFracIRratio[["Arres.UpNuc"]][,"IRratio"], sigFracIRratio[["Arres.DownNuc"]][,"IRratio"], alternative = "greater")
#data:  sigFracIRratio[["Arres.UpNuc"]][, "IRratio"] and sigFracIRratio[["Arres.DownNuc"]][, "IRratio"]
#t = 14.163, df = 1406.2, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.0281229       Inf
#sample estimates:
#  mean of x  mean of y 
#0.04186320 0.01004225

# What does the distribution of IR Ratios by Sig group look like?
gr = c("Adult:PolyA:Cytosolic","Adult:PolyA:Nuclear","Prenatal:PolyA:Cytosolic","Prenatal:PolyA:Nuclear","Adult:RiboZero:Cytosolic","Adult:RiboZero:Nuclear",
               "Prenatal:RiboZero:Cytosolic","Prenatal:RiboZero:Nuclear","Prenatal:PolyA:Cytosolic\n(Downsampled)","Prenatal:PolyA:Nuclear\n(Downsampled)")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_DEG_byFraction.pdf")
for (i in 1:length(sigFracIRratio)){
  x = ggplot(sigFracIRratio[[i]], aes(x=IRratio)) +
    geom_density(aes(group=Group, colour=Group)) +
    ylab("") + 
    xlim(0,0.15) +
    xlab("IR Ratio") +
    ggtitle(paste0("Intron Retention: ", gr[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = c(0.8, 0.55)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(x)
}
dev.off() 

## look at these comparisons in genes DE by age in individual groups
# Get the IR ratio for developmental genes
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres),
               Crres = data.frame(Crres),Nrres = data.frame(Nrres), 
               Cpres.down = data.frame(Cpres.down))
AgeList = Map(cbind, AgeList, lapply(AgeList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigAgeList)
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpPrenatal", "DownPrenatal"))
SigAgeList = Map(cbind, SigAgeList, Sign = Sign)
SigAgeList = lapply(SigAgeList, function(x) split(x, x$Sign))
SigList = unlist(SigAgeList, recursive = F) 
lapply(SigList, head)
sigAgeIRratio = lapply(SigList, function(x) max[which(max$genes %in% x$ensemblID),])
elementNROWS(sigAgeIRratio)

# genes >50% retained in at least one sample
perc = lapply(sigAgeIRratio, function(x) list(sig.more=as.character(unique(x[which(x$IRratio>=0.50),colnames(x)=="genes"])),
                                               sig.less=as.character(unique(x[which(x$IRratio<0.50),colnames(x)=="genes"])),
                                               all.more=as.character(unique(max[which(max$IRratio>=0.50),colnames(max)=="genes"])),
                                               all.less=as.character(unique(max[which(max$IRratio<0.50),colnames(max)=="genes"]))))
names = c("CP_Increasing","CP_Decreasing","NP_Increasing","NP_Decreasing","CR_Increasing","CR_Decreasing",
          "NR_Increasing","NR_Decreasing","CP_Downsampled_Increasing","CP_Downsampled_Decreasing")
for (i in 1:length(perc)){
  venn.diagram(perc[[i]], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",names[i], ".50.percent.IRratio.jpeg"), 
               main=names[i], col = "transparent", 
               fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
               label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                             "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
               fontfamily = "Arial", fontface = "bold", 
               cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
               cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 50.percent.IRratio.overlap_DEG_byGroup_age.pdf

#test whether genes with >50% intron retention in any sample are preferentially significantly DEG in the pattern listed
CP_Increasing = data.frame(c(10,3+49), c(1788,7467))
fisher.test(CP_Increasing)
#data:  CP_Increasing
#p-value = 0.6292
CP_Decreasing = data.frame(c(1+17,2+42), c(2238,7017))
fisher.test(CP_Decreasing)
#data:  CP_Decreasing
#p-value = 0.3735
NP_Increasing =  data.frame(c(13,3+46), c(1624,7631))
fisher.test(NP_Increasing)
#data:  NP_Increasing
#p-value = 0.502
NP_Decreasing =  data.frame(c(10,3+49), c(1918,7337))
fisher.test(NP_Decreasing)
#data:  NP_Decreasing
#p-value = 0.4342
CP_Downsampled_Increasing =  data.frame(c(10,3+49), c(1863,7392))
fisher.test(CP_Downsampled_Increasing)
#data:  CP_Downsampled_Increasing
#p-value = 0.5257
CP_Downsampled_Decreasing =  data.frame(c(1+17,2+42), c(2149,7106))
fisher.test(CP_Downsampled_Decreasing)
#data:  CP_Downsampled_Decreasing
#p-value = 0.2912
CR_Increasing =  data.frame(c(11,3+48), c(1711,7544))
fisher.test(CR_Increasing)
#data:  CR_Increasing
#p-value = 1
CR_Decreasing =  data.frame(c(13,3+46), c(1581,7674))
fisher.test(CR_Decreasing)
#data:  CR_Decreasing
#p-value = 0.3991
NR_Increasing = data.frame(c(9,3+50), c(1252,8003))
fisher.test(NR_Increasing)
#data:  NR_Increasing
#p-value = 0.8516
NR_Decreasing = data.frame(c(7,3+52), c(1531,7724))
fisher.test(NR_Decreasing)
#data:  NR_Decreasing
#p-value = 0.3074

# Intron retention between retained and exported genes
t.test(sigAgeIRratio[["Cpres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#data:  sigAgeIRratio[["Cpres.UpPrenatal"]][, "IRratio"] and sigAgeIRratio[["Cpres.DownPrenatal"]][, "IRratio"]
#t = 9.5703, df = 29991, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.004292059         Inf
#sample estimates:
#  mean of x  mean of y 
#0.02396868 0.01878580
t.test(sigAgeIRratio[["Npres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Npres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#data:  sigAgeIRratio[["Npres.UpPrenatal"]][, "IRratio"] and sigAgeIRratio[["Npres.DownPrenatal"]][, "IRratio"]
#t = -4.1706, df = 23669, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.003316672          Inf
#sample estimates:
#  mean of x  mean of y 
#0.01959545 0.02197400
t.test(sigAgeIRratio[["Crres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Crres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#data:  sigAgeIRratio[["Crres.UpPrenatal"]][, "IRratio"] and sigAgeIRratio[["Crres.DownPrenatal"]][, "IRratio"]
#t = 3.5496, df = 23800, p-value = 0.0001933
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.001206793         Inf
#sample estimates:
#  mean of x  mean of y 
#0.02314145 0.02089245
t.test(sigAgeIRratio[["Nrres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Nrres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#data:  sigAgeIRratio[["Nrres.UpPrenatal"]][, "IRratio"] and sigAgeIRratio[["Nrres.DownPrenatal"]][, "IRratio"]
#t = -6.726, df = 17589, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.005533654          Inf
#sample estimates:
#  mean of x  mean of y 
#0.01898629 0.02343255 
t.test(sigAgeIRratio[["Cpres.down.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.down.DownPrenatal"]][,"IRratio"], alternative = "greater")
#data:  sigAgeIRratio[["Cpres.down.UpPrenatal"]][, "IRratio"] and sigAgeIRratio[["Cpres.down.DownPrenatal"]][, "IRratio"]
#t = 9.5899, df = 30259, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.004286707         Inf
#sample estimates:
#  mean of x  mean of y 
#0.02391771 0.01874349

# What does the distribution of IR Ratios by Sig group look like?
gr = c("Cytosol:PolyA:Increasing","Cytosol:PolyA:Decreasing","Nucleus:PolyA:Increasing","Nucleus:PolyA:Decreasing","Cytosol:RiboZero:Increasing","Cytosol:RiboZero:Decreasing",
       "Nucleus:RiboZero:Increasing","Nucleus:RiboZero:Decreasing","Cytosol:PolyA:Increasing\n(Downsampled)","Cytosol:PolyA:Decreasing\n(Downsampled)")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/IR_sig_group_density_DEG_byAge.pdf")
for (i in 1:length(sigAgeIRratio)){
  x = ggplot(sigAgeIRratio[[i]], aes(x=IRratio)) +
    geom_density(aes(group=Group, colour=Group)) +
    ylab("") + 
    xlim(0,0.15) +
    xlab("IR Ratio") +
    ggtitle(paste0("Intron Retention: ", gr[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    theme(legend.position = c(0.8, 0.55)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(x)
}
dev.off()