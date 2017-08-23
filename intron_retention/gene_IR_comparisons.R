library(GenomicRanges)
library(data.table)
library(VennDiagram)
library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Load IRFinder Results
names = scan("./Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames

# Filter introns
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings=="-" | x$Warnings=="MinorIsoform"),])
IRfiltered = lapply(IRfiltered, function(x) x[grep("clean", x$GeneIntronDetails, fixed=T),])
string = lapply(IRfiltered, function(x) strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE))
string = lapply(string, function(x) unlist(x, recursive = F))
IRfiltered = Map(cbind, IRfiltered, genes = lapply(string, function(x) x[grep("ENSG", x)]), 
                 intronID = lapply(IRfiltered, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
head(IRfiltered[[1]]) # all introns passing QC per sample
irbyGene = lapply(IRfiltered, function(x) data.table(x, key="genes"))
irbyGene = lapply(irbyGene,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="genes"])) # maximum IR for each gene with at least one intron passing QC

#Number of genes with introns passing filtering steps (not overlapping anything, enough coverage, etc.)
elementNROWS(irbyGene)
#Number of Genes in sig gene set
elementNROWS(sig)

#Get IR info for genes significantly regulated by fraction 
max = do.call(rbind, irbyGene)
max$SampleID = gsub("\\..*","", rownames(max))
max$rownum = c(1:nrow(max))
max$Fraction = ifelse(max$rownum %in% grep("C", max$SampleID), "Cytosol", "Nucleus")
max$Age = ifelse(max$rownum %in% grep("Br53", max$SampleID), "Prenatal", "Adult")
max$Group = paste(max$Age, max$Fraction, sep=":")
sigIR = lapply(sig, function(x) max[which(max$genes %in% x$ensID),])
unlist(lapply(sigIR, function(x) length(unique(x$genes))))

# genes >50% retained in at least one sample
perc = lapply(sigIR, function(x) list(sig.more=as.character(unique(x[which(x$IRratio>=0.50),"genes"])),
                                      sig.less=as.character(unique(x[which(x$IRratio<0.50),"genes"])),
                                      all.more=as.character(unique(max[which(max$IRratio>=0.50),"genes"])),
                                      all.less=as.character(unique(max[which(max$IRratio<0.50),"genes"]))))
for (i in 1:length(perc)){
venn.diagram(perc[[i]], paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/",
                               names(perc)[i], ".50.percent.IRratio.jpeg"), 
             main=names(perc)[i], col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"), alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
} # Saved as pdfs together as 50.percent.IRratio.overlap.pdf

#test whether genes with >50% intron retention in any sample are preferentially significantly DEG in the pattern listed
both_exported = data.frame(c(0,4+88), c(786,8587))
fisher.test(both_exported)
#data:  both_exported
#p-value = 0.0008095
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000000 0.4481183
#sample estimates:
#  odds ratio 
#0
both_retained = data.frame(c(22,4+66), c(472,8901))
fisher.test(both_retained)
#data:  both_retained
#p-value = 1.17e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.461073 9.783295
#sample estimates:
#  odds ratio 
#5.924431
Ad_exported =  data.frame(c(8,4+80), c(613,8760))
fisher.test(Ad_exported)
#data:  Ad_exported
#p-value = 0.3935
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5663101 2.8241444
#sample estimates:
#  odds ratio 
#1.360929
Ad_retained =  data.frame(c(2+39,2+49), c(1652,7721))
fisher.test(Ad_retained)
#data:  Ad_retained
#p-value = 2.549e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.419184 5.801455
#sample estimates:
#  odds ratio 
#3.756583
Fet_exported =  data.frame(c(1,4+87), c(86,9287))
fisher.test(Fet_exported)
#data:  Fet_exported
#p-value = 0.5742
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.02938441 6.95788584
#sample estimates:
#  odds ratio 
#1.18666
Fet_retained =  data.frame(c(4,4+84), c(111,9262))
fisher.test(Fet_retained)
#data:  Fet_retained
#p-value = 0.02552
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9934781 10.3332631
#sample estimates:
#  odds ratio 
#3.791818
ret_Ad_exp_Fet = data.frame(c(0,4+88), c(10,9363))
fisher.test(ret_Ad_exp_Fet)
#data:  ret_Ad_exp_Fet
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00000 46.13797
#sample estimates:
#  odds ratio 
#0
ret_Fet_exp_Ad = data.frame(c(0,4+88), c(2,9371))
fisher.test(ret_Fet_exp_Ad)
#data:  ret_Fet_exp_Ad
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 540.6667
#sample estimates:
#  odds ratio 
#0
interacting = data.frame(c(1+11,3+77), c(870,8503))
fisher.test(interacting)
#data:  interacting
#p-value = 0.2072
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7242188 2.7183980
#sample estimates:
#  odds ratio 
#1.465956

# Intron retention between retained and exported genes
t.test(sigIR[["both_retained"]][,"IRratio"], sigIR[["both_exported"]][,"IRratio"], alternative = "greater")
#t = 31.185, df = 4210.3, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.04638024        Inf
#sample estimates:
#  mean of x   mean of y 
#0.058140973 0.009177554
t.test(sigIR[["Ad_retained"]][,"IRratio"], sigIR[["Ad_exported"]][,"IRratio"], alternative = "greater")
#t = 18.93, df = 14659, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.01632906        Inf
#sample estimates:
#  mean of x  mean of y 
#0.04133407 0.02345103
t.test(sigIR[["Fet_retained"]][,"IRratio"], sigIR[["Fet_exported"]][,"IRratio"], alternative = "greater")
#t = 5.812, df = 1513.6, p-value = 3.758e-09
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.01204279        Inf
#sample estimates:
#  mean of x  mean of y 
#0.04109963 0.02429922

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
SigFracList = Map(cbind, SigFracList, Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc")))
SigFracList = lapply(SigFracList, function(x) split(x, x$Sign))
SigList = unlist(SigFracList, recursive = F) 
lapply(SigList, head)
sigFracIRratio = lapply(SigList, function(x) max[which(max$genes %in% x$ensemblID),])
elementNROWS(sigFracIRratio)

# genes >50% retained in at least one sample
perc = lapply(sigFracIRratio, function(x) list(sig.more=as.character(unique(x[which(x$IRratio>=0.50),"genes"])),
                                      sig.less=as.character(unique(x[which(x$IRratio<0.50),"genes"])),
                                      all.more=as.character(unique(max[which(max$IRratio>=0.50),"genes"])),
                                      all.less=as.character(unique(max[which(max$IRratio<0.50),"genes"]))))
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
AP_Cytosolic = data.frame(c(1,4+87), c(536,8837))
fisher.test(AP_Cytosolic)
#data:  AP_Cytosolic
#p-value = 0.06488
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.004532804 1.039777991
#sample estimates:
#  odds ratio 
#0.1811917
AP_Nuclear = data.frame(c(1+28,3+60), c(394,8979))
fisher.test(AP_Nuclear)
#data:  AP_Nuclear
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  6.431021 16.737862
#sample estimates:
#  odds ratio 
#10.48313
FP_Cytosolic_down =  data.frame(c(0,4+88), c(0,9373))
fisher.test(FP_Cytosolic)
#data:  FP_Cytosolic_down
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
FP_Nuclear_down =  data.frame(c(6,4+82), c(16,9357))
fisher.test(FP_Nuclear)
#data:  FP_Nuclear_down
#p-value = 4.709e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  12.72757 112.57973
#sample estimates:
#  odds ratio 
#40.70501
AR_Cytosolic =  data.frame(c(1,4+87), c(580,8793))
fisher.test(AR_Cytosolic)
#data:  AR_Cytosolic
#p-value = 0.04505
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.004168995 0.955798954
#sample estimates:
#  odds ratio 
#0.1666128 
AR_Nuclear =  data.frame(c(8,4+80), c(182,9191))
fisher.test(AR_Nuclear)
#data:  AR_Nuclear
#p-value = 0.0005019
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.980513 10.112449
#sample estimates:
#  odds ratio 
#4.808445 
FR_Cytosolic = data.frame(c(0,4+88), c(0,9373))
fisher.test(FR_Cytosolic)
#data:  FR_Cytosolic
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0 Inf
#sample estimates:
#  odds ratio 
#0
FR_Nuclear = data.frame(c(0,4+88), c(7,9366))
fisher.test(FR_Nuclear)
#data:  FR_Nuclear
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00000 71.59552
#sample estimates:
#  odds ratio 
#0

# Intron retention between retained and exported genes
t.test(sigFracIRratio[["Apres.UpNuc"]][,"IRratio"], sigFracIRratio[["Apres.DownNuc"]][,"IRratio"], alternative = "greater")
#t = 31.222, df = 3666.8, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.06045412        Inf
#sample estimates:
#  mean of x   mean of y 
#0.072864367 0.009047333
t.test(sigFracIRratio[["Arres.UpNuc"]][,"IRratio"], sigFracIRratio[["Arres.DownNuc"]][,"IRratio"], alternative = "greater")
#t = 15.337, df = 1408.4, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.03620151        Inf
#sample estimates:
#  mean of x  mean of y 
#0.05187191 0.01131835

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
SigAgeList = Map(cbind, SigAgeList, Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpPrenatal", "DownPrenatal")))
SigAgeList = lapply(SigAgeList, function(x) split(x, x$Sign))
SigList = unlist(SigAgeList, recursive = F) 
lapply(SigList, head)
sigAgeIRratio = lapply(SigList, function(x) max[which(max$genes %in% x$ensemblID),])
elementNROWS(sigAgeIRratio)

# genes >50% retained in at least one sample
perc = lapply(sigAgeIRratio, function(x) list(sig.more=as.character(unique(x[which(x$IRratio>=0.50),"genes"])),
                                              sig.less=as.character(unique(x[which(x$IRratio<0.50),"genes"])),
                                              all.more=as.character(unique(max[which(max$IRratio>=0.50),"genes"])),
                                              all.less=as.character(unique(max[which(max$IRratio<0.50),"genes"]))))
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
NP_Increasing =  data.frame(c(20,4+68), c(1641,7732))
fisher.test(NP_Increasing)
#data:  NP_Increasing
#p-value = 0.2728
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7529227 2.1798240
#sample estimates:
#  odds ratio 
#1.308781
NP_Decreasing =  data.frame(c(17,4+71), c(1955,7418))
fisher.test(NP_Decreasing)
#data:  NP_Decreasing
#p-value = 0.6986
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4749613 1.4745719
#sample estimates:
#  odds ratio 
#0.860072
CP_Downsampled_Increasing =  data.frame(c(17,4+71), c(1883,7487))
fisher.test(CP_Downsampled_Increasing)
#data:  CP_Downsampled_Increasing
#p-value = 0.7941
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4977058 1.5452708
#sample estimates:
#  odds ratio 
#0.9012597
CP_Downsampled_Decreasing =  data.frame(c(1+27,3+61), c(2180,7193))
fisher.test(CP_Downsampled_Decreasing)
#data:  CP_Downsampled_Decreasing
#p-value = 0.1081
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8888387 2.2904252
#sample estimates:
#  odds ratio 
#1.443487
CR_Increasing =  data.frame(c(20,4+68), c(1722,7651))
fisher.test(CR_Increasing)
#data:  CR_Increasing
#p-value = 0.4169
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7100811 2.0553009
#sample estimates:
#  odds ratio 
#1.234162
CR_Decreasing =  data.frame(c(21,4+67), c(1616,7757))
fisher.test(CR_Decreasing)
#data:  CR_Decreasing
#p-value = 0.1656
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8259503 2.3456162
#sample estimates:
#  odds ratio 
#1.419694
NR_Increasing = data.frame(c(18,4+70), c(1260,8113))
fisher.test(NR_Increasing)
#p-value = 0.09176
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8772005 2.6594069
#sample estimates:
#  odds ratio 
#1.566124
NR_Decreasing = data.frame(c(13,4+75), c(1562,7811))
fisher.test(NR_Decreasing)
#data:  NR_Decreasing
#p-value = 0.5764
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4186881 1.4942415
#sample estimates:
#  odds ratio 
#0.8229066

## Intron retention between retained and exported genes
t.test(sigAgeIRratio[["Cpres.down.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.down.DownPrenatal"]][,"IRratio"], alternative = "greater")
#t = 9.5913, df = 31015, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.004920644         Inf
#sample estimates:
#  mean of x  mean of y 
#0.02753185 0.02159264
t.test(sigAgeIRratio[["Crres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Crres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#data:  sigAgeIRratio[["Crres.UpPrenatal"]][, "IRratio"] and sigAgeIRratio[["Crres.DownPrenatal"]][, "IRratio"]
#t = 3.16, df = 24650, p-value = 0.0007898
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.001087942         Inf
#sample estimates:
#  mean of x  mean of y 
#0.02634000 0.02407089 
t.test(sigAgeIRratio[["Npres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Npres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#t = -5.0201, df = 24383, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.004271234          Inf
#sample estimates:
#  mean of x  mean of y 
#0.02202647 0.02524357
t.test(sigAgeIRratio[["Nrres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Nrres.DownPrenatal"]][,"IRratio"], alternative = "greater")
#t = -7.2555, df = 17789, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.006739706          Inf
#sample estimates:
#  mean of x  mean of y 
#0.02154327 0.02703737
t.test(sigAgeIRratio[["Npres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Npres.DownPrenatal"]][,"IRratio"], alternative = "less")
#t = -5.0201, df = 24383, p-value = 2.6e-07
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.002162971
#sample estimates:
#  mean of x  mean of y 
#0.02202647 0.02524357
t.test(sigAgeIRratio[["Nrres.UpPrenatal"]][,"IRratio"], sigAgeIRratio[["Nrres.DownPrenatal"]][,"IRratio"], alternative = "less")
#t = -7.2555, df = 17789, p-value = 2.084e-13
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.004248494
#sample estimates:
#  mean of x  mean of y 
#0.02154327 0.02703737



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


# Plot LFC by Fraction and IR
FracByAge = list(Ares = data.frame(Ares), Fres.down = data.frame(Fres.down),Ires.down = data.frame(Ires.down))
FracByAge = Map(cbind, FracByAge,lapply(FracByAge, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Adult", "Prenatal", "Interaction"))
r = lapply(irbyGene, function(x) x[which(as.character(x$genes) %in% as.character(FracByAge[["Ares"]][,"ensemblID"])),])
FracByAgeIR = list(Map(cbind, r, lapply(r, function(x) FracByAge[["Ares"]][match(x$genes, FracByAge[["Ares"]][,"ensemblID"]),])),
                   Map(cbind, r, lapply(r, function(x) FracByAge[["Fres.down"]][match(x$genes, FracByAge[["Fres.down"]][,"ensemblID"]),])))
FracByAgeIR = lapply(FracByAgeIR, function(x) do.call(rbind, x))
FracByAgeIR = do.call(rbind, FracByAgeIR)
FracByAgeIR$SampleID = gsub("\\..*","", rownames(FracByAgeIR))
FracByAgeIR$rownum = c(1:nrow(FracByAgeIR))
FracByAgeIR$Fraction = ifelse((FracByAgeIR$rownum %in% grep("C", FracByAgeIR$SampleID)), "Cytosol", "Nucleus")
FracByAgeIR$Age = ifelse(FracByAgeIR$rownum %in% grep("53", FracByAgeIR$SampleID), "Prenatal", "Adult")
FracByAgeIR$Group = paste(FracByAgeIR$Age, FracByAgeIR$Fraction, sep=":")
FracByAgeIR$FDR = ifelse(FracByAgeIR$padj<=0.05, "FDR<0.05", "FDR>0.05")
FracByAgeIR$IR = factor(ifelse(FracByAgeIR$IRratio>=0.5, ">0.5", "<0.5"))
FracByAgeIR = FracByAgeIR[which(FracByAgeIR$padj!="NA"),]
dim(FracByAgeIR[which(FracByAgeIR$genes!=FracByAgeIR$ensemblID),])

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/RNA_localization_byIRratio.pdf",width=7,height=5)
ggplot(FracByAgeIR, aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # saved as RNA_localization_byIRratio.pdf
dev.off()

# Plot LFC by Age and IR
AgebyFrac = list(Cres.down = data.frame(Cres.down), Nres = data.frame(Nres))
AgebyFrac = Map(cbind, AgebyFrac,lapply(AgebyFrac, function(x) geneMap[match(rownames(x),rownames(geneMap)),]),
                Comparison = list("Cytosol", "Nucleus"))
r = lapply(irbyGene, function(x) x[which((as.character(x$genes) %in% as.character(AgebyFrac[["Cres.down"]][,"ensemblID"])) &
                                           (as.character(x$genes) %in% as.character(AgebyFrac[["Nres"]][,"ensemblID"]))),])
AgebyFracIR = list(Map(cbind, r, lapply(r, function(x) AgebyFrac[["Cres.down"]][match(x$genes, AgebyFrac[["Cres.down"]][,"ensemblID"]),])),
                   Map(cbind, r, lapply(r, function(x) AgebyFrac[["Nres"]][match(x$genes, AgebyFrac[["Nres"]][,"ensemblID"]),])))
AgebyFracIR = lapply(AgebyFracIR, function(x) do.call(rbind, x))
AgebyFracIR = do.call(rbind, AgebyFracIR)
AgebyFracIR$SampleID = gsub("\\..*","", rownames(AgebyFracIR))
AgebyFracIR$rownum = c(1:nrow(AgebyFracIR))
AgebyFracIR$Fraction = ifelse((AgebyFracIR$rownum %in% grep("C", AgebyFracIR$SampleID)), "Cytosol", "Nucleus")
AgebyFracIR$Age = ifelse(AgebyFracIR$rownum %in% grep("53", AgebyFracIR$SampleID), "Prenatal", "Adult")
AgebyFracIR$Group = paste(AgebyFracIR$Age, AgebyFracIR$Fraction, sep=":")
AgebyFracIR$FDR = ifelse(AgebyFracIR$padj<=0.05, "FDR<0.05", "FDR>0.05")
AgebyFracIR$IR = factor(ifelse(AgebyFracIR$IRratio>=0.5, ">0.5", "<0.5"))
AgebyFracIR = AgebyFracIR[which(AgebyFracIR$padj!="NA"),]
dim(AgebyFracIR[which(AgebyFracIR$genes!=AgebyFracIR$ensemblID),])

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/Devel_expression_trajectory_byIRratio.pdf",width=7,height=5)
ggplot(AgebyFracIR, aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("Developmental Expression\nTrajectory by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent")) # Save as Devel_expression_trajectory_byIRratio
dev.off()

fracDevel = lapply(sig, function(x) AgebyFracIR[which(AgebyFracIR$genes %in% x$ensID),])
gr = c("Both Retained", "Both Exported", "Retained in Prenatal", "Retained in Adult",
       "Exported in Prenatal", "Exported in Adult", "Retained in Adult/\nExported in Prenatal",
       "Retained in Prenatal/\nExported in Adult", "Interaction")
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/gene_IR_comparisons/RetainedbyAge_LFCxFDRxIR.pdf", width=8, height=6)
for (i in 1:length(fracDevel)){
  g = ggplot(fracDevel[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
    ylab("Log2 Fold Change") + 
    xlab("IR Ratio") +
    ggtitle(paste0("Age Expression Changes by IR Ratio:\n",gr[i])) + 
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
  print(g)
}
dev.off()

### Is there a relationship between direction of expression of significantly DE genes by age and IR being > or < 0.50?
# In cytosol:
c = lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 &
                                                                             x$padj<=0.05 & x$IR=="<0.5"),"genes"])),
                                                       length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 &
                                                                               x$padj<=0.05 & x$IR=="<0.5"),"genes"]))),
                                                     c(length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange>0 &
                                                                               x$padj<=0.05 & x$IR==">0.5"),"genes"])),
                                                       length(unique(x[which(x$Comparison=="Cytosol" & x$log2FoldChange<0 &
                                                                               x$padj<=0.05 & x$IR==">0.5"),"genes"]))))))
c = unlist(lapply(c,function(x) as.character(x)[1]))

# In nucleus
n = lapply(fracDevel, function(x) fisher.test(data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                               x$padj<=0.05 & x$IR=="<0.5"),"genes"])),
                                                       length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                               x$padj<=0.05 & x$IR=="<0.5"),"genes"]))),
                                                     c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                               x$padj<=0.05 & x$IR==">0.5"),"genes"])),
                                                       length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                               x$padj<=0.05 & x$IR==">0.5"),"genes"]))))))
n = unlist(lapply(n,function(x) as.character(x)[1]))
data.frame(Cytosol = c, Nucleus = n, FracGroup = names(c), row.names = NULL)

#Cytosol           Nucleus      FracGroup
#0.801558648580367 0.607847289114744  both_retained
#                1                 1  both_exported
#0.551721724524076 0.360028348688873   Fet_retained
#0.423864414617098 0.133782661254863    Ad_retained
#                1 0.491525423728814   Fet_exported
#                1                 1    Ad_exported
#                1                 1 ret_Ad_exp_Fet
#                1                 1 ret_Fet_exp_Ad
# 0.12361142068079 0.182098693245707    interacting

lapply(fracDevel, function(x) data.frame(c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                               x$padj<=0.05 & x$IR=="<0.5"),"genes"])),
                                                       length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                               x$padj<=0.05 & x$IR=="<0.5"),"genes"]))),
                                                     c(length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange>0 &
                                                                               x$padj<=0.05 & x$IR==">0.5"),"genes"])),
                                                       length(unique(x[which(x$Comparison=="Nucleus" & x$log2FoldChange<0 &
                                                                               x$padj<=0.05 & x$IR==">0.5"),"genes"])))))