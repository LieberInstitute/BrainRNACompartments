library(ggplot2)
library(reshape2)
install.packages("VennDiagram")
library(VennDiagram)
library(data.table)
library("GenomicRanges")

load("./Dropbox/sorted_figures/new/retained.byAge.rda")

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRresP = IRresR = list()
for (i in 1:length(shortenedNames)){
  IRresP[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
  IRresR[[i]] = read.table(paste0(path,"Ribozero/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
}
names(IRresP) = paste0(shortenedNames, "_polyA")
names(IRresR) = paste0(shortenedNames, "_RiboZero")
IRres = c(IRresP, IRresR)

# Filter introns
IRfiltered = lapply(IRresP, function(x) x[which(x$Warnings=="-"),])
string = lapply(IRfiltered, function(x) as.character(x$GeneIntronDetails))
string = lapply(string, function(x) strsplit(x, "/", fixed = TRUE))
x = lapply(string, function(x) unlist(x, recursive = FALSE))
y = lapply(x, function(x) grep("ENSG", x))
c = lapply(x, function(x) seq.int(from = 3, to = length(x), by =3))
comments = genes = list()
for (i in 1:length(string)){
  tmp = x[[i]]
  genes[[i]] = tmp[y[[i]]]
  comments[[i]] = tmp[c[[i]]]
}
names(genes) = names(comments) = names(x)
IRfiltered = Map(cbind, IRfiltered, ensID = genes, comments = comments)
IRfiltered = lapply(IRfiltered, function(x) x[which(x$comments=="clean"),])
ir = lapply(IRfiltered, function(x) data.table(x, key="ensID"))
ir = lapply(ir,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
#code not working
IR = list()
for (i in 1:23){
  tmp = IRfiltered[[i]]
  tmp2 = ir[[i]]
  tmp3 = IR[[i]] 
  for j in (1:nrow(tmp2)){
    tmp3 = data.frame(tmp2[j,], Chr = tmp[which(tmp2$IRratio==tmp$IRratio),colnames(tmp=="Chr")],
                      Start = tmp[which(tmp2$IRratio==tmp$IRratio),colnames(tmp=="Start")],
                      End = tmp[which(tmp2$IRratio==tmp$IRratio),colnames(tmp=="End")]) 
  }}

#Number of introns passing filtering steps (not overlapping anything, enough coverage, etc.)
elementNROWS(ir)
#Number of Genes in sig gene set
elementNROWS(sig)

all = do.call(rbind, ir)
all$SampleID = gsub("\\..*","", rownames(all))
all$Fraction = all$Age = "NA"
all[grep("C", all$SampleID), colnames(all)=="Fraction"] = "Cytosol"
all[grep("N", all$SampleID), colnames(all)=="Fraction"] = "Nucleus"
all[c(grep("1113", all$SampleID),grep("2046", all$SampleID),grep("2074", all$SampleID)),
    colnames(all)=="Age"] = "Adult"
all[c(grep("5339", all$SampleID),grep("5340", all$SampleID),grep("5341", all$SampleID)),
    colnames(all)=="Age"] = "Fetal"
all$Group = paste(all$Age, all$Fraction, sep=":")

sigIR.bySample = list()
for (i in 1:length(sig)){
  tmp = sig[[i]]
  sigIR.bySample[[i]] = lapply(ir, function(x) x[which(x$ensID %in% tmp$EnsID),])}
names(sigIR.bySample) = names(sig)

both_retained = sigIR.bySample[["both_retained"]]
both_exported = sigIR.bySample[["both_exported"]]
Fet_retained = sigIR.bySample[["Fet_retained"]]
Ad_retained = sigIR.bySample[["Ad_retained"]]
ret_Ad_exp_Fet = sigIR.bySample[["ret_Ad_exp_Fet"]]
ret_Fet_exp_Ad = sigIR.bySample[["ret_Fet_exp_Ad"]]
interacting = sigIR.bySample[["interacting"]]
sigIR = list(both_retained = do.call(rbind, both_retained),
             both_exported = do.call(rbind, both_exported),
             Fet_retained = do.call(rbind, Fet_retained),
             Ad_retained = do.call(rbind, Ad_retained),
             ret_Ad_exp_Fet = do.call(rbind, ret_Ad_exp_Fet),
             ret_Fet_exp_Ad = do.call(rbind, ret_Fet_exp_Ad),
             interacting = do.call(rbind, interacting))
both_retained = sigIR[["both_retained"]]
both_exported = sigIR[["both_exported"]]
Fet_retained = sigIR[["Fet_retained"]]
Ad_retained = sigIR[["Ad_retained"]]
ret_Ad_exp_Fet = sigIR[["ret_Ad_exp_Fet"]]
ret_Fet_exp_Ad = sigIR[["ret_Fet_exp_Ad"]]
interacting = sigIR[["interacting"]]

for (i in 1:length(sigIR)){
  tmp = sigIR[[i]]
  tmp$SampleID = gsub("\\..*","", rownames(tmp))
  tmp$Fraction = tmp$Age = "NA"
  tmp[grep("C", tmp$SampleID), colnames(tmp)=="Fraction"] = "Cytosol"
  tmp[grep("N", tmp$SampleID), colnames(tmp)=="Fraction"] = "Nucleus"
  tmp[c(grep("1113", tmp$SampleID),grep("2046", tmp$SampleID),grep("2074", tmp$SampleID)),
      colnames(tmp)=="Age"] = "Adult"
  tmp[c(grep("5339", tmp$SampleID),grep("5340", tmp$SampleID),grep("5341", tmp$SampleID)),
      colnames(tmp)=="Age"] = "Fetal"
  tmp$Group = paste(tmp$Age, tmp$Fraction, sep=":")
  sigIR[[i]] = tmp
}
sigGenes = lapply(sigIR, function(x) data.table(x, key="ensID"))
sigGenes = lapply(sigGenes,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
allGenes = data.table(all, key="ensID")
allGenes = data.frame(allGenes[, list(IRratio=max(IRratio)), by="ensID"])

# genes >50% retained in at least one sample
sigIR.50 = lapply(sigGenes, function(x) x[which(x$IRratio>=.50),colnames(x)=="ensID"])
sigIR.u50 = lapply(sigGenes, function(x) x[which(x$IRratio<.50),colnames(x)=="ensID"])
allgenes.50 = allGenes[which(allGenes$IRratio>=.50),colnames(all)=="ensID"]
allgenes.u50 = allGenes[which(allGenes$IRratio<.50),colnames(all)=="ensID"]

both_retained = list(sig.more = sigIR.50[["both_retained"]], sig.less = sigIR.u50[["both_retained"]],
                     all.more = allgenes.50, all.less = allgenes.u50)
both_exported =  list(sig.more = sigIR.50[["both_exported"]], sig.less = sigIR.u50[["both_exported"]],
                      all.more = allgenes.50, all.less = allgenes.u50)
Fet_retained =  list(sig.more = sigIR.50[["Fet_retained"]], sig.less = sigIR.u50[["Fet_retained"]],
                     all.more = allgenes.50, all.less = allgenes.u50)
Ad_retained =  list(sig.more = sigIR.50[["Ad_retained"]], sig.less = sigIR.u50[["Ad_retained"]],
                    all.more = allgenes.50, all.less = allgenes.u50)
ret_Ad_exp_Fet =  list(sig.more = sigIR.50[["ret_Ad_exp_Fet"]], sig.less = sigIR.u50[["ret_Ad_exp_Fet"]],
                       all.more = allgenes.50, all.less = allgenes.u50)
ret_Fet_exp_Ad =  list(sig.more = sigIR.50[["ret_Fet_exp_Ad"]], sig.less = sigIR.u50[["ret_Fet_exp_Ad"]],
                       all.more = allgenes.50, all.less = allgenes.u50)
interacting =  list(sig.more = sigIR.50[["interacting"]], sig.less = sigIR.u50[["interacting"]],
                    all.more = allgenes.50, all.less = allgenes.u50)

both_retained = lapply(both_retained, function(x) as.character(unique(x)))
both_exported = lapply(both_exported, function(x) as.character(unique(x))) 
Fet_retained =  lapply(Fet_retained, function(x) as.character(unique(x)))
Ad_retained =  lapply(Ad_retained, function(x) as.character(unique(x)))
ret_Ad_exp_Fet = lapply(ret_Ad_exp_Fet, function(x) as.character(unique(x)))
ret_Fet_exp_Ad = lapply(ret_Fet_exp_Ad, function(x) as.character(unique(x)))
interacting = lapply(interacting, function(x) as.character(unique(x)))

venn.diagram(both_retained, "/Users/amandaprice/Dropbox/sorted_figures/new/both_retained.50.jpeg", 
             main="both_retained.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(both_exported, "/Users/amandaprice/Dropbox/sorted_figures/new/both_exported.50.jpeg", 
             main="both_exported.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(Fet_retained, "/Users/amandaprice/Dropbox/sorted_figures/new/Fet_retained.50.jpeg", 
             main="Fet_retained.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(Ad_retained, "/Users/amandaprice/Dropbox/sorted_figures/new/Ad_retained.50.jpeg", 
             main="Ad_retained.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(ret_Ad_exp_Fet, "/Users/amandaprice/Dropbox/sorted_figures/new/ret_Ad_exp_Fet.50.jpeg", 
             main="ret_Ad_exp_Fet.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(ret_Fet_exp_Ad, "/Users/amandaprice/Dropbox/sorted_figures/new/ret_Fet_exp_Ad.50.jpeg", 
             main="ret_Fet_exp_Ad.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)
venn.diagram(interacting, "/Users/amandaprice/Dropbox/sorted_figures/new/interacting.50.jpeg", 
             main="interacting.50", col = "transparent", 
             fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),alpha = 0.50,
             label.col = c("olivedrab4", "white", "darkorchid4", "white", "white", "white", "white", 
                           "white", "palevioletred4", "white", "white", "white", "white", "darkblue", "white"),
             fontfamily = "Arial", fontface = "bold", 
             cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
             cat.fontfamily = "Arial", margin=0.2)

both_retained = data.frame(c(20,469), c(59,9666))
fisher.test(both_retained)
#data:  both_retained
#p-value = 5.193e-10
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  3.947864 11.886203
#sample estimates:
#  odds ratio 
#6.983946
both_exported = data.frame(c(0,906), c(79,9229))
#fisher.test(both_exported)
#data:  both_exported
#p-value = 0.001091
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000000 0.4880419
#sample estimates:
#  odds ratio 
#0
Fet_retained =  data.frame(c(6,139), c(73,9996))
#fisher.test(Fet_retained)
#data:  Fet_retained
#p-value = 0.0008756
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.064977 13.814472
#sample estimates:
#  odds ratio 
#5.908085
Ad_retained =  data.frame(c(28,1641), c(51,8494))
fisher.test(Ad_retained)
#data:  Ad_retained
#p-value = 3.928e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.719687 4.606926
#sample estimates:
#  odds ratio 
#2.841634
ret_Ad_exp_Fet = data.frame(c(0,11), c(79,10124))
#fisher.test(ret_Ad_exp_Fet)
#data:  ret_Ad_exp_Fet
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 52.0122
#sample estimates:
#  odds ratio 
#0
ret_Fet_exp_Ad = data.frame(c(0,4), c(79,10131))
fisher.test(ret_Fet_exp_Ad)
#data:  ret_Fet_exp_Ad
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 196.9234
#sample estimates:
#  odds ratio 
#0
interacting = data.frame(c(13,916), c(66,9219))
fisher.test(interacting)
##data:  interacting
#p-value = 0.02981
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.999046 3.644737
#sample estimates:
#  odds ratio 
#1.982242

# Intron retention between retained and exported genes
ret = sigGenes[["both_retained"]]
exp = sigGenes[["both_exported"]]
t.test(ret$IRratio, exp$IRratio, alternative = "greater")
#data:  ret$IRratio and exp$IRratio
#t = 13.297, df = 521.58, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.07837013        Inf
#sample estimates:
#  mean of x mean of y 
#0.1166033 0.0271476

# What does the distribution of IR Ratios by Sig group look like?
gr = c("Both Retained", "Both Exported", "Retained in Fetal",
       "Retained in Adult", "Retained in Adult/\nExported in Fetal",
       "Retained in Fetal/\nExported in Adult", "Interaction")
pdf("/Users/amandaprice/Dropbox/sorted_figures/new/IR_sig_group_density.pdf")
for (i in 1:7){
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

 



## Differences by intron
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRcomp = list(Adult_PolyA_Zone = read.table(paste0(path, "Adult_PolyA_Zone.tab"), header = TRUE, comment.char="#"),
              Fetal_PolyA_Zone = read.table(paste0(path, "Fetal_PolyA_Zone.tab"), header = TRUE, comment.char="#"),
              Cytosol_PolyA_Age = read.table(paste0(path, "Cytosol_PolyA_Age.tab"), header = TRUE, comment.char="#"),
              Nuclear_PolyA_Age = read.table(paste0(path, "Nuclear_PolyA_Age.tab"), header = TRUE, comment.char="#"),
              PolyA_Zone = read.table(paste0(path, "PolyA_Zone.tab"), header = TRUE, comment.char="#"),
              PolyA_Age = read.table(paste0(path, "PolyA_Age.tab"), header = TRUE, comment.char="#"),
              Adult_Ribozero_Zone = read.table(paste0(path, "Adult_Ribozero_Zone.tab"), header = TRUE, comment.char="#"),
              Fetal_Ribozero_Zone = read.table(paste0(path, "Fetal_Ribozero_Zone.tab"), header = TRUE, comment.char="#"),
              Cytosol_Ribozero_Age = read.table(paste0(path, "Cytosol_Ribozero_Age.tab"), header = TRUE, comment.char="#"),
              Nuclear_Ribozero_Age = read.table(paste0(path, "Nuclear_Ribozero_Age.tab"), header = TRUE, comment.char="#"),
              Ribozero_Zone = read.table(paste0(path, "Ribozero_Zone.tab"), header = TRUE, comment.char="#"),
              Ribozero_Age = read.table(paste0(path, "Ribozero_Age.tab"), header = TRUE, comment.char="#"))
elementNROWS(IRcomp)
string = lapply(IRcomp, function(x) as.character(x$Intron.GeneName.GeneID))
string = lapply(string, function(x) strsplit(x, "/", fixed = TRUE))
x = lapply(string, function(x) unlist(x, recursive = FALSE))
y = lapply(x, function(x) grep("ENSG", x))
c = lapply(x, function(x) seq.int(from = 3, to=length(x), by=3))
comments = genes = list()
for (i in 1:length(string)){
  tmp = x[[i]]
  genes[[i]] = tmp[y[[i]]]
  comments[[i]] = tmp[c[[i]]]
}
names(genes) = names(comments) = names(x)
IR.diff = lapply(IRcomp, function(x) x$A.IRratio - x$B.IRratio)
Sign = lapply(IR.diff, function(x) ifelse(x < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult"))
IRcomp = Map(cbind, IRcomp, ensID = genes, comments = comments, IR.diff = IR.diff, Sign = Sign)

polya = IRcomp[1:6]
IRclean = lapply(polya, function(x) 
  x[which(x$A.warnings=="-" & x$B.warnings=="-" & x$comments=="clean"),])

Are the DEGs DE *because* of the retained introns, which is a developmentally influenced process?
Prove:
#nucleus has more IR than cytosol (despite age): yes
#adult has more IR than fetal (despite fraction): true overall and in nucleus but not cytosol
#that this is true for highly retained introns and not just small pre-mRNA influenced changes: 
  #adult does not have more >50/75% retained introns than fetal overall or in the cytosol but borderline does in nucleus
#that >50% introns are preferentially retained in nucleus: yes? 
#that >50% introns are preferentially present in adult
#that DIR genes are more likely to be DE by fraction
#that genes DE in nucleus are more likely to have >50% IR in at least one intron
#that genes DE by fraction:age are more likely to have >50% IR in at least one intron
Intron comparisons (.tab) between fraction and age are properly enriched in DEGs by the two variables

# Plot LFC by Fraction and IR
Ares$Comparison = "Adult"
Fres$Comparison = "Fetal"

r = lapply(ir, function(x) x[which(x$ensID %in% Ares$EnsID),])
xA = lapply(r, function(x) Ares[which(Ares$EnsID %in% x$ensID),])
xF = lapply(r, function(x) Fres[which(Fres$EnsID %in% x$ensID),])
xA = lapply(xA, function(x) x[order(x$EnsID),])
xF = lapply(xF, function(x) x[order(x$EnsID),])
r = lapply(r, function(x) x[order(x$ensID),])
xA = Map(cbind, xA, r)
xF = Map(cbind, xF, r)
xA = do.call(rbind, xA)
xF = do.call(rbind, xF)

xA$SampleID = gsub("\\..*","", rownames(xA))
xA$Fraction = xA$Age = "NA"
xA[grep("C", xA$SampleID), colnames(xA)=="Fraction"] = "Cytosol"
xA[grep("N", xA$SampleID), colnames(xA)=="Fraction"] = "Nucleus"
xA[c(grep("1113", xA$SampleID),grep("2046", xA$SampleID),grep("2074", xA$SampleID)),
    colnames(xA)=="Age"] = "Adult"
xA[c(grep("5339", xA$SampleID),grep("5340", xA$SampleID),grep("5341", xA$SampleID)),
    colnames(xA)=="Age"] = "Fetal"
xA$Group = paste(xA$Age, xA$Fraction, sep=":")

xF$SampleID = gsub("\\..*","", rownames(xF))
xF$Fraction = xF$Age = "NA"
xF[grep("C", xF$SampleID), colnames(xF)=="Fraction"] = "Cytosol"
xF[grep("N", xF$SampleID), colnames(xF)=="Fraction"] = "Nucleus"
xF[c(grep("1113", xF$SampleID),grep("2046", xF$SampleID),grep("2074", xF$SampleID)),
    colnames(xF)=="Age"] = "Adult"
xF[c(grep("5339", xF$SampleID),grep("5340", xF$SampleID),grep("5341", xF$SampleID)),
    colnames(xF)=="Age"] = "Fetal"
xF$Group = paste(xF$Age, xF$Fraction, sep=":")

total = rbind(xA,xF)
total$FDR = ifelse(total$padj<=0.05, "FDR<0.05", "FDR>0.05")
total$IR = factor(ifelse(total$IRratio>=0.5, ">0.5", "<0.5"))
total = total[which(total$padj!="NA"),]

ggplot(total, aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(total, aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_boxplot() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Make the same plots for Age differences
Cres = read.csv("/Users/amanda/Dropbox/sorted_figures/new/Cres.csv")
Nres = read.csv("/Users/amanda/Dropbox/sorted_figures/new/Nres.csv")
Cres$Comparison = "Cytosol"
Nres$Comparison = "Nucleus"

r = lapply(ir, function(x) x[which(x$ensID %in% Cres$EnsID),])
xC = lapply(r, function(x) Cres[which(Cres$EnsID %in% x$ensID),])
xN = lapply(r, function(x) Nres[which(Nres$EnsID %in% x$ensID),])
xC = lapply(xC, function(x) x[order(x$EnsID),])
xN = lapply(xN, function(x) x[order(x$EnsID),])
r = lapply(r, function(x) x[order(x$ensID),])
xC = Map(cbind, xC, r)
xN = Map(cbind, xN, r)
xC = do.call(rbind, xC)
xN = do.call(rbind, xN)

xC$SampleID = gsub("\\..*","", rownames(xC))
xC$Fraction = xC$Age = "NA"
xC[grep("C", xC$SampleID), colnames(xC)=="Fraction"] = "Cytosol"
xC[grep("N", xC$SampleID), colnames(xC)=="Fraction"] = "Nucleus"
xC[c(grep("1113", xC$SampleID),grep("2046", xC$SampleID),grep("2074", xC$SampleID)),
   colnames(xC)=="Age"] = "Adult"
xC[c(grep("5339", xC$SampleID),grep("5340", xC$SampleID),grep("5341", xC$SampleID)),
   colnames(xC)=="Age"] = "Fetal"
xC$Group = paste(xC$Age, xC$Fraction, sep=":")

xN$SampleID = gsub("\\..*","", rownames(xN))
xN$Fraction = xN$Age = "NA"
xN[grep("C", xN$SampleID), colnames(xN)=="Fraction"] = "Cytosol"
xN[grep("N", xN$SampleID), colnames(xN)=="Fraction"] = "Nucleus"
xN[c(grep("1113", xN$SampleID),grep("2046", xN$SampleID),grep("2074", xN$SampleID)),
   colnames(xN)=="Age"] = "Adult"
xN[c(grep("5339", xN$SampleID),grep("5340", xN$SampleID),grep("5341", xN$SampleID)),
   colnames(xN)=="Age"] = "Fetal"
xN$Group = paste(xN$Age, xN$Fraction, sep=":")

total = rbind(xC,xN)
total$FDR = ifelse(total$padj<=0.05, "FDR<0.05", "FDR>0.05")
total$IR = factor(ifelse(total$IRratio>=0.5, ">0.5", "<0.5"))
total = total[which(total$padj!="NA"),]

ggplot(total, aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
  ylab("Log2 Fold Change") + 
  xlab("IR Ratio") +
  ggtitle("RNA Localization by IR Ratio") + 
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

fracDevel = lapply(sigIR, function(x) total[which(total$ensID %in% x$ensID),])
pdf("/Users/amandaprice/Dropbox/sorted_figures/new/FracDevel_plots.pdf", width=8, height=8)
names = c("Both Nuclear","Both Cytosolic", "Nuclear in Fetal Only",
          "Nuclear in Adult Only", "Nuclear in Adult/Cytosolic in Fetal",
          "Nuclear in Fetal/Cytosolic in Adult", "Interaction Effect")
for (i in 1:length(fracDevel)){
g = ggplot(fracDevel[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
  geom_violin() +
  facet_grid(. ~ Comparison) +
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

# IR regulated genes 
DirList = lapply(IRclean, function(x) split(x, x$Sign))
IRlist = list("Adult:Nucleus"=DirList[["Adult_PolyA_Zone"]][["MoreIRInNuc.Fetal"]],
              "Adult:Cytosol"=DirList[["Adult_PolyA_Zone"]][["MoreIRInCyt.Adult"]],
              "Fetal:Nucleus"=DirList[["Fetal_PolyA_Zone"]][["MoreIRInNuc.Fetal"]],
              "Fetal:Cytosol"=DirList[["Fetal_PolyA_Zone"]][["MoreIRInCyt.Adult"]],
              "Cytosol:Fetal"=DirList[["Cytosol_PolyA_Age"]][["MoreIRInNuc.Fetal"]],
              "Cytosol:Adult"=DirList[["Cytosol_PolyA_Age"]][["MoreIRInCyt.Adult"]],
              "Nucleus:Fetal"=DirList[["Nuclear_PolyA_Age"]][["MoreIRInNuc.Fetal"]],
              "Nucleus:Adult"=DirList[["Nuclear_PolyA_Age"]][["MoreIRInCyt.Adult"]],
              "Nuclear-Enriched"=DirList[["PolyA_Zone"]][["MoreIRInNuc.Fetal"]],
              "Fetal-Enriched"=DirList[["PolyA_Age"]][["MoreIRInNuc.Fetal"]],
              "Adult-Enriched"=DirList[["PolyA_Age"]][["MoreIRInCyt.Adult"]])

totalF = rbind(xA,xF)
totalF$FDR = ifelse(totalF$padj<=0.05, "FDR<0.05", "FDR>0.05")
totalF$IR = factor(ifelse(totalF$IRratio>=0.5, ">0.5", "<0.5"))
totalF = totalF[which(totalF$padj!="NA"),]
combined.IRlist = list(nuc.A.F = rbind(IRlist[["Adult:Nucleus"]], IRlist[["Fetal:Nucleus"]]),
                       cyt.A.F = rbind(IRlist[["Adult:Cytosol"]], IRlist[["Fetal:Cytosol"]]),
                       adult.C.N = rbind(IRlist[["Cytosol:Adult"]], IRlist[["Nucleus:Adult"]]),
                       fetal.C.N = rbind(IRlist[["Cytosol:Fetal"]], IRlist[["Nucleus:Fetal"]]))
IRfrac = lapply(IRlist, function(x) totalF[which(totalF$EnsID %in% x$ensID),])
IRfrac2 = lapply(combined.IRlist, function(x) totalF[which(totalF$EnsID %in% x$ensID),])

pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byFrac_plots.pdf", width=8, height=8)
names = names(IRfrac)
for (i in 1:length(IRfrac)){
  g = ggplot(IRfrac[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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
pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byFrac_combined.pdf", width=8, height=8)
names = names(IRfrac2)
for (i in 1:length(IRfrac2)){
  g = ggplot(IRfrac2[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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


# Make the same plots for Age differences
totalA = rbind(xC,xN)
totalA$FDR = ifelse(totalA$padj<=0.05, "FDR<0.05", "FDR>0.05")
totalA$IR = factor(ifelse(totalA$IRratio>=0.5, ">0.5", "<0.5"))
totalA = totalA[which(totalA$padj!="NA"),]
IRAge = lapply(IRlist, function(x) totalA[which(totalA$ensID %in% x$ensID),])
IRAge2 = lapply(combined.IRlist, function(x) totalA[which(totalA$EnsID %in% x$ensID),])

pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byAge_plots.pdf", width=8, height=8)
names = names(IRlist)
for (i in 1:length(IRAge)){
  g = ggplot(IRAge[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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
pdf("/Users/amanda/Dropbox/sorted_figures/new/IRgenes_byAge_combined.pdf", width=8, height=8)
names = names(IRAge2)
for (i in 1:length(IRAge2)){
  g = ggplot(IRAge2[[i]], aes(x=IR, y=log2FoldChange, fill=FDR), color=FDR) + 
    geom_violin() +
    facet_grid(. ~ Comparison) +
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


# Are the retained introns the last introns?