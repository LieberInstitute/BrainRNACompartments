library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Make list of length of significant genes in a list
FracList = list(Apres = data.frame(Apres), Fpres = data.frame(Fpres), 
                Arres = data.frame(Arres), Frres = data.frame(Frres),
                Apres.down = data.frame(Apres.down), Fpres.down = data.frame(Fpres.down), 
                Arres.down = data.frame(Arres.down), Frres.down = data.frame(Frres.down))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracBySign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(sigFracBySign, function(x) split(x, x$Sign))
SigList = unlist(DirList, recursive=F)
elementNROWS(SigList)
SigList = Map(cbind, SigList,Group1=list("Adult:PolyA:Cytosolic","Adult:PolyA:Nuclear","Prenatal:PolyA:Cytosolic","Prenatal:PolyA:Nuclear",
                                         "Adult:RiboZero:Cytosolic","Adult:RiboZero:Nuclear","Prenatal:RiboZero:Cytosolic","Prenatal:RiboZero:Nuclear",
                                         "Adult:PolyA:Cytosolic","Adult:PolyA:Nuclear","Prenatal:PolyA:Cytosolic","Prenatal:PolyA:Nuclear",
                                         "Adult:RiboZero:Cytosolic","Adult:RiboZero:Nuclear","Prenatal:RiboZero:Cytosolic","Prenatal:RiboZero:Nuclear"),
              Group=list("Adult:Cytosolic","Adult:Nuclear","Prenatal:Cytosolic","Prenatal:Nuclear","Adult:Cytosolic","Adult:Nuclear","Prenatal:Cytosolic","Prenatal:Nuclear",
                         "Adult:Cytosolic","Adult:Nuclear","Prenatal:Cytosolic","Prenatal:Nuclear","Adult:Cytosolic","Adult:Nuclear","Prenatal:Cytosolic","Prenatal:Nuclear"),
              Age=list("Adult","Adult","Prenatal","Prenatal","Adult","Adult","Prenatal","Prenatal","Adult","Adult","Prenatal","Prenatal","Adult","Adult","Prenatal","Prenatal"),
              Library=list("PolyA","PolyA","PolyA","PolyA","RiboZero","RiboZero","RiboZero","RiboZero","PolyA","PolyA","PolyA","PolyA","RiboZero","RiboZero","RiboZero","RiboZero"),
              Fraction=list("Cytosolic","Nuclear","Cytosolic","Nuclear","Cytosolic","Nuclear","Cytosolic","Nuclear","Cytosolic","Nuclear","Cytosolic","Nuclear","Cytosolic","Nuclear","Cytosolic","Nuclear"))
lapply(SigList, function(x) head(x))
SigList = lapply(SigList, function(x) data.frame(x, geneMap[match(rownames(x), geneMap$gencodeID),]))
allgenes = data.frame(Sign = NA, Group1 = "All Genes", Group = "All Genes", 
                      Age = NA, Library = "None", Fraction = NA, geneMap)
length = rbind(do.call(rbind, lapply(SigList[1:8], function(x) x[,c(8:12,17:18)])),allgenes[,c(2:6,11:12)])
length.down = rbind(do.call(rbind, lapply(SigList[9:16], function(x) x[,c(8:12,17:18)])),allgenes[,c(2:6,11:12)])

# All 8 groups
ggplot(length.down, aes(x=Length/1000)) + geom_density(aes(group=Group1, colour=Group1)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Group") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(length.down[which(length.down$Library=="PolyA" | length.down$Library=="None"),], aes(x=Length/1000)) + 
  geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Group (PolyA)") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(length.down[which(length.down$Library=="RiboZero" | length.down$Library=="None"),], aes(x=Length/1000)) + 
  geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Group (RiboZero)") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# By Fraction where genes are enriched
ggplot(length.down, aes(x=Length/1000)) + geom_density(aes(group=Fraction, colour=Fraction)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Localization") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# By Age
ggplot(length.down, aes(x=Length/1000)) + geom_density(aes(group=Age, colour=Age)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Age") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# By Library
ggplot(length.down, aes(x=Length/1000)) + geom_density(aes(group=Library, colour=Library)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Library") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))


stats = lapply(SigList[9:16], function(x) 
  c(Number = nrow(x), Min = min(x$Length), Max = max(x$Length), 
     Mean = mean(x$Length), Median = median(x$Length), Std = sd(x$Length)))
Genelength = do.call(rbind, stats)
Genelength = data.frame(Genelength, do.call(rbind, lapply(SigList[9:16], function(x) 
  data.frame(Group=unique(x$Group),Fraction=unique(x$Fraction),
             Age=unique(x$Age),Library=unique(x$Library)))))

dodge <- position_dodge(width=0.9)
limits <- aes(ymax = (Mean/1000 + Std/1000), ymin = (Mean/1000 - Std/1000))
ggplot(Genelength[which(Genelength$Library=="PolyA"),], aes(x=Age, y=Mean/1000, fill=Fraction), color=Fraction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  geom_errorbar(mapping = limits, position = dodge, width=0.25) +
  ylab("Gene Length (Kb)") +
  xlab("") +
  ggtitle("Mean Gene Length (PolyA)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.76, 0.86)) + 
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# t test of length difference between nuclear- and cytosolic-enriched genes

t.test(length.down[which(length.down$Library=="PolyA" & length.down$Fraction=="Nuclear"),"Length"],
       length.down[which(length.down$Library=="PolyA" & length.down$Fraction=="Cytosolic"),"Length"], 
       alternative = "greater")
#data:  length.down[which(length.down$Library == "PolyA" & length.down$Fraction ==  and length.down[which(length.down$Library == "PolyA" & length.down$Fraction ==     "Nuclear"), "Length"] and     "Cytosolic"), "Length"]
#t = 15.061, df = 1204.5, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2748.695      Inf
#sample estimates:
#  mean of x mean of y 
#5869.323  2783.350 
t.test(length.down[which(length.down$Library=="RiboZero" & length.down$Fraction=="Nuclear"),"Length"],
       length.down[which(length.down$Library=="RiboZero" & length.down$Fraction=="Cytosolic"),"Length"], 
       alternative = "greater")
#data:  length.down[which(length.down$Library == "RiboZero" & length.down$Fraction ==  and length.down[which(length.down$Library == "RiboZero" & length.down$Fraction ==     "Nuclear"), "Length"] and     "Cytosolic"), "Length"]
#t = 8.7688, df = 990.26, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  1846.253      Inf
#sample estimates:
#  mean of x mean of y 
#5318.323  3045.298
t.test(length.down[which(length.down$Group1=="Adult:PolyA:Nuclear"),"Length"],
       length.down[which(length.down$Group1=="Adult:PolyA:Cytosolic"),"Length"], 
       alternative = "greater")
#data:  length.down[which(length.down$Group1 == "Adult:PolyA:Nuclear"),  and length.down[which(length.down$Group1 == "Adult:PolyA:Cytosolic"),     "Length"] and     "Length"]
#t = 14.701, df = 1144.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2763.33     Inf
#sample estimates:
#  mean of x mean of y 
#5897.791  2785.998
t.test(length.down[which(length.down$Group1=="Adult:RiboZero:Nuclear"),"Length"],
       length.down[which(length.down$Group1=="Adult:RiboZero:Cytosolic"),"Length"], 
       alternative = "greater")
#data:  length.down[which(length.down$Group1 == "Adult:RiboZero:Nuclear"),  and length.down[which(length.down$Group1 == "Adult:RiboZero:Cytosolic"),     "Length"] and     "Length"]
#t = 8.6815, df = 977.19, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  1722.442      Inf
#sample estimates:
#  mean of x mean of y 
#5189.689  3064.147
t.test(length.down[which(length.down$Group1=="Prenatal:PolyA:Nuclear"),"Length"],
       length.down[which(length.down$Group1=="Prenatal:PolyA:Cytosolic"),"Length"], 
       alternative = "greater")
#Not enough samples
t.test(length.down[which(length.down$Group1=="Prenatal:RiboZero:Nuclear"),"Length"],
       length.down[which(length.down$Group1=="Prenatal:RiboZero:Cytosolic"),"Length"], 
       alternative = "greater")
#data:  length.down[which(length.down$Group1 == "Prenatal:RiboZero:Nuclear"),  and length.down[which(length.down$Group1 == "Prenatal:RiboZero:Cytosolic"),     "Length"] and     "Length"]
#t = 2.5724, df = 22.013, p-value = 0.008685
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  3286.644      Inf
#sample estimates:
#  mean of x  mean of y 
#10172.8696   287.8571

# Compare length of fraction-enriched genes to all genes 

t.test(length.down[length.down$Group1=="Adult:PolyA:Nuclear","Length"], allgenes$Length, alternative = "greater")
#data:  length.down[length.down$Group1 == "Adult:PolyA:Nuclear", "Length"] and allgenes$Length
#t = 18.034, df = 963.26, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  3314.445      Inf
#sample estimates:
#  mean of x mean of y 
#5897.791  2250.348
t.test(length.down[length.down$Group1=="Adult:RiboZero:Nuclear","Length"], allgenes$Length, alternative = "greater")
#data:  length.down[length.down$Group1 == "Adult:RiboZero:Nuclear", "Length"] and allgenes$Length
#t = 12.361, df = 872.41, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2547.8    Inf
#sample estimates:
#  mean of x mean of y 
#5189.689  2250.348
t.test(length.down[length.down$Group1=="Prenatal:PolyA:Nuclear","Length"], allgenes$Length, alternative = "greater")
#data:  length.down[length.down$Group1 == "Prenatal:PolyA:Nuclear", "Length"] and allgenes$Length
#t = 6.6766, df = 38.07, p-value = 3.35e-08
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2183.534      Inf
#sample estimates:
#  mean of x mean of y 
#5171.487  2250.348
t.test(length.down[length.down$Group1=="Prenatal:RiboZero:Nuclear","Length"], allgenes$Length, alternative = "greater")
#data:  length.down[length.down$Group1 == "Prenatal:RiboZero:Nuclear",  and allgenes$Length    "Length"] and allgenes$Length
#t = 2.062, df = 22.001, p-value = 0.02561
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  1324.948      Inf
#sample estimates:
#  mean of x mean of y 
#10172.870  2250.348
t.test(allgenes$Length, length.down[length.down$Group1=="Adult:PolyA:Cytosolic","Length"], alternative = "greater")
#data:  allgenes$Length and length.down[length.down$Group1 == "Adult:PolyA:Cytosolic", "Length"]
#t = -8.2139, df = 1019.5, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -643.0118       Inf
#sample estimates:
#  mean of x mean of y 
#2250.348  2785.998
t.test(allgenes$Length, length.down[length.down$Group1=="Adult:RiboZero:Cytosolic","Length"], alternative = "greater")
#data:  allgenes$Length and length.down[length.down$Group1 == "Adult:RiboZero:Cytosolic", allgenes$Length and     "Length"]
#t = -13.285, df = 1125.9, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -914.6443       Inf
#sample estimates:
#  mean of x mean of y 
#2250.348  3064.147
t.test(allgenes$Length, length.down[length.down$Group1=="Prenatal:PolyA:Cytosolic","Length"], alternative = "greater")
# not enough observations
t.test(allgenes$Length, length.down[length.down$Group1=="Prenatal:RiboZero:Cytosolic","Length"], alternative = "greater")
#data:  allgenes$Length and length.down[length.down$Group1 == "Prenatal:RiboZero:Cytosolic", allgenes$Length and     "Length"]
#t = 28.713, df = 6.4781, p-value = 2.247e-08
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  1831.412      Inf
#sample estimates:
#  mean of x mean of y 
#2250.3481  287.8571