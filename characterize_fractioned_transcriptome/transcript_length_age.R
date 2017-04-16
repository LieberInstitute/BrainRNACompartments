library(ggplot2)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

# Make list of length of significant genes in a list
AgeList = list(Cpres = data.frame(Cpres), Npres = data.frame(Npres), 
                Crres = data.frame(Crres), Nrres = data.frame(Nrres),
                Cpres.down = data.frame(Cpres.down), Npres.down = data.frame(Npres.down), 
                Crres.down = data.frame(Crres.down), Nrres.down = data.frame(Nrres.down))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigAgeList)
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "UpAdult"))
sigAgeBySign = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(sigAgeBySign, function(x) split(x, x$Sign))
SigList = unlist(DirList, recursive=F)
elementNROWS(SigList)
SigList = Map(cbind, SigList,Group1=list("Cytosol:PolyA:Increasing","Cytosol:PolyA:Decreasing","Nucleus:PolyA:Increasing","Nucleus:PolyA:Decreasing",
                                         "Cytosol:RiboZero:Increasing","Cytosol:RiboZero:Decreasing","Nucleus:RiboZero:Increasing","Nucleus:RiboZero:Decreasing",
                                         "Cytosol:PolyA:Increasing","Cytosol:PolyA:Decreasing","Nucleus:PolyA:Increasing","Nucleus:PolyA:Decreasing",
                                         "Cytosol:RiboZero:Increasing","Cytosol:RiboZero:Decreasing","Nucleus:RiboZero:Increasing","Nucleus:RiboZero:Decreasing"),
              Group=list("Cytosol:Increasing","Cytosol:Decreasing","Nucleus:Increasing","Nucleus:Decreasing","Cytosol:Increasing","Cytosol:Decreasing","Nucleus:Increasing","Nucleus:Decreasing",
                         "Cytosol:Increasing","Cytosol:Decreasing","Nucleus:Increasing","Nucleus:Decreasing","Cytosol:Increasing","Cytosol:Decreasing","Nucleus:Increasing","Nucleus:Decreasing"),
              Development=list("Cytosol","Cytosol","Nucleus","Nucleus","Cytosol","Cytosol","Nucleus","Nucleus","Cytosol","Cytosol","Nucleus","Nucleus","Cytosol","Cytosol","Nucleus","Nucleus"),
              Library=list("PolyA","PolyA","PolyA","PolyA","RiboZero","RiboZero","RiboZero","RiboZero","PolyA","PolyA","PolyA","PolyA","RiboZero","RiboZero","RiboZero","RiboZero"),
              Development=list("Increasing","Decreasing","Increasing","Decreasing","Increasing","Decreasing","Increasing","Decreasing","Increasing","Decreasing","Increasing","Decreasing","Increasing","Decreasing","Increasing","Decreasing"))
lapply(SigList, function(x) head(x))
SigList = lapply(SigList, function(x) data.frame(x, geneMap[match(rownames(x), geneMap$gencodeID),]))
allgenes = data.frame(Sign = NA, Group1 = "All Genes", Group = "All Genes", 
                      Development = NA, Library = "None", Development = NA, geneMap)
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

# By Age where genes are enriched
ggplot(length.down, aes(x=Length/1000)) + geom_density(aes(group=Development, colour=Development)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Developmental Pattern") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# By RNA Development where Age differences were measured
ggplot(length.down, aes(x=Length/1000)) + geom_density(aes(group=Development, colour=Development)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Development") +
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
  data.frame(Group=unique(x$Group),Development=unique(x$Development),
             Development=unique(x$Development),Library=unique(x$Library)))))

dodge <- position_dodge(width=0.9)
limits <- aes(ymax = (Mean/1000 + Std/1000), ymin = (Mean/1000 - Std/1000))
ggplot(Genelength[which(Genelength$Library=="PolyA"),], aes(x=Development, y=Mean/1000, fill=Development), color=Development) + 
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

# t test of length difference between prenatal- and adult-enriched genes

t.test(length.down[which(length.down$Library=="PolyA" & length.down$Development=="Decreasing"),"Length"],
       length.down[which(length.down$Library=="PolyA" & length.down$Development=="Increasing"),"Length"], 
       alternative = "two.sided")
#data:  length.down[which(length.down$Library == "PolyA" & length.down$Development ==  and length.down[which(length.down$Library == "PolyA" & length.down$Development ==     "Decreasing"), "Length"] and     "Increasing"), "Length"]
#t = 20.606, df = 13911, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1155.098 1397.958
#sample estimates:
#  mean of x mean of y 
#5887.831  4611.303 
t.test(length.down[which(length.down$Library=="RiboZero" & length.down$Development=="Decreasing"),"Length"],
       length.down[which(length.down$Library=="RiboZero" & length.down$Development=="Increasing"),"Length"], 
       alternative = "two.sided")
#data:  length.down[which(length.down$Library == "RiboZero" & length.down$Development ==  and length.down[which(length.down$Library == "RiboZero" & length.down$Development ==     "Decreasing"), "Length"] and     "Increasing"), "Length"]
#t = -6.4931, df = 12592, p-value = 8.726e-11
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -619.6619 -332.2850
#sample estimates:
#  mean of x mean of y 
#5064.824  5540.797
t.test(length.down[which(length.down$Group1=="Cytosol:PolyA:Decreasing"),"Length"],
       length.down[which(length.down$Group1=="Cytosol:PolyA:Increasing"),"Length"], 
       alternative = "two.sided")
#data:  length.down[which(length.down$Group1 == "Cytosol:PolyA:Decreasing"),  and length.down[which(length.down$Group1 == "Cytosol:PolyA:Increasing"),     "Length"] and     "Length"]
#t = 19.136, df = 7190, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1463.178 1797.169
#sample estimates:
#  mean of x mean of y 
#6086.762  4456.589
t.test(length.down[which(length.down$Group1=="Cytosol:RiboZero:Decreasing"),"Length"],
       length.down[which(length.down$Group1=="Cytosol:RiboZero:Increasing"),"Length"], 
       alternative = "two.sided")
#data:  length.down[which(length.down$Group1 == "Cytosol:RiboZero:Decreasing"),  and length.down[which(length.down$Group1 == "Cytosol:RiboZero:Increasing"),     "Length"] and     "Length"]
#t = -1.3911, df = 6835.5, p-value = 0.1642
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -321.19637   54.55573
#sample estimates:
#  mean of x mean of y 
#5149.520  5282.841
t.test(length.down[which(length.down$Group1=="Nucleus:PolyA:Decreasing"),"Length"],
       length.down[which(length.down$Group1=="Nucleus:PolyA:Increasing"),"Length"], 
       alternative = "two.sided")
#data:  length.down[which(length.down$Group1 == "Nucleus:PolyA:Decreasing"),  and length.down[which(length.down$Group1 == "Nucleus:PolyA:Increasing"),     "Length"] and     "Length"]
#t = 9.8798, df = 6713.9, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  712.6059 1065.3896
#sample estimates:
#  mean of x mean of y 
#5666.138  4777.140
t.test(length.down[which(length.down$Group1=="Nucleus:RiboZero:Decreasing"),"Length"],
       length.down[which(length.down$Group1=="Nucleus:RiboZero:Increasing"),"Length"], 
       alternative = "two.sided")
#data:  length.down[which(length.down$Group1 == "Nucleus:RiboZero:Decreasing"),  and length.down[which(length.down$Group1 == "Nucleus:RiboZero:Increasing"),     "Length"] and     "Length"]
#t = -8.2479, df = 5194, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1175.452  -723.979
#sample estimates:
#  mean of x mean of y 
#4965.697  5915.412

# Compare length of Development-enriched genes to all genes 

t.test(length.down[length.down$Group1=="Cytosol:PolyA:Decreasing","Length"], allgenes$Length, alternative = "two.sided")
#data:  length.down[length.down$Group1 == "Cytosol:PolyA:Decreasing",  and allgenes$Length    "Length"] and allgenes$Length
#t = 54.6, df = 4092.8, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  3698.658 3974.170
#sample estimates:
#  mean of x mean of y 
#6086.762  2250.348
t.test(length.down[length.down$Group1=="Cytosol:RiboZero:Decreasing","Length"], allgenes$Length, alternative = "two.sided")
#data:  length.down[length.down$Group1 == "Cytosol:RiboZero:Decreasing",  and allgenes$Length    "Length"] and allgenes$Length
#t = 40.723, df = 3264.2, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  2759.584 3038.760
#sample estimates:
#  mean of x mean of y 
#5149.520  2250.348
t.test(length.down[length.down$Group1=="Nucleus:PolyA:Decreasing","Length"], allgenes$Length, alternative = "two.sided")
#data:  length.down[length.down$Group1 == "Nucleus:PolyA:Decreasing",  and allgenes$Length    "Length"] and allgenes$Length
#t = 46.935, df = 3654.2, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  3273.100 3558.479
#sample estimates:
#  mean of x mean of y 
#5666.138  2250.348
t.test(length.down[length.down$Group1=="Nucleus:RiboZero:Decreasing","Length"], allgenes$Length, alternative = "two.sided")
#data:  length.down[length.down$Group1 == "Nucleus:RiboZero:Decreasing",  and allgenes$Length    "Length"] and allgenes$Length
#t = 38.063, df = 2788, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  2575.467 2855.231
#sample estimates:
#  mean of x mean of y 
#4965.697  2250.348
t.test(allgenes$Length, length.down[length.down$Group1=="Cytosol:PolyA:Increasing","Length"], alternative = "two.sided")
#data:  allgenes$Length and length.down[length.down$Group1 == "Cytosol:PolyA:Increasing", allgenes$Length and     "Length"]
#t = -42.685, df = 5207.8, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2307.567 -2104.914
#sample estimates:
#  mean of x mean of y 
#2250.348  4456.589
t.test(allgenes$Length, length.down[length.down$Group1=="Cytosol:RiboZero:Increasing","Length"], alternative = "two.sided")
#data:  allgenes$Length and length.down[length.down$Group1 == "Cytosol:RiboZero:Increasing", allgenes$Length and     "Length"]
#t = -45.366, df = 4462.6, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3163.541 -2901.444
#sample estimates:
#  mean of x mean of y 
#2250.348  5282.841
t.test(allgenes$Length, length.down[length.down$Group1=="Nucleus:PolyA:Increasing","Length"], alternative = "two.sided")
#data:  allgenes$Length and length.down[length.down$Group1 == "Nucleus:PolyA:Increasing", allgenes$Length and     "Length"]
#t = -45.011, df = 4756.5, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2636.847 -2416.736
#sample estimates:
#  mean of x mean of y 
#2250.348  4777.140
t.test(allgenes$Length, length.down[length.down$Group1=="Nucleus:RiboZero:Increasing","Length"], alternative = "two.sided")
#data:  allgenes$Length and length.down[length.down$Group1 == "Nucleus:RiboZero:Increasing", allgenes$Length and     "Length"]
#t = -39.704, df = 2956.6, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3846.061 -3484.067
#sample estimates:
#  mean of x mean of y 
#2250.348  5915.412