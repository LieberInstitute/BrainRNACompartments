library(ggplot2)

# Make list of length of significant genes in a list
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
sigFracBySign = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(sigFracBySign, function(x) split(x, x$Sign))
SigList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]]) 
elementLengths(SigList)
APu = data.frame(SigList[["Apres.Up"]], Group1 = "Adult:PolyA:Nucleus", Group = "Adult:Nucleus",Age = "Adult", Library = "PolyA", Fraction = "Nucleus") 
APd = data.frame(SigList[["Apres.Down"]], Group1 = "Adult:PolyA:Cytosol", Group = "Adult:Cytosol",Age = "Adult", Library = "PolyA", Fraction = "Cytosol") 
FPu = data.frame(SigList[["Fpres.Up"]], Group1 = "Fetal:PolyA:Nucleus", Group = "Fetal:Nucleus",Age = "Fetal", Library = "PolyA", Fraction = "Nucleus") 
FPd = data.frame(SigList[["Fpres.Down"]], Group1 = "Fetal:PolyA:Cytosol", Group = "Fetal:Cytosol",Age = "Fetal", Library = "PolyA", Fraction = "Cytosol") 
ARu = data.frame(SigList[["Arres.Up"]], Group1 = "Adult:Ribozero:Nucleus", Group = "Adult:Nucleus",Age = "Adult", Library = "Ribozero", Fraction = "Nucleus") 
ARd = data.frame(SigList[["Arres.Down"]], Group1 = "Adult:Ribozero:Cytosol", Group = "Adult:Cytosol",Age = "Adult", Library = "Ribozero", Fraction = "Cytosol") 
FRu = data.frame(SigList[["Frres.Up"]], Group1 = "Fetal:Ribozero:Nucleus", Group = "Fetal:Nucleus",Age = "Fetal", Library = "Ribozero", Fraction = "Nucleus")
FRd = data.frame(SigList[["Frres.Down"]], Group1 = "Fetal:Ribozero:Cytosol", Group = "Fetal:Cytosol",Age = "Fetal", Library = "Ribozero", Fraction = "Cytosol")
allgenes = data.frame(FracList[[1]], Sign = NA, Group1 = "All Genes", Group = "All Genes", Age = NA, Library = "None", Fraction = NA)
lengthList = list(APu,APd,FPu,FPd,ARu,ARd,FRu,FRd)
names(lengthList) = c("Adult:PolyA:Nucleus", "Adult:PolyA:Cytosol", "Fetal:PolyA:Nucleus", "Fetal:PolyA:Cytosol",
                      "Adult:Ribozero:Nucleus", "Adult:Ribozero:Cytosol", "Fetal:Ribozero:Nucleus", "Fetal:Ribozero:Cytosol")
length = rbind(APu,FPu,ARu,FRu,APd,FPd,ARd,FRd, allgenes)

# All 8 groups
ggplot(length, aes(x=Length/1000)) + geom_density(aes(group=Group1, colour=Group1)) +
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

ggplot(length[which(length$Library=="PolyA" | length$Library=="None"),], aes(x=Length/1000)) + geom_density(aes(group=Group, colour=Group)) +
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

ggplot(length[which(length$Library=="Ribozero"| length$Library=="None"),], aes(x=Length/1000)) + geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("Gene Length (Kb)") +
  ggtitle("Gene Length By Group (Ribozero)") +
  xlim(0,20) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# By Fraction
length = rbind(APu,FPu,ARu,FRu,APd,FPd,ARd,FRd)
ggplot(length, aes(x=Length/1000)) + geom_density(aes(group=Fraction, colour=Fraction)) +
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
ggplot(length, aes(x=Length/1000)) + geom_density(aes(group=Age, colour=Age)) +
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
ggplot(length, aes(x=Length/1000)) + geom_density(aes(group=Library, colour=Library)) +
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


stats = lapply(lengthList, function(x) {c(Number = nrow(x), Min = min(x$Length), Max = max(x$Length), 
                                           Mean = mean(x$Length), Median = median(x$Length), Std = sd(x$Length))})
Genelength = do.call(rbind, stats)
Genelength = data.frame(Genelength, Group = c("Adult:PolyA","Adult:PolyA","Fetal:PolyA","Fetal:PolyA","Adult:Ribozero",
                                              "Adult:Ribozero","Fetal:Ribozero","Fetal:Ribozero"), 
                        Fraction = c("Nucleus", "Cytosol", "Nucleus", "Cytosol", "Nucleus", "Cytosol", "Nucleus", "Cytosol"))

dodge <- position_dodge(width=0.9)
limits <- aes(ymax = (Mean + Std), ymin = (Mean - Std))
ggplot(Genelength, aes(x=Group, y=Mean, fill=Fraction), color=Fraction) + 
  stat_summary(position=position_dodge(),geom="bar") +
  geom_errorbar(mapping = limits, position = dodge, width=0.25) +
  ylab("Gene Length (BP)") +
  xlab("") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.76, 0.86)) + 
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# t test of Up vs Down length
Up.polya = rbind(lengthList[["Adult:PolyA:Nucleus"]], lengthList[["Fetal:PolyA:Nucleus"]])
Up.ribo = rbind(lengthList[["Adult:Ribozero:Nucleus"]], lengthList[["Fetal:Ribozero:Nucleus"]])
D.polya = rbind(lengthList[["Adult:PolyA:Cytosol"]], lengthList[["Fetal:PolyA:Cytosol"]])
D.ribo = rbind(lengthList[["Adult:Ribozero:Cytosol"]], lengthList[["Fetal:Ribozero:Cytosol"]])

t.test(Up.polya$Length, D.polya$Length, alternative = "greater")
#Welch Two Sample t-test

#data:  Up.polya and D.polya
#t = 15.211, df = 1394.6, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2357.217      Inf
#sample estimates:
#  mean of x mean of y 
#5553.771  2910.527 

t.test(Up.ribo, D.ribo, alternative = "greater")
#data:  Up.ribo and D.ribo
#t = 1.3322, df = 2348.9, p-value = 0.09147
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -45.76775       Inf
#sample estimates:
#  mean of x mean of y 
#3390.511  3195.936

length = lapply(lengthList, function(x) x$Length)

t.test(length[["Adult:PolyA:Nucleus"]], length[["Adult:PolyA:Cytosol"]], alternative = "greater")
#data:  length[["Adult:PolyA:Nucleus"]] and length[["Adult:PolyA:Cytosol"]]
#t = 14.744, df = 1262, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2404.53     Inf
#sample estimates:
#  mean of x mean of y 
#5620.999  2914.288

t.test(length[["Adult:Ribozero:Nucleus"]], length[["Adult:Ribozero:Cytosol"]], alternative = "greater")
#data:  length[["Adult:Ribozero:Nucleus"]] and length[["Adult:Ribozero:Cytosol"]]
#t = 0.95917, df = 1958.9, p-value = 0.1688
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -96.81393       Inf
#sample estimates:
#  mean of x mean of y 
#3351.997  3216.721

t.test(length[["Fetal:PolyA:Nucleus"]], length[["Fetal:PolyA:Cytosol"]], alternative = "greater")
#Not enough samples

t.test(length[["Fetal:Ribozero:Nucleus"]], length[["Fetal:Ribozero:Cytosol"]], alternative = "greater")
#data:  length[["Fetal:Ribozero:Nucleus"]] and length[["Fetal:Ribozero:Cytosol"]]
#t = 6.4975, df = 258.76, p-value = 2.084e-10
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2384.553      Inf
#sample estimates:
#  mean of x mean of y 
#3561.72    365.00 

t.test(length[["Adult:PolyA:Nucleus"]], allgenes$Length, alternative = "greater")
t.test(length[["Adult:Ribozero:Nucleus"]], allgenes$Length, alternative = "greater")
t.test(length[["Fetal:PolyA:Nucleus"]], allgenes$Length, alternative = "greater")
t.test(length[["Fetal:Ribozero:Nucleus"]], allgenes$Length, alternative = "greater")
t.test(allgenes$Length, length[["Adult:PolyA:Cytosol"]], alternative = "greater")
t.test(allgenes$Length, length[["Adult:Ribozero:Cytosol"]], alternative = "greater")
t.test(allgenes$Length, length[["Fetal:PolyA:Cytosol"]], alternative = "greater")
t.test(allgenes$Length, length[["Fetal:Ribozero:Cytosol"]], alternative = "greater")

t.test(length[["Adult:PolyA:Nucleus"]], allgenes$Length, alternative = "greater")
#data:  length[["Adult:PolyA:Nucleus"]] and allgenes$Length
#t = 17.38, df = 968.3, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  2675.554      Inf
#sample estimates:
#  mean of x mean of y 
#5620.999  2665.459 

t.test(length[["Adult:Ribozero:Nucleus"]], allgenes$Length, alternative = "greater")
#data:  length[["Adult:Ribozero:Nucleus"]] and allgenes$Length
#t = 5.1986, df = 1546.4, p-value = 1.138e-07
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  469.1836      Inf
#sample estimates:
#  mean of x mean of y 
#3351.997  2665.459 

t.test(length[["Fetal:PolyA:Nucleus"]], allgenes$Length, alternative = "greater")
#data:  length[["Fetal:PolyA:Nucleus"]] and allgenes$Length
#t = 7.0056, df = 73.404, p-value = 4.968e-10
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  1540.894      Inf
#sample estimates:
#  mean of x mean of y 
#4687.081  2665.459 

t.test(length[["Fetal:Ribozero:Nucleus"]], allgenes$Length, alternative = "greater")
#data:  length[["Fetal:Ribozero:Nucleus"]] and allgenes$Length
#t = 1.9395, df = 338.73, p-value = 0.02663
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  134.0781      Inf
#sample estimates:
#  mean of x mean of y 
#3561.720  2665.459 

t.test(allgenes$Length, length[["Adult:PolyA:Cytosol"]], alternative = "greater")
#data:  allgenes$Length and length[["Adult:PolyA:Cytosol"]]
#t = -3.4373, df = 757.88, p-value = 0.9997
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -368.0469       Inf
#sample estimates:
#  mean of x mean of y 
#2665.459  2914.288 

t.test(allgenes$Length, length[["Adult:Ribozero:Cytosol"]], alternative = "greater")
#data:  allgenes$Length and length[["Adult:Ribozero:Cytosol"]]
#t = -10.221, df = 1603.7, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -640.0254       Inf
#sample estimates:
#  mean of x mean of y 
#2665.459  3216.721 

t.test(allgenes$Length, length[["Fetal:PolyA:Cytosol"]], alternative = "greater")
#not enough 'y' observations
t.test(allgenes$Length, length[["Fetal:Ribozero:Cytosol"]], alternative = "greater")
#data:  allgenes$Length and length[["Fetal:Ribozero:Cytosol"]]
#t = 13.515, df = 9.1442, p-value = 1.187e-07
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#1988.986      Inf
#sample estimates:
#mean of x mean of y 
#2665.459   365.000

t.test(allgenes$Length, length[["Adult:PolyA:Cytosol"]])
t.test(allgenes$Length, length[["Adult:Ribozero:Cytosol"]])

t.test(allgenes$Length, length[["Adult:PolyA:Cytosol"]])
#data:  allgenes$Length and length[["Adult:PolyA:Cytosol"]]
#t = -3.4373, df = 757.88, p-value = 0.0006196
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -390.9392 -106.7190
#sample estimates:
#  mean of x mean of y 
#2665.459  2914.288 

t.test(allgenes$Length, length[["Adult:Ribozero:Cytosol"]])
#data:  allgenes$Length and length[["Adult:Ribozero:Cytosol"]]
#t = -10.221, df = 1603.7, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -657.0489 -445.4749
#sample estimates:
#  mean of x mean of y 
#2665.459  3216.721