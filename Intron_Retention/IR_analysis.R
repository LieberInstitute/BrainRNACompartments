library(ggplot2)
library(reshape2)
library(VennDiagram)
library(data.table)

load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
}
names(IRres) = shortenedNames

# Filter introns
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings=="-"),])       
clean = lapply(IRfiltered, function(x) grep("clean", x$GeneIntronDetails, fixed=T))
IRfiltered2 = list()
for (i in 1:length(IRfiltered)){
  cl = clean[[i]]
  irf = IRfiltered[[i]]
  IRfiltered2[[i]] = irf[cl,]
}
names(IRfiltered2) = names(IRfiltered)
string = lapply(IRfiltered2, function(x) strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE))
string = lapply(string, function(x) unlist(x, recursive = FALSE))
y = lapply(string, function(x) grep("ENSG", x))
genes = list()
for (i in 1:length(string)){
  tmp = string[[i]]
  genes[[i]] = tmp[y[[i]]]
}
names(genes) = names(string)
intronID = lapply(IRfiltered2, function(x) paste0(x$Chr,":",x$Start,"-",x$End))
IRfiltered2 = Map(cbind, IRfiltered2, genes = genes, intronID = intronID)
head(IRfiltered2[[1]])

## Look at distribution of filtered introns by group
introns = data.frame(num.introns = elementNROWS(IRfiltered2))
introns$SampleID = as.factor(rownames(introns))
introns$Group = factor(x=c("Adult:Cytosol","Adult:Nucleus","Adult:Cytosol","Adult:Nucleus",
                           "Adult:Cytosol","Adult:Nucleus","Prenatal:Cytosol","Prenatal:Nucleus",
                           "Prenatal:Cytosol","Prenatal:Nucleus","Prenatal:Cytosol","Prenatal:Nucleus"),
                       levels = c("Adult:Cytosol","Prenatal:Cytosol",
                                  "Adult:Nucleus","Prenatal:Nucleus"))
introns$num.genes = as.numeric(lapply(genes, function(x) length(unique(x))))
introns$MeanIntronDepth = as.numeric(lapply(IRfiltered2, function(x) mean(x$IntronDepth)))
write.csv(introns, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/global_IR_comparisons/filtered.intron.stats.csv")

ggplot(introns, aes(x=Group, y=num.introns)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + 
  xlab("") +
  ggtitle("Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,60000) +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(introns, aes(x=Group, y=num.genes)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + 
  xlab("") +
  ggtitle("Unique Genes Containing Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,8000)

ggplot(introns, aes(x=Group, y=MeanIntronDepth)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + 
  xlab("") +
  ggtitle("Mean Intron Depth (Passing QC)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,.35)

allintrons = Reduce(intersect, intronID)
length(allintrons) # 3621 of the introns pass QC in all four groups


### Global IR Comparisons:
# What does the distribution of IR Ratios look like?

IRratios = lapply(IRfiltered2, function(x) data.frame(IRratio=x$IRratio,intronID=x$intronID))
IRratios = Map(cbind, IRratios, Age=NA, Fraction=NA)
for (i in 1:6){IRratios[[i]][,"Age"] = "Adult"}
for (i in 7:12){IRratios[[i]][,"Age"] = "Prenatal"}
for (i in c(2,4,6,8,10,12)){IRratios[[i]][,"Fraction"] = "Nucleus"}
for (i in c(1,3,5,7,9,11)){IRratios[[i]][,"Fraction"] = "Cytosol"}
IRratios = do.call(rbind, IRratios)
IRratios$Group = paste(IRratios$Age, IRratios$Fraction, sep=":")
IRratios$Group = factor(IRratios$Group, levels = c("Adult:Cytosol", "Prenatal:Cytosol",
                                                   "Adult:Nucleus", "Prenatal:Nucleus"))

ggplot(IRratios, aes(x=IRratio)) + geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("IR Ratio") +
  ggtitle("IR Ratios By Group") +
  xlim(0,.06) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Just in shared introns
ratio.overlaps = IRratios[which(IRratios$intronID %in% allintrons),]
ggplot(ratio.overlaps, aes(x=IRratio)) + geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("IR Ratio") +
  ggtitle("IR Ratios By Group\n(Introns Passing QC In All Groups)") +
  xlim(0,.06) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))


# How many introns are retained at different thresholds?
PercentIRs = data.frame(Constitutively.Spliced = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio==0),]))),
                 ">0%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>0),]))),
                 "less.five" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.05),]))),
                 ">10%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.10),]))),
                 ">20%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.20),]))),
                 ">30%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.30),]))),
                 ">40%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.40),]))),
                 "fifty" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.50),]))),
                 ">75%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.75),]))),
                 ">90%" = unlist(lapply(IRfiltered2, function(x) nrow(x[which(x$IRratio>=.9),]))))
Total = unlist(lapply(IRfiltered2, function(x) nrow(x)))
PercentIRs = round(PercentIRs / Total *100, digits = 5)
PercentIRs = data.frame(melt(PercentIRs), SampleID = rep.int(rownames(PercentIRs), 10))
PercentIRs$Age = c(rep.int("Adult", 6), rep.int("Prenatal", 6))
PercentIRs$Fraction = rep.int(c("Cytosol", "Nucleus"), 6)
PercentIRs$Group = factor(paste(PercentIRs$Age, PercentIRs$Fraction, sep=":"), 
                          levels = c("Adult:Cytosol", "Prenatal:Cytosol", "Adult:Nucleus", "Prenatal:Nucleus"))
PercentIRs$variable = gsub("Constitutively.Spliced", "Constitutively\nSpliced", PercentIRs$variable)
PercentIRs$variable = gsub("X.0.", ">0%", PercentIRs$variable)
PercentIRs$variable = gsub("less.five",">5%", PercentIRs$variable)
PercentIRs$variable = gsub("X.10.",">10%", PercentIRs$variable)
PercentIRs$variable = gsub("X.20.",">20%", PercentIRs$variable)
PercentIRs$variable = gsub("X.30.",">30%", PercentIRs$variable)
PercentIRs$variable = gsub("X.40.",">40%", PercentIRs$variable)
PercentIRs$variable = gsub("fifty",">50%", PercentIRs$variable)
PercentIRs$variable = gsub("X.75.",">75%", PercentIRs$variable)
PercentIRs$variable = gsub("X.90.",">90%", PercentIRs$variable)
PercentIRs$variable = factor(PercentIRs$variable, levels=c("Constitutively\nSpliced", ">0%", ">5%", 
                                                           ">10%", ">20%", ">30%", ">40%", ">50%", ">75%", ">90%"))
ggplot(PercentIRs, aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(PercentIRs[which(PercentIRs$variable != "Constitutively\nSpliced" & PercentIRs$variable != ">0%"),],
       aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# In only the introns in all four groups
overlaps = lapply(IRfiltered2, function(x) x[which(x$intronID %in% allintrons),])
PercentIRsOverlap = data.frame(Constitutively.Spliced = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio==0),]))),
                        ">0%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>0),]))),
                        "less.five" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.05),]))),
                        ">10%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.10),]))),
                        ">20%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.20),]))),
                        ">30%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.30),]))),
                        ">40%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.40),]))),
                        "fifty" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.50),]))),
                        ">75%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.75),]))),
                        ">90%" = unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.9),]))))
PercentIRsOverlap = round(PercentIRsOverlap / length(allintrons) *100, digits = 5)
PercentIRsOverlap = data.frame(melt(PercentIRsOverlap), SampleID = rep.int(rownames(PercentIRsOverlap), 10))
PercentIRsOverlap$Age = c(rep.int("Adult", 6), rep.int("Prenatal", 6))
PercentIRsOverlap$Fraction = rep.int(c("Cytosol", "Nucleus"), 6)
PercentIRsOverlap$Group = factor(paste(PercentIRsOverlap$Age, PercentIRsOverlap$Fraction, sep=":"), 
                          levels = c("Adult:Cytosol", "Prenatal:Cytosol", "Adult:Nucleus", "Prenatal:Nucleus"))
PercentIRsOverlap$variable = gsub("Constitutively.Spliced", "Constitutively\nSpliced", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.0.", ">0%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("less.five",">5%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.10.",">10%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.20.",">20%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.30.",">30%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.40.",">40%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("fifty",">50%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.75.",">75%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = gsub("X.90.",">90%", PercentIRsOverlap$variable)
PercentIRsOverlap$variable = factor(PercentIRsOverlap$variable, levels=c("Constitutively\nSpliced", ">0%", ">5%", 
                                                           ">10%", ">20%", ">30%", ">40%", ">50%", ">75%", ">90%"))
ggplot(PercentIRsOverlap, aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio\n(Introns Passing QC In All Groups)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(PercentIRsOverlap[which(PercentIRsOverlap$variable != "Constitutively\nSpliced" & 
                                 PercentIRsOverlap$variable != ">0%"),],
       aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio\n(Introns Passing QC In All Groups)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Is IR increased in the nucleus?
# Main effect: All filtered introns
t.test(x = IRratios[which(IRratios$Fraction=="Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Fraction=="Cytosol"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Fraction == "Nucleus"), "IRratio"] and IRratios[which(IRratios$Fraction == "Cytosol"), "IRratio"]
#t = 57.47, df = 373260, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.003913777         Inf
#sample estimates:
#  mean of x   mean of y 
#0.008557266 0.004528172

# in adult polyA
t.test(x = IRratios[which(IRratios$Group=="Adult:Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Adult:Cytosol"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Group == "Adult:Nucleus"), "IRratio"] and IRratios[which(IRratios$Group == "Adult:Cytosol"), "IRratio"]
#t = 49.278, df = 130080, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.005493063         Inf
#sample estimates:
#  mean of x   mean of y 
#0.009611667 0.003928916

# in Prenatal polyA
t.test(x = IRratios[which(IRratios$Group=="Prenatal:Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Cytosol"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Group == "Prenatal:Nucleus"), "IRratio"] and IRratios[which(IRratios$Group == "Prenatal:Cytosol"), "IRratio"]
#t = 34.448, df = 242280, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.002878125         Inf
#sample estimates:
#  mean of x   mean of y 
#0.007889297 0.004866854

# Main effect: introns >=75% retained
head(PercentIRs)
# in all polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$variable == ">75%"),"value"],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$variable == ">75%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$variable ==     ">75%"), "value"] and     ">75%"), "value"]
#t = 0.4709, df = 8.9399, p-value = 0.3245
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.0026255        Inf
#sample estimates:
#  mean of x   mean of y 
#0.003896667 0.002990000
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$variable == ">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$variable == ">50%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = 3.0601, df = 5.5816, p-value = 0.01219
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.01136551        Inf
#sample estimates:
#  mean of x  mean of y 
#0.04321833 0.01132500

# in adult polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">75%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">75%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Adult:Cytosol" & PercentIRs$variable ==     ">75%"), "value"] and     ">75%"), "value"]
#t = 0.93677, df = 3.8054, p-value = 0.2022
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.004404569          Inf
#sample estimates:
#  mean of x   mean of y 
#0.005596667 0.002233333
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Adult:Cytosol" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = 2.8918, df = 2.1153, p-value = 0.04764
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.0012081       Inf
#sample estimates:
#  mean of x  mean of y 
#0.05553333 0.01050333

# in Prenatal polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">75%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">75%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Prenatal:Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Cytosol" & PercentIRs$variable ==     ">75%"), "value"] and     ">75%"), "value"]
#t = -1.0492, df = 3.1303, p-value = 0.8159
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.004968579          Inf
#sample estimates:
#  mean of x   mean of y 
#0.002196667 0.003746667
t.test(x = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Prenatal:Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Cytosol" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = 1.5401, df = 2.6922, p-value = 0.1156
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.01129211         Inf
#sample estimates:
#  mean of x  mean of y 
#0.03090333 0.01214667 


## Is IR increased in adults?
# Main effect: All filtered introns
t.test(x = IRratios[which(IRratios$Age=="Adult"),"IRratio"],
       y = IRratios[which(IRratios$Age=="Prenatal"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Age == "Adult"), "IRratio"] and IRratios[which(IRratios$Age == "Prenatal"), "IRratio"]
#t = 7.572, df = 323310, p-value = 1.843e-14
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.0004424013          Inf
#sample estimates:
#  mean of x   mean of y 
#0.006955816 0.006390642

# in Nuclear polyA
t.test(x = IRratios[which(IRratios$Group=="Adult:Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Nucleus"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Group == "Adult:Nucleus"), "IRratio"] and IRratios[which(IRratios$Group == "Prenatal:Nucleus"), "IRratio"]
#t = 13.864, df = 166000, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.001518024         Inf
#sample estimates:
#  mean of x   mean of y 
#0.009611667 0.007889297 

# in Cytosol polyA
t.test(x = IRratios[which(IRratios$Group=="Adult:Cytosol"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Cytosol"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Group == "Adult:Cytosol"), "IRratio"] and IRratios[which(IRratios$Group == "Prenatal:Cytosol"), "IRratio"]
#t = -12.575, df = 182850, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.001060624          Inf
#sample estimates:
#  mean of x   mean of y 
#0.003928916 0.004866854
t.test(x = IRratios[which(IRratios$Group=="Adult:Cytosol"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Cytosol"),"IRratio"],
       alternative = "less")
#data:  IRratios[which(IRratios$Group == "Adult:Cytosol"), "IRratio"] and IRratios[which(IRratios$Group == "Prenatal:Cytosol"), "IRratio"]
#t = -12.575, df = 182850, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.0008152517
#sample estimates:
#  mean of x   mean of y 
#0.003928916 0.004866854

# Main effect: introns >=75% retained
head(PercentIRs)
# in all polyA
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">75%"),"value"],
       y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">75%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Age == "Adult" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Age == "Prenatal" & PercentIRs$variable ==     ">75%"), "value"] and     ">75%"), "value"]
#t = 0.49039, df = 6.7169, p-value = 0.3197
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.002724383          Inf
#sample estimates:
#  mean of x   mean of y 
#0.003915000 0.002971667
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Age == "Adult" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Age == "Prenatal" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = 0.81858, df = 7.8679, p-value = 0.2186
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.01467262         Inf
#sample estimates:
#  mean of x  mean of y 
#0.03301833 0.02152500
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">50%"),"value"],
       alternative = "less")
#data:  PercentIRs[which(PercentIRs$Age == "Adult" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Age == "Prenatal" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = 0.81858, df = 7.8679, p-value = 0.7814
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.03765929
#sample estimates:
#  mean of x  mean of y 
#0.03301833 0.02152500

# in Nuclear polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">75%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">75%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Nucleus" & PercentIRs$variable ==     ">75%"), "value"] and     ">75%"), "value"]
#t = 1.0991, df = 2.8076, p-value = 0.1785
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.004087901          Inf
#sample estimates:
#  mean of x   mean of y 
#0.005596667 0.002196667
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Nucleus" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Nucleus" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = 1.2954, df = 3.6621, p-value = 0.1354
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.01700336         Inf
#sample estimates:
#  mean of x  mean of y 
#0.05553333 0.03090333

# in Cytosol polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">75%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">75%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Cytosol" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Cytosol" & PercentIRs$variable ==     ">75%"), "value"] and     ">75%"), "value"]
#t = -0.64506, df = 2.4095, p-value = 0.7125
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.007632187          Inf
#sample estimates:
#  mean of x   mean of y 
#0.002233333 0.003746667
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Cytosol" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Cytosol" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = -0.30372, df = 3.109, p-value = 0.6097
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.01419662         Inf
#sample estimates:
#  mean of x  mean of y 
#0.01050333 0.01214667
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "less")
#data:  PercentIRs[which(PercentIRs$Group == "Adult:Cytosol" & PercentIRs$variable ==  and PercentIRs[which(PercentIRs$Group == "Prenatal:Cytosol" & PercentIRs$variable ==     ">50%"), "value"] and     ">50%"), "value"]
#t = -0.30372, df = 3.109, p-value = 0.3903
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.01090996
#sample estimates:
#  mean of x  mean of y 
#0.01050333 0.01214667

### Differential IR by group
# read in results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age",
          "Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
IRcomp = rat.1 = nonconst.66warn = nonconst.nowarn = nonconst = list()
for (i in 1:length(comps)){
  IRcomp[[i]] = read.table(paste0(path, comps[i], ".tab"), header = TRUE, comment.char="#")
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
  rat.1[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_0.1.tab"), header = TRUE, comment.char="#")
  nonconst.nowarn[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.nowarn.tab"), header = TRUE, comment.char="#")
  nonconst.66warn[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.66warn.tab"), header = TRUE, comment.char="#")
}
names(IRcomp) = names(nonconst) = names(rat.1) = names(nonconst.nowarn) = names(nonconst.66warn) = comps
dIR = list(IRcomp=IRcomp, nonconst=nonconst, rat.1=rat.1, nonconst.nowarn=nonconst.nowarn, nonconst.66warn=nonconst.66warn)
              
elementNROWS(IRcomp)
elementNROWS(nonconst)
elementNROWS(rat.1)
elementNROWS(nonconst.nowarn)
elementNROWS(nonconst.66warn)

string = lapply(dIR, function(x) lapply(x, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),
                                                                   "/", fixed = TRUE),recursive = FALSE)))
genes = lapply(string, function(x) lapply(x, function(y) y[grep("ENSG", y)]))
comments = lapply(string, function(x) lapply(x, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])))
IR.diff = lapply(dIR, function(x) lapply(x, function(y) y$A.IRratio - y$B.IRratio))
Sign = lapply(IR.diff, function(x) lapply(x, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")))
IR = list()
for (i in 1:length(dIR)){
  IR[[i]] = Map(cbind, dIR[[i]], ensID = genes[[i]], comments = comments[[i]], IR.diff = IR.diff[[i]], Sign = Sign[[i]])
}
names(IR) = names(dIR)
dIR = IR
IRclean = lapply(dIR, function(x) lapply(x, function(y) 
                  y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & 
                            y$B.warnings!="NonUniformIntronCover" & y$B.warnings!="LowCover" & 
                            y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                            y$comments=="clean"),]))


# Explore the results
lapply(IRclean, function(x) elementNROWS(x))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1097              1001               689              1214               136               168 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1097              1001               689              1214               136               168 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#785               661               421               877               131               155 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#184               330                94               237                 5                 3 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#646               710               409               694                56                57
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1087               975               576               595               136               104 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#1087               975               576               595               136               104 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#776               641               334               380               131                95 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#182               322                80               137                 5                 3 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#638               698               344               355                56                38 
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$p.diff<=0.05),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#314               213               171               351                99                95 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#314               213               171               351                99                95 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#268               192               141               323                97                93 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#57                30                18                48                 3                 2 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#203               134                83               191                42                36
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$p.diff<=0.05 & y$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#312               212               140               165                99                62 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#312               212               140               165                99                62 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#266               191               112               145                97                60 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#56                30                14                23                 3                 2 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#201               134                67                91                42                25
lapply(IRclean, function(x) elementNROWS(lapply(x, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 3                16                 3                 7 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 3                16                 3                 7 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 3                16                 3                 7 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#0                 0                 0                 1                 0                 0 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2                 1                 0                 9                 0                 1      
lapply(IRclean, function(x) elementNROWS(lapply(lapply(x, function(y) y[which(y$A.IRratio>=0.5 | y$B.IRratio>=0.5),]), 
                                                function(z) z[which(z$Sign=="MoreIRInNuc.Fetal"),])))
#IRcomp
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 1                 6                 3                 5 
#nonconst
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 1                 6                 3                 5 
#rat.1
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#8                 5                 1                 6                 3                 5 
#nonconst.nowarn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#0                 0                 0                 1                 0                 0 
#nonconst.66warn
#Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
#2                 1                 0                 4                 0                 1

#moving forward with the nonconstitutively spliced results because they capture more introns, 
  #and coverage issues are filtered after differential retention calculation

# Comparison of significantly vs nonsignificantly retained introns by zone/age
fisher.test(data.frame(c(312,2),c(775,8))) 
# adult zone
#p-value = 0.7332
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3188043 15.6435734
#sample estimates:
#  odds ratio 
#1.609729
fisher.test(data.frame(c(212,1),c(763,25)))
# fetal zone
#p-value = 0.02666
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.122308 286.209119
#sample estimates:
#  odds ratio 
#6.938855
fisher.test(data.frame(c(141,30),c(435,83)))
# cytosol age
#p-value = 0.6352
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5572585 1.4726430
#sample estimates:
#  odds ratio 
#0.8969333 
fisher.test(data.frame(c(165,186),c(430,433)))
# nucleus age
#p-value = 0.3762
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6912077 1.1541152
#sample estimates:
#  odds ratio 
#0.8933697

## Do DEG by fraction contain introns with higher IR ratios than non-DEG by fraction?
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
allgenes = as.character(rownames(FracList[[1]]))
sigFracIRratio = list(list(), list(), list(),list(),list(), list(),list(), list(),list(), list()) 
for (i in 1:length(SigList)){
  for (j in 1:length(IRfiltered2)){
    filt = IRfiltered2[[j]][,c(22,23,20)]
    sigFracIRratio[[i]][[j]] = filt[which(as.character(SigList[[i]][,"ensemblID"]) %in% as.character(filt$genes)),]
  }}
names(sigFracIRratio) = names(SigList)
for (i in 1:length(sigFracIRratio)){names(sigFracIRratio[[i]]) = names(IRfiltered2)}
lapply(sigFracIRratio, elementNROWS)
sigFracIRratio = lapply(sigFracIRratio, function(x) lapply(x, function(y) data.table(y, key="genes")))
sigFracIRratio = lapply(sigFracIRratio, function(x) lapply(x, function(y) data.frame(y[, list(IRratio=max(IRratio)), by="genes"]))) # limit to intron with highest IR ratio per gene per sample
sigFracIRratio = lapply(sigFracIRratio, function(x) do.call(rbind,x))
elementNROWS(sigFracIRratio)

# Compare the intron with the greatest retention per gene per sample in different sets of DEGs
t.test(sigFracIRratio[["Apres.DownNuc"]][,"IRratio"], sigFracIRratio[["Apres.UpNuc"]][,"IRratio"], alternative = "less")
#data:  sigFracIRratio[["Apres.DownNuc"]][, "IRratio"] and sigFracIRratio[["Apres.UpNuc"]][, "IRratio"]
#t = -0.77008, df = 2768, p-value = 0.2207
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.00141993
#sample estimates:
#  mean of x  mean of y 
#0.01539956 0.01664877 
t.test(sigFracIRratio[["Arres.DownNuc"]][,"IRratio"], sigFracIRratio[["Arres.UpNuc"]][,"IRratio"], alternative = "two.sided")
#data:  sigFracIRratio[["Arres.DownNuc"]][, "IRratio"] and sigFracIRratio[["Arres.UpNuc"]][, "IRratio"]
#t = 3.6311, df = 1971.5, p-value = 0.0002894
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.002722987 0.009118719
#sample estimates:
#  mean of x  mean of y 
#0.01730208 0.01138122  

# Get the DEG Fraction p-value and LFC sign for differentially retained introns
nonconst = IRclean[["nonconst"]]
elementNROWS(lapply(nonconst, function(x) x[which(x$p.diff<=0.05),]))
nonconst = Map(cbind, nonconst, 
               AP.sig = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"padj"]),
               AP.LFC = lapply(nonconst, function(x) FracList[["Apres"]][match(x$ensID, FracList[["Apres"]][,"ensemblID"]),"log2FoldChange"]),
               FP.sig = lapply(nonconst, function(x) FracList[["Fpres"]][match(x$ensID, FracList[["Fpres"]][,"ensemblID"]),"padj"]),
               FP.LFC = lapply(nonconst, function(x) FracList[["Fpres"]][match(x$ensID, FracList[["Fpres"]][,"ensemblID"]),"log2FoldChange"]),
               AR.sig = lapply(nonconst, function(x) FracList[["Arres"]][match(x$ensID, FracList[["Arres"]][,"ensemblID"]),"padj"]),
               AR.LFC = lapply(nonconst, function(x) FracList[["Arres"]][match(x$ensID, FracList[["Arres"]][,"ensemblID"]),"log2FoldChange"]),
               FR.sig = lapply(nonconst, function(x) FracList[["Frres"]][match(x$ensID, FracList[["Frres"]][,"ensemblID"]),"padj"]),
               FR.LFC = lapply(nonconst, function(x) FracList[["Frres"]][match(x$ensID, FracList[["Frres"]][,"ensemblID"]),"log2FoldChange"]),
               FP.down.sig = lapply(nonconst, function(x) FracList[["Fpres.down"]][match(x$ensID, FracList[["Fpres.down"]][,"ensemblID"]),"padj"]),
               FP.down.LFC = lapply(nonconst, function(x) FracList[["Fpres.down"]][match(x$ensID, FracList[["Fpres.down"]][,"ensemblID"]),"log2FoldChange"]))
elementNROWS(nonconst)
lapply(nonconst, head)

# Are significantly dIR genes more likely to be significantly DEG by fraction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]>=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]<=0.05),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]>=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]<=0.05),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.down.sig"]>=0.05),])
fisher.test(data.frame(c(164,148), c(286,489))) #adult
#data:  data.frame(c(164, 148), c(286, 489))
#p-value = 2.539e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.439242 2.493470
#sample estimates:
#  odds ratio 
#1.893514
fisher.test(data.frame(c(42,167), c(47,717))) #prenatal
#data:  data.frame(c(42, 167), c(47, 717))
#p-value = 1.507e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.379715 6.151949
#sample estimates:
#  odds ratio 
#3.829787

# Of dIR Fraction genes, does the LFC go in the same direction?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]>0),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"AP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]<0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]>0),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"FP.LFC"]<0),])
fisher.test(data.frame(c(143,69), c(447,337))) #prenatal
#data:  data.frame(c(143, 69), c(447, 337))
#p-value = 0.007306
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.122700 2.186824
#sample estimates:
#  odds ratio 
#1.561769
fisher.test(data.frame(c(151,162), c(439,341))) #adult
#data:  data.frame(c(151, 162), c(439, 341))
#p-value = 0.01873
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5517484 0.9501778
#sample estimates:
#  odds ratio 
#0.7242496

# Of dIR Fraction genes, does the direction of dIR more in nucleus?
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Adult_PolyA_Zone"]][which(nonconst[["Adult_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Adult_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]<=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Fetal_PolyA_Zone"]][which(nonconst[["Fetal_PolyA_Zone"]][,"p.diff"]>=0.05 & nonconst[["Fetal_PolyA_Zone"]][,"Sign"]=="MoreIRInCyt.Adult"),])
fisher.test(data.frame(c(312,2), c(775,8))) #adult
#data:  data.frame(c(312, 2), c(775, 8))
#p-value = 0.7332
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3188043 15.6435734
#sample estimates:
#  odds ratio 
#1.609729
fisher.test(data.frame(c(212,1), c(763,25))) #prenatal
#data:  data.frame(c(212, 1), c(763, 25))
#p-value = 0.02666
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.122308 286.209119
#sample estimates:
#  odds ratio 
#6.938855


## Do DEG by age contain introns with higher IR ratios than non-DEG by age?
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
allgenes = as.character(rownames(AgeList[[1]]))
sigAgeIRratio = list(list(), list(), list(),list(),list(), list(),list(), list(),list(), list()) 
for (i in 1:length(SigList)){
  for (j in 1:length(IRfiltered2)){
    filt = IRfiltered2[[j]][,c(22,23,20)]
    sigAgeIRratio[[i]][[j]] = filt[which(as.character(SigList[[i]][,"ensemblID"]) %in% as.character(filt$genes)),]
  }}
names(sigAgeIRratio) = names(SigList)
for (i in 1:length(sigAgeIRratio)){names(sigAgeIRratio[[i]]) = names(IRfiltered2)}
lapply(sigAgeIRratio, elementNROWS)
sigAgeIRratio = lapply(sigAgeIRratio, function(x) lapply(x, function(y) data.table(y, key="genes")))
sigAgeIRratio = lapply(sigAgeIRratio, function(x) lapply(x, function(y) data.frame(y[, list(IRratio=max(IRratio)), by="genes"]))) # limit to intron with highest IR ratio per gene per sample
sigAgeIRratio = lapply(sigAgeIRratio, function(x) do.call(rbind,x))
elementNROWS(sigAgeIRratio)

# Compare the intron with the greatest retention per gene per sample in different sets of DEGs
t.test(sigAgeIRratio[["Cpres.DownPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.UpPrenatal"]][,"IRratio"], alternative = "less")
#data:  sigAgeIRratio[["Cpres.DownPrenatal"]][, "IRratio"] and sigAgeIRratio[["Cpres.UpPrenatal"]][, "IRratio"]
#t = -4.6249, df = 11983, p-value = 1.893e-06
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.001788791
#sample estimates:
#  mean of x  mean of y 
#0.01217757 0.01495380 
t.test(sigAgeIRratio[["Npres.DownPrenatal"]][,"IRratio"], sigAgeIRratio[["Npres.UpPrenatal"]][,"IRratio"], alternative = "less")
#data:  sigAgeIRratio[["Npres.DownPrenatal"]][, "IRratio"] and sigAgeIRratio[["Npres.UpPrenatal"]][, "IRratio"]
#t = -5.3515, df = 10065, p-value = 4.457e-08
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.002540159
#sample estimates:
#  mean of x  mean of y 
#0.01248776 0.01615528
t.test(sigAgeIRratio[["Cpres.down.DownPrenatal"]][,"IRratio"], sigAgeIRratio[["Cpres.down.UpPrenatal"]][,"IRratio"], alternative = "less")
#data:  sigAgeIRratio[["Cpres.down.DownPrenatal"]][, "IRratio"] and sigAgeIRratio[["Cpres.down.UpPrenatal"]][, "IRratio"]
#t = -5.3569, df = 10730, p-value = 4.319e-08
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.002534245
#sample estimates:
#  mean of x  mean of y 
#0.01280146 0.01645878

# Get the DEG Age p-value and LFC sign for differentially retained introns
nonconst = Map(cbind, nonconst, 
               CP.sig = lapply(nonconst, function(x) AgeList[["Cpres"]][match(x$ensID, AgeList[["Cpres"]][,"ensemblID"]),"padj"]),
               CP.LFC = lapply(nonconst, function(x) AgeList[["Cpres"]][match(x$ensID, AgeList[["Cpres"]][,"ensemblID"]),"log2FoldChange"]),
               NP.sig = lapply(nonconst, function(x) AgeList[["Npres"]][match(x$ensID, AgeList[["Npres"]][,"ensemblID"]),"padj"]),
               NP.LFC = lapply(nonconst, function(x) AgeList[["Npres"]][match(x$ensID, AgeList[["Npres"]][,"ensemblID"]),"log2FoldChange"]),
               CR.sig = lapply(nonconst, function(x) AgeList[["Crres"]][match(x$ensID, AgeList[["Crres"]][,"ensemblID"]),"padj"]),
               CR.LFC = lapply(nonconst, function(x) AgeList[["Crres"]][match(x$ensID, AgeList[["Crres"]][,"ensemblID"]),"log2FoldChange"]),
               NR.sig = lapply(nonconst, function(x) AgeList[["Nrres"]][match(x$ensID, AgeList[["Nrres"]][,"ensemblID"]),"padj"]),
               NR.LFC = lapply(nonconst, function(x) AgeList[["Nrres"]][match(x$ensID, AgeList[["Nrres"]][,"ensemblID"]),"log2FoldChange"]),
               CP.down.sig = lapply(nonconst, function(x) AgeList[["Cpres.down"]][match(x$ensID, AgeList[["Cpres.down"]][,"ensemblID"]),"padj"]),
               CP.down.LFC = lapply(nonconst, function(x) AgeList[["Cpres.down"]][match(x$ensID, AgeList[["Cpres.down"]][,"ensemblID"]),"log2FoldChange"]))
elementNROWS(nonconst)
lapply(nonconst, head)

# Are significantly dIR genes more likely to be significantly DEG by age?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]>=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]<=0.05),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]>=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]<=0.05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.sig"]>=0.05),])
fisher.test(data.frame(c(118,50), c(302,207))) #cytosol
#data:  data.frame(c(164, 148), c(286, 489))
#data:  data.frame(c(118, 50), c(302, 207))
#p-value = 0.0132
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.096525 2.406682
#sample estimates:
#  odds ratio 
#1.61649
fisher.test(data.frame(c(235,112), c(513,339))) #nucleus
#data:  data.frame(c(235, 112), c(513, 339))
#p-value = 0.01509
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.057376 1.823272
#sample estimates:
#  odds ratio 
#1.386161

# Of dIR Age genes, does the LFC go in the same direction?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]>0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]<0),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"CP.down.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]<0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]>0),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]<05),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"NP.LFC"]>0),])
fisher.test(data.frame(c(109,59), c(248,261))) #cytosol
#data:  data.frame(c(109, 59), c(248, 261))
#p-value = 0.0003479
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.335597 2.844027
#sample estimates:
#  odds ratio 
#1.942424
fisher.test(data.frame(c(200,148), c(853,392))) #nucleus
#data:  data.frame(c(200, 148), c(853, 392))
#p-value = 0.0001536
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4833510 0.7992836
#sample estimates:
#  odds ratio 
#0.6211912

# Of dIR Age genes, does the direction of dIR more in adult?
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Cytosol_PolyA_Age"]][which(nonconst[["Cytosol_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Cytosol_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]<=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInNuc.Fetal"),])
dim(nonconst[["Nuclear_PolyA_Age"]][which(nonconst[["Nuclear_PolyA_Age"]][,"p.diff"]>=0.05 & nonconst[["Nuclear_PolyA_Age"]][,"Sign"]=="MoreIRInCyt.Adult"),])
fisher.test(data.frame(c(140,31), c(436,82))) #cytosol
#data:  data.frame(c(140, 31), c(436, 82))
#p-value = 0.4768
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5296991 1.3881400
#sample estimates:
#  odds ratio 
#0.8495898
fisher.test(data.frame(c(165,186), c(430,433))) #nucleus
#data:  data.frame(c(165, 186), c(430, 433))
#p-value = 0.3762
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6912077 1.1541152
#sample estimates:
#  odds ratio 
#0.8933697