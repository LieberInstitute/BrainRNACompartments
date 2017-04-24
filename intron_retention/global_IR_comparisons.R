library(ggplot2)
library(reshape2)

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