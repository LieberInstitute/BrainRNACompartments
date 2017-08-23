library(ggplot2)
library(reshape2)
library(GenomicRanges)

# Load IRFinder Results
names = scan("./Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
}
names(IRres) = shortenedNames

# Filter introns
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings=="-" | x$Warnings=="MinorIsoform"),])
IRfiltered = lapply(IRfiltered, function(x) x[grep("clean", x$GeneIntronDetails, fixed=T),])
string = lapply(IRfiltered, function(x) strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE))
string = lapply(string, function(x) unlist(x, recursive = F))
IRfiltered = Map(cbind, IRfiltered, genes = lapply(string, function(x) x[grep("ENSG", x)]), 
                 intronID = lapply(IRfiltered, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
head(IRfiltered[[1]])

## Look at distribution of filtered introns by group
introns = data.frame(num.introns = elementNROWS(IRfiltered))
introns$SampleID = as.factor(rownames(introns))
introns$rownum = c(1:nrow(introns))
introns$Fraction = ifelse(introns$rownum %in% grep("N",introns$SampleID), "Nucleus","Cytosol")
introns$Age = ifelse(introns$rownum %in% grep("Br53",introns$SampleID), "Prenatal","Adult")
introns$Group = factor(paste0(introns$Age,":",introns$Fraction), levels = c("Adult:Cytosol","Prenatal:Cytosol", "Adult:Nucleus","Prenatal:Nucleus"))
introns$num.genes = as.numeric(lapply(lapply(string, function(x) x[grep("ENSG", x)]), function(x) length(unique(x))))
introns$MeanIntronDepth = as.numeric(lapply(IRfiltered, function(x) mean(x$IntronDepth)))
write.csv(introns, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/global_IR_comparisons/filtered.intron.stats.csv")

pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/introns_passingQC.pdf", width =8,height = 5)
ggplot(introns, aes(x=Group, y=num.introns)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + xlab("") +
  ggtitle("Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,60000)
dev.off()
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/unique_genes_introns_passingQC.pdf", width =8,height = 5)
ggplot(introns, aes(x=Group, y=num.genes)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + xlab("") +
  ggtitle("Unique Genes Containing Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,8000)
dev.off()
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/meanIntronDepth_introns_passingQC.pdf", width =8,height = 5)
ggplot(introns, aes(x=Group, y=MeanIntronDepth)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + 
  xlab("") +
  ggtitle("Mean Intron Depth (Passing QC)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,.35)
dev.off()

allintrons = Reduce(intersect, lapply(IRfiltered, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
length(allintrons) # 3889 of the introns pass QC in all four groups


### Global IR Comparisons:
# What does the distribution of IR Ratios look like?

IRratios = lapply(IRfiltered, function(x) data.frame(IRratio=x$IRratio,intronID=x$intronID))
IRratios = Map(cbind, IRratios, Age=NA, Fraction=NA)
for (i in 1:6){IRratios[[i]][,"Age"] = "Adult"}
for (i in 7:12){IRratios[[i]][,"Age"] = "Prenatal"}
for (i in c(2,4,6,8,10,12)){IRratios[[i]][,"Fraction"] = "Nucleus"}
for (i in c(1,3,5,7,9,11)){IRratios[[i]][,"Fraction"] = "Cytosol"}
IRratios = do.call(rbind, IRratios)
IRratios$Group = factor(paste(IRratios$Age, IRratios$Fraction, sep=":"), levels = c("Adult:Cytosol", "Prenatal:Cytosol",
                                                                                    "Adult:Nucleus", "Prenatal:Nucleus"))

pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/IRratio_byGroup.pdf", width =6,height = 5)
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
dev.off()

# Just in shared introns
ratio.overlaps = IRratios[which(IRratios$intronID %in% allintrons),]
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/IRratio_byGroup_QC_introns_only.pdf", width =6.25,height = 5)
ggplot(ratio.overlaps, aes(x=IRratio)) + geom_density(aes(group=Group, colour=Group)) +
  ylab("") + xlab("IR Ratio") +
  ggtitle("IR Ratios By Group\n(Introns Passing QC In All Groups)") +
  xlim(0,.06) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

# How many introns are retained at different thresholds?
PercentIRs = data.frame(Constitutively.Spliced = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio==0),]))),
                 ">0%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>0),]))),
                 "less.five" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.05),]))),
                 ">10%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.10),]))),
                 ">20%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.20),]))),
                 ">30%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.30),]))),
                 ">40%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.40),]))),
                 "fifty" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.50),]))),
                 ">75%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.75),]))),
                 ">90%" = unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.9),]))))
Total = unlist(lapply(IRfiltered, function(x) nrow(x)))
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
(PercentIRs[which(PercentIRs$variable==">0%"),"value"]-PercentIRs[which(PercentIRs$variable==">5%"),"value"])[order(PercentIRs[which(PercentIRs$variable==">0%"),"value"]-PercentIRs[which(PercentIRs$variable==">5%"),"value"])]
# range of 34.42-56.11% have IR ratio between 0-5%
match(PercentIRs[which(PercentIRs$variable==">0%"),"SampleID"],PercentIRs[which(PercentIRs$variable==">5%"),"SampleID"])

pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio.pdf", 
    width =11,height = 5)
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
dev.off()
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio_greaterthan5perc.pdf", 
    width =8,height = 5)
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
dev.off()

# In only the introns in all four groups
overlaps = lapply(IRfiltered, function(x) x[which(x$intronID %in% allintrons),])
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
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio_QC_intronsOnly.pdf", 
    width =11,height = 5)
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
dev.off() 
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio_greaterthan5perc_QC_introns.pdf", 
    width =8,height = 5)
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
dev.off()

# Is IR increased in the nucleus?
# Main effect: All filtered introns
t.test(x = IRratios[which(IRratios$Fraction=="Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Fraction=="Cytosol"),"IRratio"],
       alternative = "greater")
#t = 59.932, df = 391890, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.004549564         Inf
#sample estimates:
#  mean of x   mean of y 
#0.009639673 0.004961720

# in adult polyA
t.test(x = IRratios[which(IRratios$Group=="Adult:Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Adult:Cytosol"),"IRratio"],
       alternative = "greater")
#data:  IRratios[which(IRratios$Group == "Adult:Nucleus"), "IRratio"] and IRratios[which(IRratios$Group == "Adult:Cytosol"), "IRratio"]
#t = 51.34, df = 134720, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.006244051         Inf
#sample estimates:
#  mean of x   mean of y 
#0.010631792 0.004181069 

# in Prenatal polyA
t.test(x = IRratios[which(IRratios$Group=="Prenatal:Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Cytosol"),"IRratio"],
       alternative = "greater")
#t = 36.518, df = 254640, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.003453338         Inf
#sample estimates:
#  mean of x   mean of y 
#0.009015836 0.005399616

# Main effect: introns >=50% retained
head(PercentIRs)
# in all polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$variable == ">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$variable == ">50%"),"value"],
       alternative = "greater")
#t = 3.6862, df = 5.5616, p-value = 0.005872
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.02269485        Inf
#sample estimates:
#  mean of x  mean of y 
#0.06533000 0.01655667

# in adult polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#t = 3.2264, df = 2.0066, p-value = 0.04188
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.006040032         Inf
#sample estimates:
#  mean of x  mean of y 
#0.07532667 0.01305333

# in Prenatal polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#t = 1.7847, df = 2.3756, p-value = 0.09795
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.01667134         Inf
#sample estimates:
#  mean of x  mean of y 
#0.05533333 0.02006000 


## Is IR increased in adults?
# Main effect: All filtered introns
t.test(x = IRratios[which(IRratios$Age=="Adult"),"IRratio"],
       y = IRratios[which(IRratios$Age=="Prenatal"),"IRratio"],
       alternative = "greater")
#t = 4.8165, df = 349850, p-value = 7.308e-07
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.00026177        Inf
#sample estimates:
#  mean of x   mean of y 
#0.007624388 0.007226860 

# in Nuclear polyA
t.test(x = IRratios[which(IRratios$Group=="Adult:Nucleus"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Nucleus"),"IRratio"],
       alternative = "greater")
#t = 11.717, df = 179990, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.001389112         Inf
#sample estimates:
#  mean of x   mean of y 
#0.010631792 0.009015836  

# in Cytosol polyA
t.test(x = IRratios[which(IRratios$Group=="Adult:Cytosol"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Cytosol"),"IRratio"],
       alternative = "greater")
#t = -15.029, df = 198370, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.001351909          Inf
#sample estimates:
#  mean of x   mean of y 
#0.004181069 0.005399616

t.test(x = IRratios[which(IRratios$Group=="Adult:Cytosol"),"IRratio"],
       y = IRratios[which(IRratios$Group=="Prenatal:Cytosol"),"IRratio"],
       alternative = "less")
#t = -15.029, df = 198370, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.001085186
#sample estimates:
#  mean of x   mean of y 
#0.004181069 0.005399616

# Main effect: introns >=50% retained
head(PercentIRs)
# in all polyA
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#t = 0.32118, df = 9.106, p-value = 0.3777
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.03051811         Inf
#sample estimates:
#  mean of x  mean of y 
#0.04419000 0.03769667
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">50%"),"value"],
       alternative = "less")
#t = 0.32118, df = 9.106, p-value = 0.6223
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.04350477
#sample estimates:
#  mean of x  mean of y 
#0.04419000 0.03769667

# in Nuclear polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#t = 0.74063, df = 3.9983, p-value = 0.25
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.03756283         Inf
#sample estimates:
#  mean of x  mean of y 
#0.07532667 0.05533333

# in Cytosol polyA
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "greater")
#t = -1.1944, df = 2.0723, p-value = 0.8244
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.02373551         Inf
#sample estimates:
#  mean of x  mean of y 
#0.01305333 0.02006000
t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosol" & PercentIRs$variable==">50%"),"value"],
       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosol" & PercentIRs$variable==">50%"),"value"],
       alternative = "less")
#t = -1.1944, df = 2.0723, p-value = 0.1756
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf 0.009722174
#sample estimates:
#  mean of x  mean of y 
#0.01305333 0.02006000