library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(data.table)

# Load IRFinder Results
names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
}
names(IRres) = shortenedNames

# Filter introns
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings!="NonUniformIntronCover"),])
IRfiltered = lapply(IRfiltered, function(x) x[grep("clean", x$GeneIntronDetails, fixed=T),])
IRfiltered = lapply(IRfiltered, function(x) x[max(x$SplicesExact,x$SplicesRight,x$SplicesLeft)>4 | (x$ExonToIntronReadsLeft>4 & x$ExonToIntronReadsRight>4),])
IRfiltered = Map(cbind, IRfiltered, 
                 genes = lapply(lapply(IRfiltered, function(x) unlist(strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE), recursive = F)), function(x) x[grep("ENSG", x)]), 
                 intronID = lapply(IRfiltered, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
elementNROWS(IRfiltered)
head(IRfiltered[[1]])


## Look at distribution of filtered introns by group

introns = data.frame(num.introns = elementNROWS(IRfiltered))
introns$SampleID = as.factor(rownames(introns))
introns$rownum = c(1:nrow(introns))
introns$Fraction = ifelse(introns$rownum %in% grep("N",introns$SampleID), "Nucleus","Cytosplasm")
introns$Age = ifelse(introns$rownum %in% grep("Br53",introns$SampleID), "Prenatal","Adult")
introns$Group = factor(paste0(introns$Age,"\n",introns$Fraction), levels = c("Adult\nCytosplasm","Prenatal\nCytosplasm", "Adult\nNucleus","Prenatal\nNucleus"))
introns$num.genes = as.numeric(lapply(lapply(lapply(IRfiltered, function(x) unlist(strsplit(as.character(x$GeneIntronDetails), "/", fixed = TRUE), recursive = F)), 
                                             function(x) x[grep("ENSG", x)]), function(x) length(unique(x))))
introns$MeanIntronDepth = as.numeric(lapply(IRfiltered, function(x) mean(x$IntronDepth)))
write.csv(introns, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/global_IR_comparisons/filtered.intron.stats.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/global_IR_comparisons/filtered.intron.stats.csv")
range(df$num.introns) # 166661 173125
range(df$num.genes) # 15345 15389


pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/introns_passingQC.pdf", width =7,height = 5)
ggplot(introns, aes(x=Group, y=num.introns)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + xlab("") +
  ggtitle("Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")
dev.off()
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/unique_genes_introns_passingQC.pdf", width =7,height = 5)
ggplot(introns, aes(x=Group, y=num.genes)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + xlab("") +
  ggtitle("Unique Genes Containing Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="")  
dev.off()
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/meanIntronDepth_introns_passingQC.pdf", width =7,height = 5)
ggplot(introns, aes(x=Group, y=MeanIntronDepth)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + 
  xlab("") +
  ggtitle("Mean Intron Depth (Passing QC)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") + ylim(0,.25)
dev.off()

allintrons = Reduce(intersect, lapply(IRfiltered, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End,"(",x$Direction,")")))
length(allintrons) # 152432 of the introns pass QC in all four groups



### Global IR Comparisons:
# What does the distribution of IR Ratios look like?

IRratios = do.call(rbind, Map(cbind, lapply(IRfiltered, function(x) data.frame(IRratio=x$IRratio,intronID=x$intronID)), SampleID = as.list(names(IRfiltered))))
IRratios$rnum = 1:nrow(IRratios)
IRratios$Age = ifelse(IRratios$rnum %in% grep("53", IRratios$SampleID), "Prenatal","Adult")
IRratios$Fraction = ifelse(IRratios$rnum %in% grep("C", IRratios$SampleID), "Cytosplasm","Nucleus")
IRratios$Group = factor(paste(IRratios$Age, IRratios$Fraction, sep=":"), levels = c("Adult:Cytosplasm","Prenatal:Cytosplasm","Adult:Nucleus","Prenatal:Nucleus"))

pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/IRratio_byGroup.pdf", width =6,height = 5)
ggplot(IRratios, aes(x=IRratio)) + geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("IR Ratio") +
  ggtitle("IR Ratios By Group") +
  xlim(0,.06) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.75, 0.55)) +
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
  theme(legend.position = c(0.75, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


# How many introns are retained at different thresholds?
PercentIRs = data.frame(unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio==0),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>0),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.05),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.10),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.20),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.30),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.40),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.50),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.75),]))),
                        unlist(lapply(IRfiltered, function(x) nrow(x[which(x$IRratio>=.9),]))))
Total = elementNROWS(IRfiltered)
PercentIRs = cbind(round(PercentIRs / Total *100, digits = 2), SampleID = rownames(PercentIRs), 
                   Age = c(rep.int("Adult", 6), rep.int("Prenatal", 6)), Fraction = rep.int(c("Cytosplasm", "Nucleus"), 6))
colnames(PercentIRs) = c("Constitutively\nSpliced",">0%",">5%",">10%",">20%",">30%",">40%",">50%",">75%",">90%","SampleID","Age","Fraction")
PercentIRs = melt(PercentIRs)
PercentIRs$Group = factor(paste(PercentIRs$Age, PercentIRs$Fraction, sep=":"), 
                          levels = c("Adult:Cytosplasm", "Prenatal:Cytosplasm", "Adult:Nucleus", "Prenatal:Nucleus"))
PercentIRs$variable = factor(PercentIRs$variable, levels=c("Constitutively\nSpliced",">0%",">5%",">10%",">20%",">30%",">40%",">50%",">75%",">90%"))
(PercentIRs[PercentIRs$variable=="Constitutively\nSpliced","value"])[order(PercentIRs[PercentIRs$variable=="Constitutively\nSpliced","value"])]
# range of 58.68-85.33% are constitutively spliced
(PercentIRs[PercentIRs$variable==">0%","value"]-PercentIRs[PercentIRs$variable==">5%","value"])[order(PercentIRs[PercentIRs$variable==">0%","value"]-PercentIRs[PercentIRs$variable==">5%","value"])]
# range of 12.20-34.63% have 0% < IR ratio < 5%


pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio.pdf",width =11,height = 5)
ggplot(PercentIRs, aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.8, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio_greaterthan5perc.pdf",width =8,height = 5)
ggplot(PercentIRs[which(PercentIRs$variable != "Constitutively\nSpliced" & PercentIRs$variable != ">0%"),],
       aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.8, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


# In only the introns in all four groups
overlaps = lapply(IRfiltered, function(x) x[which(x$intronID %in% allintrons),])
PercentIRsOverlap = data.frame(unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio==0),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>0),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.05),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.10),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.20),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.30),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.40),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.50),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.75),]))),
                               unlist(lapply(overlaps, function(x) nrow(x[which(x$IRratio>=.9),]))))
PercentIRsOverlap = cbind(round(PercentIRsOverlap / length(allintrons) * 100, digits = 2), SampleID = rownames(PercentIRsOverlap), 
                   Age = c(rep.int("Adult", 6), rep.int("Prenatal", 6)), Fraction = rep.int(c("Cytosplasm", "Nucleus"), 6))
colnames(PercentIRsOverlap) = c("Constitutively\nSpliced",">0%",">5%",">10%",">20%",">30%",">40%",">50%",">75%",">90%","SampleID","Age","Fraction")
PercentIRsOverlap = melt(PercentIRsOverlap)
PercentIRsOverlap$Group = factor(paste(PercentIRsOverlap$Age, PercentIRsOverlap$Fraction, sep=":"), 
                          levels = c("Adult:Cytosplasm", "Prenatal:Cytosplasm", "Adult:Nucleus", "Prenatal:Nucleus"))
PercentIRsOverlap$variable = factor(PercentIRsOverlap$variable, levels=c("Constitutively\nSpliced",">0%",">5%",">10%",">20%",">30%",">40%",">50%",">75%",">90%"))
(PercentIRsOverlap[PercentIRsOverlap$variable=="Constitutively\nSpliced","value"])[order(PercentIRsOverlap[PercentIRsOverlap$variable=="Constitutively\nSpliced","value"])]
# range of 61.33-87.02% are constitutively spliced
(PercentIRsOverlap[PercentIRsOverlap$variable==">0%","value"]-PercentIRsOverlap[PercentIRsOverlap$variable==">5%","value"])[order(PercentIRsOverlap[PercentIRsOverlap$variable==">0%","value"]-PercentIRsOverlap[PercentIRsOverlap$variable==">5%","value"])]
# range of 10.76-32.48% have 0% < IR ratio < 5%


pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio_QC_intronsOnly.pdf",width =11,height = 5)
ggplot(PercentIRsOverlap, aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio\n(Introns Passing QC In All Groups)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.8, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off() 
pdf(file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/global_IR_comparisons/percent_introns_byIRRatio_greaterthan5perc_QC_introns.pdf",width =8,height = 5)
ggplot(PercentIRsOverlap[which(PercentIRsOverlap$variable != "Constitutively\nSpliced" & PercentIRsOverlap$variable != ">0%"),],
       aes(x=variable, y=value, fill=Group), color=Group) + 
  geom_boxplot() +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio\n(Introns Passing QC In All Groups)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.8, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()


# Is IR increased in the nucleus or in adults?

ttest = list(byFraction = t.test(x = IRratios[which(IRratios$Fraction=="Nucleus"),"IRratio"],y = IRratios[which(IRratios$Fraction=="Cytosplasm"),"IRratio"]),
             byFrac.Adult = t.test(x = IRratios[which(IRratios$Group=="Adult:Nucleus"),"IRratio"],y = IRratios[which(IRratios$Group=="Adult:Cytosplasm"),"IRratio"]),
             byFrac.Prenatal = t.test(x = IRratios[which(IRratios$Group=="Prenatal:Nucleus"),"IRratio"],y = IRratios[which(IRratios$Group=="Prenatal:Cytosplasm"),"IRratio"]),
             byFrac.50 = t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$variable == ">50%"),"value"],
                                y = PercentIRs[which(PercentIRs$Fraction=="Cytosplasm" & PercentIRs$variable == ">50%"),"value"]),
             byFrac.50.Adult = t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">50%"),"value"],
                                      y = PercentIRs[which(PercentIRs$Group=="Adult:Cytosplasm" & PercentIRs$variable==">50%"),"value"]),
             byFrac.50.Prenatal = t.test(x = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">50%"),"value"],
                                         y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosplasm" & PercentIRs$variable==">50%"),"value"]),
             byFrac.10 = t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$variable == ">10%"),"value"],
                                y = PercentIRs[which(PercentIRs$Fraction=="Cytosplasm" & PercentIRs$variable == ">10%"),"value"]),
             byFrac.10.Adult = t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">10%"),"value"],
                                      y = PercentIRs[which(PercentIRs$Group=="Adult:Cytosplasm" & PercentIRs$variable==">10%"),"value"]),
             byFrac.10.Prenatal = t.test(x = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">10%"),"value"],
                                         y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosplasm" & PercentIRs$variable==">10%"),"value"]),
             byAge = t.test(x = IRratios[which(IRratios$Age=="Adult"),"IRratio"],y = IRratios[which(IRratios$Age=="Prenatal"),"IRratio"]),
             byAge.Nucleus = t.test(x = IRratios[which(IRratios$Group=="Adult:Nucleus"),"IRratio"],
                                    y = IRratios[which(IRratios$Group=="Prenatal:Nucleus"),"IRratio"]),
             byAge.Cytosplasm = t.test(x = IRratios[which(IRratios$Group=="Adult:Cytosplasm"),"IRratio"],
                                    y = IRratios[which(IRratios$Group=="Prenatal:Cytosplasm"),"IRratio"]),
             byAge.50 = t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">50%"),"value"],
                               y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">50%"),"value"]),
             byAge.50.Nucleus = t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">50%"),"value"],
                                       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">50%"),"value"]),
             byAge.50.Cytosplasm = t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosplasm" & PercentIRs$variable==">50%"),"value"],
                                       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosplasm" & PercentIRs$variable==">50%"),"value"]),
             byAge.10 = t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$variable==">10%"),"value"],
                               y = PercentIRs[which(PercentIRs$Age=="Prenatal" & PercentIRs$variable==">10%"),"value"]),
             byAge.10.Nucleus = t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Nucleus" & PercentIRs$variable==">10%"),"value"],
                                       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Nucleus" & PercentIRs$variable==">10%"),"value"]),
             byAge.10.Cytosplasm = t.test(x = PercentIRs[which(PercentIRs$Group=="Adult:Cytosplasm" & PercentIRs$variable==">10%"),"value"],
                                       y = PercentIRs[which(PercentIRs$Group=="Prenatal:Cytosplasm" & PercentIRs$variable==">10%"),"value"]))
x = data.frame(Comp = names(ttest), Tstat = unlist(lapply(ttest, function(x) round(x$statistic,3))), pval = unlist(lapply(ttest, function(x) x$p.value)),
               confInt1 = unlist(lapply(ttest, function(x) round(x$conf.int,3)[1])), confInt2 = unlist(lapply(ttest, function(x) round(x$conf.int,3)[2])),
               estMeans1 = unlist(lapply(ttest, function(x) round(x$estimate,3)[1])), estMeans2 = unlist(lapply(ttest, function(x) round(x$estimate,3)[2])),
               row.names = NULL)
x$FDR = p.adjust(x$pval, method="fdr")
write.csv(x, quote=F, file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/global_IR_comparisons/Ttest_results_IRratio_byFraction_byAge.csv")
df = read.csv("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/global_IR_comparisons/Ttest_results_IRratio_byFraction_byAge.csv")
df

