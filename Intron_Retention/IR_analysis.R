library(ggplot2)
library(reshape2)
library(VennDiagram)
library(data.table)

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
IRfiltered = lapply(IRres, function(x) x[which(x$Warnings=="-"),])       
introns = data.frame(elementLengths(IRfiltered))
introns$SampleID = as.factor(rownames(introns))
introns$SampleID = gsub("_RiboZero", "\nRibozero", introns$SampleID)
introns$SampleID = gsub("_polyA", "\npolyA", introns$SampleID)
introns$Group = as.factor(c("Adult:Cytosol:PolyA", "Adult:Nucleus:PolyA","Adult:Cytosol:PolyA", "Adult:Nucleus:PolyA","Adult:Cytosol:PolyA", "Adult:Nucleus:PolyA",
                            "Fetal:Cytosol:PolyA",  "Fetal:Nucleus:PolyA","Fetal:Cytosol:PolyA",  "Fetal:Nucleus:PolyA","Fetal:Cytosol:PolyA",  "Fetal:Nucleus:PolyA",
                            "Adult:Cytosol:RiboZero","Adult:Nucleus:RiboZero","Adult:Cytosol:RiboZero","Adult:Nucleus:RiboZero","Adult:Cytosol:RiboZero","Adult:Nucleus:RiboZero",
                            "Fetal:Cytosol:RiboZero","Fetal:Nucleus:RiboZero","Fetal:Cytosol:RiboZero","Fetal:Nucleus:RiboZero","Fetal:Cytosol:RiboZero","Fetal:Nucleus:RiboZero"))
levels(introns$Group) = c("Adult:Cytosol:PolyA","Fetal:Cytosol:PolyA","Adult:Nucleus:PolyA","Fetal:Nucleus:PolyA",
                          "Adult:Cytosol:RiboZero","Fetal:Cytosol:RiboZero","Adult:Nucleus:RiboZero","Fetal:Nucleus:RiboZero")
string = lapply(IRfiltered, function(x) as.character(x$GeneIntronDetails))
string = lapply(string, function(x) strsplit(x, "/", fixed = TRUE))
x = lapply(string, function(x) unlist(x, recursive = FALSE))
y = lapply(x, function(x) grep("ENSG", x))
c = lapply(x, function(x) seq.int(3, length(x), 3))
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


ggplot(introns, aes(x=Group, y=elementLengths.IRfiltered.)) + 
  geom_boxplot() + geom_jitter() +
  ylab("Number") + 
  xlab("") +
  ggtitle("Introns Passing QC") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# What does the distribution of IR Ratios by individual intron look like?

IRratios = lapply(IRfiltered, function(x) x$IRratio) 
elementLengths(IRratios)
IRratios = Map(cbind, IRratios, Age=NA, Fraction=NA, Library=NA)
IRratios = lapply(IRratios, function(x) as.data.frame(x))
pol = IRratios[1:12]
rib = IRratios[13:24]
pol = lapply(pol, function(x) data.frame(IRratios = x[,1], Age = x[,2], Fraction = x[,3], Library = "polyA"))
rib = lapply(rib, function(x) data.frame(IRratios = x[,1], Age = x[,2], Fraction = x[,3], Library = "Ribozero"))
IRratios = c(pol, rib)
ad = IRratios[c(1:6, 13:18)]
fet = IRratios[c(7:12, 19:24)]
ad = lapply(ad, function(x) data.frame(IRratios = x[,1], Age = "Adult", Fraction = x[,3], Library = x$Library))
fet = lapply(fet, function(x) data.frame(IRratios = x[,1], Age = "Fetal", Fraction = x[,3], Library = x$Library))
IRratios = c(ad,fet)
cyt = IRratios[seq.int(1,24,2)]
nuc = IRratios[seq.int(2,24,2)]
cyt = lapply(cyt, function(x) data.frame(IRratios = x[,1], Age = x$Age, Fraction = "Cytosol", Library = x$Library))
nuc = lapply(nuc, function(x) data.frame(IRratios = x[,1], Age = x$Age, Fraction = "Nucleus", Library = x$Library))
IRratios = c(cyt, nuc)
IRratios = IRratios[c(1:15,17:24)]
IRratios = do.call(rbind, IRratios)
IRratios$Group = paste(IRratios$Age, IRratios$Fraction, IRratios$Library, sep=":")
IRratios$Group = factor(IRratios$Group) 
levels(IRratios$Group) = c("Adult:Cytosol:PolyA", "Fetal:Cytosol:PolyA", "Adult:Nucleus:PolyA", "Fetal:Nucleus:PolyA",
                           "Adult:Cytosol:RiboZero", "Fetal:Cytosol:RiboZero", "Adult:Nucleus:RiboZero", "Fetal:Nucleus:RiboZero")

ggplot(IRratios[which(IRratios$Library=="polyA"),], aes(x=IRratios)) + geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("IR Ratio") +
  ggtitle("IR Ratios By Group (PolyA)") +
  xlim(0,.06) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(IRratios[which(IRratios$Library=="Ribozero"),], aes(x=IRratios)) + 
  geom_density(aes(group=Group, colour=Group)) +
  ylab("") + 
  xlab("IR Ratio") +
  ggtitle("IR Ratios By Group (RiboZero)") + 
  xlim(0,.06) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
  
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
PercentIRs = PercentIRs[which(rownames(PercentIRs)!="Br1113N1_RiboZero"),]
PercentIRs = data.frame(melt(PercentIRs), SampleID = rep.int(rownames(PercentIRs), 10))
PercentIRs$Age = c(rep.int("Adult", 6), rep.int("Fetal", 6), rep.int("Adult", 5), rep.int("Fetal", 6))
PercentIRs$Library = rep.int(c(rep.int("PolyA", 12), rep.int("RiboZero", 11)), 10)
PercentIRs$Fraction = c(rep.int(c("Cytosol", "Nucleus"), 6), "Cytosol", rep.int(c("Cytosol", "Nucleus"), 5))
PercentIRs$Group = factor(paste(PercentIRs$Age, PercentIRs$Fraction, PercentIRs$Library, sep=":"), 
                          levels = c("Adult:Cytosol:PolyA", "Fetal:Cytosol:PolyA", "Adult:Nucleus:PolyA", "Fetal:Nucleus:PolyA",
                                     "Adult:Cytosol:RiboZero", "Fetal:Cytosol:RiboZero", "Adult:Nucleus:RiboZero", "Fetal:Nucleus:RiboZero"))
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
dodge <- position_dodge(width=0.9)
ggplot(PercentIRs, aes(x=variable, y=value, fill=Group), color=Group) + 
  stat_summary(position=position_dodge(),geom="bar") +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(PercentIRs[which(PercentIRs$Library=="PolyA"),],
       aes(x=variable, y=value, fill=Group), color=Group) + 
  stat_summary(position=position_dodge(),geom="bar") +
  ylab("Percent") + 
  xlab("Intron Retention") +
  ggtitle("Percent Introns By IR Ratio (PolyA)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

ggplot(PercentIRs[which(PercentIRs$Library=="PolyA" & PercentIRs$variable != "Constitutively\nSpliced" & PercentIRs$variable != ">0%"),],
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

ggplot(PercentIRs[which(PercentIRs$Library=="RiboZero" & PercentIRs$variable != "Constitutively\nSpliced" & PercentIRs$variable != ">0%"),],
       aes(x=variable, y=value, fill=Group), color=Group) + 
  stat_summary(position=position_dodge(),geom="bar") +
  ylab("Percent") + 
  xlab("") +
  ggtitle("Percent Introns By IR Ratio (RiboZero)") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.position = c(.85, 0.6)) +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))


# Is IR increased in the nucleus?
# Main effect: All filtered introns
elementLengths(IRratios)
t.test(x = IRratios[which(IRratios$Fraction=="Nucleus" & IRratios$Library=="polyA"),which(colnames(IRratios)=="IRratios")],
       y = IRratios[which(IRratios$Fraction=="Cytosol" & IRratios$Library=="polyA"),which(colnames(IRratios)=="IRratios")],
       alternative = "greater")
#data:  IRratios[which(IRratios$Fraction == "Nucleus" & IRratios$Library ==  and IRratios[which(IRratios$Fraction == "Cytosol" & IRratios$Library ==     "polyA"), which(colnames(IRratios) == "IRratios")] and     "polyA"), which(colnames(IRratios) == "IRratios")]
#t = 41.942, df = 555970, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.005518557         Inf
#sample estimates:
#  mean of x  mean of y 
#0.01578728 0.01004347

# in adult polyA
t.test(x = IRratios[which(IRratios$Fraction=="Nucleus" & IRratios$Library=="polyA" & IRratios$Age=="Adult"),which(colnames(IRratios)=="IRratios")],
       y = IRratios[which(IRratios$Fraction=="Cytosol" & IRratios$Library=="polyA" & IRratios$Age=="Adult"),which(colnames(IRratios)=="IRratios")],
       alternative = "greater")
#data:  IRratios[which(IRratios$Fraction == "Nucleus" & IRratios$Library ==  and IRratios[which(IRratios$Fraction == "Cytosol" & IRratios$Library ==     "polyA" & IRratios$Age == "Adult"), which(colnames(IRratios) ==  and     "polyA" & IRratios$Age == "Adult"), which(colnames(IRratios) ==     "IRratios")] and     "IRratios")]
#t = 51.234, df = 187730, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.01060462        Inf
#sample estimates:
#  mean of x   mean of y 
#0.018589335 0.007632962 

# in fetal polyA
t.test(x = IRratios[which(IRratios$Fraction=="Nucleus" & IRratios$Library=="polyA" & IRratios$Age=="Fetal"),which(colnames(IRratios)=="IRratios")],
       y = IRratios[which(IRratios$Fraction=="Cytosol" & IRratios$Library=="polyA" & IRratios$Age=="Fetal"),which(colnames(IRratios)=="IRratios")],
       alternative = "greater")
#data:  IRratios[which(IRratios$Fraction == "Nucleus" & IRratios$Library ==  and IRratios[which(IRratios$Fraction == "Cytosol" & IRratios$Library ==     "polyA" & IRratios$Age == "Fetal"), which(colnames(IRratios) ==  and     "polyA" & IRratios$Age == "Fetal"), which(colnames(IRratios) ==     "IRratios")] and     "IRratios")]
#t = 14.141, df = 350770, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.002187468         Inf
#sample estimates:
#  mean of x  mean of y 
#0.01389595 0.01142055


# Main effect: introns >=75% retained
head(PercentIRs)
# in all polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$variable == ">75%"), which(colnames(PercentIRs) ==  and     "PolyA" & PercentIRs$variable == ">75%"), which(colnames(PercentIRs) ==     "value")] and     "value")]
#t = 0.11417, df = 7.6003, p-value = 0.456
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.056107       Inf
#sample estimates:
#  mean of x mean of y 
#0.1046333 0.1009900
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$variable == ">50%"), which(colnames(PercentIRs) ==  and     "PolyA" & PercentIRs$variable == ">50%"), which(colnames(PercentIRs) ==     "value")] and     "value")]
#t = 1.7437, df = 9.9664, p-value = 0.05596
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.005334226          Inf
#sample estimates:
#  mean of x mean of y 
#0.3478683 0.2138433

# in adult polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Adult" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Adult" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Adult" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Adult" & PercentIRs$variable ==     ">75%"), which(colnames(PercentIRs) == "value")] and     ">75%"), which(colnames(PercentIRs) == "value")]
#t = 2.644, df = 3.1771, p-value = 0.03644
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.008993908         Inf
#sample estimates:
#  mean of x  mean of y 
#0.12810667 0.05880333

# in fetal polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==     ">75%"), which(colnames(PercentIRs) == "value")] and     ">75%"), which(colnames(PercentIRs) == "value")]
#t = -1.3617, df = 2.0908, p-value = 0.8493
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.1911508        Inf
#sample estimates:
#  mean of x mean of y 
#0.0811600 0.1431767
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==     ">50%"), which(colnames(PercentIRs) == "value")] and     ">50%"), which(colnames(PercentIRs) == "value")]
#t = -0.40062, df = 2.7842, p-value = 0.6413
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.2564962        Inf
#sample estimates:
#  mean of x mean of y 
#0.2604233 0.2967300 


## Is IR increased in adults?
# Main effect: All filtered introns
elementLengths(IRratios)
t.test(x = IRratios[which(IRratios$Age=="Adult" & IRratios$Library=="polyA"),which(colnames(IRratios)=="IRratios")],
       y = IRratios[which(IRratios$Age=="Fetal" & IRratios$Library=="polyA"),which(colnames(IRratios)=="IRratios")],
       alternative = "greater")
#data:  IRratios[which(IRratios$Age == "Adult" & IRratios$Library ==  and IRratios[which(IRratios$Age == "Fetal" & IRratios$Library == "polyA"), which(colnames(IRratios) == "IRratios")] and     "polyA"), which(colnames(IRratios) == "IRratios")]
#t = 5.6143, df = 464860, p-value = 9.874e-09
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.0005593861          Inf
#sample estimates:
#  mean of x  mean of y 
#0.01341493 0.01262374

# in Nuclear polyA
t.test(x = IRratios[which(IRratios$Age=="Adult" & IRratios$Library=="polyA" & IRratios$Fraction=="Nucleus"),which(colnames(IRratios)=="IRratios")],
       y = IRratios[which(IRratios$Age=="Fetal" & IRratios$Library=="polyA" & IRratios$Fraction=="Nucleus"),which(colnames(IRratios)=="IRratios")],
       alternative = "greater")
#data:  IRratios[which(IRratios$Age == "Adult" & IRratios$Library ==  and IRratios[which(IRratios$Age == "Fetal" & IRratios$Library ==     "polyA" & IRratios$Fraction == "Nucleus"), which(colnames(IRratios) ==  and     "polyA" & IRratios$Fraction == "Nucleus"), which(colnames(IRratios) ==     "IRratios")] and     "IRratios")]
#t = 21.301, df = 214160, p-value < 2.2e-16
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  0.004330958         Inf
#sample estimates:
#  mean of x  mean of y 
#0.01858933 0.01389595 

# in Cytosol polyA
t.test(x = IRratios[which(IRratios$Age=="Adult" & IRratios$Library=="polyA" & IRratios$Fraction=="Cytosol"),which(colnames(IRratios)=="IRratios")],
       y = IRratios[which(IRratios$Age=="Fetal" & IRratios$Library=="polyA" & IRratios$Fraction=="Cytosol"),which(colnames(IRratios)=="IRratios")],
       alternative = "greater")
#data:  IRratios[which(IRratios$Age == "Adult" & IRratios$Library ==  and IRratios[which(IRratios$Age == "Fetal" & IRratios$Library ==     "polyA" & IRratios$Fraction == "Cytosol"), which(colnames(IRratios) ==  and     "polyA" & IRratios$Fraction == "Cytosol"), which(colnames(IRratios) ==     "IRratios")] and     "IRratios")]
#t = -22.705, df = 275350, p-value = 1
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.004061975          Inf
#sample estimates:
#  mean of x   mean of y 
#0.007632962 0.011420549


# Main effect: introns >=75% retained
head(PercentIRs)
# in all polyA
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Age=="Fetal" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Age == "Adult" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Age == "Fetal" & PercentIRs$Library ==     "PolyA" & PercentIRs$variable == ">75%"), which(colnames(PercentIRs) ==  and     "PolyA" & PercentIRs$variable == ">75%"), which(colnames(PercentIRs) ==     "value")] and     "value")]
#t = -0.59634, df = 9.4844, p-value = 0.7175
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.07590452         Inf
#sample estimates:
#  mean of x mean of y 
#0.0934550 0.1121683
t.test(x = PercentIRs[which(PercentIRs$Age=="Adult" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Age=="Fetal" & PercentIRs$Library=="PolyA" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Age == "Adult" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Age == "Fetal" & PercentIRs$Library ==     "PolyA" & PercentIRs$variable == ">50%"), which(colnames(PercentIRs) ==  and     "PolyA" & PercentIRs$variable == ">50%"), which(colnames(PercentIRs) ==     "value")] and     "value")]
#t = 0.05194, df = 7.6364, p-value = 0.48
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.1596521        Inf
#sample estimates:
#  mean of x mean of y 
#0.2831350 0.2785767

# in Nuclear polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Adult" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Adult" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==     ">75%"), which(colnames(PercentIRs) == "value")] and     ">75%"), which(colnames(PercentIRs) == "value")]
#t = 1.9761, df = 2.3527, p-value = 0.08367
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.01582913         Inf
#sample estimates:
#  mean of x mean of y 
#0.1281067 0.0811600
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Adult" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Nucleus" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Nucleus" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Adult" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==     ">50%"), which(colnames(PercentIRs) == "value")] and     ">50%"), which(colnames(PercentIRs) == "value")]
#t = 1.9576, df = 2.8095, p-value = 0.07569
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.04129456         Inf
#sample estimates:
#  mean of x mean of y 
#0.4353133 0.2604233

# in cytosol polyA
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Adult" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">75%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater")
#data:  PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Adult" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==     ">75%"), which(colnames(PercentIRs) == "value")] and     ">75%"), which(colnames(PercentIRs) == "value")]
#t = -1.8001, df = 2.3304, p-value = 0.9023
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.2088809        Inf
#sample estimates:
#  mean of x  mean of y 
#0.05880333 0.14317667
t.test(x = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Adult" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       y = PercentIRs[which(PercentIRs$Fraction=="Cytosol" & PercentIRs$Library=="PolyA" & 
                              PercentIRs$Age=="Fetal" & PercentIRs$variable == ">50%"),which(colnames(PercentIRs)=="value")],
       alternative = "greater") 
#data:  PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==  and PercentIRs[which(PercentIRs$Fraction == "Cytosol" & PercentIRs$Library ==     "PolyA" & PercentIRs$Age == "Adult" & PercentIRs$variable ==  and     "PolyA" & PercentIRs$Age == "Fetal" & PercentIRs$variable ==     ">50%"), which(colnames(PercentIRs) == "value")] and     ">50%"), which(colnames(PercentIRs) == "value")]
#t = -1.9755, df = 2.1297, p-value = 0.9105
#alternative hypothesis: true difference in means is greater than 0
#95 percent confidence interval:
#  -0.4009251        Inf
#sample estimates:
#  mean of x mean of y 
#0.1309567 0.2967300

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
elementLengths(IRcomp)
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
IRclean = lapply(polya, function(x) x[which(x$A.warnings=="-" & x$B.warnings=="-" & x$comments=="clean"),])

# Explore the results
elementLengths(IRclean)
      #Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
      #761               481               517               662                53                82
elementLengths(lapply(IRclean, function(x) x[which(x$Sign=="MoreIRInNuc.Fetal"),]))
      #Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
      #752               445               461               343                53                67
elementLengths(lapply(IRclean, function(x) x[which(x$p.diff<=0.05),]))
      #Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
      #161               166                97               156                34                36 
elementLengths(lapply(IRclean, function(x) x[which(x$p.diff<=0.05 & x$Sign=="MoreIRInNuc.Fetal"),]))
      #Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
      #159               164                77                64                34                26
elementLengths(lapply(IRclean, function(x) x[which(x$A.IRratio>=0.5 | x$B.IRratio>=0.5),]))
      #Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
      #4                 3                 4                10                 2                 1
elementLengths(lapply(lapply(IRclean, function(x) x[which(x$A.IRratio>=0.5 | x$B.IRratio>=0.5),]), function(x) 
  x[which(x$Sign=="MoreIRInNuc.Fetal"),]))
      #Adult_PolyA_Zone  Fetal_PolyA_Zone Cytosol_PolyA_Age Nuclear_PolyA_Age        PolyA_Zone         PolyA_Age 
      #4                 2                 4                 3                 2                 1
fisher.test(data.frame(c(159,2),c(593,7))) 
#data:  polyA adult
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.176366 9.346762
#sample estimates:
#  odds ratio 
#0.9385205 
fisher.test(data.frame(c(164,2),c(281,36)))
#data:  polyA fetal
#p-value = 2.19e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.636957 90.956349
#sample estimates:
#  odds ratio 
#10.47494
fisher.test(data.frame(c(77,20),c(384,36)))
#data:  data.frame(c(77, 20), c(384, 36))
#p-value = 0.001621
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1918414 0.6970453
#sample estimates:
#  odds ratio 
#0.3618495 
fisher.test(data.frame(c(64,92),c(279,227)))
#data:  data.frame(c(64, 92), c(279, 227))
#p-value = 0.002444
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3862584 0.8272246
#sample estimates:
#  odds ratio 
#0.5665001


# Are IR genes significantly regulated by fraction?
FracList = list(Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv"),
                Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv"),
                Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv"),
                Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv"))
SigFracList = lapply(FracList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigFracList)
Sign = lapply(SigFracList, function(x) ifelse(x$log2FoldChange > 0,"UpNuc", "DownNuc"))
SigFracList = Map(cbind, SigFracList, Sign = Sign)
DirList = lapply(SigFracList, function(x) split(x, x$Sign))
SigList = list(Apres.Up = DirList[["Apres"]][["UpNuc"]], Apres.Down = DirList[["Apres"]][["DownNuc"]],
               Fpres.Up = DirList[["Fpres"]][["UpNuc"]], Fpres.Down = DirList[["Fpres"]][["DownNuc"]],
               Arres.Up = DirList[["Arres"]][["UpNuc"]], Arres.Down = DirList[["Arres"]][["DownNuc"]],
               Frres.Up = DirList[["Frres"]][["UpNuc"]], Frres.Down = DirList[["Frres"]][["DownNuc"]]) 
SigList = lapply(SigList, function(x) as.character(x$X))
names(SigList) = c("Adult\nPolyA\nNucleus", "Adult\nPolyA\nCytosol", "Fetal\nPolyA\nNucleus", "Fetal\nPolyA\nCytosol", 
                   "Adult\nRibozero\nNucleus", "Adult\nRibozero\nCytosol", "Fetal\nRibozero\nNucleus", "Fetal\nRibozero\nCytosol")
elementLengths(SigList)
allgenes = FracList[[1]]
allgenes = as.character(allgenes$X)
string = lapply(IRfiltered, function(x) as.character(x$GeneIntronDetails))
string = lapply(string, function(x) strsplit(x, "/", fixed = TRUE))
x = lapply(string, function(x) unlist(x, recursive = FALSE))
x = x[1:12]
y = lapply(x, function(x) grep("ENSG", x))
c = lapply(x, function(x) seq.int(3, length(x), 3))
comments = genes = list()
for (i in 1:length(x)){
  tmp = x[[i]]
  genes[[i]] = tmp[y[[i]]]
  comments[[i]] = tmp[c[[i]]]}
names(genes) = names(comments) = names(x)
polyaFilt = Map(cbind, IRfiltered[1:12], ensID = genes, comments = comments)

# Get the IR ratio for Fraction genes
AP = lapply(SigFracList[1], function(x) as.character(x$X))
FP = lapply(SigFracList[2], function(x) as.character(x$X))
APC = SigList[["Adult\nPolyA\nCytosol"]]
APN = SigList[["Adult\nPolyA\nNucleus"]]
FPC = SigList[["Fetal\nPolyA\nCytosol"]]
FPN = SigList[["Fetal\nPolyA\nNucleus"]]

AP = lapply(polyaFilt, function(x) x[which(x$ensID%in%AP),])
AP = lapply(AP, function(x) data.table(x, key="ensID"))
AP = lapply(AP,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
AP = do.call(rbind,AP)
FP = lapply(polyaFilt, function(x) x[which(x$ensID%in%FP),])
FP = lapply(FP, function(x) data.table(x, key="ensID"))
FP = lapply(FP,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
FP = do.call(rbind,FP)
APC = lapply(polyaFilt, function(x) x[which(x$ensID%in%APC),])
APC = lapply(APC, function(x) data.table(x, key="ensID"))
APC = lapply(APC,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
APC = do.call(rbind,APC)
APN = lapply(polyaFilt, function(x) x[which(x$ensID%in%APN),])
APN = lapply(APN, function(x) data.table(x, key="ensID"))
APN = lapply(APN,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
APN = do.call(rbind,APN)
FPC = lapply(polyaFilt, function(x) x[which(x$ensID%in%FPC),])
FPC = lapply(FPC, function(x) data.table(x, key="ensID"))
FPC = lapply(FPC,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
FPC = do.call(rbind,FPC)
FPN = lapply(polyaFilt, function(x) x[which(x$ensID%in%FPN),])
FPN = lapply(FPN, function(x) data.table(x, key="ensID"))
FPN = lapply(FPN,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
FPN = do.call(rbind,FPN)

t.test(APC$IRratio, APN$IRratio, alternative = "less")
#For LFC>1: data:  APC$IRratio and APN$IRratio
#t = -38.132, df = 5636.5, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.111235
#sample estimates:
#  mean of x  mean of y 
#0.02157524 0.13782572

#NO LFC Cutoff: data:  APC$IRratio and APN$IRratio
#t = -50.035, df = 17027, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.0717663
#sample estimates:
#  mean of x  mean of y 
#0.02447499 0.09868088
t.test(FPC$IRratio, FPN$IRratio, alternative = "less")
#data:  FPC$IRratio and FPN$IRratio
#t = -33.876, df = 2177, p-value < 2.2e-16
#alternative hypothesis: true difference in means is less than 0
#95 percent confidence interval:
#  -Inf -0.1452936
#sample estimates:
#  mean of x   mean of y 
#0.004986751 0.157698449 

AP = lapply(SigFracList[1], function(x) as.character(x$X))
FP = lapply(SigFracList[2], function(x) as.character(x$X))
APC = SigList[["Adult\nPolyA\nCytosol"]]
APN = SigList[["Adult\nPolyA\nNucleus"]]
FPC = SigList[["Fetal\nPolyA\nCytosol"]]
FPN = SigList[["Fetal\nPolyA\nNucleus"]]
elementLengths(IRcomp)
head(IRcomp[[1]])
IR.APC = IRclean[["Adult_PolyA_Zone"]]
IR.APN = IRclean[["Adult_PolyA_Zone"]]
IR.FPC = IRclean[["Fetal_PolyA_Zone"]]
IR.FPN = IRclean[["Fetal_PolyA_Zone"]]
IR.APC.sig = IR.APC[which(IR.APC$p.diff<=0.05 & IR.APC$Sign=="MoreIRInCyt.Adult"),]
IR.APN.sig = IR.APN[which(IR.APN$p.diff<=0.05 & IR.APN$Sign=="MoreIRInNuc.Fetal"),]
IR.FPC.sig = IR.FPC[which(IR.FPC$p.diff<=0.05 & IR.FPC$Sign=="MoreIRInCyt.Adult"),]
IR.FPN.sig = IR.FPN[which(IR.FPN$p.diff<=0.05 & IR.FPN$Sign=="MoreIRInNuc.Fetal"),]
IR.APC.perc = IR.APC[which((IR.APC$A.IRratio>=0.5 | IR.APC$B.IRratio>=0.5) & IR.APC$Sign=="MoreIRInCyt.Adult"),]
IR.APN.perc = IR.APN[which((IR.APN$A.IRratio>=0.5 | IR.APN$B.IRratio>=0.5) & IR.APN$Sign=="MoreIRInNuc.Fetal"),]
IR.FPC.perc = IR.FPC[which((IR.FPC$A.IRratio>=0.5 | IR.FPC$B.IRratio>=0.5) & IR.FPC$Sign=="MoreIRInCyt.Adult"),]
IR.FPN.perc = IR.FPN[which((IR.FPN$A.IRratio>=0.5 | IR.FPN$B.IRratio>=0.5) & IR.FPN$Sign=="MoreIRInNuc.Fetal"),]
IR.APC.sig = as.character(IR.APC.sig$ensID)
IR.APN.sig = as.character(IR.APN.sig$ensID)
IR.FPC.sig = as.character(IR.FPC.sig$ensID)
IR.FPN.sig = as.character(IR.FPN.sig$ensID)
IR.APC.perc = as.character(IR.APC.perc$ensID)
IR.APN.perc = as.character(IR.APN.perc$ensID)
IR.FPC.perc = as.character(IR.FPC.perc$ensID)
IR.FPN.perc = as.character(IR.FPN.perc$ensID)

AP = FracList[["Apres"]]
APC = AP[which(AP$padj<=0.05),]
APN = AP[which(AP$padj>0.05),]
AdultfracIR.sig = list(IR.APC = IR.APC.sig, IR.APN = IR.APN.sig, YesFrac = as.character(APC$X), NoFrac = as.character(APN$X))
FP = FracList[["Apres"]]
FPC = FP[which(FP$padj<=0.05),]
FPN = FP[which(FP$padj>0.05),]
FetalfracIR.sig = list(IR.FPC = IR.FPC.sig, IR.FPN = IR.FPN.sig, YesFrac = as.character(FPC$X), NoFrac = as.character(FPN$X))
venn.Age <- venn.diagram(AdultfracIR.sig, "/Users/amanda/Desktop/AdultfracIR.sig.jpeg", main="AdultfracIR.sig",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)
elementLengths(FetalfracIR.sig)
venn.Age <- venn.diagram(FetalfracIR.sig, "/Users/amanda/Desktop/FetalfracIR.sig.jpeg", main="FetalfracIR.sig",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)

adult = data.frame(c(69,4348), c(70,12568))
adult = as.matrix(adult)
fetal = data.frame(c(69,4348), c(71,12567))
fetal = as.matrix(fetal)

fisher.test(adult)
#data:  adult
#p-value = 2.074e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.00934 4.03866
#sample estimates:
#  odds ratio 
#2.849022
fisher.test(fetal)
#data:  fetal
#p-value = 2.687e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.983151 3.976034
#sample estimates:
#  odds ratio 
#2.808675

AdultfracIR.perc = list(IR.APC = IR.APC.perc, IR.APN = IR.APN.perc, YesFrac = as.character(FPC$X), NoFrac = as.character(FPN$X))
FetalfracIR.perc = list(IR.FPC = IR.FPC.perc, IR.FPN = IR.FPN.perc, YesFrac = as.character(FPC$X), NoFrac = as.character(FPN$X))
elementLengths(AdultfracIR.perc)
venn.Age <- venn.diagram(AdultfracIR.perc, "/Users/amanda/Desktop/AdultfracIR.perc.jpeg", main="AdultfracIR.perc",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)
elementLengths(FetalfracIR.perc)
venn.Age <- venn.diagram(FetalfracIR.perc, "/Users/amanda/Desktop/FetalfracIR.perc.jpeg", main="FetalfracIR.perc",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)

# Get the IR ratio for developmental genes
AgeList = list(Cpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Cpres.csv"),
               Npres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Npres.csv"))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementLengths(SigAgeList)
Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpFetal", "UpAdult"))
SigAgeList = Map(cbind, SigAgeList, Sign = Sign)
DirList = lapply(SigAgeList, function(x) split(x, x$Sign))
SigList = list(Apres.Up = DirList[["Cpres"]][["UpFetal"]], Apres.Down = DirList[["Cpres"]][["UpAdult"]],
               Fpres.Up = DirList[["Npres"]][["UpFetal"]], Fpres.Down = DirList[["Npres"]][["UpAdult"]]) 
SigList = lapply(SigList, function(x) as.character(x$X))
names(SigList) = c("Cytosol\nPolyA\nFetal", "Cytosol\nPolyA\nAdult", "Nucleus\nPolyA\nFetal", "Nucleus\nPolyA\nAdult")
elementLengths(SigList)
allgenes = AgeList[[1]]
allgenes = as.character(allgenes$X)

CP = lapply(SigAgeList[1], function(x) as.character(x$X))
NP = lapply(SigAgeList[2], function(x) as.character(x$X))
CPA = SigList[["Cytosol\nPolyA\nAdult"]]
CPF = SigList[["Cytosol\nPolyA\nFetal"]]
NPA = SigList[["Nucleus\nPolyA\nAdult"]]
NPF = SigList[["Nucleus\nPolyA\nFetal"]]

CP = lapply(polyaFilt, function(x) x[which(x$ensID%in%CP),])
CP = lapply(CP, function(x) data.table(x, key="ensID"))
CP = lapply(CP,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
CP = do.call(rbind,CP)
NP = lapply(polyaFilt, function(x) x[which(x$ensID%in%NP),])
NP = lapply(NP, function(x) data.table(x, key="ensID"))
NP = lapply(NP,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
NP = do.call(rbind,NP)
CPA = lapply(polyaFilt, function(x) x[which(x$ensID%in%CPA),])
CPA = lapply(CPA, function(x) data.table(x, key="ensID"))
CPA = lapply(CPA,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
CPA = do.call(rbind,CPA)
CPF = lapply(polyaFilt, function(x) x[which(x$ensID%in%CPF),])
CPF = lapply(CPF, function(x) data.table(x, key="ensID"))
CPF = lapply(CPF,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
CPF = do.call(rbind,CPF)
NPA = lapply(polyaFilt, function(x) x[which(x$ensID%in%NPA),])
NPA = lapply(NPA, function(x) data.table(x, key="ensID"))
NPA = lapply(NPA,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
NPA = do.call(rbind,NPA)
NPF = lapply(polyaFilt, function(x) x[which(x$ensID%in%NPF),])
NPF = lapply(NPF, function(x) data.table(x, key="ensID"))
NPF = lapply(NPF,function(x) data.frame(x[, list(IRratio=max(IRratio)), by="ensID"]))
NPF = do.call(rbind,NPF)

t.test(CPA$IRratio, CPF$IRratio)
#data:  CPA$IRratio and CPF$IRratio
#t = -7.7772, df = 39783, p-value = 7.593e-15
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.010248666 -0.006122721
#sample estimates:
#  mean of x mean of y 
#0.0398155 0.0480012 
t.test(NPA$IRratio, NPF$IRratio)
#data:  NPA$IRratio and NPF$IRratio
#t = 6.7818, df = 33380, p-value = 1.206e-11
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.005641063 0.010227192
#sample estimates:
#  mean of x  mean of y 
#0.04770549 0.03977136

IR.CPA = IRclean[["Cytosol_PolyA_Age"]]
IR.CPF = IRclean[["Cytosol_PolyA_Age"]]
IR.NPA = IRclean[["Nuclear_PolyA_Age"]]
IR.NPF = IRclean[["Nuclear_PolyA_Age"]]

IR.CPA.sig = IR.CPA[which(IR.CPA$p.diff<=0.05 & IR.CPA$Sign=="MoreIRInCyt.Adult"),]
IR.CPF.sig = IR.CPF[which(IR.CPF$p.diff<=0.05 & IR.CPF$Sign=="MoreIRInNuc.Fetal"),]
IR.NPA.sig = IR.NPA[which(IR.NPA$p.diff<=0.05 & IR.NPA$Sign=="MoreIRInCyt.Adult"),]
IR.NPF.sig = IR.NPF[which(IR.NPF$p.diff<=0.05 & IR.NPF$Sign=="MoreIRInNuc.Fetal"),]
IR.CPA.perc = IR.CPA[which((IR.CPA$A.IRratio>=0.5 | IR.CPA$B.IRratio>=0.5) & IR.CPA$Sign=="MoreIRInCyt.Adult" & IR.CPA$p.diff<=0.05),]
IR.CPF.perc = IR.CPF[which((IR.CPF$A.IRratio>=0.5 | IR.CPF$B.IRratio>=0.5) & IR.CPF$Sign=="MoreIRInNuc.Fetal" & IR.CPF$p.diff<=0.05),]
IR.NPA.perc = IR.NPA[which((IR.NPA$A.IRratio>=0.5 | IR.NPA$B.IRratio>=0.5) & IR.NPA$Sign=="MoreIRInCyt.Adult" & IR.NPA$p.diff<=0.05),]
IR.NPF.perc = IR.NPF[which((IR.NPF$A.IRratio>=0.5 | IR.NPF$B.IRratio>=0.5) & IR.NPF$Sign=="MoreIRInNuc.Fetal" & IR.NPF$p.diff<=0.05),]
IR.CPA.sig = as.character(IR.CPA.sig$ensID)
IR.CPF.sig = as.character(IR.CPF.sig$ensID)
IR.NPA.sig = as.character(IR.NPA.sig$ensID)
IR.NPF.sig = as.character(IR.NPF.sig$ensID)
IR.CPA.perc = as.character(IR.CPA.perc$ensID)
IR.CPF.perc = as.character(IR.CPF.perc$ensID)
IR.NPA.perc = as.character(IR.NPA.perc$ensID)
IR.NPF.perc = as.character(IR.NPF.perc$ensID)

CP = AgeList[["Cpres"]]
CPA = AP[which(CP$padj<=0.05),]
CPF = AP[which(CP$padj>0.05),]
CytAgeIR.sig = list(IR.CPA = IR.CPA.sig, IR.CPF = IR.CPF.sig, YesAge = as.character(CPA$X), NoAge = as.character(CPF$X))
NP = AgeList[["Npres"]]
NPA = CP[which(CP$padj<=0.05),]
NPF = CP[which(CP$padj>0.05),]
NucAgeIR.sig = list(IR.NPA = IR.NPA.sig, IR.NPF = IR.NPF.sig, YesAge = as.character(NPA$X), NoAge = as.character(NPF$X))
venn.Age <- venn.diagram(CytAgeIR.sig, "/Users/amanda/Desktop/CytAgeIR.sig.jpeg", main="CytAgeIR.sig",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)
elementLengths(NucAgeIR.sig)
venn.Age <- venn.diagram(NucAgeIR.sig, "/Users/amanda/Desktop/NucAgeIR.sig.jpeg", main="NucAgeIR.sig",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)

cyt = data.frame(c(14,60), c(6,11))
cyt = as.matrix(cyt)
nuc = data.frame(c(62,42), c(25,17))
nuc = as.matrix(nuc)

fisher.test(cyt)
#data:  cyt
#p-value = 0.1922
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1198702 1.6747812
#sample estimates:
#  odds ratio 
#0.4324056 
fisher.test(nuc)
#data:  nuc
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4495005 2.2078870
#sample estimates:
#  odds ratio 
#1.003784

CytAgeIR.perc = list(IR.CPA = IR.CPA.perc, IR.CPF = IR.CPF.perc, YesAge = as.character(NPA$X), NoAge = as.character(NPF$X))
NucAgeIR.perc = list(IR.NPA = IR.NPA.perc, IR.NPF = IR.NPF.perc, YesAge = as.character(NPA$X), NoAge = as.character(NPF$X))
elementLengths(CytAgeIR.perc)
venn.Age <- venn.diagram(CytAgeIR.perc, "/Users/amanda/Desktop/CytAgeIR.perc.jpeg", main="CytAgeIR.perc",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)
elementLengths(NucAgeIR.perc)
venn.Age <- venn.diagram(NucAgeIR.perc, "/Users/amanda/Desktop/NucAgeIR.perc.jpeg", main="NucAgeIR.perc",
                         col = "transparent",
                         fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                         alpha = 0.50,
                         label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                       "white", "white", "white", "white", "palevioletred4", "white",
                                       "white", "white", "white", "darkblue", "white"),
                         fontfamily = "Arial",
                         fontface = "bold",
                         cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                         cat.fontfamily = "Arial", margin=0.2)

# Are the retained introns the last introns?













####### From before

# Do the gene sets overlap?
geneIDs = lapply(lapply(IRclean, function(x) x[which(x$p.diff<=0.05),]), function(x) unique(x$GeneID))
IRZone = list("Adult:PolyA"=geneIDs[["AdultPolyA"]], "Fetal:PolyA"=geneIDs[["FetalPolyA"]], 
                "Adult:RiboZero"=geneIDs[["AdultRibo"]], "Fetal:RiboZero"=geneIDs[["FetalRibo"]])

venn.IRZone <- venn.diagram(IRZone, "/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/Figure_PDFs/Venn.IRbyZone.jpeg", 
                            main="Genes with Significantly Differentially Retained Introns by Fraction",
                            col = "transparent",
                            fill = c("lightpink2","cornflowerblue", "olivedrab2", "darkorchid1"),
                            alpha = 0.50,
                            label.col = c("olivedrab4", "white", "darkorchid4", "white",
                                          "white", "white", "white", "white", "palevioletred4", "white",
                                          "white", "white", "white", "darkblue", "white"),
                            fontfamily = "Arial",
                            fontface = "bold",
                            cat.col = c("palevioletred4", "darkblue", "olivedrab4", "darkorchid4"),
                            cat.fontfamily = "Arial", margin=0.2)

# How many are differentially expressed by fraction?
Apres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Apres.csv")
Fpres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Fpres.csv")
Arres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Arres.csv")
Frres = read.csv("/Users/amanda/Dropbox/NucVsCytosol/Manuscript_Materials/RDAs/Frres.csv")

DE.IR = list(AP.IR = Apres[which(Apres$X %in% geneIDs[["AdultPolyA"]]),], FP.IR = Fpres[which(Fpres$X %in% geneIDs[["FetalPolyA"]]),], 
                 AR.IR = Arres[which(Arres$X %in% geneIDs[["AdultRibo"]]),], FR.IR = Frres[which(Frres$X %in% geneIDs[["FetalRibo"]]),])
    elementLengths(DE.IR)
      #AP.IR FP.IR  AR.IR  FR.IR 
      #144   149  1015   881
    elementLengths(lapply(DE.IR, function(x) x[which(x$padj<=0.05),])) 
      #AP.IR FP.IR AR.IR FR.IR 
      #71    19   543   105 
    elementLengths(lapply(FracList, function(x) x[which(x$padj<=0.05),]))
      #Apres Fpres Arres Frres 
      #4417   297  5605  1515
fisher.test(data.frame(c(71,4346), c(73,39951)))
    #data:  data.frame(c(71, 4346), c(73, 39951))
    #p-value < 2.2e-16
    #alternative hypothesis: true odds ratio is not equal to 1
    #95 percent confidence interval:
    #  6.346241 12.585631
    #sample estimates:
    #  odds ratio 
    #8.939848
fisher.test(data.frame(c(19,278),c(130,44014)))
    #data:  data.frame(c(19, 278), c(130, 44014))
    #p-value < 2.2e-16
    #alternative hypothesis: true odds ratio is not equal to 1
    #95 percent confidence interval:
    #  13.29770 38.24939
    #sample estimates:
    #  odds ratio 
    #23.12528

# Are the genes also differentially polyadenylated?
DE.dapars = list(AP.Dapars = Apres[which(Apres$X %in% DaparsGenes[["AdultPolyA"]]),], FP.Dapars = Fpres[which(Fpres$X %in% DaparsGenes[["FetalPolyA"]]),], 
                 AR.Dapars = Arres[which(Arres$X %in% DaparsGenes[["AdultRiboZero"]]),], FR.Dapars = Frres[which(Frres$X %in% DaparsGenes[["FetalRiboZero"]]),])
    elementLengths(DE.dapars)
      #AP.Dapars FP.Dapars AR.Dapars FR.Dapars
      #422   609   132   162 
    elementLengths(lapply(DE.dapars, function(x) x[which(x$padj<=0.05),])) 
      #AP.Dapars FP.Dapars AR.Dapars FR.Dapars 
      #181    23   103    38 
    elementLengths(lapply(DE.dapars, function(x) x[which(x$padj<=0.05 & x$log2FoldChange >= 1),])) #greater exp in nuc
      #AP.Dapars FP.Dapars AR.Dapars FR.Dapars
      #12         3         4         3 
Daparsgenes = list(AdultPolyA=DaparsGenes[["AdultPolyA"]], FetalPolyA=DaparsGenes[["FetalPolyA"]], AdultRiboZero=DaparsGenes[["AdultRiboZero"]], FetalRiboZero=DaparsGenes[["FetalRiboZero"]])
Dapars.IR = list()
    for (i in 1:4){
    tmp=IRclean[[i]]
    Dapars.IR[[i]]=tmp[which(tmp$GeneID %in% Daparsgenes[[i]] & tmp$p.diff <=0.05),]}
names(Dapars.IR) = c("AdultPolyA", "FetalPolyA", "AdultRibo", "FetalRibo")  
    elementLengths(Dapars.IR)
      #AdultPolyA FetalPolyA  AdultRibo  FetalRibo 
      #11         11         44         79  
    elementLengths(lapply(Dapars.IR, function(x) x[which(x$Sign=="MoreIRInNuc"),]))
      #AdultPolyA FetalPolyA  AdultRibo  FetalRibo 
      #11         11         44         79

