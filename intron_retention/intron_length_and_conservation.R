library(GenomicRanges)
library(ggplot2)
library(plyr)

### Prepare Intron Lists

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames
lapply(IRres, head)
# Filter introns
allIntrons = IRres[[1]]
allIntrons = allIntrons[grep("clean", allIntrons$GeneIntronDetails, fixed=T),]
string = unlist(strsplit(as.character(allIntrons$GeneIntronDetails), "/", fixed = TRUE), recursive = FALSE)
allIntrons = data.frame(allIntrons[,c(1:4,6)], genes = string[grep("ENSG", string)],
                        intronID = paste0("chr",allIntrons$Chr,":",allIntrons$Start,"-",allIntrons$End))
# read in Differential IR results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
nonconst = list()
for (i in 1:length(comps)){
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
}
names(nonconst) = comps
elementNROWS(nonconst)
string = lapply(nonconst, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE))
IR.diff = lapply(nonconst, function(y) y$A.IRratio - y$B.IRratio)
dIR = Map(cbind, nonconst, ensID = lapply(string, function(y) y[grep("ENSG", y)]), 
          comments = lapply(string, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])), 
          IR.diff = IR.diff, Sign = lapply(IR.diff, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")), 
          intronID = lapply(nonconst, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End)))
dIRclean = lapply(dIR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" &
                                             y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                                             y$comments=="clean"),])
sigdIR = lapply(dIRclean, function(x) x[which(x$p.diff<=0.05),])
sigdIR = lapply(sigdIR, function(x) split(x, x$Sign))
sigdIR = unlist(sigdIR, recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                  "Cytosol-Increased", "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)]))
introns = c("All Introns" = list(allIntrons), "Introns (Fraction)" = list(do.call(rbind, lapply(dIRclean[1:2], function(x) x[,c(1:4,6,8,30,32,33,34)]))), 
            "Introns (Age)" = list(do.call(rbind, lapply(dIRclean[3:4], function(x) x[,c(1:4,6,8,30,32,33,34)]))), sigdIR)
introns[["Introns (Fraction)"]] = introns[["Introns (Fraction)"]][unique(introns[["Introns (Fraction)"]][,"intronID"]),]
introns[["Introns (Age)"]] = introns[["Introns (Age)"]][unique(introns[["Introns (Age)"]][,"intronID"]),]

# get conservation information from the UCSC table browser
elementNROWS(dIRclean)
coord = do.call(rbind, lapply(dIRclean, function(x) x[,c(1:3,6)]))
coord = coord[!duplicated(coord), ]
coord$Chr = paste0("chr",coord$Chr)
coord$ID = paste0(coord$Chr,":",coord$Start,"-",coord$End)
dim(coord) # 2234    3
write.table(coord[1:1000,1:3], file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.1000.txt",
            quote = F, sep = "\t", row.names = F, col.names = F) # The table browser can handle up to 1000 entries at a time
write.table(coord[1001:2000,1:3], file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.2000.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(coord[2001:2234,1:3], file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.2234.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)
# Go to the table browser, and under regions click "define regions" and upload the first 1000 regions in the list written above. Click submit.
# To get GERP scores for hg19, go to the "Comparative Genomics" group and "GERP" track. Change the filter 10 million lines. Label the output "repeatmasker.introns.1000.txt" and click "Get Output."
# Repeat for the remaining 1,234 introns

#################################################################
######################### intron length #########################
#################################################################

### Check length distribution of differentially retained introns

introns = Map(cbind, introns, length = lapply(introns, function(x) abs(x$End - x$Start)))
elementNROWS(introns)
introns = Map(cbind, introns[c(1:11,13:15)], 
              Comparison = list("All Introns", "Introns (Fraction)", "Introns (Age)","Adult:Cytosol-Increased","Adult:Nucleus-Increased",
                                "Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased","Cytosol:Adult-Increased","Cytosol:Prenatal-Increased",
                                "Nucleus:Adult-Increased","Nucleus:Prenatal-Increased","Nucleus-Increased", "Adult-Increased", "Prenatal-Increased"))
length = do.call(rbind, lapply(introns, function(x) data.frame(length = x$length, Comparison=x$Comparison)))
head(length)
levels(length$Comparison) = c("All Introns", "Introns (Fraction)", "Introns (Age)","Adult:Cytosol-Increased","Prenatal:Cytosol-Increased","Adult:Nucleus-Increased",
                              "Prenatal:Nucleus-Increased","Cytosol:Adult-Increased","Nucleus:Adult-Increased","Cytosol:Prenatal-Increased",
                              "Nucleus:Prenatal-Increased","Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")

# Plot density distribution of intron lengths using all "clean" introns as background
ggplot(length[which(length$Comparison=="All Introns" | length$Comparison=="Adult:Cytosol-Increased" | 
                      length$Comparison=="Adult:Nucleus-Increased" |
                      length$Comparison=="Prenatal:Cytosol-Increased" | length$Comparison=="Prenatal:Nucleus-Increased"),],
       aes(x=length/1000)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("") + 
  xlab("Intron Length (Kb)") +
  ggtitle("Intron Length By Group") +
  xlim(0,10) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(length[which(length$Comparison=="All Introns" | length$Comparison=="Cytosol:Adult-Increased" | 
                      length$Comparison=="Cytosol:Prenatal-Increased" |
                      length$Comparison=="Nucleus:Adult-Increased" | length$Comparison=="Nucleus:Prenatal-Increased"),],
       aes(x=length/1000)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("") + 
  xlab("Intron Length (Kb)") +
  ggtitle("Intron Length By Group") +
  xlim(0,10) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Plot density distribution of intron lengths using all introns reported in dIR output as background
ggplot(length[which(length$Comparison=="Introns (Fraction)" | length$Comparison=="Adult:Cytosol-Increased" | 
                      length$Comparison=="Adult:Nucleus-Increased" |
                      length$Comparison=="Prenatal:Cytosol-Increased" | length$Comparison=="Prenatal:Nucleus-Increased"),],
       aes(x=length/1000)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("") + 
  xlab("Intron Length (Kb)") +
  ggtitle("Intron Length By Group") +
  xlim(0,10) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(length[which(length$Comparison=="Introns (Age)" | length$Comparison=="Cytosol:Adult-Increased" | 
                      length$Comparison=="Cytosol:Prenatal-Increased" |
                      length$Comparison=="Nucleus:Adult-Increased" | length$Comparison=="Nucleus:Prenatal-Increased"),],
       aes(x=length/1000)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("") + 
  xlab("Intron Length (Kb)") +
  ggtitle("Intron Length By Group") +
  xlim(0,10) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.55)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Measure length differences in groups of introns
t.test(c(length[length$Comparison=="Adult:Cytosol-Increased","length"],length[length$Comparison=="Prenatal:Cytosol-Increased","length"]),
       c(length[length$Comparison=="Adult:Nucleus-Increased","length"],length[length$Comparison=="Prenatal:Nucleus-Increased","length"]))
#data: pooled cytosolic vs pooled nuclear introns
#t = 0.62686, df = 183.33, p-value = 0.5315
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -125.6413  242.6588
#sample estimates:
#  mean of x mean of y 
#556.7346  498.2258
t.test(length[length$Comparison=="Prenatal:Cytosol-Increased","length"],length[length$Comparison=="Prenatal:Nucleus-Increased","length"])
#data: Prenatal cytosolic vs Prenatal nuclear introns
#t = 0.65007, df = 181.32, p-value = 0.5165
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -124.8337  247.5038
#sample estimates:
#  mean of x mean of y 
#560.4438  499.1087
t.test(c(length[length$Comparison=="Adult:Cytosol-Increased","length"],length[length$Comparison=="Prenatal:Cytosol-Increased","length"],
         length[length$Comparison=="Adult:Nucleus-Increased","length"],length[length$Comparison=="Prenatal:Nucleus-Increased","length"]),
       length[length$Comparison=="Introns (Fraction)","length"])
#data: pooled significantly differentially retained introns by fraction vs all introns reported in comparison
#t = -9.7567, df = 830.88, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1039.3288  -691.1897
#sample estimates:
#  mean of x mean of y 
#535.3961 1400.6553
t.test(c(length[length$Comparison=="Cytosol:Adult-Increased","length"],length[length$Comparison=="Nucleus:Adult-Increased","length"]),
       c(length[length$Comparison=="Cytosol:Prenatal-Increased","length"],length[length$Comparison=="Nucleus:Prenatal-Increased","length"]))
#data: pooled adult vs pooled prenatal introns
#t = -1.3803, df = 256.93, p-value = 0.1687
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -543.37728   95.53985
#sample estimates:
#  mean of x mean of y 
#752.4500  976.3687
t.test(length[length$Comparison=="Cytosol:Adult-Increased","length"],length[length$Comparison=="Cytosol:Prenatal-Increased","length"])
#data: Cytosol:adult vs Cytosol:prenatal introns
#t = 0.12642, df = 33.723, p-value = 0.9002
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -634.4213  718.5585
#sample estimates:
#  mean of x mean of y 
#1072.412  1030.343
t.test(length[length$Comparison=="Nucleus:Adult-Increased","length"],length[length$Comparison=="Nucleus:Prenatal-Increased","length"])
#data: Nucleus:adult vs Nucleus:prenatal introns
#t = -1.1565, df = 106.55, p-value = 0.2501
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -648.0555  170.5374
#sample estimates:
#  mean of x mean of y 
#666.1111  904.8701
t.test(c(length[length$Comparison=="Cytosol:Adult-Increased","length"],length[length$Comparison=="Nucleus:Adult-Increased","length"],
         length[length$Comparison=="Cytosol:Prenatal-Increased","length"],length[length$Comparison=="Nucleus:Prenatal-Increased","length"]),
       length[length$Comparison=="Introns (Age)","length"])
#data: pooled significantly differentially retained introns by age vs all introns reported in comparison
#t = -5.1573, df = 681.04, p-value = 3.289e-07
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1023.1356  -458.9016
#sample estimates:
#  mean of x mean of y 
#907.2046 1648.2232

#############################################################################
######################### evolutionary conservation #########################
#############################################################################

### Check evolutionary conservation via GERP score

gerp1000 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP.introns.1000.txt", comment.char = "t", 
                      col.names = c("chromosome", "Start", "End", "Score"))
gerp2000 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP.introns.2000.txt", comment.char = "t", 
                      col.names = c("chromosome", "Start", "End", "Score"))
gerp2234 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP.introns.2234.txt", comment.char = "t", 
                      col.names = c("chromosome", "Start", "End", "Score"))
gerp = rbind(gerp1000,gerp2000,gerp2234)
gerp = makeGRangesFromDataFrame(gerp, keep.extra.columns=TRUE)
coord = makeGRangesFromDataFrame(coord, strand.field="Direction", keep.extra.columns=T)
hits = as.data.frame(findOverlaps(coord, gerp))
cons = list()
for (i in 1:length(coord)){
  cons[[i]] = gerp[hits[which(hits$queryHits==i),"subjectHits"]]$Score
}
names(cons) = coord$ID
meancons = lapply(cons, function(x) mean(x))
meancons = data.frame(intronID = names(meancons), mean.GERP = unlist(meancons))
head(meancons)

# Plot mean GERP scores using all "clean" introns as background
introns = Map(cbind, introns, lapply(introns, function(x) meancons[match(x$intronID, meancons$intronID),]))
gerp = do.call(rbind, lapply(introns, function(x) data.frame(mean.GERP = x$mean.GERP, Comparison=x$Comparison)))
gerp$Samples = gerp$Dir = "NA"
gerp[grep("Adult:", gerp$Comparison), "Samples"] = "In Adult"
gerp[grep("Prenatal:", gerp$Comparison), "Samples"] = "In Prenatal"
gerp[grep("Nucleus:", gerp$Comparison), "Samples"] = "In Nucleus"
gerp[grep("Cytosol:", gerp$Comparison), "Samples"] = "In Cytosol"
gerp[grep("Adult-", gerp$Comparison), "Dir"] = "Increasing\nRetention"
gerp[grep("Prenatal-", gerp$Comparison), "Dir"] = "Decreasing\nRetention"
gerp[grep("Nucleus-", gerp$Comparison), "Dir"] = "Nuclear\nRetention"
gerp[grep("Cytosol-", gerp$Comparison), "Dir"] = "Cytosolic\nRetention"

ggplot(gerp[which(gerp$Samples=="In Adult" | gerp$Samples=="In Prenatal"),],
       aes(x=Dir, y=mean.GERP, fill=Samples), color=Samples) + geom_boxplot() +
  xlab("") + 
  ylab("Mean GERP") +
  ggtitle("Mean Base Conservation (GERP) in Introns\nDifferentially Retained By Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(gerp[which(gerp$Samples=="In Cytosol" | gerp$Samples=="In Nucleus"),],
       aes(x=Dir, y=mean.GERP, fill=Samples), color=Samples) + geom_boxplot() +
  xlab("") + 
  ylab("Mean GERP") +
  ggtitle("Mean Base Conservation (GERP) in Introns\nDifferentially Retained Over Development") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))

# Measure gerp differences in groups of introns
t.test(c(gerp[gerp$Comparison=="Adult:Cytosol-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Cytosol-Increased","mean.GERP"]),
       c(gerp[gerp$Comparison=="Adult:Nucleus-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Nucleus-Increased","mean.GERP"]))
#data: pooled cytosolic vs pooled nuclear introns
#t = 0.38464, df = 146.15, p-value = 0.7011
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2020140  0.2996507
#sample estimates:
#  mean of x  mean of y 
#-0.5531988 -0.6020171
t.test(gerp[gerp$Comparison=="Prenatal:Cytosol-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Nucleus-Increased","mean.GERP"])
#data: Prenatal cytosolic vs Prenatal nuclear introns
#t = 0.41074, df = 144.77, p-value = 0.6819
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2008175  0.3061764
#sample estimates:
#  mean of x  mean of y 
#-0.5539564 -0.6066359
t.test(c(gerp[gerp$Comparison=="Adult:Cytosol-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Cytosol-Increased","mean.GERP"],
         gerp[gerp$Comparison=="Adult:Nucleus-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Nucleus-Increased","mean.GERP"]),
       gerp[gerp$Comparison=="Introns (Fraction)","mean.GERP"])
#data: pooled significantly differentially retained introns by fraction vs all introns reported in comparison
#t = -0.66658, df = 389.49, p-value = 0.5054
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.16454257  0.08121887
#sample estimates:
#  mean of x  mean of y 
#-0.5710031 -0.5293412
t.test(c(gerp[gerp$Comparison=="Cytosol:Adult-Increased","mean.GERP"],gerp[gerp$Comparison=="Nucleus:Adult-Increased","mean.GERP"]),
       c(gerp[gerp$Comparison=="Cytosol:Prenatal-Increased","mean.GERP"],gerp[gerp$Comparison=="Nucleus:Prenatal-Increased","mean.GERP"]))
#data: pooled adult vs pooled prenatal introns
#t = 0.89377, df = 161.03, p-value = 0.3728
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1301288  0.3453030
#sample estimates:
#  mean of x  mean of y 
#-0.3517381 -0.4593253
t.test(gerp[gerp$Comparison=="Cytosol:Adult-Increased","mean.GERP"],gerp[gerp$Comparison=="Cytosol:Prenatal-Increased","mean.GERP"])
#data: Cytosol:adult vs Cytosol:prenatal introns
#t = 0.22215, df = 20.084, p-value = 0.8264
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5150967  0.6379259
#sample estimates:
#  mean of x  mean of y 
#-0.4640919 -0.5255065
t.test(gerp[gerp$Comparison=="Nucleus:Adult-Increased","mean.GERP"],gerp[gerp$Comparison=="Nucleus:Prenatal-Increased","mean.GERP"])
#data: Nucleus:adult vs Nucleus:prenatal introns
#t = 0.33485, df = 137.56, p-value = 0.7382
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2464155  0.3468878
#sample estimates:
#  mean of x  mean of y 
#-0.3214205 -0.3716566
t.test(c(gerp[gerp$Comparison=="Cytosol:Adult-Increased","mean.GERP"],gerp[gerp$Comparison=="Nucleus:Adult-Increased","mean.GERP"],
         gerp[gerp$Comparison=="Cytosol:Prenatal-Increased","mean.GERP"],gerp[gerp$Comparison=="Nucleus:Prenatal-Increased","mean.GERP"]),
       gerp[gerp$Comparison=="Introns (Age)","mean.GERP"])
#data: pooled significantly differentially retained introns by age vs all introns reported in comparison
#t = -1.0988, df = 505.41, p-value = 0.2724
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.21335355  0.06029901
#sample estimates:
#  mean of x  mean of y 
#-0.4260937 -0.3495665

# More info about GERP can be found at http://mendel.stanford.edu/SidowLab/downloads/gerp/
