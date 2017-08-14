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

# get repeat information from the UCSC table browser
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
# To get repeat information for hg19, change the group to "Repeats" and the "RepeatMasker" track. Label the output "repeatmasker.introns.1000.txt" and click "Get Output."
# Repeat for the remaining 1,234 introns


#######################################################################
######################### Repetitive Elements #########################
#######################################################################

### Check repetitive elements present in introns

rpmsk1000 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/repeatmasker.introns.1000.txt", header=T)
rpmsk2000 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/repeatmasker.introns.2000.txt", header=T)
rpmsk2234 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/repeatmasker.introns.2234.txt", header=T)
rpmsk = rbind(rpmsk1000,rpmsk2000,rpmsk2234)
# Track info can be found at http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema
simrp1000 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/simplerepeats.introns.1000.txt", header=T)
simrp2000 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/simplerepeats.introns.2000.txt", header=T)
simrp2234 = read.table("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/simplerepeats.introns.2234.txt", header=T)
simrp = rbind(simrp1000,simrp2000,simrp2234)
# track info can be found at http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=215014537&c=chrX&g=simpleRepeat
rpmsk = rpmsk[,which(colnames(rpmsk)!="strand")]
rpmsk = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE)
simrp = makeGRangesFromDataFrame(simrp, seqnames.field="chrom",start.field="chromStart",end.field="chromEnd",keep.extra.columns=TRUE)
coord = makeGRangesFromDataFrame(coord, keep.extra.columns = T)
hits = as.data.frame(findOverlaps(coord, rpmsk))
simhits = as.data.frame(findOverlaps(coord, simrp))
rps = srps = list()
for (i in 1:length(coord)){
  rps[[i]] = as.data.frame(rpmsk[hits[which(hits$queryHits==i),"subjectHits"]][,2:13])
  srps[[i]] = as.data.frame(simrp[simhits[which(simhits$queryHits==i),"subjectHits"]][,2:13])
  if (nrow(rps[[i]])>0){ rps[[i]][,"intronID"] = coord[i]$ID }
  if (nrow(srps[[i]])>0){srps[[i]][,"intronID"] = coord[i]$ID }
}
names(rps) = names(srps) = coord$ID
rps = do.call(rbind, rps)
srps = do.call(rbind,srps)

repeatintrons = Map(cbind, introns, lapply(introns, function(x) rps[match(x$intronID, rps$intronID),]))
length(repeatintrons)
for (i in c(1:11,13:15)){
  print(i)
  tmp = repeatintrons[[i]]
  tmp$repName = as.character(tmp$repName)
  tmp[is.na(tmp$repName),"repName"] = "No.Repeats"
  tmp$repClass = as.character(tmp$repClass)
  tmp[is.na(tmp$repClass),"repClass"] = "No.Repeats"
  tmp$repFamily = as.character(tmp$repFamily)
  tmp[is.na(tmp$repFamily),"repFamily"] = "No.Repeats"
  repeatintrons[[i]] = tmp
}

# Count the number of repeats
repNameFreq = lapply(repeatintrons[2:11], function(x) count(x, vars = "repName"))
repClassFreq = lapply(repeatintrons[2:11], function(x) count(x, vars = "repClass"))
repFamilyFreq = lapply(repeatintrons[2:11], function(x) count(x, vars = "repFamily"))
uniquenames = unique(unlist(lapply(repNameFreq, function(x) as.character(x$repName))))
uniqueclasses = unique(unlist(lapply(repClassFreq, function(x) as.character(x$repClass))))
uniqueFamily = unique(unlist(lapply(repFamilyFreq, function(x) as.character(x$repFamily))))

repName = data.frame(Names = rep.int(uniquenames,10), Count = NA, 
                     Group = factor(x=c(rep.int("Introns (Fraction)", length(uniquenames)),rep.int("Introns (Age)", length(uniquenames)),
                                        rep.int("Adult:Cytosol-Increased", length(uniquenames)),rep.int("Prenatal:Cytosol-Increased", length(uniquenames)),
                                        rep.int("Adult:Nucleus-Increased", length(uniquenames)),rep.int("Prenatal:Nucleus-Increased", length(uniquenames)),
                                        rep.int("Cytosol:Adult-Increased", length(uniquenames)),rep.int("Nucleus:Adult-Increased", length(uniquenames)),
                                        rep.int("Cytosol:Prenatal-Increased", length(uniquenames)),rep.int("Nucleus:Prenatal-Increased", length(uniquenames)))), 
                     Sample = factor(x=c(rep.int("NA", (length(uniquenames)*2)),rep.int("In Adult", length(uniquenames)),rep.int("In Prenatal", length(uniquenames)),
                                         rep.int("In Adult", length(uniquenames)),rep.int("In Prenatal", length(uniquenames)),rep.int("In Cytosol", length(uniquenames)),
                                         rep.int("In Nucleus", length(uniquenames)),rep.int("In Cytosol", length(uniquenames)),rep.int("In Nucleus", length(uniquenames)))), 
                     Dir = factor(x=c(rep.int("NA", (length(uniquenames)*2)),rep.int("Cytosolic\nRetention", (length(uniquenames)*2)),rep.int("Nuclear\nRetention", (length(uniquenames)*2)),
                                      rep.int("Increasing\nRetention", (length(uniquenames)*2)),rep.int("Decreasing\nRetention", (length(uniquenames)*2)))))
repClass = data.frame(Class = rep.int(uniqueclasses,10), Count = NA, 
                      Group = factor(x=c(rep.int("Introns (Fraction)", length(uniqueclasses)),rep.int("Introns (Age)", length(uniqueclasses)),
                                         rep.int("Adult:Cytosol-Increased", length(uniqueclasses)),rep.int("Prenatal:Cytosol-Increased", length(uniqueclasses)),
                                         rep.int("Adult:Nucleus-Increased", length(uniqueclasses)),rep.int("Prenatal:Nucleus-Increased", length(uniqueclasses)),
                                         rep.int("Cytosol:Adult-Increased", length(uniqueclasses)),rep.int("Nucleus:Adult-Increased", length(uniqueclasses)),
                                         rep.int("Cytosol:Prenatal-Increased", length(uniqueclasses)),rep.int("Nucleus:Prenatal-Increased", length(uniqueclasses)))),
                      Sample = factor(x=c(rep.int("NA", (length(uniqueclasses)*2)),rep.int("In Adult", length(uniqueclasses)),rep.int("In Prenatal", length(uniqueclasses)),
                                          rep.int("In Adult", length(uniqueclasses)),rep.int("In Prenatal", length(uniqueclasses)),rep.int("In Cytosol", length(uniqueclasses)),
                                          rep.int("In Nucleus", length(uniqueclasses)),rep.int("In Cytosol", length(uniqueclasses)),rep.int("In Nucleus", length(uniqueclasses)))), 
                      Dir = factor(x=c(rep.int("NA", (length(uniqueclasses)*2)),rep.int("Cytosolic\nRetention", (length(uniqueclasses)*2)),rep.int("Nuclear\nRetention", (length(uniqueclasses)*2)),
                                       rep.int("Increasing\nRetention", (length(uniqueclasses)*2)),rep.int("Decreasing\nRetention", (length(uniqueclasses)*2)))))
repFamily = data.frame(Family = rep.int(uniqueFamily,10), Count = NA, 
                       Group = factor(x=c(rep.int("Introns (Fraction)", length(uniqueFamily)),rep.int("Introns (Age)", length(uniqueFamily)),
                                          rep.int("Adult:Cytosol-Increased", length(uniqueFamily)),rep.int("Prenatal:Cytosol-Increased", length(uniqueFamily)),
                                          rep.int("Adult:Nucleus-Increased", length(uniqueFamily)),rep.int("Prenatal:Nucleus-Increased", length(uniqueFamily)),
                                          rep.int("Cytosol:Adult-Increased", length(uniqueFamily)),rep.int("Nucleus:Adult-Increased", length(uniqueFamily)),
                                          rep.int("Cytosol:Prenatal-Increased", length(uniqueFamily)),rep.int("Nucleus:Prenatal-Increased", length(uniqueFamily)))),
                       Sample = factor(x=c(rep.int("NA", (length(uniqueFamily)*2)),rep.int("In Adult", length(uniqueFamily)),rep.int("In Prenatal", length(uniqueFamily)),
                                           rep.int("In Adult", length(uniqueFamily)),rep.int("In Prenatal", length(uniqueFamily)),rep.int("In Cytosol", length(uniqueFamily)),
                                           rep.int("In Nucleus", length(uniqueFamily)),rep.int("In Cytosol", length(uniqueFamily)),rep.int("In Nucleus", length(uniqueFamily)))), 
                       Dir = factor(x=c(rep.int("NA", (length(uniqueFamily)*2)),rep.int("Cytosolic\nRetention", (length(uniqueFamily)*2)),rep.int("Nuclear\nRetention", (length(uniqueFamily)*2)),
                                        rep.int("Increasing\nRetention", (length(uniqueFamily)*2)),rep.int("Decreasing\nRetention", (length(uniqueFamily)*2)))))
comparisons = names(repNameFreq)
for (i in 1:length(uniquenames)){
  for (j in 1:length(repNameFreq)){
    repName[which(repName$Names==uniquenames[i] & repName$Group==comparisons[j]),2] = 
      repNameFreq[[j]][match(uniquenames[i], repNameFreq[[j]][,1]),"freq"]
  }}
for (i in 1:length(uniqueclasses)){
  for (j in 1:length(repClassFreq)){
    repClass[which(repClass$Class==uniqueclasses[i] & repClass$Group==comparisons[j]),2] = 
      repClassFreq[[j]][match(uniqueclasses[i], repClassFreq[[j]][,1]),"freq"]
  }}
for (i in 1:length(uniqueFamily)){
  for (j in 1:length(repFamilyFreq)){
    repFamily[which(repFamily$Family==uniqueFamily[i] & repFamily$Group==comparisons[j]),2] = 
      repFamilyFreq[[j]][match(uniqueFamily[i], repFamilyFreq[[j]][,1]),"freq"]
  }}

# Plot frequencies of different types of repetitive elements (including "no repeats")
ggplot(repName[which(repName$Sample=="In Adult" | repName$Sample=="In Prenatal"),], aes(x = Dir, y = Count, fill = Names)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repClass[which(repClass$Sample=="In Adult" | repClass$Sample=="In Prenatal"),], aes(x = Dir, y = Count, fill = Class)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repFamily[which(repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal"),], 
       aes(x = Dir, y = Count, fill = Family)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repName[which(repName$Sample=="In Cytosol" | repName$Sample=="In Nucleus"),], 
       aes(x = Dir, y = Count, fill = Names)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repClass[which(repClass$Sample=="In Cytosol" | repClass$Sample=="In Nucleus"),], 
       aes(x = Dir, y = Count, fill = Class)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repFamily[which(repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus"),], 
       aes(x = Dir, y = Count, fill = Family)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

# Plot frequencies of different types of repetitive elements (excluding "no repeats")
ggplot(repName[which((repName$Sample=="In Adult" | repName$Sample=="In Prenatal") & repName$Names!="No.Repeats"),], 
       aes(x = Dir, y = Count, fill = Names)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repClass[which((repClass$Sample=="In Adult" | repClass$Sample=="In Prenatal") & repClass$Class!="No.Repeats"),],
       aes(x = Dir, y = Count, fill = Class)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family!="No.Repeats"),], 
       aes(x = Dir, y = Count, fill = Family)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repName[which((repName$Sample=="In Cytosol" | repName$Sample=="In Nucleus") & repName$Names!="No.Repeats"),], 
       aes(x = Dir, y = Count, fill = Names)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repClass[which((repClass$Sample=="In Cytosol" | repClass$Sample=="In Nucleus") & repClass$Class!="No.Repeats"),],
       aes(x = Dir, y = Count, fill = Class)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
ggplot(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family!="No.Repeats"),], 
       aes(x = Dir, y = Count, fill = Family)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Sample) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle("Repetitive Elements Present in\nDifferentially Retained Introns by Age") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))

## Is there a difference between the proportion of repeat-containing introns in significantly/nonsignificantly retained fraction introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family!="No.Repeats"), "Count"])) # 68
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family=="No.Repeats"), "Count"])) # 187
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family!="No.Repeats"), "Count"])) # 320
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family=="No.Repeats"), "Count"])) # 1087
fisher.test(data.frame(c(68,187),c(320-68,1087-187)))
#p-value = 0.09947
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9363112 1.7869159
#sample estimates:
#  odds ratio 
#1.298459
# adult only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family!="No.Repeats"), "Count"])) # 47
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family=="No.Repeats"), "Count"])) # 115
fisher.test(data.frame(c(47,115),c(320-47,1087-115)))
#p-value = 0.04651
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9867554 2.1190699
#sample estimates:
#  odds ratio 
#1.454718
# prenatal only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family!="No.Repeats"), "Count"])) # 21
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family=="No.Repeats"), "Count"])) # 72
fisher.test(data.frame(c(21,72),c(320-21,1087-72)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5682908 1.6613859
#sample estimates:
#  odds ratio 
#0.9901112
# all Nuclear\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 67
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 185
fisher.test(data.frame(c(67,185),c(320-67,1087-185)))
#p-value = 0.1152
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9292028 1.7794860
#sample estimates:
#  odds ratio 
#1.290954
# all Cytosolic\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 2
fisher.test(data.frame(c(1,2),c(320-1,1087-2)))
#p-value = 0.5392
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0287364 32.7781600
#sample estimates:
#  odds ratio 
#1.699895

## between the proportion of repeat-containing introns in significantly retained introns in cytosol vs nucleus?
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 1 
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 67
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 2
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 185
fisher.test(data.frame(c(1,2),c(67,185)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.02306073 26.88253947
#sample estimates:
#  odds ratio 
#1.378729

## between the proportion of repeat-containing introns in significantly/nonsignificantly retained age introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family!="No.Repeats"), "Count"])) # 88
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family=="No.Repeats"), "Count"])) # 171
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family!="No.Repeats"), "Count"])) # 248
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family=="No.Repeats"), "Count"])) # 897
fisher.test(data.frame(c(88,171),c(248-88,897-171)))
#p-value = 1.346e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.690218 3.212252
#sample estimates:
#  odds ratio 
#2.333124
# cytosol only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family!="No.Repeats"), "Count"])) # 27
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family=="No.Repeats"), "Count"])) # 53
fisher.test(data.frame(c(27,53),c(248-27,897-53)))
#p-value = 0.01061
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.147616 3.230467
#sample estimates:
#  odds ratio 
#1.944298
# nucleus only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family!="No.Repeats"), "Count"])) # 61
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family=="No.Repeats"), "Count"])) # 118
fisher.test(data.frame(c(61,118),c(248-61,897-118)))
#p-value = 2.84e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.491657 3.085126
#sample estimates:
#  odds ratio 
#2.151842
# all Nuclear Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 39
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 101
fisher.test(data.frame(c(39,101),c(248-39,897-101)))
#p-value = 0.06278
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9583922 2.2216494
#sample estimates:
#  odds ratio 
#1.470116
# all Adult Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 49 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 70
fisher.test(data.frame(c(49,70),c(248-49,897-70)))
#p-value = 4.251e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.910338 4.396355
#sample estimates:
#  odds ratio 
#2.905605

## between the proportion of repeat-containing introns in significantly retained introns in adult vs prenatal?
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 49 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="No.Repeats"), "Count"])) # 39
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 70
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="No.Repeats"), "Count"])) # 101
fisher.test(data.frame(c(49,70),c(39,101)))
#p-value = 0.02595
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.043702 3.152920
#sample estimates:
#  odds ratio 
#1.808618

## In repeat-containing introns, what is the proportion of ALU in significantly/nonsignificantly retained fraction introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family=="Alu"), "Count"])) # 13
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 55
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family=="Alu"), "Count"])) # 81
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 239
fisher.test(data.frame(c(13,55),c(81-13,239-55)))
#p-value = 0.2109
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3014435 1.2805666
#sample estimates:
#  odds ratio 
#0.6404092
# adult only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family=="Alu"), "Count"])) # 10
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 37
fisher.test(data.frame(c(10,37),c(81-10,239-37)))
#p-value = 0.5875
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3239917 1.6818272
#sample estimates:
#  odds ratio 
#0.7695173
# prenatal only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family=="Alu"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 18
fisher.test(data.frame(c(3,18),c(81-3,239-18)))
#p-value = 0.3038
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.08693487 1.68558516
#sample estimates:
#  odds ratio 
#0.473148
# all Nuclear\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="Alu"), "Count"])) # 13
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 54
fisher.test(data.frame(c(13,54),c(81-13,239-54)))
#p-value = 0.2686
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3084188 1.3128428
#sample estimates:
#  odds ratio 
#0.6557703
# all Cytosolic\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="Alu"), "Count"])) # 0
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 1
fisher.test(data.frame(c(0,1),c(81-0,239-1)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 114.8693
#sample estimates:
#  odds ratio 
#0

## In repeat-containing introns, what is the proportion of ALU in in significantly retained introns in cytosol vs nucleus?
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="Alu"), "Count"])) # 0 
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="Alu"), "Count"])) # 13
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 54
fisher.test(data.frame(c(0,1),c(13,54)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 164.5783
#sample estimates:
#  odds ratio 
#0 

## In repeat-containing introns, what is the proportion of ALU in significantly/nonsignificantly retained age introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family=="Alu"), "Count"])) # 16
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 72
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family=="Alu"), "Count"])) # 71
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 177
fisher.test(data.frame(c(16,72),c(71-16,177-72)))
#p-value = 0.008057
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2103679 0.8253307
#sample estimates:
#  odds ratio 
#0.4256395
# cytosol only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family=="Alu"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 24
fisher.test(data.frame(c(3,24),c(71-3,177-24)))
#p-value = 0.04077
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0526617 0.9777675
#sample estimates:
#  odds ratio 
#0.2823631
# nucleus only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family=="Alu"), "Count"])) # 13
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 48
fisher.test(data.frame(c(13,48),c(71-13,177-48)))
#p-value = 0.1916
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2778936 1.2396240
#sample estimates:
#  odds ratio 
#0.6035362
# all Nuclear Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="Alu"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 36
fisher.test(data.frame(c(3,36),c(71-3,177-36)))
#p-value = 0.0009313
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.03304419 0.58036870
#sample estimates:
#  odds ratio 
#0.1736778
# all Adult Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="Alu"), "Count"])) # 13 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 36
fisher.test(data.frame(c(13,36),c(71-13,177-36)))
#p-value = 0.8602
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3974408 1.8464278
#sample estimates:
#  odds ratio 
#0.8783052

## In repeat-containing introns, what is the proportion of ALU in in significantly retained introns in adult vs prenatal?
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="Alu"), "Count"])) # 13 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="Alu"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 36
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="Alu" & repFamily$Family!="No.Repeats"), "Count"])) # 36
fisher.test(data.frame(c(13,36),c(3,36)))
#p-value = 0.02737
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.048261 25.315645
#sample estimates:
#  odds ratio 
#4.267787

## In repeat-containing introns, what is the proportion of L1 in significantly/nonsignificantly retained fraction introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family=="L1"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 65
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family=="L1"), "Count"])) # 22
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 298
fisher.test(data.frame(c(3,65),c(22-3,298-65)))
#p-value = 0.5882
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1042259 2.0146724
#sample estimates:
#  odds ratio 
#0.5668281
# adult only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family=="L1"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 44
fisher.test(data.frame(c(3,44),c(22-3,298-44)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1659044 3.2934220
#sample estimates:
#  odds ratio 
#0.9117564
# prenatal only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family=="L1"), "Count"])) # 0
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 21
fisher.test(data.frame(c(0,21),c(22-0,298-21)))
#p-value = 0.3796
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00000 2.62893
#sample estimates:
#  odds ratio 
#0
# all Nuclear\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="L1"), "Count"])) # 13
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 54
fisher.test(data.frame(c(3,65),c(22-3,298-65)))
#p-value = 0.3796
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00000 2.62893
#sample estimates:
#  odds ratio 
#0
# all Cytosolic Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="L1"), "Count"])) # 0
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 1
fisher.test(data.frame(c(0,1),c(22-0,298-1)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.000 523.956
#sample estimates:
#  odds ratio 
#0

## In repeat-containing introns, what is the proportion of L1 in in significantly retained introns in cytosol vs nucleus?
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="L1"), "Count"])) # 0 
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="L1"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 64
fisher.test(data.frame(c(0,1),c(3,64)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 834.0016
#sample estimates:
#  odds ratio 
#0 

## In repeat-containing introns, what is the proportion of L1 in significantly/nonsignificantly retained age introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family=="L1"), "Count"])) # 4
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 84
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family=="L1"), "Count"])) # 16
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 232
fisher.test(data.frame(c(4,84),c(16-4,232-84)))
#p-value = 0.4299
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1340304 2.0214308
#sample estimates:
#  odds ratio 
#0.5884661
# cytosol only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family=="L1"), "Count"])) # 0
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 27
fisher.test(data.frame(c(0,27),c(16-0,232-27)))
#p-value = 0.2299
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.00000 2.11177
#sample estimates:
#  odds ratio 
#0
# nucleus only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family=="L1"), "Count"])) # 4
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 57
fisher.test(data.frame(c(4,57),c(16-4,232-57)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2313622 3.5513235
#sample estimates:
#  odds ratio 
#1.023297
# all Prenatal Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="L1"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 38
fisher.test(data.frame(c(1,38),c(16-1,232-38)))
#p-value = 0.4793
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.007883896 2.351199005
#sample estimates:
#  odds ratio 
#0.3414002
# all Adult Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="L1"), "Count"])) # 3 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 46
fisher.test(data.frame(c(3,46),c(16-3,232-46)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1638889 3.5970086
#sample estimates:
#  odds ratio 
#0.9333779

## In repeat-containing introns, what is the proportion of L1 in in significantly retained introns in adult vs prenatal?
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="L1"), "Count"])) # 3 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="L1"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 46
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="L1" & repFamily$Family!="No.Repeats"), "Count"])) # 38
fisher.test(data.frame(c(3,46),c(1,38)))
#p-value = 0.6263
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1882939 133.3736704
#sample estimates:
#  odds ratio 
#2.455071

## In repeat-containing introns, what is the proportion of L2 in significantly/nonsignificantly retained fraction introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family=="L2"), "Count"])) # 9
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 59
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family=="L2"), "Count"])) # 44
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 276
fisher.test(data.frame(c(9,59),c(44-9,276-59)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3781307 2.1521339
#sample estimates:
#  odds ratio 
#0.9459259
# adult only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family=="L2"), "Count"])) # 6
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 41
fisher.test(data.frame(c(6,41),c(44-6,276-41)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2939733 2.3503259
#sample estimates:
#  odds ratio 
#0.9052822
# prenatal only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family=="L2"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 18
fisher.test(data.frame(c(3,18),c(44-3,276-18)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1895543 3.8313121
#sample estimates:
#  odds ratio 
#1.04863
# all Nuclear\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="L2"), "Count"])) # 8
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 59
fisher.test(data.frame(c(8,59),c(44-8,276-59)))
#p-value = 0.6952
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3113465 1.9151389
#sample estimates:
#  odds ratio 
#0.8178292
# all Cytosolic Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="L2"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 0
fisher.test(data.frame(c(1,0),c(44-1,276-0)))
#p-value = 0.1375
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1608391       Inf
#sample estimates:
#  odds ratio 
#Inf

## In repeat-containing introns, what is the proportion of L2 in in significantly retained introns in cytosol vs nucleus?
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="L2"), "Count"])) # 1 
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="L2"), "Count"])) # 8
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 0
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 59
fisher.test(data.frame(c(1,0),c(8,59)))
#p-value = 0.1324
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1680911       Inf
#sample estimates:
#  odds ratio 
#Inf 

## In repeat-containing introns, what is the proportion of L2 in significantly/nonsignificantly retained age introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family=="L2"), "Count"])) # 10
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 78
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family=="L2"), "Count"])) # 22
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 226
fisher.test(data.frame(c(10,78),c(22-10,226-78)))
#p-value = 0.3528
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5824394 4.1888503
#sample estimates:
#  odds ratio 
#1.578201
# cytosol only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family=="L2"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 24
fisher.test(data.frame(c(3,24),c(22-3,226-24)))
#p-value = 0.717
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2345641 5.0233143
#sample estimates:
#  odds ratio 
#1.327279
# nucleus only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family=="L2"), "Count"])) # 7
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 54
fisher.test(data.frame(c(7,54),c(22-7,226-54)))
#p-value = 0.4389
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4855653 4.1146998
#sample estimates:
#  odds ratio 
#1.483853
# all Prenatal Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="L2"), "Count"])) # 6
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 33
fisher.test(data.frame(c(6,33),c(22-6,226-33)))
#p-value = 0.1284
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6518612 6.4332697
#sample estimates:
#  odds ratio 
#2.184656
# all Adult Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="L2"), "Count"])) # 4 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 45
fisher.test(data.frame(c(4,45),c(22-4,226-45)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2098217 2.9014020
#sample estimates:
#  odds ratio 
#0.8942214

## In repeat-containing introns, what is the proportion of L2 in in significantly retained introns in adult vs prenatal?
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="L2"), "Count"])) # 4 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="L2"), "Count"])) # 6
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 45
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="L2" & repFamily$Family!="No.Repeats"), "Count"])) # 33
fisher.test(data.frame(c(4,45),c(6,33)))
#p-value = 0.328
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.09443793 2.26968605
#sample estimates:
#  odds ratio 
#0.4929316

## In repeat-containing introns, what is the proportion of simple repeats in significantly/nonsignificantly retained fraction introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family=="Simple_repeat"), "Count"])) # 11
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult" | repFamily$Sample=="In Prenatal") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 57
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family=="Simple_repeat"), "Count"])) # 24
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Fraction)" & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 296
fisher.test(data.frame(c(11,57),c(24-11,296-57)))
#p-value = 0.006974
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.355914 9.046724
#sample estimates:
#  odds ratio 
#3.52981
# adult only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family=="Simple_repeat"), "Count"])) # 8
sum(na.omit(repFamily[which((repFamily$Sample=="In Adult") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 39
fisher.test(data.frame(c(8,39),c(24-8,296-39)))
#p-value = 0.01385
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.135907 8.787496
#sample estimates:
#  odds ratio 
#3.277978
# prenatal only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family=="Simple_repeat"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Sample=="In Prenatal") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 18
fisher.test(data.frame(c(3,18),c(24-3,296-18)))
#p-value = 0.2011
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3845129 8.4777160
#sample estimates:
#  odds ratio 
#2.199099
# all Nuclear\nRetention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 11
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 56
fisher.test(data.frame(c(11,56),c(24-11,296-56)))
#p-value = 0.006696
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.384785 9.252297
#sample estimates:
#  odds ratio 
#3.607416
# all Cytosolic Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 0
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 1
fisher.test(data.frame(c(0,1),c(24-0,296-1)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 477.4194
#sample estimates:
#  odds ratio 
#0

## In repeat-containing introns, what is the proportion of simple repeats in in significantly retained introns in cytosol vs nucleus?
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 0 
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 11
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Cytosol-Increased" | repFamily$Group=="Prenatal:Cytosol-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 1
sum(na.omit(repFamily[which((repFamily$Group=="Adult:Nucleus-Increased" | repFamily$Group=="Prenatal:Nucleus-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 56
fisher.test(data.frame(c(0,1),c(11,56)))
#p-value = 1
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.0000 201.4579
#sample estimates:
#  odds ratio 
#0 

## In repeat-containing introns, what is the proportion of simple repeats in significantly/nonsignificantly retained age introns?
# all sig introns vs all ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family=="Simple_repeat"), "Count"])) # 12
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol" | repFamily$Sample=="In Nucleus") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 76
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family=="Simple_repeat"), "Count"])) # 14
sum(na.omit(repFamily[which(repFamily$Group=="Introns (Age)" & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 234
fisher.test(data.frame(c(12,76),c(14-12,234-76)))
#p-value = 0.000107
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.650573 116.444542
#sample estimates:
#  odds ratio 
#12.35338
# cytosol only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family=="Simple_repeat"), "Count"])) # 5
sum(na.omit(repFamily[which((repFamily$Sample=="In Cytosol") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 22
fisher.test(data.frame(c(5,22),c(14-5,234-22)))
#p-value = 0.0106
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.278135 19.536780
#sample estimates:
#  odds ratio 
#5.292931
# nucleus only sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family=="Simple_repeat"), "Count"])) # 7
sum(na.omit(repFamily[which((repFamily$Sample=="In Nucleus") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 54
fisher.test(data.frame(c(7,54),c(14-7,234-54)))
#p-value = 0.04798
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9462269 11.6107667
#sample estimates:
#  odds ratio 
#3.313398 
# all Prenatal Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 36
fisher.test(data.frame(c(3,36),c(14-3,234-36)))
#p-value = 0.4671
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2557606 6.0583376
#sample estimates:
#  odds ratio 
#1.497236
# all Adult Retention sig introns vs ns introns
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 9 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 40
fisher.test(data.frame(c(9,40),c(14-9,234-40)))
#p-value = 0.0002146
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.442843 34.558576
#sample estimates:
#  odds ratio 
#8.616687

## In repeat-containing introns, what is the proportion of simple repeats in in significantly retained introns in adult vs prenatal?
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 9 
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family=="Simple_repeat"), "Count"])) # 3
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Adult-Increased" | repFamily$Group=="Nucleus:Adult-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 40
sum(na.omit(repFamily[which((repFamily$Group=="Cytosol:Prenatal-Increased" | repFamily$Group=="Nucleus:Prenatal-Increased") & repFamily$Family!="Simple_repeat" & repFamily$Family!="No.Repeats"), "Count"])) # 36
fisher.test(data.frame(c(9,40),c(3,36)))
#p-value = 0.2136
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6048291 16.5253128
#sample estimates:
#  odds ratio 
#2.671732