library(GenomicRanges)
library(ggplot2)
library(data.table)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/introns_object.rda")

### Check repetitive elements present in introns

intronsdf = unique(introns[[1]][,colnames(introns[[1]]) %in% c("Chr","Start","End","Direction","intronID","ensID","ID")])
rpmsk = read.table("./Dropbox/sorted_figures/new/github_controlled/RepeatMasker_genomewide_HG19.txt", header=T)
rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE)
intronsgr = makeGRangesFromDataFrame(intronsdf,seqnames.field="Chr",start.field="Start",end.field="End",
                                     strand.field="Direction",keep.extra.columns = T)
overlaps = findOverlaps(rpmskgr, intronsgr)

intronMap = rbind(data.frame(intronsdf[subjectHits(overlaps),], ContainsRepeat = "Yes", rpmsk[queryHits(overlaps),c("repName","repClass","repFamily")], AgeFrac = "Both", Dir = "All"),
                  data.frame(intronsdf[-unique(subjectHits(overlaps)),], ContainsRepeat = "No", "repName"="No Repeats","repClass"="No Repeats","repFamily"="No Repeats", AgeFrac = "Both", Dir = "All"))
intronMap = rbind(intronMap, do.call(rbind, Map(cbind, lapply(introns[elementNROWS(introns)>0 & elementNROWS(introns)<1000], 
                                                              function(x) unique(intronMap[intronMap$intronID %in% x$intronID,!colnames(intronMap) %in% c("AgeFrac","Dir")])),
                                                AgeFrac = list("In Adult", "In Prenatal", "In Cytoplasm", "In Cytoplasm", "In Nucleus", "In Nucleus"),
                                                Dir = list("Nuclear\nRetention","Nuclear\nRetention","Increasing\nRetention","Decreasing\nRetention","Increasing\nRetention","Decreasing\nRetention"))))
intronMap$Group = paste(intronMap$AgeFrac, intronMap$Dir, sep = ":")
b = intronMap[intronMap$Group=="Both:All",]
b = b[!b$intronID %in% intronMap[intronMap$Group!="Both:All", "intronID"],]
intronMap = rbind(b, intronMap[intronMap$Group!="Both:All",])

head(intronMap)
intronMap$AgeFrac = factor(intronMap$AgeFrac, levels = c("Both","In Adult","In Prenatal","In Cytoplasm","In Nucleus"))
intronMap$Dir = factor(intronMap$Dir, levels = c("All","Nuclear\nRetention","Increasing\nRetention","Decreasing\nRetention"))


## Count the number of repeats
intronMapdt = data.table(intronMap)
Freq = list(Name = intronMapdt[,list(Count = length(unique(intronID))), by = c("repName", "Dir", "AgeFrac", "Group")],
            Class = intronMapdt[,list(Count = length(unique(intronID))), by = c("repClass", "Dir", "AgeFrac", "Group")],
            Family = intronMapdt[,list(Count = length(unique(intronID))), by = c("repFamily", "Dir", "AgeFrac", "Group")],
            Repeat = intronMapdt[,list(Count = length(unique(intronID))), by = c("ContainsRepeat", "Dir", "AgeFrac", "Group")])
Freq = lapply(Freq, data.frame)
for (i in 1:length(Freq)){
  for (j in 1:length(unique(Freq[[i]][,"Group"]))){ 
  Freq[[i]][Freq[[i]][,"Group"]==unique(Freq[[i]][,"Group"])[j],"Total"] = sum(Freq[[i]][Freq[[i]][,"Group"]==unique(Freq[[i]][,"Group"])[j],"Count"])
  Freq[[i]][,"Percent"] = round(Freq[[i]][,"Count"]/Freq[[i]][,"Total"]*100, 2)
  }
}


## Plot frequencies of different types of repetitive elements

path = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/"

pdf(paste0(path,"repetitive_elements_bySigIntrons_incl_NoRepeats.pdf"), width = 8.5, height = 5.5)
for (i in 1:length(Freq)) {
g = ggplot(Freq[[i]][which(Freq[[i]][,"AgeFrac"] %in% c("In Adult","In Prenatal")),], 
       aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  facet_grid(. ~ AgeFrac) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
g = ggplot(Freq[[i]][which(Freq[[i]][,"AgeFrac"] %in% c("In Cytoplasm","In Nucleus")),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  facet_grid(. ~ AgeFrac) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Age: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
}
dev.off()

pdf(paste0(path,"repetitive_elements_AllIntrons_incl_NoRepeats.pdf"), width = 9, height = 6.5)
for (i in 2:length(Freq)) {
g = ggplot(Freq[[i]][which(Freq[[i]][,"Group"]=="Both:All"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in All Introns: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
}
dev.off()

pdf(paste0(path,"repetitive_elements_sigIntrons.pdf"), width = 8, height = 6)
for (i in 1:length(Freq)) {
g = ggplot(Freq[[i]][which((Freq[[i]][,"AgeFrac"] %in% c("In Adult", "In Prenatal")) & 
                             Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  facet_grid(. ~ AgeFrac) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
print(g)
g = ggplot(Freq[[i]][which((Freq[[i]][,"AgeFrac"] %in% c("In Cytoplasm", "In Nucleus")) & 
                             Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  facet_grid(. ~ AgeFrac) +
  labs(fill="") + ylab("Count") + xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Age: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
print(g)
}
dev.off()

pdf(paste0(path,"repetitive_elements_AllIntrons.pdf"), width = 9, height = 6.5)
for (i in 2:length(Freq)) {
  g = ggplot(Freq[[i]][which(Freq[[i]][,"Group"]=="Both:All" & Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
             aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
    labs(fill="") + ylab("Count") + xlab("") +
    ggtitle(paste0("Repetitive Elements Present in All Introns: ", names(Freq)[i])) +
    theme(title = element_text(size = 20)) + theme(text = element_text(size = 20))
  print(g)  
}
dev.off()



### Is there a difference between the proportion of repeat-containing introns in significantly/nonsignificantly retained introns?

intronMap = rbind(data.frame(intronsdf[subjectHits(overlaps),], ContainsRepeat = "Yes", rpmsk[queryHits(overlaps),c("repName","repClass","repFamily")]),
                  data.frame(intronsdf[-unique(subjectHits(overlaps)),], ContainsRepeat = "No", "repName"="No Repeats","repClass"="No Repeats","repFamily"="No Repeats"))
head(intronMap)

intronMap$Nuclear = ifelse(as.character(intronMap$intronID) %in% c(as.character(introns$`Adult:Nucleus-Increased`$intronID), as.character(introns$`Prenatal:Nucleus-Increased`$intronID)), "YES","NO")
intronMap$Cytoplasmic = ifelse(as.character(intronMap$intronID) %in% c(as.character(introns$`Adult:Cytoplasm-Increased`$intronID), as.character(introns$`Prenatal:Cytoplasm-Increased`$intronID)), "YES","NO")
intronMap$NuclearinAd = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Adult:Nucleus-Increased`$intronID), "YES","NO")
intronMap$CytoplasmicinAd = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Adult:Cytoplasm-Increased`$intronID), "YES","NO")
intronMap$NuclearinPren = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Prenatal:Nucleus-Increased`$intronID), "YES","NO")
intronMap$CytoplasmicinPren = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Prenatal:Cytoplasm-Increased`$intronID), "YES","NO")
intronMap$byFrac = ifelse(as.character(intronMap$intronID) %in% c(as.character(introns$`Adult:Nucleus-Increased`$intronID), as.character(introns$`Prenatal:Nucleus-Increased`$intronID),
                                                                  as.character(introns$`Adult:Cytoplasm-Increased`$intronID), as.character(introns$`Prenatal:Cytoplasm-Increased`$intronID)), "YES","NO")

intronMap$Increasing = ifelse(as.character(intronMap$intronID) %in% c(as.character(introns$`Cytoplasm:Adult-Increased`$intronID), as.character(introns$`Nucleus:Adult-Increased`)), "YES","NO")
intronMap$Decreasing = ifelse(as.character(intronMap$intronID) %in% c(as.character(introns$`Cytoplasm:Prenatal-Increased`$intronID), as.character(introns$`Nucleus:Prenatal-Increased`$intronID)), "YES","NO")
intronMap$IncreasinginCyt = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Cytoplasm:Adult-Increased`$intronID), "YES","NO")
intronMap$DecreasinginCyt = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Cytoplasm:Prenatal-Increased`$intronID), "YES","NO")
intronMap$IncreasinginNuc = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Nucleus:Adult-Increased`), "YES","NO")
intronMap$DecreasinginNuc = ifelse(as.character(intronMap$intronID) %in% as.character(introns$`Nucleus:Prenatal-Increased`$intronID), "YES","NO")
intronMap$byAge = ifelse(as.character(intronMap$intronID) %in% c(as.character(introns$`Cytoplasm:Adult-Increased`$intronID), as.character(introns$`Nucleus:Adult-Increased`),
                                                                 as.character(introns$`Cytoplasm:Prenatal-Increased`$intronID), as.character(introns$`Nucleus:Prenatal-Increased`$intronID)), "YES","NO")
cols = c("Nuclear","Cytoplasmic","NuclearinAd","CytoplasmicinAd","NuclearinPren","CytoplasmicinPren","Increasing","Decreasing",
         "IncreasinginCyt","DecreasinginCyt","IncreasinginNuc","DecreasinginNuc","byFrac","byAge")
fam = c("Alu","L1","L2","Simple_repeat")

tables = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(cols)) {
  for (j in 1:length(fam)) {
    tables[[i]][[j]] = data.frame(rep = c(length(unique(intronMap[intronMap$repFamily==fam[j] & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"])),
                                          length(unique(intronMap[intronMap$repFamily==fam[j] & 
                                                                    !intronMap$intronID %in% intronMap[intronMap$repFamily==fam[j] & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"],"intronID"]))),
                                  norep = c(length(unique(intronMap[intronMap[,colnames(intronMap)==cols[i]]=="YES" &
                                                                      !intronMap$intronID %in% intronMap[intronMap$repFamily==fam[j],"intronID"],"intronID"])),
                                            length(unique(intronMap[!intronMap$intronID %in% intronMap[intronMap$repFamily==fam[j],"intronID"] & 
                                                                      !intronMap$intronID %in% intronMap[intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"],"intronID"]))))
  }
  tables[[i]][[5]] = data.frame(rep = c(length(unique(intronMap[intronMap$ContainsRepeat=="Yes" & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"])),
                                        length(unique(intronMap[intronMap$ContainsRepeat=="Yes" & 
                                                                  !intronMap$intronID %in% intronMap[intronMap$ContainsRepeat=="Yes" & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"],"intronID"]))),
                                norep = c(length(unique(intronMap[intronMap$ContainsRepeat=="No" & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"])),
                                          length(unique(intronMap[intronMap$ContainsRepeat=="No" & 
                                                                    !intronMap$intronID %in% intronMap[intronMap$ContainsRepeat=="No" & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"],"intronID"]))))
  names(tables[[i]]) = c(fam, "ContainsRepeats")
}
names(tables) = cols
lapply(tables, lapply, sum)
tables = lapply(tables, lapply, fisher.test)

df = do.call(rbind, Map(cbind, Comparison = as.list(names(tables)), 
                        lapply(tables, function(t) do.call(rbind, Map(cbind, Repeat = as.list(names(t)),lapply(t, function(x) data.frame(pval = x$p.value, OddsRatio = x$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonsig_introns_ContainsRepeat.csv")
df[df$FDR<=0.05,]
#        Comparison         pval  OddsRatio          FDR
#1          Nuclear 1.480288e-16 0.03400510 1.036202e-15
#3      NuclearinAd 1.508800e-14 0.03825958 7.041067e-14
#5    NuclearinPren 5.946167e-04 0.00000000 1.189233e-03
#8       Decreasing 1.423377e-06 0.11341137 4.981820e-06
#10 DecreasinginCyt 3.463295e-03 0.11668319 6.060767e-03
#12 DecreasinginNuc 9.416230e-06 0.10888234 2.197120e-05
#13          byFrac 1.480288e-16 0.03400510 1.036202e-15
#14           byAge 4.138591e-06 0.13609814 1.158805e-05


## In repeat-containing introns, what is the proportion of ALU, L1,L2, and simple repeats in significantly/nonsignificantly retained fraction introns?

fam = c("Alu","L1","L2","Simple_repeat")

tables = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(cols)) {
  for (j in 1:length(fam)) {
    tables[[i]][[j]] = data.frame(rep = c(length(unique(intronMap[intronMap$repFamily==fam[j] & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"])),
                                          length(unique(intronMap[intronMap$repFamily==fam[j] & 
                                                                    !intronMap$intronID %in% intronMap[intronMap$repFamily==fam[j] & intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"],"intronID"]))),
                                  norep = c(length(unique(intronMap[intronMap$ContainsRepeat=="Yes" & intronMap[,colnames(intronMap)==cols[i]]=="YES" &
                                                                      !intronMap$intronID %in% intronMap[intronMap$repFamily==fam[j],"intronID"],"intronID"])),
                                            length(unique(intronMap[intronMap$ContainsRepeat=="Yes" & !intronMap$intronID %in% intronMap[intronMap$repFamily==fam[j],"intronID"] & 
                                                                      !intronMap$intronID %in% intronMap[intronMap[,colnames(intronMap)==cols[i]]=="YES","intronID"],"intronID"]))))
  }
  names(tables[[i]]) = fam
}
names(tables) = cols
lapply(tables, lapply, sum)
tables = lapply(tables, lapply, fisher.test)

df = do.call(rbind, Map(cbind, Comparison = as.list(names(tables)), 
                        lapply(tables, function(t) do.call(rbind, Map(cbind, Repeat = as.list(names(t)),lapply(t, function(x) data.frame(pval = x$p.value, OddsRatio = x$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonsig_introns_isAluL1L2simplerepeats.csv")
df[df$FDR<=0.05,]