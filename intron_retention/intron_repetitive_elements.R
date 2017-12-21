library(GenomicRanges)
library(ggplot2)
library(data.table)
library(plyr)

### Prepare Intron Lists

comps = c("Adult_PolyA_Zone_cleanIntrons_adultShared","Fetal_PolyA_Zone_cleanIntrons_prenatalShared",
          "Cytosol_PolyA_Age_cleanIntrons_cytosolShared","Nuclear_PolyA_Age_cleanIntrons_nucleusShared")
IRcomp = list()
for (i in 1:length(comps)){
  IRcomp[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/", comps[i], ".tab"), header = TRUE, comment.char="#")
}
names(IRcomp) = c("Adult_byFraction","Fetal_byFraction","Cytosol_byAge","Nuclear_byAge")
IRcomp = Map(cbind, IRcomp,
             intronID = lapply(IRcomp, function(x) paste0(x$Intron.GeneName.GeneID,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)),
             ensID = lapply(lapply(IRcomp, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE)), function(y) y[grep("ENSG", y)]),
             IR.diff = lapply(IRcomp, function(y) y$A.IRratio - y$B.IRratio),
             Sign = lapply(IRcomp, function(y) ifelse((y$A.IRratio - y$B.IRratio) < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")))
full = list(adult = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-adultShared.txt", header = TRUE),
            prenatal = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-prenatalShared.txt", header = TRUE),
            cytosol = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br5339C1/IRFinder-IR-nondir-cleanIntrons-cytosolShared.txt", header = TRUE),
            nucleus = read.table("/Users/amanda/Dropbox/sorted_figures/IRfinder/PolyA/Br2046N/IRFinder-IR-nondir-cleanIntrons-nucleusShared.txt", header = TRUE))
total = as.list(elementNROWS(full))
names(total) = names(IRcomp)
IRcomp = Map(cbind, IRcomp, padj = mapply(function(p,t) p.adjust(p, method = "fdr", n = t), lapply(IRcomp, function(x) x$p.diff), total))

sigdIR = lapply(IRcomp, function(x) x[which(x$padj<=0.05),])
sigdIR = unlist(lapply(sigdIR, function(x) split(x, x$Sign)), recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased")

introns = c("Introns (Fraction)" = list(do.call(rbind, IRcomp[1:2])),"Introns (Age)" = list(do.call(rbind, IRcomp[3:4])), sigdIR)
for (i in 1:length(introns)) { if (nrow(introns[[i]]) > 0) { introns[[i]][,"Chr"] = paste0("chr", introns[[i]][,"Chr"]) } }
elementNROWS(introns)
intronsdf = do.call(rbind, Map(cbind, introns[elementNROWS(introns)>0], Group = as.list(names(introns)[elementNROWS(introns)>0]))) 


### Check repetitive elements present in introns

rpmsk = read.table("./Dropbox/sorted_figures/new/github_controlled/RepeatMasker_genomewide_HG19.txt", header=T)
rpmskgr = makeGRangesFromDataFrame(rpmsk, seqnames.field="genoName",start.field="genoStart",end.field="genoEnd",keep.extra.columns=TRUE)
intronsgr = makeGRangesFromDataFrame(intronsdf,seqnames.field="Chr",start.field="Start",end.field="End",
                                     strand.field="Direction",keep.extra.columns = T)
overlaps = findOverlaps(rpmskgr, intronsgr)

intronMap = rbind(data.frame(intronsdf[subjectHits(overlaps),], ContainsRepeat = "Yes", rpmsk[queryHits(overlaps),c("repName","repClass","repFamily")]),
                  data.frame(intronsdf[-unique(subjectHits(overlaps)),], ContainsRepeat = "No", "repName"="No Repeats","repClass"="No Repeats","repFamily"="No Repeats"))

intronMap[c(grep("Fraction", intronMap$Group),grep("Adult:", intronMap$Group),grep("Prenatal:", intronMap$Group)),"Comparison"] = "Fraction"
intronMap[c(grep("Age", intronMap$Group),grep("Cytosol:", intronMap$Group),grep("Nucleus:", intronMap$Group)),"Comparison"] = "Age"
intronMap[grep("Introns", intronMap$Group),"AgeFrac"] = "Both"
intronMap[grep("Adult:", intronMap$Group),"AgeFrac"] = "Adult"
intronMap[grep("Prenatal:", intronMap$Group),"AgeFrac"] = "Prenatal"
intronMap[grep("Cytosol:", intronMap$Group),"AgeFrac"] = "Cytosol"
intronMap[grep("Nucleus:", intronMap$Group),"AgeFrac"] = "Nucleus"
intronMap[grep("Cytosol-Increased", intronMap$Group),"Dir"] = "Cytosolic\nRetention"
intronMap[grep("Nucleus-Increased", intronMap$Group),"Dir"] = "Nuclear\nRetention"
intronMap[grep("Adult-Increased", intronMap$Group),"Dir"] = "Increasing\nRetention"
intronMap[grep("Prenatal-Increased", intronMap$Group),"Dir"] = "Decreasing\nRetention"
intronMap[grep("Introns", intronMap$Group),"Dir"] = "All"
head(intronMap)
intronMap$AgeFrac = factor(intronMap$AgeFrac, levels = c("Both","Adult","Prenatal","Cytosol","Nucleus"))
intronMap$Dir = factor(intronMap$Dir, levels = c("All","Nuclear\nRetention","Increasing\nRetention","Decreasing\nRetention"))



## Count the number of repeats
intronMapdt = data.table(intronMap)
Freq = list(Name = intronMapdt[,list(Count = length(unique(intronID))), by = c("repName", "Group", "Comparison", "Dir", "AgeFrac")],
            Class = intronMapdt[,list(Count = length(unique(intronID))), by = c("repClass", "Group", "Comparison", "Dir", "AgeFrac")],
            Family = intronMapdt[,list(Count = length(unique(intronID))), by = c("repFamily", "Group", "Comparison", "Dir", "AgeFrac")],
            Repeat = intronMapdt[,list(Count = length(unique(intronID))), by = c("ContainsRepeat", "Group", "Comparison", "Dir", "AgeFrac")])
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
g = ggplot(Freq[[i]][which(Freq[[i]][,"AgeFrac"]=="Adult" | Freq[[i]][,"AgeFrac"]=="Prenatal"),], 
       aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  facet_grid(. ~ AgeFrac) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
g = ggplot(Freq[[i]][which(Freq[[i]][,"AgeFrac"]=="Cytosol" | Freq[[i]][,"AgeFrac"]=="Nucleus"),], 
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
g = ggplot(Freq[[i]][which(Freq[[i]][,"Group"]=="Introns (Fraction)"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
g = ggplot(Freq[[i]][which(Freq[[i]][,"Group"]=="Introns (Age)"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Age: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
}
dev.off()

pdf(paste0(path,"repetitive_elements_sigIntrons.pdf"), width = 8, height = 6)
for (i in 1:length(Freq)) {
g = ggplot(Freq[[i]][which((Freq[[i]][,"AgeFrac"]=="Adult" | Freq[[i]][,"AgeFrac"]=="Prenatal") & 
                             Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
  facet_grid(. ~ AgeFrac) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction: ", names(Freq)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
g = ggplot(Freq[[i]][which((Freq[[i]][,"AgeFrac"]=="Cytosol" | Freq[[i]][,"AgeFrac"]=="Nucleus") & 
                             Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
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

pdf(paste0(path,"repetitive_elements_AllIntrons.pdf"), width = 9, height = 6.5)
for (i in 2:length(Freq)) {
  g = ggplot(Freq[[i]][which(Freq[[i]][,"Group"]=="Introns (Fraction)" & Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
             aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
    labs(fill="") +
    ylab("Count") + 
    xlab("") +
    ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Fraction: ", names(Freq)[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)  
  g = ggplot(Freq[[i]][which(Freq[[i]][,"Group"]=="Introns (Age)" & Freq[[i]][,1]!="No Repeats" & Freq[[i]][,1]!="No"),], 
           aes(x = Dir, y = Count, fill = get(colnames(Freq[[i]])[1]))) + geom_bar(stat = "identity") +
    labs(fill="") + 
    ylab("Count") + 
    xlab("") +
    ggtitle(paste0("Repetitive Elements Present in\nDifferentially Retained Introns by Age: ", names(Freq)[i])) +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20))
  print(g)
}
dev.off()



### Is there a difference between the proportion of repeat-containing introns in significantly/nonsignificantly retained introns?

intronMapdt$Sig = ifelse(intronMapdt$padj<=0.05, "Yes", "No")
nodirF = as.data.frame(intronMapdt[Comparison=="Fraction", length(unique(intronID)), by = c("Sig", "ContainsRepeat", "AgeFrac")])
dirF = as.data.frame(intronMapdt[Comparison=="Fraction", length(unique(intronID)), by = c("Sig", "ContainsRepeat", "Dir")])
nodirA = as.data.frame(intronMapdt[Comparison=="Age", length(unique(intronID)), by = c("Sig", "ContainsRepeat", "AgeFrac")])
dirA = as.data.frame(intronMapdt[Comparison=="Age", length(unique(intronID)), by = c("Sig", "ContainsRepeat", "Dir")])

tables = list(bothAges = data.frame(rep = c(nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="Yes" & nodirF$AgeFrac=="Both","V1"],
                                        nodirF[nodirF$Sig=="No" & nodirF$ContainsRepeat=="Yes" & nodirF$AgeFrac=="Both","V1"]),
                                norep = c(nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Both","V1"],
                                          nodirF[nodirF$Sig=="No" & nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Both","V1"]), row.names = c("sig","NS")),
              adult = data.frame(rep = c(nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="Yes" & nodirF$AgeFrac=="Adult","V1"],
                                         sum(nodirF[nodirF$ContainsRepeat=="Yes" & nodirF$AgeFrac=="Both","V1"])-
                                           nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="Yes" & nodirF$AgeFrac=="Adult","V1"]),
                                 norep = c(nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Adult","V1"],
                                           sum(nodirF[nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Both","V1"])-
                                             nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Adult","V1"]), row.names = c("sig","NS")),
              prenatal = data.frame(rep = c(0, sum(nodirF[nodirF$ContainsRepeat=="Yes" & nodirF$AgeFrac=="Both","V1"])),
                                    norep = c(nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Prenatal","V1"],
                                              sum(nodirF[nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Both","V1"])-
                                                nodirF[nodirF$Sig=="Yes" & nodirF$ContainsRepeat=="No" & nodirF$AgeFrac=="Prenatal","V1"]), row.names = c("sig","NS")),
              bothFracs = data.frame(rep = c(nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Both","V1"],
                                            nodirA[nodirA$Sig=="No" & nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Both","V1"]),
                                    norep = c(nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Both","V1"],
                                              nodirA[nodirA$Sig=="No" & nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Both","V1"]), row.names = c("sig","NS")),
              Cytosol = data.frame(rep = c(nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Cytosol","V1"],
                                         sum(nodirA[nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Both","V1"])-
                                           nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Cytosol","V1"]),
                                 norep = c(nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Cytosol","V1"],
                                           sum(nodirA[nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Both","V1"])-
                                             nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Cytosol","V1"]), row.names = c("sig","NS")),
              Nucleus = data.frame(rep = c(nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Nucleus","V1"],
                                           sum(nodirA[nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Both","V1"])-
                                             nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="Yes" & nodirA$AgeFrac=="Nucleus","V1"]),
                                   norep = c(nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Nucleus","V1"],
                                             sum(nodirA[nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Both","V1"])-
                                               nodirA[nodirA$Sig=="Yes" & nodirA$ContainsRepeat=="No" & nodirA$AgeFrac=="Nucleus","V1"]), row.names = c("sig","NS")))
fisher = lapply(tables, fisher.test)
write.csv(data.frame(rbind(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value)),unlist(lapply(lapply(tables, fisher.test), function(x) x$estimate))),
           row.names = c("pval", "OR")), quote = F, 
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonsig_introns_ContainsRepeat.csv")
names(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value))[unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value))<=(0.05/length(tables))])
#"bothAges" "adult" repeats are depleted in significantly nuclear-retained introns in adult



## In repeat-containing introns, what is the proportion of ALU in significantly/nonsignificantly retained fraction introns?

intronMapdt$isAlu = ifelse(intronMapdt$repFamily=="Alu", "Yes", "No")
Fr = as.data.frame(intronMapdt[Comparison=="Fraction", length(unique(intronID)), by = c("Sig", "isAlu", "AgeFrac")])
Ag = as.data.frame(intronMapdt[Comparison=="Age", length(unique(intronID)), by = c("Sig", "isAlu", "AgeFrac")])

tables = list(bothAges = data.frame(rep = c(Fr[Fr$Sig=="Yes" & Fr$isAlu=="Yes" & Fr$AgeFrac=="Both","V1"],
                                            Fr[Fr$Sig=="No" & Fr$isAlu=="Yes" & Fr$AgeFrac=="Both","V1"]),
                                    norep = c(Fr[Fr$Sig=="Yes" & Fr$isAlu=="No" & Fr$AgeFrac=="Both","V1"],
                                              Fr[Fr$Sig=="No" & Fr$isAlu=="No" & Fr$AgeFrac=="Both","V1"]), row.names = c("sig","NS")),
              adult = data.frame(rep = c(Fr[Fr$Sig=="Yes" & Fr$isAlu=="Yes" & Fr$AgeFrac=="Adult","V1"],
                                         sum(Fr[Fr$isAlu=="Yes" & Fr$AgeFrac=="Both","V1"])-
                                           Fr[Fr$Sig=="Yes" & Fr$isAlu=="Yes" & Fr$AgeFrac=="Adult","V1"]),
                                 norep = c(Fr[Fr$Sig=="Yes" & Fr$isAlu=="No" & Fr$AgeFrac=="Adult","V1"],
                                           sum(Fr[Fr$isAlu=="No" & Fr$AgeFrac=="Both","V1"])-
                                             Fr[Fr$Sig=="Yes" & Fr$isAlu=="No" & Fr$AgeFrac=="Adult","V1"]), row.names = c("sig","NS")),
              prenatal = data.frame(rep = c(0, sum(Fr[Fr$isAlu=="Yes" & Fr$AgeFrac=="Both","V1"])),
                                    norep = c(Fr[Fr$Sig=="Yes" & Fr$isAlu=="No" & Fr$AgeFrac=="Prenatal","V1"],
                                              sum(Fr[Fr$isAlu=="No" & Fr$AgeFrac=="Both","V1"])-
                                                Fr[Fr$Sig=="Yes" & Fr$isAlu=="No" & Fr$AgeFrac=="Prenatal","V1"]), row.names = c("sig","NS")),
              bothFracs = data.frame(rep = c(Ag[Ag$Sig=="Yes" & Ag$isAlu=="Yes" & Ag$AgeFrac=="Both","V1"],
                                             Ag[Ag$Sig=="No" & Ag$isAlu=="Yes" & Ag$AgeFrac=="Both","V1"]),
                                     norep = c(Ag[Ag$Sig=="Yes" & Ag$isAlu=="No" & Ag$AgeFrac=="Both","V1"],
                                               Ag[Ag$Sig=="No" & Ag$isAlu=="No" & Ag$AgeFrac=="Both","V1"]), row.names = c("sig","NS")),
              Cytosol = data.frame(rep = c(Ag[Ag$Sig=="Yes" & Ag$isAlu=="Yes" & Ag$AgeFrac=="Cytosol","V1"],
                                           sum(Ag[Ag$isAlu=="Yes" & Ag$AgeFrac=="Both","V1"])-
                                             Ag[Ag$Sig=="Yes" & Ag$isAlu=="Yes" & Ag$AgeFrac=="Cytosol","V1"]),
                                   norep = c(Ag[Ag$Sig=="Yes" & Ag$isAlu=="No" & Ag$AgeFrac=="Cytosol","V1"],
                                             sum(Ag[Ag$isAlu=="No" & Ag$AgeFrac=="Both","V1"])-
                                               Ag[Ag$Sig=="Yes" & Ag$isAlu=="No" & Ag$AgeFrac=="Cytosol","V1"]), row.names = c("sig","NS")),
              Nucleus = data.frame(rep = c(Ag[Ag$Sig=="Yes" & Ag$isAlu=="Yes" & Ag$AgeFrac=="Nucleus","V1"],
                                           sum(Ag[Ag$isAlu=="Yes" & Ag$AgeFrac=="Both","V1"])-
                                             Ag[Ag$Sig=="Yes" & Ag$isAlu=="Yes" & Ag$AgeFrac=="Nucleus","V1"]),
                                   norep = c(Ag[Ag$Sig=="Yes" & Ag$isAlu=="No" & Ag$AgeFrac=="Nucleus","V1"],
                                             sum(Ag[Ag$isAlu=="No" & Ag$AgeFrac=="Both","V1"])-
                                               Ag[Ag$Sig=="Yes" & Ag$isAlu=="No" & Ag$AgeFrac=="Nucleus","V1"]), row.names = c("sig","NS")))
fisher = lapply(tables, fisher.test)
write.csv(data.frame(rbind(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value)),unlist(lapply(lapply(tables, fisher.test), function(x) x$estimate))),
                     row.names = c("pval", "OR")), quote = F, 
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonsig_introns_isAlu.csv")
names(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value))[unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value))<=(0.05/length(tables))])



## There are no significant L1-containing introns
## 2 prenatal-increasing introns in nuclear RNA contain L2, including one in SHANK3


## In repeat-containing introns, what is the proportion of simple repeats in significantly/nonsignificantly retained fraction introns?

intronMapdt$isSimple_repeat = ifelse(intronMapdt$repFamily=="Simple_repeat", "Yes", "No")
Fr = as.data.frame(intronMapdt[Comparison=="Fraction", length(unique(intronID)), by = c("Sig", "isSimple_repeat", "AgeFrac")])
Ag = as.data.frame(intronMapdt[Comparison=="Age", length(unique(intronID)), by = c("Sig", "isSimple_repeat", "AgeFrac")])

tables = list(bothAges = data.frame(rep = c(Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="Yes" & Fr$AgeFrac=="Both","V1"],
                                            Fr[Fr$Sig=="No" & Fr$isSimple_repeat=="Yes" & Fr$AgeFrac=="Both","V1"]),
                                    norep = c(Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Both","V1"],
                                              Fr[Fr$Sig=="No" & Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Both","V1"]), row.names = c("sig","NS")),
              adult = data.frame(rep = c(Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="Yes" & Fr$AgeFrac=="Adult","V1"],
                                         sum(Fr[Fr$isSimple_repeat=="Yes" & Fr$AgeFrac=="Both","V1"])-
                                           Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="Yes" & Fr$AgeFrac=="Adult","V1"]),
                                 norep = c(Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Adult","V1"],
                                           sum(Fr[Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Both","V1"])-
                                             Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Adult","V1"]), row.names = c("sig","NS")),
              prenatal = data.frame(rep = c(0, sum(Fr[Fr$isSimple_repeat=="Yes" & Fr$AgeFrac=="Both","V1"])),
                                    norep = c(Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Prenatal","V1"],
                                              sum(Fr[Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Both","V1"])-
                                                Fr[Fr$Sig=="Yes" & Fr$isSimple_repeat=="No" & Fr$AgeFrac=="Prenatal","V1"]), row.names = c("sig","NS")),
              bothFracs = data.frame(rep = c(Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Both","V1"],
                                             Ag[Ag$Sig=="No" & Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Both","V1"]),
                                     norep = c(Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Both","V1"],
                                               Ag[Ag$Sig=="No" & Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Both","V1"]), row.names = c("sig","NS")),
              Cytosol = data.frame(rep = c(Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Cytosol","V1"],
                                           sum(Ag[Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Both","V1"])-
                                             Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Cytosol","V1"]),
                                   norep = c(Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Cytosol","V1"],
                                             sum(Ag[Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Both","V1"])-
                                               Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Cytosol","V1"]), row.names = c("sig","NS")),
              Nucleus = data.frame(rep = c(Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Nucleus","V1"],
                                           sum(Ag[Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Both","V1"])-
                                             Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="Yes" & Ag$AgeFrac=="Nucleus","V1"]),
                                   norep = c(Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Nucleus","V1"],
                                             sum(Ag[Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Both","V1"])-
                                               Ag[Ag$Sig=="Yes" & Ag$isSimple_repeat=="No" & Ag$AgeFrac=="Nucleus","V1"]), row.names = c("sig","NS")))
fisher = lapply(tables, fisher.test)
write.csv(data.frame(rbind(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value)),unlist(lapply(lapply(tables, fisher.test), function(x) x$estimate))),
                     row.names = c("pval", "OR")), quote = F, 
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/fisher_sig_vs_nonsig_introns_isSimple_repeat.csv")
names(unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value))[unlist(lapply(lapply(tables, fisher.test), function(x) x$p.value))<=(0.05/length(tables))])