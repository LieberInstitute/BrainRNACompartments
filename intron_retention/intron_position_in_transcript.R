library(GenomicRanges)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")

## read in Differential IR results files
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


## Map introns to transcripts

head(exonMap)
length(exonMap[grep(";", exonMap$gencodeTx),])
exonMap$exonID = rownames(exonMap)
x = strsplit(exonMap$gencodeTx, ";", fixed=T)
names(x) = paste0(rownames(exonMap),".")
x = unlist(x)
x = data.frame(exonID = names(x), TxID = x)
x$exonID = gsub("\\..*","",x$exonID)
exonMap = cbind(x, exonMap[match(x$exonID, exonMap$exonID),])
lastinTx = ddply(exonMap, .(TxID), summarise, LastStart=max(Start), FirstStart=min(Start), TxEnd=max(End))
lastinTx$matchID = paste0(lastinTx$TxID, ":", lastinTx$LastStart)
exonMap$matchID = paste0(exonMap$TxID, ":", exonMap$Start)
exonMap$lastinTx = ifelse(exonMap$matchID %in% lastinTx$matchID, "Last", "NotLast")
exonMap = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
introns = lapply(introns, function(x) makeGRangesFromDataFrame(na.omit(x), 
                                                               seqnames.field="Chr",
                                                               start.field="Start",
                                                               end.field="End",
                                                               strand.field="Direction",
                                                               keep.extra.columns = T))
introns = c(lapply(introns[c("Introns (Fraction)","Introns (Age)")],unique), 
            introns[c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                      "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased")])
overlaps = lapply(introns, function(x) findOverlaps(exonMap, x))
introns = lapply(introns, function(x) as.data.frame(x))
exonMap = as.data.frame(exonMap)
colnames(exonMap) = c("seqnames","exonStart","exonEnd","width","strand","exonID","TxID","Length","gencodeID","ensemblID",
                      "gene_type","Symbol","EntrezID","Class","meanExprs","NumTx","gencodeTx","exonID.1","matchID","lastinTx")

intronMap = intronTx = Tx = list()
for (i in 1:length(introns)){
  ov = overlaps[[i]]
  int = introns[[i]]
  intronMap[[i]] = cbind(int[subjectHits(ov),],exonMap[queryHits(ov),])
  map = intronMap[[i]]
  map$afterExon = ifelse((map$start==map$exonEnd), "YES","NO")
  map$beforeExon = ifelse((map$end==map$exonStart), "YES","NO")
  map = map[map$afterExon=="YES" | map$beforeExon=="YES",]
  intronTx[[i]] = lastinTx[match(map$TxID, lastinTx$TxID),]
  tx = intronTx[[i]]
  Tx[[i]] = cbind(tx, map[match(tx$TxID, map$TxID),])
  intronMap[[i]] = map
}
names(intronMap) = names(intronTx) = names(Tx) = names(introns)



## Check intron position in the transcript
# calculate the distance of midpoint of intron from the end of the transcript divided by the length of the transcript
position = lapply(Tx, function(x) (x$TxEnd-(x$start + x$width/2))/(x$TxEnd-x$FirstStart)) # the smaller the value, the closer to the 3' end
lapply(position, head)

# summary statistics on position
stats = data.frame(mean = unlist(lapply(position, mean)), median = unlist(lapply(position, median)), SD = unlist(lapply(position, sd)))
stats$group = c("All", "All", "Cytosolic Retention\nIn Adult","Nuclear Retention\nIn Adult","Cytosolic Retention\nIn Prenatal","Nuclear Retention\nIn Prenatal",
                "Adult Retention\nIn Cytosol","Prenatal Retention\nIn Cytosol","Adult Retention\nIn Nucleus","Prenatal Retention\nIn Nucleus")
stats$retention = c("NA","NA","Cytosolic","Nuclear","Cytosolic","Nuclear","Adult","Prenatal","Adult","Prenatal")


# distribution of introns measured
posDF = Map(cbind, position = lapply(position, function(x) if (length(x)>0) { data.frame(x) } else { data.frame("no") }), Comparison = as.list(names(position)))
for (i in 1:length(posDF)) { colnames(posDF[[i]]) = c("Position", "Comparison") }
posDF = do.call(rbind, posDF)
posDF = posDF[posDF$Position!="no",]
posDF$Comparison = factor(posDF$Comparison, 
                          levels = c("Introns (Fraction)","Introns (Age)","Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased",
                                     "Prenatal:Nucleus-Increased","Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased"))
posDF$FracAge = NA
posDF[c(grep("Nucleus-Increased", posDF$Comparison),grep("Cytosol-Increased", posDF$Comparison),grep("Fraction", posDF$Comparison)),"FracAge"] = "Fraction"
posDF[c(grep("Adult-Increased", posDF$Comparison),grep("Prenatal-Increased", posDF$Comparison),grep("Age", posDF$Comparison)),"FracAge"] = "Age"
posDF$Position = as.numeric(posDF$Position)

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_proportionalDistance_fromIntron_toTxEnd_fraction.pdf", width=8.5,height=5)
ggplot(posDF[which(posDF$FracAge=="Fraction"),], aes(x=Position)) +
    geom_density(aes(group=Comparison, colour=Comparison)) +
    ylab("") +
    xlim(0,1) +
    xlab("Proportion of Transcript") +
    ggtitle("Distance from Transcript End") +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_proportionalDistance_fromIntron_toTxEnd_age.pdf", width=8.5,height=5)
ggplot(posDF[which(posDF$FracAge=="Age"),], aes(x=Position)) +
  geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("") +
  xlim(0,1) +
  xlab("Proportion of Transcript") +
  ggtitle("Distance from Transcript End") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

## Is there a difference by fraction

byFraction = t.test(x=c(position$'Adult:Cytosol-Increased', position$'Prenatal:Cytosol-Increased'), 
                    y=c(position$'Adult:Nucleus-Increased', position$'Prenatal:Nucleus-Increased'))
byFrac.InAdults = t.test(x=position$'Adult:Cytosol-Increased', y=position$'Adult:Nucleus-Increased')
byFrac.InPrenatal = t.test(x=position$'Prenatal:Cytosol-Increased', y=position$'Prenatal:Nucleus-Increased')
upCyt.VS.allFrac = t.test(x=c(position$'Adult:Cytosol-Increased', position$'Prenatal:Cytosol-Increased'), 
                          y=position$'Introns (Fraction)')
upCyt.inAdult.Vs.allFrac = t.test(x=position$'Adult:Cytosol-Increased', y=position$'Introns (Fraction)')
upCyt.inPrenatal.Vs.allFrac = t.test(x=position$'Prenatal:Cytosol-Increased', y=position$'Introns (Fraction)')

ttests = list(IRfrac.VS.allFrac = t.test(x=c(position$'Adult:Cytosol-Increased', position$'Prenatal:Cytosol-Increased',
                                             position$'Adult:Nucleus-Increased', position$'Prenatal:Nucleus-Increased'), 
                                         y=position$'Introns (Fraction)'),
              upNuc.VS.allFrac = t.test(x=c(position$'Adult:Nucleus-Increased', position$'Prenatal:Nucleus-Increased'), 
                                        y=position$'Introns (Fraction)'),
              upNuc.inAdult.Vs.allFrac = t.test(x=position$'Adult:Nucleus-Increased', y=position$'Introns (Fraction)'),
              upNuc.inPrenatal.Vs.allFrac = t.test(x=position$'Prenatal:Nucleus-Increased', y=position$'Introns (Fraction)'),
              byAge = t.test(x=c(position$'Cytosol:Adult-Increased', position$'Nucleus:Adult-Increased'), 
                             y=c(position$'Cytosol:Prenatal-Increased', position$'Nucleus:Prenatal-Increased')),
              byAge.inCytosol = t.test(x=position$'Cytosol:Adult-Increased', y=position$'Cytosol:Prenatal-Increased'),
              byAge.inNucleus = t.test(x=position$'Nucleus:Adult-Increased', y=position$'Nucleus:Prenatal-Increased'),
              IRage.VS.allAge = t.test(x=c(position$'Cytosol:Adult-Increased', position$'Nucleus:Adult-Increased',
                                           position$'Cytosol:Prenatal-Increased', position$'Nucleus:Prenatal-Increased'), 
                                       y=position$'Introns (Age)'),
              upAd.inCytosol.Vs.allAge = t.test(x=position$'Cytosol:Adult-Increased', y=position$'Introns (Age)'),
              upAd.inNucleus.Vs.allAge = t.test(x=position$'Nucleus:Adult-Increased', y=position$'Introns (Age)'),
              upPren.inCytosol.Vs.allAge = t.test(x=position$'Cytosol:Prenatal-Increased', y=position$'Introns (Age)'),
              upPren.inNucleus.Vs.allAge = t.test(x=position$'Nucleus:Prenatal-Increased', y=position$'Introns (Age)'))

ttests = data.frame(rbind(unlist(lapply(ttests, function(x) x$statistic)),unlist(lapply(ttests, function(x) x$p.value)),data.frame(lapply(ttests, function(x) x$conf.int)),
                          data.frame(lapply(ttests, function(x) x$estimate))), row.names=c("Tstat","p.value","conf.int1","conf.int2","mean of x","mean of y"))
write.csv(ttests, quote = F,
          file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ttest_intronPosition_inTx.csv")
colnames(ttests[,which(ttests["p.value",]<=(0.05/18))])
#"byAge" "byAge.inNucleus" "upAd.inNucleus.Vs.allAge"
# adult-increased are more 3' than prenatal-increased, and nuclear Adult-enriched introns are more 3' than all the general introns.

