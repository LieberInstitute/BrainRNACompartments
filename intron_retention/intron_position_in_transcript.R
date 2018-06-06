library(GenomicRanges)
library(plyr)
library(ggplot2)

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

for (i in 1:length(full)) { colnames(full[[i]])[1:7] = c("Chr","Start","End","Intron.GeneName.GeneID","X.","Direction","ExcludedBases") }
full = lapply(full, function(x) data.frame(x[,1:7], "p.diff"=NA,"p.increased"=NA,"p.decreased"=NA,"A.IRratio"=NA,"A.warnings"=NA,"A.IntronCover"=NA,"A.IntronDepth"=NA,"A.SplicesMax"=NA,
                                           "A.SplicesExact"=NA,"B.IRratio"=NA,"B.warnings"=NA,"B.IntronCover"=NA,"B.IntronDepth"=NA,"B.SplicesMax"=NA,"B.SplicesExact"=NA,"replicates"=NA,
                                           "A1.IRratio"=NA,"A2.IRratio"=NA,"A3.IRratio"=NA,"B1.IRratio"=NA,"B2.IRratio"=NA,"B3.IRratio"=NA))
full = Map(cbind, full, intronID = lapply(full, function(x) paste0(x$Intron.GeneName.GeneID,"/",x$Chr,":",x$Start,"-",x$End,":",x$Direction)),
           ensID = lapply(lapply(full, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE)), function(y) y[grep("ENSG", y)]),
           IR.diff = NA, Sign = NA, padj=1)
IRcomp$Adult_byFraction = rbind(IRcomp$Adult_byFraction, full$adult[-which(full$adult$intronID %in% IRcomp$Adult_byFraction$intronID),])
IRcomp$Fetal_byFraction = rbind(IRcomp$Fetal_byFraction, full$prenatal[-which(full$prenatal$intronID %in% IRcomp$Fetal_byFraction$intronID),])
IRcomp$Cytosol_byAge = rbind(IRcomp$Cytosol_byAge, full$cytosol[-which(full$cytosol$intronID %in% IRcomp$Cytosol_byAge$intronID),])
IRcomp$Nuclear_byAge = rbind(IRcomp$Nuclear_byAge, full$nucleus[-which(full$nucleus$intronID %in% IRcomp$Nuclear_byAge$intronID),])


sigdIR = lapply(IRcomp, function(x) x[which(x$padj<=0.05),])
sigdIR = unlist(lapply(sigdIR, function(x) split(x, x$Sign)), recursive=F)
names(sigdIR) = c("Adult:Cytoplasm-Increased","Adult:Nucleus-Increased","Prenatal:Cytoplasm-Increased","Prenatal:Nucleus-Increased",
                  "Cytoplasm:Adult-Increased","Cytoplasm:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased")
introns = c(list("All Introns" = do.call(rbind, IRcomp)), sigdIR)
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
introns = lapply(introns, function(x) reduce(makeGRangesFromDataFrame(x, seqnames.field="Chr", start.field="Start", end.field="End", strand.field="Direction")))
overlaps = lapply(introns, function(x) findOverlaps(exonMap, x))
introns = lapply(introns, as.data.frame)
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
stats$group = c("All", "Cytoplasmic\nIn Adult","Nuclear\nIn Adult","Cytoplasmic\nIn Prenatal","Nuclear\nIn Prenatal",
                "Adult\nIn Cytoplasm","Prenatal\nIn Cytoplasm","Adult\nIn Nucleus","Prenatal\nIn Nucleus")
stats$retention = c("NA","Cytoplasmic","Nuclear","Cytoplasmic","Nuclear","Adult","Prenatal","Adult","Prenatal")


# distribution of introns measured
posDF = Map(cbind, position = lapply(position, function(x) if (length(x)>0) { data.frame(x) } else { data.frame("no") }), Comparison = as.list(names(position)))
for (i in 1:length(posDF)) { colnames(posDF[[i]]) = c("Position", "Comparison") }
posDF = do.call(rbind, posDF)
posDF = posDF[posDF$Position!="no",]
posDF$Comparison = factor(posDF$Comparison, 
                          levels = c("All Introns","Adult:Cytoplasm-Increased","Adult:Nucleus-Increased","Prenatal:Cytoplasm-Increased","Prenatal:Nucleus-Increased",
                                     "Cytoplasm:Adult-Increased","Cytoplasm:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased"))
posDF$FracAge = NA
posDF[c(grep("Nucleus-Increased", posDF$Comparison),grep("Cytoplasm-Increased", posDF$Comparison),grep("Fraction", posDF$Comparison)),"FracAge"] = "Fraction"
posDF[c(grep("Adult-Increased", posDF$Comparison),grep("Prenatal-Increased", posDF$Comparison),grep("Age", posDF$Comparison)),"FracAge"] = "Age"
posDF$Position = as.numeric(posDF$Position)

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_proportionalDistance_fromIntron_toTxEnd_fraction.pdf", width=7,height=5)
ggplot(posDF[which(posDF$FracAge=="Fraction"),], aes(x=Position)) +
    geom_density(aes(group=Comparison, colour=Comparison)) +
    ylab("") +
    xlim(0,1) +
    xlab("Proportion of Transcript") + ylab("Density") +
    ggtitle("Distance from Transcript End") +
    theme(title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    labs(fill="") +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", color = "transparent"), legend.position = "bottom", legend.title = element_blank())
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_proportionalDistance_fromIntron_toTxEnd_age.pdf", width=8.5,height=5)
ggplot(posDF[which(posDF$FracAge=="Age"),], aes(x=Position)) +
  geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("") +
  xlim(0,1) +
  xlab("Proportion of Transcript") + ylab("Density") +
  ggtitle("Distance from Transcript End") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", color = "transparent"),legend.title = element_blank())
dev.off()

## Is there a difference by fraction

byFraction = t.test(x=c(position$'Adult:Cytoplasm-Increased', position$'Prenatal:Cytoplasm-Increased'), 
                    y=c(position$'Adult:Nucleus-Increased', position$'Prenatal:Nucleus-Increased'))
byFrac.InAdults = t.test(x=position$'Adult:Cytoplasm-Increased', y=position$'Adult:Nucleus-Increased')
byFrac.InPrenatal = t.test(x=position$'Prenatal:Cytoplasm-Increased', y=position$'Prenatal:Nucleus-Increased')
upCyt.VS.allFrac = t.test(x=c(position$'Adult:Cytoplasm-Increased', position$'Prenatal:Cytoplasm-Increased'), 
                          y=position$'All Introns')
upCyt.inAdult.Vs.allFrac = t.test(x=position$'Adult:Cytoplasm-Increased', y=position$'All Introns')
upCyt.inPrenatal.Vs.allFrac = t.test(x=position$'Prenatal:Cytoplasm-Increased', y=position$'All Introns')

ttests = list(IRfrac.VS.all = t.test(x=c(position$'Adult:Cytoplasm-Increased', position$'Prenatal:Cytoplasm-Increased',
                                         position$'Adult:Nucleus-Increased', position$'Prenatal:Nucleus-Increased'), 
                                     y=position$'All Introns'),
              upNuc.VS.all = t.test(x=c(position$'Adult:Nucleus-Increased', position$'Prenatal:Nucleus-Increased'), 
                                    y=position$'All Introns'),
              upNuc.inAdult.Vs.all = t.test(x=position$'Adult:Nucleus-Increased', y=position$'All Introns'),
              upNuc.inPrenatal.Vs.all = t.test(x=position$'Prenatal:Nucleus-Increased', y=position$'All Introns'),
              byAge = t.test(x=c(position$'Cytoplasm:Adult-Increased', position$'Nucleus:Adult-Increased'), 
                             y=c(position$'Cytoplasm:Prenatal-Increased', position$'Nucleus:Prenatal-Increased')),
              byAge.inCytoplasm = t.test(x=position$'Cytoplasm:Adult-Increased', y=position$'Cytoplasm:Prenatal-Increased'),
              byAge.inNucleus = t.test(x=position$'Nucleus:Adult-Increased', y=position$'Nucleus:Prenatal-Increased'),
              IRage.VS.all = t.test(x=c(position$'Cytoplasm:Adult-Increased', position$'Nucleus:Adult-Increased',
                                        position$'Cytoplasm:Prenatal-Increased', position$'Nucleus:Prenatal-Increased'), 
                                    y=position$'All Introns'),
              upAd.inCytoplasm.Vs.all = t.test(x=position$'Cytoplasm:Adult-Increased', y=position$'All Introns'),
              upAd.inNucleus.Vs.all = t.test(x=position$'Nucleus:Adult-Increased', y=position$'All Introns'),
              upPren.inCytoplasm.Vs.all = t.test(x=position$'Cytoplasm:Prenatal-Increased', y=position$'All Introns'),
              upPren.inNucleus.Vs.all = t.test(x=position$'Nucleus:Prenatal-Increased', y=position$'All Introns'))

df = data.frame(Comp = names(ttests), tstat = unlist(lapply(ttests, function(x) x$statistic)),pval = unlist(lapply(ttests, function(x) x$p.value)),
                ConfInt1 = unlist(lapply(ttests, function(x) x$conf.int[1])), ConfInt2 = unlist(lapply(ttests, function(x) x$conf.int[2])),
                Mean1 = unlist(lapply(ttests, function(x) x$estimate[1])), Mean2 = unlist(lapply(ttests, function(x) x$estimate[2])), row.names = NULL)
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote=F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/ttest_intronPosition_inTx.csv")
df[df$FDR<=0.05,]


