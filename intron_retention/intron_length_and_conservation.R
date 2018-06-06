library(GenomicRanges)
library(ggplot2)
library(plyr)

### Prepare Intron Lists

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
introns[elementNROWS(introns)>0] = Map(cbind, introns[elementNROWS(introns)>0], ID = lapply(introns[elementNROWS(introns)>0], 
                                                                                            function(x) paste0(x$Chr,":",x$Start,"-",x$End,":",x$Direction)))
save(introns, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/introns_object.rda")

# get conservation information from the UCSC table browser

coord = lapply(introns, function(x) x[which(x$Chr %in% c(paste0("chr", c(1:22, "X","Y","M")))),])
coord = lapply(coord, function(x) reduce(makeGRangesFromDataFrame(x, seqnames.field = "Chr",start.field="Start", end.field="End", strand.field="Direction")))
all = as.data.frame(reduce(coord[[1]], min.gapwidth = 1))
all$ID = paste0(all$seqnames,":",all$start,"-",all$end)
all = all[,!colnames(all) %in% c("width","strand")]
dim(all) # 81367     4

for ( i in 1:102) {
  write.table(all[((i-1)*800+1):(i*800),1:3], file=paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.",i,".txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
for ( i in 1:102) {
  if (i==102) { gerp[[i]] = makeGRangesFromDataFrame(all[((i-1)*800+1):81367,1:3]) } else {
    gerp[[i]] = makeGRangesFromDataFrame(all[((i-1)*800+1):(i*800),1:3])
  }
}
names(gerp) = as.character(1:length(gerp))
lapply(gerp, function(x) sum(as.numeric(width(x)))/10000000)[lapply(gerp, function(x) sum(as.numeric(width(x)))/10000000)>1]

num = c(13,19,20,24,30,33,37,43,46,47,56,86)
for (i in 1:length(num)) {
  x = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.",num[i],".txt"))
  write.table(x[1:400,], file=paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.",num[i],"a.txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(x[401:800,], file=paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.",num[i],"b.txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}

# Go to the table browser, and under regions click "define regions" and upload the first 1000 regions in the list written above. Click submit.
# To get GERP scores for hg19, go to the "Comparative Genomics" group and "GERP" track. Change the filter 10 million lines. Label the output "repeatmasker.introns.1000.txt" and click "Get Output."
# Repeat for the remaining introns

gerp = list()
for (i in 1:102) { if (i %in% num) {
  gerp[[i]] = list(a = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,"a.txt"),
                                  comment.char = "t", col.names = c("chromosome", "start", "end", "score")),
                   b = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,"b.txt"),
                                  comment.char = "t", col.names = c("chromosome", "start", "end", "score"))) } 
  else {
    gerp[[i]] = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,".txt"),
                           comment.char = "t", col.names = c("chromosome", "start", "end", "score"))
  } }

elementNROWS(gerp)
x = unlist(gerp[elementNROWS(gerp)==2], recursive = F)
gerp = c(gerp[elementNROWS(gerp)!=2], x)
gerp = do.call(getMethod(c, "GenomicRanges"), GRangesList(lapply(gerp, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))))
rm(x)

save(gerp, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/gerp.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/gerp.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp.rda")

coord = lapply(introns, function(x) x[which(x$Chr %in% c(paste0("chr", c(1:22, "X","Y","M")))),])
coord = lapply(coord, function(x) reduce(makeGRangesFromDataFrame(x, seqnames.field = "Chr",start.field="Start", end.field="End", strand.field="Direction")))
all = as.data.frame(reduce(coord[[1]]))
allgr = makeGRangesFromDataFrame(all)
names(allgr) = paste0(all$seqnames,":",all$start,"-",all$end,":",all$strand)
save(allgr, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/gerp2.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp2.rda", verbose = T)

hits = findOverlaps(allgr, gerp)
id = split(subjectHits(hits), queryHits(hits))
y = width(allgr[unique(queryHits(hits))])
names(y) = as.character(unique(queryHits(hits)))
num = as.integer(names(y[(y-1)!=elementNROWS(id)]))

## apparently skipped some
missing = c(allgr[-unique(queryHits(hits))],allgr[num])
missing = as.data.frame(missing)
for ( i in 1:44) {
  write.table(missing[((i-1)*800+1):(i*800),1:3], file=paste0("/dcl01/lieber/ajaffe/Amanda/intron.coordinates.",i,".txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
x=list()
for ( i in 1:44) {
  x[[i]] = makeGRangesFromDataFrame(missing[((i-1)*800+1):(i*800),1:3])
}
names(x) = as.character(1:length(x))
lapply(x, function(x) sum(as.numeric(width(x)))/10000000)[lapply(x, function(x) sum(as.numeric(width(x)))/10000000)>1]

num = c(2,6,12,29,32,34,36)
for (i in 1:length(num)) {
  x = read.table(paste0("/dcl01/lieber/ajaffe/Amanda/intron.coordinates.",num[i],".txt"))
  write.table(x[1:400,], file=paste0("/dcl01/lieber/ajaffe/Amanda/intron.coordinates.",num[i],"a.txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(x[401:800,], file=paste0("/dcl01/lieber/ajaffe/Amanda/intron.coordinates.",num[i],"b.txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
# scp aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/intron.coordinates.* ../sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/

gerp = list()
for (i in 1:44) { if (i %in% num) {
  gerp[[i]] = list(a = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,"a.txt"),
                                  comment.char = "t", col.names = c("chromosome", "start", "end", "score")),
                   b = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,"b.txt"),
                                  comment.char = "t", col.names = c("chromosome", "start", "end", "score"))) } 
  else {
    gerp[[i]] = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,".txt"),
                           comment.char = "t", col.names = c("chromosome", "start", "end", "score"))
  } }

elementNROWS(gerp)
x = unlist(gerp[elementNROWS(gerp)==2], recursive = F)
gerp3 = c(gerp[elementNROWS(gerp)!=2], x)
gerp3 = do.call(getMethod(c, "GenomicRanges"), GRangesList(lapply(gerp3, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))))
rm(x)
save(gerp3, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/gerp3.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp3.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp.rda", verbose = T)
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp2.rda", verbose = T)

gerp = unique(c(gerp,gerp3))

hits = findOverlaps(allgr, gerp)
id = split(subjectHits(hits), queryHits(hits))
y = width(allgr[unique(queryHits(hits))])
names(y) = as.character(unique(queryHits(hits)))
num = as.integer(names(y[(y-1)!=elementNROWS(id)]))

## apparently skipped some
missing = c(allgr[-unique(queryHits(hits))],allgr[num])
missing = as.data.frame(missing)
for ( i in 1:5) {
  write.table(missing[((i-1)*800+1):(i*800),1:3], file=paste0("/dcl01/lieber/ajaffe/Amanda/intron.coordinates.",i,".txt"),
              quote = F, sep = "\t", row.names = F, col.names = F)
}
x=list()
for ( i in 1:5) {
  x[[i]] = makeGRangesFromDataFrame(missing[((i-1)*800+1):(i*800),1:3])
}
names(x) = as.character(1:length(x))
lapply(x, function(x) sum(as.numeric(width(x)))/10000000)[lapply(x, function(x) sum(as.numeric(width(x)))/10000000)>1]

x = read.table("/dcl01/lieber/ajaffe/Amanda/intron.coordinates.2.txt")
write.table(x[1:400,], file="/dcl01/lieber/ajaffe/Amanda/intron.coordinates.2a.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x[401:800,], file="/dcl01/lieber/ajaffe/Amanda/intron.coordinates.2b.txt", quote = F, sep = "\t", row.names = F, col.names = F)
# scp aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/intron.coordinates.* ../sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/
gerp = list()
for (i in 1:5) { if (i==2) {
  gerp[[i]] = list(a = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,"a.txt"),
                                  comment.char = "t", col.names = c("chromosome", "start", "end", "score")),
                   b = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,"b.txt"),
                                  comment.char = "t", col.names = c("chromosome", "start", "end", "score"))) } 
  else {
    gerp[[i]] = read.table(paste0("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/GERP_",i,".txt"),
                           comment.char = "t", col.names = c("chromosome", "start", "end", "score"))
  } }

elementNROWS(gerp)
x = unlist(gerp[elementNROWS(gerp)==2], recursive = F)
gerp4 = c(gerp[elementNROWS(gerp)!=2], x)
gerp4 = do.call(getMethod(c, "GenomicRanges"), GRangesList(lapply(gerp4, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))))
save(gerp4, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP/gerp4.rda")
load("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp4.rda", verbose = T)

gerp = unique(c(gerp,gerp4))

hits = findOverlaps(allgr, gerp)
id = split(subjectHits(hits), queryHits(hits))
y = width(allgr[unique(queryHits(hits))])
names(y) = as.character(unique(queryHits(hits)))
num = as.integer(names(y[(y-1)!=elementNROWS(id)]))

## apparently skipped some
missing = c(allgr[-unique(queryHits(hits))],allgr[num])
missing = as.data.frame(missing)
write.table(missing[,1:3], file="/dcl01/lieber/ajaffe/Amanda/intron.coordinates.1.txt", quote = F, sep = "\t", row.names = F, col.names = F)
# scp aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/Amanda/intron.coordinates.* Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/
# scp Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/GERP_1.txt aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/brain-epigenomics/misc/
gerp4 = read.table("/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/GERP_1.txt", comment.char = "t", col.names = c("chromosome", "start", "end", "score"))
dim(gerp4)
gerp4 = makeGRangesFromDataFrame(gerp4, keep.extra.columns = T)

gerp = unique(c(gerp,gerp4))
save(gerp, allgr, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp.rda")
load("Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/gerp.rda")

uniq = makeGRangesFromDataFrame(introns$"All Introns")
names(uniq) = introns$"All Introns"$ID
length(uniq)
uniq = unique(uniq)
u1 = uniq[1:62973]
u2 = uniq[62974:length(uniq)]

hits1 = findOverlaps(u1, gerp)
id = split(subjectHits(hits1), queryHits(hits1))
y = width(u1[unique(queryHits(hits1))])
names(y) = as.character(unique(queryHits(hits1)))
missing = u1[-unique(queryHits(hits1))]
table(abs(y-elementNROWS(id))>1)
#FALSE 
#62973
cons1 = lapply(id, function(x) gerp[x]$score)
names(cons1) = names(as.character(u1[unique(queryHits(hits1))]))
meancons = data.frame(intronID = names(cons1), mean.GERP = unlist(lapply(cons1, mean)))

hits2 = findOverlaps(u2, gerp)
id = split(subjectHits(hits2), queryHits(hits2))
y = width(u2[unique(queryHits(hits2))])
names(y) = as.character(unique(queryHits(hits2)))
missing = u2[-unique(queryHits(hits2))]
table(abs(y-elementNROWS(id))>1)
#FALSE 
#62973
cons2 = lapply(id, function(x) gerp[x]$score)
names(cons2) = names(as.character(u2[unique(queryHits(hits2))]))
meancons = rbind(meancons, data.frame(intronID = names(cons2), mean.GERP = unlist(lapply(cons2, mean))))
rownames(meancons) = NULL

save(meancons, file = "/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/meanGerp.rda")
# scp aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/meanGerp.rda Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/
# scp aprice26@jhpce01.jhsph.edu:/dcl01/lieber/ajaffe/lab/brain-epigenomics/misc/gerp.rda Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/


#################################################################
######################### intron length #########################
#################################################################

### Check length distribution of differentially retained introns

introns = Map(cbind, introns, length = lapply(introns, function(x) abs(x$End - x$Start)))
elementNROWS(introns)
length = do.call(rbind, lapply(Map(cbind, introns[elementNROWS(introns)>0], Comparison = as.list(names(introns[elementNROWS(introns)>0]))), function(x) data.frame(length = x$length, Comparison=x$Comparison)))
length$Comparison = factor(length$Comparison, levels = c("All Introns","Adult:Cytoplasm-Increased","Prenatal:Cytoplasm-Increased","Adult:Nucleus-Increased","Prenatal:Nucleus-Increased",
                              "Cytoplasm:Adult-Increased","Cytoplasm:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased"))


# Plot density distribution of intron lengths using all "clean" introns as background
pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_intron_length_byFraction_allIntrons.pdf", width = 8, height = 5)
ggplot(length[which(length$Comparison=="All Introns" | length$Comparison=="Adult:Cytosol-Increased" | 
                    length$Comparison=="Adult:Nucleus-Increased" |
                    length$Comparison=="Prenatal:Cytosol-Increased" | length$Comparison=="Prenatal:Nucleus-Increased"),],
       aes(x=length/1000)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("Density") + 
  xlab("Intron Length (Kb)") +
  ggtitle("Intron Length By Group") +
  xlim(0,10) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.7, 0.6)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_intron_length_byAge_allIntrons.pdf", width = 8, height = 5)
ggplot(length[which(length$Comparison=="All Introns" | length$Comparison=="Cytoplasm:Adult-Increased" | 
                      length$Comparison=="Cytoplasm:Prenatal-Increased" |
                      length$Comparison=="Nucleus:Adult-Increased" | length$Comparison=="Nucleus:Prenatal-Increased"),],
       aes(x=length/1000)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("Density") + 
  xlab("Intron Length (Kb)") +
  ggtitle("Intron Length By Group") +
  xlim(0,10) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.7, 0.6)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

# Measure length differences in groups of introns

t.test(c(length[length$Comparison=="Adult:Cytoplasm-Increased","length"],length[length$Comparison=="Prenatal:Cytoplasm-Increased","length"]),
       c(length[length$Comparison=="Adult:Nucleus-Increased","length"],length[length$Comparison=="Prenatal:Nucleus-Increased","length"]))
t.test(length[length$Comparison=="Adult:Cytoplasm-Increased","length"],length[length$Comparison=="Adult:Nucleus-Increased","length"])
t.test(length[length$Comparison=="Prenatal:Cytoplasm-Increased","length"],length[length$Comparison=="Prenatal:Nucleus-Increased","length"])
t.test(length[length$Comparison=="Cytoplasm:Adult-Increased","length"],length[length$Comparison=="Cytoplasm:Prenatal-Increased","length"])

ttest = list(Frac.VS.AllIntrons = t.test(c(length[length$Comparison=="Adult:Cytoplasm-Increased","length"],length[length$Comparison=="Prenatal:Cytoplasm-Increased","length"],
                                         length[length$Comparison=="Adult:Nucleus-Increased","length"],length[length$Comparison=="Prenatal:Nucleus-Increased","length"]),
                                       length[length$Comparison=="All Introns","length"]),
             byAge = t.test(c(length[length$Comparison=="Cytoplasm:Adult-Increased","length"],length[length$Comparison=="Nucleus:Adult-Increased","length"]),
                            c(length[length$Comparison=="Cytoplasm:Prenatal-Increased","length"],length[length$Comparison=="Nucleus:Prenatal-Increased","length"])),
             byAge.inNuc = t.test(length[length$Comparison=="Nucleus:Adult-Increased","length"],length[length$Comparison=="Nucleus:Prenatal-Increased","length"]),
             Age.VS.AllIntrons = t.test(c(length[length$Comparison=="Cytoplasm:Adult-Increased","length"],length[length$Comparison=="Nucleus:Adult-Increased","length"],
                                        length[length$Comparison=="Cytoplasm:Prenatal-Increased","length"],length[length$Comparison=="Nucleus:Prenatal-Increased","length"]),
                                      length[length$Comparison=="All Introns","length"]))
df = do.call(rbind, lapply(ttest, function(x) data.frame(Tstat = x$statistic, pval = x$p.value, Mean1 = x$estimate[1], Mean2 = x$estimate[2])))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote=F, file = "./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron_length_ttests.csv")


#############################################################################
######################### evolutionary conservation #########################
#############################################################################

### Check evolutionary conservation via GERP score

# Positive scores scale with the level of constraint, such that the greater the score, the greater the level of evolutionary constraint inferred to be acting on that site.

load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/meanGerp.rda")
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/introns_object.rda")

# Plot mean GERP scores using all "clean" introns as background
introns = Map(cbind, introns[elementNROWS(introns)>0], lapply(introns[elementNROWS(introns)>0], function(x) meancons[match(x$ID, meancons$intronID),]), Comparison = as.list(names(introns[elementNROWS(introns)>0])))
gerp = do.call(rbind, lapply(introns, function(x) data.frame(ID = x$ID, mean.GERP = x$mean.GERP, Comparison=x$Comparison)))
gerp$Comparison = factor(gerp$Comparison, 
                         levels = c("All Introns","Adult:Cytoplasm-Increased","Prenatal:Cytoplasm-Increased","Adult:Nucleus-Increased","Prenatal:Nucleus-Increased",
                                    "Cytoplasm:Adult-Increased","Cytoplasm:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased"))
gerp$Samples = gerp$Dir = "NA"
gerp[grep("Adult:", gerp$Comparison), "Samples"] = "In Adult"
gerp[grep("Prenatal:", gerp$Comparison), "Samples"] = "In Prenatal"
gerp[grep("Nucleus:", gerp$Comparison), "Samples"] = "In Nucleus"
gerp[grep("Cytoplasm:", gerp$Comparison), "Samples"] = "In Cytoplasm"
gerp[grep("Adult-", gerp$Comparison), "Dir"] = "Developmentally\nIncreasing"
gerp[grep("Prenatal-", gerp$Comparison), "Dir"] = "Developmentally\nDecreasing"
gerp[grep("Nucleus-", gerp$Comparison), "Dir"] = "Nuclear"
gerp[grep("Cytoplasm-", gerp$Comparison), "Dir"] = "Cytoplasmic"

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/density_GERP.pdf", width = 9, height = 5.5)
ggplot(gerp[which(gerp$Comparison %in% c("All Introns","Cytoplasm:Adult-Increased","Cytoplasm:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased")),],
       aes(x=mean.GERP)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("Density") + 
  xlab("GERP)") +
  ggtitle("GERP Scores By Group") +
  xlim(-5,5) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.8)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
ggplot(gerp[which(gerp$Comparison %in% c("All Introns","Adult:Cytoplasm-Increased","Adult:Nucleus-Increased","Prenatal:Cytoplasm-Increased","Prenatal:Nucleus-Increased")),],
       aes(x=mean.GERP)) + geom_density(aes(group=Comparison, colour=Comparison)) +
  ylab("Density") + 
  xlab("GERP)") +
  ggtitle("GERP Scores By Group") +
  xlim(-5,5) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = c(0.8, 0.8)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/gerp_introns_byFraction.pdf", width = 6, height = 5)
ggplot(gerp[which(gerp$Samples=="In Adult" | gerp$Samples=="In Prenatal"),],
       aes(x=Dir, y=mean.GERP, fill=Samples), color=Samples) + geom_boxplot() +
  xlab("") + 
  ylab("Mean GERP") +
  scale_fill_manual(values=c("red1","dodgerblue1")) +
  ggtitle("Base Conservation\nin Introns Differentially\nRetained By Fraction") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off()

pdf("./Dropbox/sorted_figures/new/github_controlled/intron_retention/figures/intron_IR_comparisons/gerp_introns_byAge.pdf", width = 7, height = 5)
ggplot(gerp[which(gerp$Samples=="In Cytoplasm" | gerp$Samples=="In Nucleus"),],
       aes(x=Dir, y=mean.GERP, fill=Samples), color=Samples) + geom_boxplot() +
  xlab("") + 
  ylab("Mean GERP") +
  scale_fill_brewer(8, palette="Dark2") +
  ggtitle("Base Conservation\nin Introns Differentially\nRetained Over Development") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  labs(fill="") +
  theme(legend.background = element_rect(fill = "transparent"), legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"))
dev.off() 

# Measure gerp differences in groups of introns

t.test(c(gerp[gerp$Comparison=="Adult:Cytoplasm-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Cytoplasm-Increased","mean.GERP"]),
       c(gerp[gerp$Comparison=="Adult:Nucleus-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Nucleus-Increased","mean.GERP"]))
t.test(gerp[gerp$Comparison=="Adult:Cytoplasm-Increased","mean.GERP"],gerp[gerp$Comparison=="Adult:Nucleus-Increased","mean.GERP"])
t.test(gerp[gerp$Comparison=="Prenatal:Cytoplasm-Increased","mean.GERP"],gerp[gerp$Comparison=="Prenatal:Nucleus-Increased","mean.GERP"])
t.test(gerp[gerp$Comparison=="Cytoplasm:Adult-Increased","mean.GERP"],gerp[gerp$Comparison=="Cytoplasm:Prenatal-Increased","mean.GERP"])

ttest = list(FracVSAllIntrons = t.test(gerp[gerp$Comparison %in% c("Adult:Cytoplasm-Increased","Prenatal:Cytoplasm-Increased","Adult:Nucleus-Increased",
                                                                   "Prenatal:Nucleus-Increased"),"mean.GERP"],gerp[gerp$Comparison=="All Introns","mean.GERP"]),
             byAge = t.test(gerp[gerp$Comparison %in% c("Cytoplasm:Adult-Increased","Nucleus:Adult-Increased"),"mean.GERP"],
                            gerp[gerp$Comparison %in% c("Cytoplasm:Prenatal-Increased","Nucleus:Prenatal-Increased"),"mean.GERP"]),
             byAge.inNuc = t.test(gerp[gerp$Comparison=="Nucleus:Adult-Increased","mean.GERP"],gerp[gerp$Comparison=="Nucleus:Prenatal-Increased","mean.GERP"]),
             AgeVSAllIntrons = t.test(gerp[gerp$Comparison %in% c("Cytoplasm:Adult-Increased","Nucleus:Adult-Increased","Cytoplasm:Prenatal-Increased",
                                                                  "Nucleus:Prenatal-Increased"),"mean.GERP"], gerp[gerp$Comparison=="All Introns","mean.GERP"]))
df = do.call(rbind, lapply(ttest, function(x) data.frame(Tstat = x$statistic, pval = x$p.value, Mean1 = x$estimate[1], Mean2 = x$estimate[2])))
df$FDR = p.adjust(df$pval, method = "fdr")
df
#                      Tstat        pval       Mean1      Mean2         FDR
#FracVSAllIntrons -3.3757377 0.001648091 -0.53314506 -0.2342608 0.006592365   # introns regulated by fraction are less conserved than the larger group
#byAge             0.3167434 0.770144281 -0.05479808 -0.1864130 0.770144281
#byAge.inNuc       0.3591365 0.767266018 -0.12014480 -0.3611805 0.770144281
#AgeVSAllIntrons   0.3322658 0.741998957 -0.17367606 -0.2342608 0.770144281



# More info about GERP can be found at http://mendel.stanford.edu/SidowLab/downloads/gerp/
