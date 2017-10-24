library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(VennDiagram)
library(plyr)
library(clusterProfiler)
require("org.Hs.eg.db")

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/rna_editing_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("/Users/amanda/Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Explore results
lapply(editingres, head)
editingres = editingres[c(1:6,8,10:14)]
names(editingres) = gsub(".downsampled", "", names(editingres))
editingres = Map(cbind, editingres, conversion = lapply(editingres, function(x) paste0(x$ref, ":", x$alt)),
                 sampleID = lapply(names(editingres), function(x) x), rnum = lapply(editingres, function(x) 1:nrow(x)),
                 editingID = lapply(editingres, function(x) paste0(x$chromosome,":",x$start,"-",x$end,":",x$ref,":",x$alt)))
for (i in 1:length(editingres)){
  tmp = editingres[[i]]
  tmp$Fraction = ifelse((tmp$rnum %in% grep("N", tmp$sampleID)), "Nucleus", "Cytosol")
  tmp$Age = ifelse((tmp$rnum %in% grep("53", tmp$sampleID)), "Prenatal", "Adult")
  tmp$Group = paste0(tmp$Age, ":", tmp$Fraction)
  tmp$collapsedconversion = NA
  tmp[which(tmp$conversion=="A:C" | tmp$conversion=="T:G"), "collapsedconversion"] = "A:C / T:G"
  tmp[which(tmp$conversion=="A:G" | tmp$conversion=="T:C"), "collapsedconversion"] = "A:G / T:C"
  tmp[which(tmp$conversion=="A:T" | tmp$conversion=="T:A"), "collapsedconversion"] = "A:T / T:A"
  tmp[which(tmp$conversion=="C:A" | tmp$conversion=="G:T"), "collapsedconversion"] = "C:A / G:T"
  tmp[which(tmp$conversion=="C:G" | tmp$conversion=="G:C"), "collapsedconversion"] = "C:G / G:C"
  tmp[which(tmp$conversion=="C:T" | tmp$conversion=="G:A"), "collapsedconversion"] = "C:T / G:A"
  editingres[[i]] = tmp
}
editingres_df = do.call(rbind, editingres)
editingres_df$rate = editingres_df$alt.count / editingres_df$valdepth
editingres_dt = data.table(editingres_df)

# Annotate editing sites to features in the genome
txdb = loadDb("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))

grediting = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grediting, y))

grediting$rnum = 1:length(grediting)
grediting$cds = ifelse(grediting$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grediting$intron = ifelse(grediting$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grediting$UTR5 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR5"]]), "UTR5", NA)
grediting$UTR3 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR3"]]), "UTR3", NA)
grediting$anno = paste0(grediting$cds,":",grediting$intron, ":", grediting$UTR5, ":", grediting$UTR3)

editing = as.data.frame(grediting)
editing[which(editing$anno == "NA:NA:NA:NA"),"annotation"] = "Other" 
editing[grep("CDS", editing$cds),"annotation"] = "CDS"
editing[which(is.na(editing$annotation) & editing$UTR3 == "UTR3"),"annotation"] = "UTR3"
editing[which(is.na(editing$annotation) & editing$UTR5 == "UTR5"),"annotation"] = "UTR5"
editing[which(is.na(editing$annotation) & editing$intron == "Intron"),"annotation"] = "Intron"

# Mapping editing sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grediting, geneMapGR)
editing$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
editing$nearestID = names(geneMapGR)[subjectHits(dA)]
editing$distToGene = mcols(dA)$distance
editing$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
editing_anno = data.table(editing)

## Limit the editing sites to those present in one group but not another

AGonly = editing_anno[collapsedconversion=="A:G / T:C",,]
cyt = AGonly[Fraction=="Cytosol",,]
nuc = AGonly[Fraction=="Nucleus",,]
ad = AGonly[Age=="Adult",,]
pren = AGonly[Age=="Prenatal",,]
AC = AGonly[Group=="Adult:Cytosol",,]
AN = AGonly[Group=="Adult:Nucleus",,]
PC = AGonly[Group=="Prenatal:Cytosol",,]
PN = AGonly[Group=="Prenatal:Nucleus",,]

unique = list(cytosolOnly = cyt[!(editingID %in% nuc$editingID),,], nucleusOnly = nuc[!(editingID %in% cyt$editingID),,], 
              adultOnly = ad[!(editingID %in% pren$editingID),,], prenatalOnly = pren[!(editingID %in% ad$editingID),,], 
              ANnotAC = AN[!(editingID %in% AC$editingID),,], ACnotAN = AC[!(editingID %in% AN$editingID),,], 
              ANnotPN = AN[!(editingID %in% PN$editingID),,], PNnotAN = PN[!(editingID %in% AN$editingID),,],
              ACnotPC = AC[!(editingID %in% PC$editingID),,], PCnotAC = PC[!(editingID %in% AC$editingID),,], 
              PCnotPN = PC[!(editingID %in% PN$editingID),,], PNnotPC = PN[!(editingID %in% PC$editingID),,])
uniqueID = lapply(unique, function(x) as.data.frame(x[,list(unique(editingID)),]))
unique = lapply(unique, function(x) as.data.frame(x))
elementNROWS(uniqueID)
all = list(cytosolAll = cyt, nucleusAll = nuc, adultAll = ad, prenatalAll = pren, allAC = AC, allAN = AN, allPC = PC, allPN = PN)

uniqueID = lapply(uniqueID, function(x) as.character(x$V1))
unique_bySamp = Map(cbind, EditingID = uniqueID, Br1113C1=NA, Br1113N1=NA, Br2046C=NA, Br2046N=NA, Br2074C=NA, Br2074N=NA, 
                                                 Br5339N1=NA, Br5340N1=NA, Br5341C1=NA, Br5341N1=NA, Br5339C1=NA, Br5340C1=NA)
for (i in 1:length(unique)){
  for (j in 1:length(uniqueID[[i]])){
  unique_bySamp[[i]][j,"Br1113C1"] = ifelse("Br1113C1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br1113C1", "no")
  unique_bySamp[[i]][j,"Br1113N1"] = ifelse("Br1113N1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br1113N1", "no")
  unique_bySamp[[i]][j,"Br2046C"] = ifelse("Br2046C" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br2046C", "no")
  unique_bySamp[[i]][j,"Br2046N"] = ifelse("Br2046N" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br2046N", "no")
  unique_bySamp[[i]][j,"Br2074C"] = ifelse("Br2074C" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br2074C", "no")
  unique_bySamp[[i]][j,"Br2074N"] = ifelse("Br2074N" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br2074N", "no")
  unique_bySamp[[i]][j,"Br5339N1"] = ifelse("Br5339N1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br5339N1", "no")
  unique_bySamp[[i]][j,"Br5340N1"] = ifelse("Br5340N1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br5340N1", "no")
  unique_bySamp[[i]][j,"Br5341C1"] = ifelse("Br5341C1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br5341C1", "no")
  unique_bySamp[[i]][j,"Br5341N1"] = ifelse("Br5341N1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br5341N1", "no")
  unique_bySamp[[i]][j,"Br5339C1"] = ifelse("Br5339C1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br5339C1", "no")
  unique_bySamp[[i]][j,"Br5340C1"] = ifelse("Br5340C1" %in% unique[[i]][which(unique[[i]][,"editingID"]==uniqueID[[i]][j]),"sampleID"], "Br5340C1", "no")
  }
}

unique_bySamp = lapply(unique_bySamp, as.data.frame)
unique_bySamp = Map(cbind, unique_bySamp, cytosolAll = lapply(unique_bySamp, function(x) paste(x$Br1113C1,x$Br2046C,x$Br2074C,x$Br5339C1,x$Br5340C1,x$Br5341C1, sep=":")),
                    nucleusAll = lapply(unique_bySamp, function(x) paste(x$Br1113N1,x$Br2046N,x$Br2074N,x$Br5339N1,x$Br5340N1,x$Br5341N1, sep=":")),
                    adultAll = lapply(unique_bySamp, function(x) paste(x$Br1113C1,x$Br2046C,x$Br2074C,x$Br1113N1,x$Br2046N,x$Br2074N, sep=":")),
                    prenatalAll = lapply(unique_bySamp, function(x) paste(x$Br5339C1,x$Br5340C1,x$Br5341C1,x$Br5339N1,x$Br5340N1,x$Br5341N1, sep=":")),
                    allAC = lapply(unique_bySamp, function(x) paste(x$Br1113C1,x$Br2046C,x$Br2074C, sep=":")),
                    allAN = lapply(unique_bySamp, function(x) paste(x$Br1113N1,x$Br2046N,x$Br2074N, sep=":")),
                    allPC = lapply(unique_bySamp, function(x) paste(x$Br5339C1,x$Br5340C1,x$Br5341C1, sep=":")),
                    allPN = lapply(unique_bySamp, function(x) paste(x$Br5339N1,x$Br5340N1,x$Br5341N1, sep=":")))
lapply(unique_bySamp, head)

unique_bySamp = Map(cbind, unique_bySamp, 
                    nos = list(cytosolOnly = rowSums(as.matrix(unique_bySamp$cytosolOnly[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1")])=="no"),
                    nucleusOnly = rowSums(as.matrix(unique_bySamp$nucleusOnly[,c("Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"),
                    adultOnly = rowSums(as.matrix(unique_bySamp$adultOnly[,c("Br1113C1","Br2046C","Br2074C","Br1113N1","Br2046N","Br2074N")])=="no"),
                    prenatalOnly = rowSums(as.matrix(unique_bySamp$prenatalOnly[,c("Br5339N1","Br5340N1","Br5341N1","Br5339C1","Br5340C1","Br5341C1")])=="no"),
                    ANnotAC = NA, ACnotAN = NA, ANnotPN = NA, PNnotAN = NA, ACnotPC = NA, PCnotAC = NA, PCnotPN = NA, PNnotPC = NA),
                    AllNos = list(cytosolOnly = rowSums(as.matrix(unique_bySamp$cytosolOnly[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                               "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"),
                                  nucleusOnly = rowSums(as.matrix(unique_bySamp$nucleusOnly[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                               "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"),
                                  adultOnly = rowSums(as.matrix(unique_bySamp$adultOnly[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                           "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"),
                                  prenatalOnly = rowSums(as.matrix(unique_bySamp$prenatalOnly[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                                 "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"),
                                  ANnotAC = rowSums(as.matrix(unique_bySamp$ANnotAC[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  ACnotAN = rowSums(as.matrix(unique_bySamp$ACnotAN[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  ANnotPN = rowSums(as.matrix(unique_bySamp$ANnotPN[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  PNnotAN = rowSums(as.matrix(unique_bySamp$PNnotAN[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  ACnotPC = rowSums(as.matrix(unique_bySamp$ACnotPC[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  PCnotAC = rowSums(as.matrix(unique_bySamp$PCnotAC[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  PCnotPN = rowSums(as.matrix(unique_bySamp$PCnotPN[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"), 
                                  PNnotPC = rowSums(as.matrix(unique_bySamp$PNnotPC[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1",
                                                                                       "Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no")))
lapply(unique_bySamp, head)

save(unique_bySamp, all, editing_anno, geneMap, file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")

elementNROWS(lapply(unique_bySamp[1:4], function(x) x[which(x$nos<3),]))
#cytosolOnly  nucleusOnly    adultOnly prenatalOnly 
#9           59          619          384
elementNROWS(lapply(unique_bySamp[1:4], function(x) x[which(x$nos<2),]))
#cytosolOnly  nucleusOnly    adultOnly prenatalOnly 
#2           11          311          208

unique_bySamp_all = list(cytosolOnly = unique_bySamp$cytosolOnly[-grep("no", unique_bySamp$cytosolOnly$cytosolAll),],
                         nucleusOnly = unique_bySamp$nucleusOnly[-grep("no", unique_bySamp$nucleusOnly$nucleusAll),], 
                         adultOnly = unique_bySamp$adultOnly[-grep("no", unique_bySamp$adultOnly$adultAll),], 
                         prenatalOnly = unique_bySamp$prenatalOnly[-grep("no", unique_bySamp$prenatalOnly$prenatalAll),], 
                         ANnotAC = unique_bySamp$ANnotAC[-grep("no", unique_bySamp$ANnotAC$allAN),], 
                         ACnotAN = unique_bySamp$ACnotAN[-grep("no", unique_bySamp$ACnotAN$allAC),], 
                         ANnotPN = unique_bySamp$ANnotPN[-grep("no", unique_bySamp$ANnotPN$allAN),], 
                         PNnotAN = unique_bySamp$PNnotAN[-grep("no", unique_bySamp$PNnotAN$allPN),],
                         ACnotPC = unique_bySamp$ACnotPC[-grep("no", unique_bySamp$ACnotPC$allAC),], 
                         PCnotAC = unique_bySamp$PCnotAC[-grep("no", unique_bySamp$PCnotAC$allPC),], 
                         PCnotPN = unique_bySamp$PCnotPN[-grep("no", unique_bySamp$PCnotPN$allPC),], 
                         PNnotPC = unique_bySamp$PNnotPC[-grep("no", unique_bySamp$PNnotPC$allPN),])
elementNROWS(unique_bySamp_all)

unique_bySamp_3more = c(lapply(unique_bySamp[1:4], function(x) as.character(x[which(x$nos<3),"EditingID"])), 
                        lapply(unique_bySamp_all[5:12], function(x) as.character(x$EditingID)))

unique_all = lapply(unique_bySamp_all, function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_3more = lapply(unique_bySamp_3more, function(x) editing_anno[which(editingID %in% x),,])
elementNROWS(unique_all)


## What is the distribution of features edited across groups?

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Genomic_features_uniqueEditingSites_inAllSamps.pdf")
for (i in 3:length(unique_all)){
tmp = unique_all[[i]]
g = ggplot(tmp[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("annotation", "Fraction", "Age")], 
       aes(x = Fraction, y = V1, fill = annotation)) + geom_bar(stat = "identity") +
  facet_grid(. ~ Age) +
  labs(fill="") +
  ylab("Count") + 
  xlab("") +
  ggtitle(paste0("Number of RNA Editing Sites\nby Feature and Group: ", names(unique_all)[i])) +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
print(g)
}
dev.off()

# Check coverage by unique editing group
for (i in 3:length(unique_all)){  unique_all[[i]] = cbind(unique_all[[i]], Status = names(unique_all)[i])  }
unique_alldt = do.call(rbind, unique_all[3:12])
unique_alldt$Status = factor(unique_alldt$Status, levels = c("cytosolOnly", "nucleusOnly", "adultOnly", "prenatalOnly", "ANnotAC",
                                                 "ACnotAN", "ANnotPN", "PNnotAN", "ACnotPC", "PCnotAC", "PCnotPN", "PNnotPC"))

pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/valdepth_unique_editing_sites_byGroup_bySample_inAllSamps.pdf", width = 28, height = 7)
ggplot(unique_alldt, aes(x = sampleID, y = valdepth, fill = Group)) + geom_boxplot() +
  facet_grid(. ~ Status, scales = "free") +
  labs(fill="") +
  ylab("High Quality Read Depth") + ylim(0,30) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  ggtitle("Validated Coverage Range By Group At Uniquely Edited Sites") +
  theme(title = element_text(size = 20)) +
  theme(text = element_text(size = 20))
dev.off()


### Are editing sites in a specific annotation more likely to be unique to a group?

x = lapply(unique_all[3:12], function(x) x[,length(unique(editingID)),by="annotation"])
x = lapply(x, function(x) as.data.frame(x))
x = lapply(x, function(x) x[order(x$annotation),])

## Are more likely to be in adult over prenatal or vice versa?
fisherres = list(CDS = fisher.test(data.frame({if("CDS" %in% x$adultOnly$annotation)
                       c(x$adultOnly[which(x$adultOnly$annotation=="CDS"),"V1"], sum(x$adultOnly$V1)-x$adultOnly[which(x$adultOnly$annotation=="CDS"),"V1"]) 
                       else c(0,sum(x$adultOnly$V1))},
                       {if("CDS" %in% x$prenatalOnly$annotation)
                       c(x$prenatalOnly[which(x$prenatalOnly$annotation=="CDS"),"V1"],sum(x$prenatalOnly$V1)-x$prenatalOnly[which(x$prenatalOnly$annotation=="CDS"),"V1"]) 
                       else c(0,sum(x$prenatalOnly$V1))})),
                 Intron = fisher.test(data.frame({if("Intron" %in% x$adultOnly$annotation)
                          c(x$adultOnly[which(x$adultOnly$annotation=="Intron"),"V1"], sum(x$adultOnly$V1)-x$adultOnly[which(x$adultOnly$annotation=="Intron"),"V1"]) 
                          else c(0,sum(x$adultOnly$V1))},
                          {if("Intron" %in% x$prenatalOnly$annotation)
                          c(x$prenatalOnly[which(x$prenatalOnly$annotation=="Intron"),"V1"],sum(x$prenatalOnly$V1)-x$prenatalOnly[which(x$prenatalOnly$annotation=="Intron"),"V1"]) 
                          else c(0,sum(x$prenatalOnly$V1))})), 
                 UTR3 = fisher.test(data.frame({if("UTR3" %in% x$adultOnly$annotation)
                        c(x$adultOnly[which(x$adultOnly$annotation=="UTR3"),"V1"], sum(x$adultOnly$V1)-x$adultOnly[which(x$adultOnly$annotation=="UTR3"),"V1"]) 
                        else c(0,sum(x$adultOnly$V1))},
                        {if("UTR3" %in% x$prenatalOnly$annotation)
                        c(x$prenatalOnly[which(x$prenatalOnly$annotation=="UTR3"),"V1"],sum(x$prenatalOnly$V1)-x$prenatalOnly[which(x$prenatalOnly$annotation=="UTR3"),"V1"]) 
                        else c(0,sum(x$prenatalOnly$V1))})), 
                 UTR5 = fisher.test(data.frame({if("UTR5" %in% x$adultOnly$annotation)
                        c(x$adultOnly[which(x$adultOnly$annotation=="UTR5"),"V1"], sum(x$adultOnly$V1)-x$adultOnly[which(x$adultOnly$annotation=="UTR5"),"V1"]) 
                        else c(0,sum(x$adultOnly$V1))},
                        {if("UTR5" %in% x$prenatalOnly$annotation)
                        c(x$prenatalOnly[which(x$prenatalOnly$annotation=="UTR5"),"V1"],sum(x$prenatalOnly$V1)-x$prenatalOnly[which(x$prenatalOnly$annotation=="UTR5"),"V1"]) 
                        else c(0,sum(x$prenatalOnly$V1))})),
                 Other = fisher.test(data.frame({if("Other" %in% x$adultOnly$annotation)
                         c(x$adultOnly[which(x$adultOnly$annotation=="Other"),"V1"], sum(x$adultOnly$V1)-x$adultOnly[which(x$adultOnly$annotation=="Other"),"V1"]) 
                         else c(0,sum(x$adultOnly$V1))},
                         {if("Other" %in% x$prenatalOnly$annotation)
                         c(x$prenatalOnly[which(x$prenatalOnly$annotation=="Other"),"V1"],sum(x$prenatalOnly$V1)-x$prenatalOnly[which(x$prenatalOnly$annotation=="Other"),"V1"]) 
                         else c(0,sum(x$prenatalOnly$V1))})))

unlist(lapply(fisherres, function(x) x$p.value))
#CDS     Intron       UTR3       UTR5      Other 
#0.08350051 0.06785168 0.68091709 1.00000000 0.04171681 

## Break the above down by age
# In Adult: Are  more likely to be in a cyt over nuc or vice versa?
fisherres.ad = list(CDS = fisher.test(data.frame({if("CDS" %in% x$ACnotAN$annotation)
                          c(x$ACnotAN[which(x$ACnotAN$annotation=="CDS"),"V1"], sum(x$ACnotAN$V1)-x$ACnotAN[which(x$ACnotAN$annotation=="CDS"),"V1"]) 
                          else c(0,sum(x$ACnotAN$V1))},
                          {if("CDS" %in% x$ANnotAC$annotation)
                          c(x$ANnotAC[which(x$ANnotAC$annotation=="CDS"),"V1"],sum(x$ANnotAC$V1)-x$ANnotAC[which(x$ANnotAC$annotation=="CDS"),"V1"]) 
                          else c(0,sum(x$ANnotAC$V1))})),
                    Intron = fisher.test(data.frame({if("Intron" %in% x$ACnotAN$annotation)
                          c(x$ACnotAN[which(x$ACnotAN$annotation=="Intron"),"V1"], sum(x$ACnotAN$V1)-x$ACnotAN[which(x$ACnotAN$annotation=="Intron"),"V1"]) 
                          else c(0,sum(x$ACnotAN$V1))},
                          {if("Intron" %in% x$ANnotAC$annotation)
                          c(x$ANnotAC[which(x$ANnotAC$annotation=="Intron"),"V1"],sum(x$ANnotAC$V1)-x$ANnotAC[which(x$ANnotAC$annotation=="Intron"),"V1"]) 
                          else c(0,sum(x$ANnotAC$V1))})), # Intron more likely to be in nucleus than cytosol in adult
                    UTR3 = fisher.test(data.frame({if("UTR3" %in% x$ACnotAN$annotation)
                          c(x$ACnotAN[which(x$ACnotAN$annotation=="UTR3"),"V1"], sum(x$ACnotAN$V1)-x$ACnotAN[which(x$ACnotAN$annotation=="UTR3"),"V1"]) 
                          else c(0,sum(x$ACnotAN$V1))},
                          {if("UTR3" %in% x$ANnotAC$annotation)
                          c(x$ANnotAC[which(x$ANnotAC$annotation=="UTR3"),"V1"],sum(x$ANnotAC$V1)-x$ANnotAC[which(x$ANnotAC$annotation=="UTR3"),"V1"]) 
                          else c(0,sum(x$ANnotAC$V1))})), # 3'UTR is more likely to be in cytosol than nucleus in adult 
                    UTR5 = fisher.test(data.frame({if("UTR5" %in% x$ACnotAN$annotation)
                          c(x$ACnotAN[which(x$ACnotAN$annotation=="UTR5"),"V1"], sum(x$ACnotAN$V1)-x$ACnotAN[which(x$ACnotAN$annotation=="UTR5"),"V1"]) 
                          else c(0,sum(x$ACnotAN$V1))},
                          {if("UTR5" %in% x$ANnotAC$annotation)
                          c(x$ANnotAC[which(x$ANnotAC$annotation=="UTR5"),"V1"],sum(x$ANnotAC$V1)-x$ANnotAC[which(x$ANnotAC$annotation=="UTR5"),"V1"]) 
                          else c(0,sum(x$ANnotAC$V1))})),
                    Other = fisher.test(data.frame({if("Other" %in% x$ACnotAN$annotation)
                          c(x$ACnotAN[which(x$ACnotAN$annotation=="Other"),"V1"], sum(x$ACnotAN$V1)-x$ACnotAN[which(x$ACnotAN$annotation=="Other"),"V1"]) 
                          else c(0,sum(x$ACnotAN$V1))},
                          {if("Other" %in% x$ANnotAC$annotation)
                          c(x$ANnotAC[which(x$ANnotAC$annotation=="Other"),"V1"],sum(x$ANnotAC$V1)-x$ANnotAC[which(x$ANnotAC$annotation=="Other"),"V1"]) 
                          else c(0,sum(x$ANnotAC$V1))})))
unlist(lapply(fisherres.ad, function(x) x$p.value))
#CDS       Intron         UTR3         UTR5        Other 
#1.000000000 0.001071672 0.001468643 1.000000000 1.000000000

# In Prenatal: Are  more likely to be in a cyt over nuc or vice versa?
fisherres.pren = list(CDS = fisher.test(data.frame({if("CDS" %in% x$PCnotPN$annotation)
                          c(x$PCnotPN[which(x$PCnotPN$annotation=="CDS"),"V1"], sum(x$PCnotPN$V1)-x$PCnotPN[which(x$PCnotPN$annotation=="CDS"),"V1"]) 
                          else c(0,sum(x$PCnotPN$V1))},
                          {if("CDS" %in% x$PNnotPC$annotation)
                          c(x$PNnotPC[which(x$PNnotPC$annotation=="CDS"),"V1"],sum(x$PNnotPC$V1)-x$PNnotPC[which(x$PNnotPC$annotation=="CDS"),"V1"]) 
                          else c(0,sum(x$PNnotPC$V1))})),
                    Intron = fisher.test(data.frame({if("Intron" %in% x$PCnotPN$annotation)
                          c(x$PCnotPN[which(x$PCnotPN$annotation=="Intron"),"V1"], sum(x$PCnotPN$V1)-x$PCnotPN[which(x$PCnotPN$annotation=="Intron"),"V1"]) 
                          else c(0,sum(x$PCnotPN$V1))},
                          {if("Intron" %in% x$PNnotPC$annotation)
                          c(x$PNnotPC[which(x$PNnotPC$annotation=="Intron"),"V1"],sum(x$PNnotPC$V1)-x$PNnotPC[which(x$PNnotPC$annotation=="Intron"),"V1"]) 
                          else c(0,sum(x$PNnotPC$V1))})), # Intron is borderline more likely to be in nucleus than cytosol in adult
                    UTR3 = fisher.test(data.frame({if("UTR3" %in% x$PCnotPN$annotation)
                          c(x$PCnotPN[which(x$PCnotPN$annotation=="UTR3"),"V1"], sum(x$PCnotPN$V1)-x$PCnotPN[which(x$PCnotPN$annotation=="UTR3"),"V1"]) 
                          else c(0,sum(x$PCnotPN$V1))},
                          {if("UTR3" %in% x$PNnotPC$annotation)
                          c(x$PNnotPC[which(x$PNnotPC$annotation=="UTR3"),"V1"],sum(x$PNnotPC$V1)-x$PNnotPC[which(x$PNnotPC$annotation=="UTR3"),"V1"]) 
                          else c(0,sum(x$PNnotPC$V1))})), # 3'UTR is borderline more likely to be in cytosol than nucleus in adult 
                    UTR5 = fisher.test(data.frame({if("UTR5" %in% x$PCnotPN$annotation)
                          c(x$PCnotPN[which(x$PCnotPN$annotation=="UTR5"),"V1"], sum(x$PCnotPN$V1)-x$PCnotPN[which(x$PCnotPN$annotation=="UTR5"),"V1"]) 
                          else c(0,sum(x$PCnotPN$V1))},
                          {if("UTR5" %in% x$PNnotPC$annotation)
                          c(x$PNnotPC[which(x$PNnotPC$annotation=="UTR5"),"V1"],sum(x$PNnotPC$V1)-x$PNnotPC[which(x$PNnotPC$annotation=="UTR5"),"V1"]) 
                          else c(0,sum(x$PNnotPC$V1))})),
                    Other = fisher.test(data.frame({if("Other" %in% x$PCnotPN$annotation)
                          c(x$PCnotPN[which(x$PCnotPN$annotation=="Other"),"V1"], sum(x$PCnotPN$V1)-x$PCnotPN[which(x$PCnotPN$annotation=="Other"),"V1"]) 
                          else c(0,sum(x$PCnotPN$V1))},
                          {if("Other" %in% x$PNnotPC$annotation)
                          c(x$PNnotPC[which(x$PNnotPC$annotation=="Other"),"V1"],sum(x$PNnotPC$V1)-x$PNnotPC[which(x$PNnotPC$annotation=="Other"),"V1"]) 
                          else c(0,sum(x$PNnotPC$V1))})))
unlist(lapply(fisherres.pren, function(x) x$p.value))
#CDS       Intron         UTR3         UTR5        Other 
#1.00000000 0.05621382 0.04105724 1.00000000 0.65031604

# In Cytosol: Are  more likely to be in a adult over prenatal or vice versa?
fisherres.cyt = list(CDS = fisher.test(data.frame({if("CDS" %in% x$ACnotPC$annotation)
                           c(x$ACnotPC[which(x$ACnotPC$annotation=="CDS"),"V1"], sum(x$ACnotPC$V1)-x$ACnotPC[which(x$ACnotPC$annotation=="CDS"),"V1"]) 
                           else c(0,sum(x$ACnotPC$V1))},
                           {if("CDS" %in% x$PCnotAC$annotation)
                           c(x$PCnotAC[which(x$PCnotAC$annotation=="CDS"),"V1"],sum(x$PCnotAC$V1)-x$PCnotAC[which(x$PCnotAC$annotation=="CDS"),"V1"]) 
                           else c(0,sum(x$PCnotAC$V1))})), # CDS more likely to be in adult than prenatal in cytosol
                     Intron = fisher.test(data.frame({if("Intron" %in% x$ACnotPC$annotation)
                              c(x$ACnotPC[which(x$ACnotPC$annotation=="Intron"),"V1"], sum(x$ACnotPC$V1)-x$ACnotPC[which(x$ACnotPC$annotation=="Intron"),"V1"]) 
                              else c(0,sum(x$ACnotPC$V1))},
                              {if("Intron" %in% x$PCnotAC$annotation)
                              c(x$PCnotAC[which(x$PCnotAC$annotation=="Intron"),"V1"],sum(x$PCnotAC$V1)-x$PCnotAC[which(x$PCnotAC$annotation=="Intron"),"V1"]) 
                              else c(0,sum(x$PCnotAC$V1))})), # Intron more likely to be in prenatal than adult in cytosol
                     UTR3 = fisher.test(data.frame({if("UTR3" %in% x$ACnotPC$annotation)
                            c(x$ACnotPC[which(x$ACnotPC$annotation=="UTR3"),"V1"], sum(x$ACnotPC$V1)-x$ACnotPC[which(x$ACnotPC$annotation=="UTR3"),"V1"]) 
                            else c(0,sum(x$ACnotPC$V1))},
                            {if("UTR3" %in% x$PCnotAC$annotation)
                            c(x$PCnotAC[which(x$PCnotAC$annotation=="UTR3"),"V1"],sum(x$PCnotAC$V1)-x$PCnotAC[which(x$PCnotAC$annotation=="UTR3"),"V1"]) 
                            else c(0,sum(x$PCnotAC$V1))})), 
                    UTR5 = fisher.test(data.frame({if("UTR5" %in% x$ACnotPC$annotation)
                           c(x$ACnotPC[which(x$ACnotPC$annotation=="UTR5"),"V1"], sum(x$ACnotPC$V1)-x$ACnotPC[which(x$ACnotPC$annotation=="UTR5"),"V1"]) 
                           else c(0,sum(x$ACnotPC$V1))},
                           {if("UTR5" %in% x$PCnotAC$annotation)
                           c(x$PCnotAC[which(x$PCnotAC$annotation=="UTR5"),"V1"],sum(x$PCnotAC$V1)-x$PCnotAC[which(x$PCnotAC$annotation=="UTR5"),"V1"]) 
                           else c(0,sum(x$PCnotAC$V1))})),
                    Other = fisher.test(data.frame({if("Other" %in% x$ACnotPC$annotation)
                            c(x$ACnotPC[which(x$ACnotPC$annotation=="Other"),"V1"], sum(x$ACnotPC$V1)-x$ACnotPC[which(x$ACnotPC$annotation=="Other"),"V1"]) 
                            else c(0,sum(x$ACnotPC$V1))},
                            {if("Other" %in% x$PCnotAC$annotation)
                            c(x$PCnotAC[which(x$PCnotAC$annotation=="Other"),"V1"],sum(x$PCnotAC$V1)-x$PCnotAC[which(x$PCnotAC$annotation=="Other"),"V1"]) 
                            else c(0,sum(x$PCnotAC$V1))})))
unlist(lapply(fisherres.cyt, function(x) x$p.value))
#         CDS       Intron         UTR3         UTR5        Other 
#0.0151781952 0.0002089156 0.1408658815 1.0000000000 0.1097467348
 

# In Nucleus: Are  more likely to be in a adult over prenatal or vice versa?
fisherres.nuc = list(CDS = fisher.test(data.frame({if("CDS" %in% x$ANnotPN$annotation)
                           c(x$ANnotPN[which(x$ANnotPN$annotation=="CDS"),"V1"], sum(x$ANnotPN$V1)-x$ANnotPN[which(x$ANnotPN$annotation=="CDS"),"V1"]) 
                           else c(0,sum(x$ANnotPN$V1))},
                           {if("CDS" %in% x$PNnotAN$annotation)
                           c(x$PNnotAN[which(x$PNnotAN$annotation=="CDS"),"V1"],sum(x$PNnotAN$V1)-x$PNnotAN[which(x$PNnotAN$annotation=="CDS"),"V1"]) 
                           else c(0,sum(x$PNnotAN$V1))})), # CDS more likely to be in adult than prenatal in nucleus
                    Intron = fisher.test(data.frame({if("Intron" %in% x$ANnotPN$annotation)
                           c(x$ANnotPN[which(x$ANnotPN$annotation=="Intron"),"V1"], sum(x$ANnotPN$V1)-x$ANnotPN[which(x$ANnotPN$annotation=="Intron"),"V1"]) 
                           else c(0,sum(x$ANnotPN$V1))},
                           {if("Intron" %in% x$PNnotAN$annotation)
                           c(x$PNnotAN[which(x$PNnotAN$annotation=="Intron"),"V1"],sum(x$PNnotAN$V1)-x$PNnotAN[which(x$PNnotAN$annotation=="Intron"),"V1"]) 
                           else c(0,sum(x$PNnotAN$V1))})), # Intron more likely to be in prenatal than adult in cytosol
                    UTR3 = fisher.test(data.frame({if("UTR3" %in% x$ANnotPN$annotation)
                           c(x$ANnotPN[which(x$ANnotPN$annotation=="UTR3"),"V1"], sum(x$ANnotPN$V1)-x$ANnotPN[which(x$ANnotPN$annotation=="UTR3"),"V1"]) 
                           else c(0,sum(x$ANnotPN$V1))},
                           {if("UTR3" %in% x$PNnotAN$annotation)
                           c(x$PNnotAN[which(x$PNnotAN$annotation=="UTR3"),"V1"],sum(x$PNnotAN$V1)-x$PNnotAN[which(x$PNnotAN$annotation=="UTR3"),"V1"]) 
                           else c(0,sum(x$PNnotAN$V1))})), 
                    UTR5 = fisher.test(data.frame({if("UTR5" %in% x$ANnotPN$annotation)
                           c(x$ANnotPN[which(x$ANnotPN$annotation=="UTR5"),"V1"], sum(x$ANnotPN$V1)-x$ANnotPN[which(x$ANnotPN$annotation=="UTR5"),"V1"]) 
                           else c(0,sum(x$ANnotPN$V1))},
                           {if("UTR5" %in% x$PNnotAN$annotation)
                           c(x$PNnotAN[which(x$PNnotAN$annotation=="UTR5"),"V1"],sum(x$PNnotAN$V1)-x$PNnotAN[which(x$PNnotAN$annotation=="UTR5"),"V1"]) 
                           else c(0,sum(x$PNnotAN$V1))})),
                    Other = fisher.test(data.frame({if("Other" %in% x$ANnotPN$annotation)
                           c(x$ANnotPN[which(x$ANnotPN$annotation=="Other"),"V1"], sum(x$ANnotPN$V1)-x$ANnotPN[which(x$ANnotPN$annotation=="Other"),"V1"]) 
                           else c(0,sum(x$ANnotPN$V1))},
                           {if("Other" %in% x$PNnotAN$annotation)
                           c(x$PNnotAN[which(x$PNnotAN$annotation=="Other"),"V1"],sum(x$PNnotAN$V1)-x$PNnotAN[which(x$PNnotAN$annotation=="Other"),"V1"]) 
                           else c(0,sum(x$PNnotAN$V1))})))
unlist(lapply(fisherres.nuc, function(x) x$p.value))
#CDS       Intron         UTR3         UTR5        Other
#0.04269830 0.01641909 0.38074743 0.26857678 0.17624242

## Gene ontology of different groups of editing sites

entrez.editing = lapply(unique_all, function(x) as.character(na.omit(x$EntrezID)))

split.anno = lapply(unique_all, function(x) split(x, f = x$annotation))
cytosolOnly = split.anno[["cytosolOnly"]]
nucleusOnly = split.anno[["nucleusOnly"]] 
adultOnly = split.anno[["adultOnly"]]   
prenatalOnly = split.anno[["prenatalOnly"]]
ANnotAC = split.anno[["ANnotAC"]]   
ACnotAN = split.anno[["ACnotAN"]]     
ANnotPN = split.anno[["ANnotPN"]]     
PNnotAN = split.anno[["PNnotAN"]]     
ACnotPC = split.anno[["ACnotPC"]]     
PCnotAC = split.anno[["PCnotAC"]]    
PCnotPN = split.anno[["PCnotPN"]]
PNnotPC = split.anno[["PNnotPC"]]

cytosolOnly = lapply(cytosolOnly, function(x) as.character(na.omit(x$EntrezID)))
nucleusOnly = lapply(nucleusOnly, function(x) as.character(na.omit(x$EntrezID)))
adultOnly = lapply(adultOnly, function(x) as.character(na.omit(x$EntrezID)))
prenatalOnly = lapply(prenatalOnly, function(x) as.character(na.omit(x$EntrezID)))
ANnotAC = lapply(ANnotAC, function(x) as.character(na.omit(x$EntrezID)))
ACnotAN = lapply(ACnotAN, function(x) as.character(na.omit(x$EntrezID)))
ANnotPN = lapply(ANnotPN, function(x) as.character(na.omit(x$EntrezID)))
PNnotAN = lapply(PNnotAN, function(x) as.character(na.omit(x$EntrezID))) 
ACnotPC = lapply(ACnotPC, function(x) as.character(na.omit(x$EntrezID)))
PCnotAC = lapply(PCnotAC, function(x) as.character(na.omit(x$EntrezID)))
PCnotPN = lapply(PCnotPN, function(x) as.character(na.omit(x$EntrezID)))
PNnotPC = lapply(PNnotPC, function(x) as.character(na.omit(x$EntrezID)))

## Compare the enriched terms between the unsplit groups
# KEGG
compareKegg = compareCluster(entrez.editing, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Biological Process
compareBP = compareCluster(entrez.editing, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Molecular Function
compareMF = compareCluster(entrez.editing, fun="enrichGO",  ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Cellular Component
compareCC = compareCluster(entrez.editing, fun="enrichGO",  ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
# Disease Ontology
compareDO = compareCluster(entrez.editing, fun="enrichDO",  ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(compareKegg, compareBP, compareMF, compareCC, compareDO, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/interaction.kegg.GO.DO.objects.RNAediting_inAllSamps.rda")
# Plot
pdf(file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/interaction.kegg.GO.DO_unsplit_by_annotation_inAllSamps.pdf", width=14, height=10)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO,colorBy="p.adjust",  showCategory = 45, title= "Disease Ontology Enrichment")
dev.off()

## Compare the enriched terms between the split groups
# KEGG
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/Kegg_split_by_annotation_inAllSamps.pdf")
Kegg.adultOnly = compareCluster(adultOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: adultOnly")
Kegg.prenatalOnly = compareCluster(prenatalOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: prenatalOnly")
Kegg.ANnotAC = compareCluster(ANnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotAC")
Kegg.ACnotAN = compareCluster(ACnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotAN")
Kegg.ANnotPN = compareCluster(ANnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotPN")
Kegg.PNnotAN = compareCluster(PNnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotAN")
Kegg.ACnotPC = compareCluster(ACnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotPC")
Kegg.PCnotAC = compareCluster(PCnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotAC")
Kegg.PNnotPC = compareCluster(PNnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(Kegg.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotPC")
dev.off()
save(Kegg.adultOnly,Kegg.prenatalOnly,Kegg.ANnotAC,Kegg.ACnotAN,
     Kegg.ANnotPN,Kegg.PNnotAN,Kegg.ACnotPC,Kegg.PCnotAC,Kegg.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/kegg.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")
# Biological Process
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/BP_split_by_annotation_inAllSamps.pdf", height = 12, width = 11)
BP.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: adultOnly")
BP.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: prenatalOnly")
BP.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotAC")
BP.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ACnotAN")
BP.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotPN")
BP.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PNnotAN")
BP.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ACnotPC")
BP.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotAC")
BP.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotPN")
BP.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(BP.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PNnotPC")
dev.off()
save(BP.adultOnly, BP.prenatalOnly, BP.ANnotAC, BP.ACnotAN, BP.ANnotPN, BP.PNnotAN, BP.ACnotPC, BP.PCnotAC, BP.PCnotPN, BP.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/BP.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")
# Molecular Function
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/MF_split_by_annotation_inAllSamps.pdf", height = 12, width = 12)
MF.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: adultOnly")
MF.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: prenatalOnly")
MF.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotAN")
MF.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotPN")
MF.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PNnotAN")
MF.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotPC")
MF.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PCnotAC")
MF.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PCnotPN")
MF.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(MF.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PNnotPC")
dev.off()
save(MF.adultOnly, MF.prenatalOnly, MF.ACnotAN, MF.ANnotPN, MF.PNnotAN, MF.ACnotPC, MF.PCnotAC, MF.PCnotPN, MF.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/MF.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")
# Cellular Component
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/CC_split_by_annotation_inAllSamps.pdf", height = 8, width = 10)
CC.cytosolOnly = compareCluster(cytosolOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.nucleusOnly = compareCluster(nucleusOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: adultOnly")
CC.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: prenatalOnly")
CC.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotAN")
CC.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ANnotPN")
CC.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotAN")
CC.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotPC")
CC.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PCnotPN")
CC.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(CC.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotPC")
dev.off()
save(CC.adultOnly, CC.prenatalOnly, CC.ACnotAN, CC.ANnotPN, CC.PNnotAN, CC.ACnotPC, CC.PCnotPN, CC.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/CC.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")
# Disease Ontology
pdf("./Dropbox/sorted_figures/new/github_controlled/rna_editing/figures/DO_split_by_annotation_inAllSamps.pdf")
DO.cytosolOnly = compareCluster(cytosolOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.nucleusOnly = compareCluster(nucleusOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.adultOnly = compareCluster(adultOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: adultOnly")
DO.prenatalOnly = compareCluster(prenatalOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: prenatalOnly")
DO.ANnotAC = compareCluster(ANnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotAC")
DO.ACnotAN = compareCluster(ACnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ANnotPN = compareCluster(ANnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotPN")
DO.PNnotAN = compareCluster(PNnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PNnotAN")
DO.ACnotPC = compareCluster(ACnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ACnotPC")
DO.PCnotAC = compareCluster(PCnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotAC")
DO.PCnotPN = compareCluster(PCnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotPN")
DO.PNnotPC = compareCluster(PNnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
plot(DO.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotAC")
dev.off()
save(DO.adultOnly, DO.prenatalOnly, DO.ANnotAC, DO.ANnotPN, DO.PNnotAN, DO.ACnotPC, DO.PCnotAC, DO.PCnotPN, DO.PNnotPC, 
     file="./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/DO.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")
