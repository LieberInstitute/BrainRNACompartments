library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(VennDiagram)
library(plyr)
library(clusterProfiler)
require("org.Hs.eg.db")

load("./Dropbox/sorted_figures/github_controlled/rna_editing/data/rna_editing_results.rda")
load("./Dropbox/sorted_figures/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

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
txdb = loadDb("./Dropbox/sorted_figures/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
features = list(genes = genes(txdb), CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))

editingres_df$strand = ifelse(editingres_df$conversion=="A:G","+","-")
editingres_df[editingres_df$collapsedconversion!="A:G / T:C","strand"] = "*"
grediting = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
editedgenes = findOverlaps(grediting, features$genes)
genes = as.data.frame(features$genes)
x = cbind(editingres_df[queryHits(editedgenes),], overlappingGene = genes[subjectHits(editedgenes),"gene_id"])
editingres_df = rbind(x, data.frame(editingres_df[-unique(queryHits(editedgenes)),], overlappingGene = "NA"))
grediting = makeGRangesFromDataFrame(editingres_df, keep.extra.columns = T)
annotation = lapply(features, function(y) findOverlaps(grediting, y))

grediting$rnum = 1:length(grediting)
grediting$cds = ifelse(grediting$rnum %in% queryHits(annotation[["CDS"]]), "CDS", NA)
grediting$intron = ifelse(grediting$rnum %in% queryHits(annotation[["Introns"]]), "Intron", NA)
grediting$UTR5 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR5"]]), "5'UTR", NA)
grediting$UTR3 = ifelse(grediting$rnum %in% queryHits(annotation[["UTR3"]]), "3'UTR", NA)
grediting$anno = paste0(grediting$cds,":",grediting$intron, ":", grediting$UTR5, ":", grediting$UTR3)

editing = as.data.frame(grediting)
editing[which(editing$anno == "NA:NA:NA:NA"),"annotation"] = "Intergenic" 
editing[grep("CDS", editing$cds),"annotation"] = "CDS"
editing[which(is.na(editing$annotation) & editing$UTR3 == "3'UTR"),"annotation"] = "3'UTR"
editing[which(is.na(editing$annotation) & editing$UTR5 == "5'UTR"),"annotation"] = "5'UTR"
editing[which(is.na(editing$annotation) & editing$intron == "Intron"),"annotation"] = "Intron"

# Mapping editing sites to nearest gene
geneMapGR = makeGRangesFromDataFrame(geneMap, start.field="Start",end.field="End",strand.field="Strand",keep=TRUE)
dA = distanceToNearest(grediting, geneMapGR)
editing$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
editing$nearestID = names(geneMapGR)[subjectHits(dA)]
editing$distToGene = mcols(dA)$distance
editing$EntrezID = ifelse(editing$overlappingGene!="NA", geneMap[match(editing$overlappingGene, geneMap$gencodeID),"EntrezID"], geneMapGR$EntrezID[subjectHits(dA)])
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
unique = lapply(unique, as.data.frame)
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

samps = c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1","Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")
unique_bySamp = Map(cbind, unique_bySamp, 
                    nos = list(cytosolOnly = rowSums(as.matrix(unique_bySamp$cytosolOnly[,c("Br1113C1","Br2046C","Br2074C","Br5339C1","Br5340C1","Br5341C1")])=="no"),
                               nucleusOnly = rowSums(as.matrix(unique_bySamp$nucleusOnly[,c("Br1113N1","Br2046N","Br2074N","Br5339N1","Br5340N1","Br5341N1")])=="no"),
                               adultOnly = rowSums(as.matrix(unique_bySamp$adultOnly[,c("Br1113C1","Br2046C","Br2074C","Br1113N1","Br2046N","Br2074N")])=="no"),
                               prenatalOnly = rowSums(as.matrix(unique_bySamp$prenatalOnly[,c("Br5339N1","Br5340N1","Br5341N1","Br5339C1","Br5340C1","Br5341C1")])=="no"),
                               ANnotAC = NA, ACnotAN = NA, ANnotPN = NA, PNnotAN = NA, ACnotPC = NA, PCnotAC = NA, PCnotPN = NA, PNnotPC = NA),
                    AllNos = list(cytosolOnly = rowSums(as.matrix(unique_bySamp$cytosolOnly[,samps])=="no"),
                                  nucleusOnly = rowSums(as.matrix(unique_bySamp$nucleusOnly[,])=="no"),
                                  adultOnly = rowSums(as.matrix(unique_bySamp$adultOnly[,samps])=="no"),
                                  prenatalOnly = rowSums(as.matrix(unique_bySamp$prenatalOnly[,samps])=="no"),
                                  ANnotAC = rowSums(as.matrix(unique_bySamp$ANnotAC[,samps])=="no"), 
                                  ACnotAN = rowSums(as.matrix(unique_bySamp$ACnotAN[,samps])=="no"), 
                                  ANnotPN = rowSums(as.matrix(unique_bySamp$ANnotPN[,samps])=="no"), 
                                  PNnotAN = rowSums(as.matrix(unique_bySamp$PNnotAN[,samps])=="no"), 
                                  ACnotPC = rowSums(as.matrix(unique_bySamp$ACnotPC[,samps])=="no"), 
                                  PCnotAC = rowSums(as.matrix(unique_bySamp$PCnotAC[,samps])=="no"), 
                                  PCnotPN = rowSums(as.matrix(unique_bySamp$PCnotPN[,samps])=="no"), 
                                  PNnotPC = rowSums(as.matrix(unique_bySamp$PNnotPC[,samps])=="no")))
lapply(unique_bySamp, head)

save(unique_bySamp, all, editing_anno, geneMap, file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")

load("./Dropbox/sorted_figures/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
elementNROWS(lapply(unique_bySamp[1:4], function(x) x[which(x$nos<3),]))
#cytosolOnly  nucleusOnly    adultOnly prenatalOnly 
#          9           59          619          384
elementNROWS(lapply(unique_bySamp[1:4], function(x) x[which(x$nos<2),]))
#cytosolOnly  nucleusOnly    adultOnly prenatalOnly 
#          2           11          311          208

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
elementNROWS(lapply(unique_all, function(x) unique(x$editingID)))
write.csv(t(data.frame(elementNROWS(lapply(unique_all, function(x) unique(x$editingID))))), quote=F,row.names = F,
          file = "./Dropbox/sorted_figures/github_controlled/rna_editing/data/number_unique_editingSite_inAll.csv")

## What is the distribution of features edited across groups?

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/Genomic_features_uniqueEditingSites_inAllSamps.pdf")
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
pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/Genomic_features_uniqueEditingSites_inAllSamps_percent.pdf")
for (i in 3:length(unique_all)){
  tmp = unique_all[[i]]
  g = ggplot(tmp[collapsedconversion=="A:G / T:C",length(unique(editingID)), by = c("annotation", "Fraction", "Age")], 
             aes(x = Fraction, y = V1, fill = annotation)) + geom_bar(position = "fill", stat = "identity") +
    facet_grid(. ~ Age) +
    labs(fill="") +
    ylab("Proportion") + 
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
unique_alldt$Status = factor(unique_alldt$Status, levels = c("cytosolOnly", "nucleusOnly", "adultOnly", "prenatalOnly", "ACnotAN", "ANnotAC", 
                                                             "PCnotPN", "PNnotPC","ACnotPC", "PCnotAC","ANnotPN", "PNnotAN"))

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/valdepth_unique_editing_sites_byGroup_bySample_inAllSamps.pdf", width = 28, height = 7)
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

x = lapply(unique_all[3:12], function(x) as.data.frame(cbind(x[,length(unique(editingID)),by="annotation"], total = length(unique(x$editingID)))))
x = Map(cbind, x, diff = lapply(x, function(y) y$total - y$V1))
comps = list(byAge = c("adultOnly","prenatalOnly"), byFracInAdult = c("ACnotAN","ANnotAC"),
             byFracinPrenatal = c("PCnotPN","PNnotPC"), byAgeinCyt = c("ACnotPC","PCnotAC"), byAgeinNuc = c("ANnotPN", "PNnotAN"))
anno = c("CDS","Intron","Intergenic","3'UTR","5'UTR")


## Are  more likely to be in cyt over nuc or vice versa?
anno.comps = list(list(),list(),list(),list(),list())
for (i in (1:length(comps))){
  for (j in (1:length(anno))){
    anno.comps[[i]][[j]] = data.frame(
      {if(anno[j] %in% x[[comps[[i]][1]]][,"annotation"])
      c(x[[comps[[i]][1]]][which(x[[comps[[i]][1]]][,"annotation"]==anno[j]),"V1"],
        x[[comps[[i]][1]]][which(x[[comps[[i]][1]]][,"annotation"]==anno[j]),"diff"]) else Cyt.Adult = c(0,x[[comps[[i]][1]]][1,"total"])},
      {if(anno[j] %in% x[[comps[[i]][2]]][,"annotation"])
      c(x[[comps[[i]][2]]][which(x[[comps[[i]][2]]][,"annotation"]==anno[j]),"V1"],
        x[[comps[[i]][2]]][which(x[[comps[[i]][2]]][,"annotation"]==anno[j]),"diff"]) else Nuc.Prenatal = c(0,x[[comps[[i]][2]]][1,"total"])}, 
      row.names = c("inAnno","notAnno"))
    colnames(anno.comps[[i]][[j]]) = c("Cyt.Adult", "Nuc.Prenatal")
  }
  names(anno.comps[[i]]) = anno
}
names(anno.comps) = names(comps)
anno.comps = lapply(anno.comps, lapply, fisher.test)
df = do.call(rbind, Map(cbind, Comparison = as.list(names(anno.comps)), lapply(anno.comps, function(x) do.call(rbind, Map(cbind, annotation = as.list(names(x)), 
                                                                                 lapply(x, function(y) data.frame(pval=y$p.value, OddsRatio=y$estimate, row.names = NULL)))))))
df$FDR = p.adjust(df$pval, method = "fdr")
write.csv(df, quote=F, file= "./Dropbox/sorted_figures/github_controlled/rna_editing/data/fisher_annotationEnrichment_inUniqueSites_inAllSamps_byGroup.csv")
df[df$FDR<=0.05,]
#    Comparison annotation         pval OddsRatio         FDR
# byFracInAdult     Intron 0.0004544449 0.1693807 0.009785001
# byFracInAdult      3'UTR 0.0007828001 3.9617189 0.009785001


## Gene ontology of different groups of editing sites

entrez.editing = lapply(unique_all, as.data.frame)
entrez.editing = lapply(entrez.editing, function(x) unique(as.character(na.omit(x[x$distToGene<500,"EntrezID"]))))
GeneUniverse = unique(as.character(na.omit(editing_anno[editing_anno$distToGene<500,]$EntrezID)))

split.anno = lapply(unique_all, function(x) split(as.data.frame(x), f = x$annotation))
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

cytosolOnly = lapply(cytosolOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
nucleusOnly = lapply(nucleusOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
adultOnly = lapply(adultOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
prenatalOnly = lapply(prenatalOnly, function(x) unique(as.character(na.omit(x$EntrezID))))
ANnotAC = lapply(ANnotAC, function(x) unique(as.character(na.omit(x$EntrezID))))
ACnotAN = lapply(ACnotAN, function(x) unique(as.character(na.omit(x$EntrezID))))
ANnotPN = lapply(ANnotPN, function(x) unique(as.character(na.omit(x$EntrezID))))
PNnotAN = lapply(PNnotAN, function(x) unique(as.character(na.omit(x$EntrezID)))) 
ACnotPC = lapply(ACnotPC, function(x) unique(as.character(na.omit(x$EntrezID))))
PCnotAC = lapply(PCnotAC, function(x) unique(as.character(na.omit(x$EntrezID))))
PCnotPN = lapply(PCnotPN, function(x) unique(as.character(na.omit(x$EntrezID))))
PNnotPC = lapply(PNnotPC, function(x) unique(as.character(na.omit(x$EntrezID))))

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
     file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/interaction.kegg.GO.DO.objects.RNAediting_inAllSamps.rda")
# Plot
pdf(file="./Dropbox/sorted_figures/github_controlled/rna_editing/figures/interaction.kegg.GO.DO_unsplit_by_annotation_inAllSamps.pdf", width=10, height=6)
plot(compareKegg,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment")
plot(compareBP,colorBy="p.adjust",  showCategory = 45, title= "Biological Process GO Enrichment")
plot(compareMF,colorBy="p.adjust",  showCategory = 45, title= "Molecular Function GO Enrichment")
plot(compareCC,colorBy="p.adjust",  showCategory = 45, title= "Cellular Compartment GO Enrichment")
plot(compareDO,colorBy="p.adjust",  showCategory = 45, title= "Disease Ontology Enrichment")
dev.off()

## Compare the enriched terms between the split groups

# KEGG
Kegg.adultOnly = compareCluster(adultOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.prenatalOnly = compareCluster(prenatalOnly, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ANnotAC = compareCluster(ANnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ACnotAN = compareCluster(ACnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ANnotPN = compareCluster(ANnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PNnotAN = compareCluster(PNnotAN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.ACnotPC = compareCluster(ACnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PCnotAC = compareCluster(PCnotAC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PNnotPC = compareCluster(PNnotPC, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
Kegg.PCnotPN = compareCluster(PCnotPN, fun="enrichKEGG", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(Kegg.adultOnly,Kegg.ANnotAC,Kegg.ACnotAN,Kegg.ANnotPN,Kegg.PNnotAN,Kegg.PCnotAC,Kegg.PCnotPN,Kegg.PNnotPC, 
     file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/kegg.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/Kegg_split_by_annotation_inAllSamps.pdf")
plot(Kegg.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: adultOnly")
plot(Kegg.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotAC")
plot(Kegg.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ACnotAN")
plot(Kegg.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: ANnotPN")
plot(Kegg.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotAN")
plot(Kegg.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotAC")
plot(Kegg.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PCnotPN")
plot(Kegg.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "KEGG Pathway Enrichment: PNnotPC")
dev.off()

# Biological Process
BP.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
BP.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
goList_BP = lapply(entrez.editing, function(x) enrichGO(entrez.editing$adultOnly, ont = "BP", OrgDb = org.Hs.eg.db, 
                                                        universe= as.character(unique(all$adultAll$EntrezID)), pAdjustMethod="BH"))

save(goList_BP, BP.adultOnly,BP.prenatalOnly,BP.ANnotAC,BP.ACnotAN,BP.ANnotPN,BP.PNnotAN,BP.ACnotPC,BP.PCnotAC,BP.PCnotPN,BP.PNnotPC, 
     file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/sorted/BP.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/BP_split_by_annotation_inAllSamps.pdf", height = 12, width = 11)
plot(BP.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: adultOnly")
plot(BP.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: prenatalOnly")
plot(BP.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotAC")
plot(BP.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ACnotAN")
plot(BP.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ANnotPN")
plot(BP.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PNnotAN")
plot(BP.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: ACnotPC")
plot(BP.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotAC")
plot(BP.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PCnotPN")
plot(BP.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "BP Pathway Enrichment: PNnotPC")
dev.off()

# Molecular Function
MF.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
MF.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "MF", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(MF.adultOnly,MF.prenatalOnly,MF.ANnotAC,MF.ACnotAN,MF.ANnotPN,MF.PNnotAN,MF.PCnotAC,MF.PCnotPN,MF.PNnotPC, 
     file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/MF.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/MF_split_by_annotation_inAllSamps.pdf", height = 12, width = 12)
plot(MF.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: adultOnly")
plot(MF.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: prenatalOnly")
plot(MF.ANnotAC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotAC")
plot(MF.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ACnotAN")
plot(MF.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: ANnotPN")
plot(MF.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PNnotAN")
plot(MF.PCnotAC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PCnotAC")
plot(MF.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PCnotPN")
plot(MF.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "MF Pathway Enrichment: PNnotPC")
dev.off()

# Cellular Component
CC.adultOnly = compareCluster(adultOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.prenatalOnly = compareCluster(prenatalOnly, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ANnotAC = compareCluster(ANnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ACnotAN = compareCluster(ACnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ANnotPN = compareCluster(ANnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PNnotAN = compareCluster(PNnotAN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.ACnotPC = compareCluster(ACnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PCnotAC = compareCluster(PCnotAC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PCnotPN = compareCluster(PCnotPN, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
CC.PNnotPC = compareCluster(PNnotPC, fun="enrichGO", ont = "CC", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(CC.adultOnly,CC.prenatalOnly,CC.ACnotAN,CC.ANnotPN,CC.PNnotAN,CC.ACnotPC,CC.PCnotPN,CC.PNnotPC,  
     file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/CC.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/CC_split_by_annotation_inAllSamps.pdf", height = 8, width = 10)
plot(CC.adultOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: adultOnly")
plot(CC.prenatalOnly,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: prenatalOnly")
plot(CC.ACnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotAN")
plot(CC.ANnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ANnotPN")
plot(CC.PNnotAN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotAN")
plot(CC.ACnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: ACnotPC")
plot(CC.PCnotPN,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PCnotPN")
plot(CC.PNnotPC,colorBy="p.adjust",  showCategory = 45, title= "CC Pathway Enrichment: PNnotPC")
dev.off()

# Disease Ontology
DO.adultOnly = compareCluster(adultOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.prenatalOnly = compareCluster(prenatalOnly, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ANnotAC = compareCluster(ANnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ACnotAN = compareCluster(ACnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ANnotPN = compareCluster(ANnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PNnotAN = compareCluster(PNnotAN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.ACnotPC = compareCluster(ACnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PCnotAC = compareCluster(PCnotAC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PCnotPN = compareCluster(PCnotPN, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
DO.PNnotPC = compareCluster(PNnotPC, fun="enrichDO", ont = "DO", qvalueCutoff = 0.05, pvalueCutoff = 0.05)
save(DO.prenatalOnly,DO.ANnotAC,DO.ANnotPN,DO.PNnotAN,DO.PCnotAC,DO.PCnotPN,DO.PNnotPC,
     file="./Dropbox/sorted_figures/github_controlled/rna_editing/data/DO.objects.RNAediting.SplitByAnnotation_inAllSamps.rda")

pdf("./Dropbox/sorted_figures/github_controlled/rna_editing/figures/DO_split_by_annotation_inAllSamps.pdf")
plot(DO.prenatalOnly, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: prenatalOnly")
plot(DO.ANnotAC, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotAC")
plot(DO.ANnotPN, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: ANnotPN")
plot(DO.PNnotAN, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PNnotAN")
plot(DO.PCnotAC, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotAC")
plot(DO.PCnotPN, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PCnotPN")
plot(DO.PNnotPC, colorBy="p.adjust",  showCategory = 45, title= "DO Pathway Enrichment: PNnotPC")
dev.off()

### Characterize the overlap with retained introns

# read in results files
load("./Dropbox/sorted_figures/github_controlled/intron_retention/data/intron_IR_comparisons/introns_object.rda")

# Find overlaps with RNA editing sites
IRranges = lapply(introns, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr", start.field="Start",end.field="End",strand.field="Direction",keep.extra.columns = T))
editing_anno = editing_anno[collapsedconversion=="A:G / T:C",,]
editing_ranges = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
IR_hits = lapply(IRranges, function(x) findOverlaps(editing_ranges, x))
editingIR = mapply(function(IR, oo) cbind(IR[subjectHits(oo),], editing_anno[queryHits(oo),]), introns[elementNROWS(IR_hits)>0], IR_hits[elementNROWS(IR_hits)>0], SIMPLIFY = F)
editingIR = do.call(rbind, Map(cbind, editingIR[elementNROWS(editingIR)>0], IRgroup = as.list(names(editingIR[elementNROWS(editingIR)>0]))))
head(editingIR)
unique(editingIR[editingIR$IRgroup!="All Introns","editingID"])
# chr1:1959229-1959229    
# chr22:51167546-51167546
editingIR[editingIR$IRgroup!="All Introns",]
