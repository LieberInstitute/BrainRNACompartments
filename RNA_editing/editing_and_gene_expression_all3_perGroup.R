library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")


### Isolate the sites present in all samples in each group

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
elementNROWS(lapply(unique_bySamp_all, function(x) x[which(x$AllNos==9),]))

unique_all = lapply(unique_bySamp_all[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
unique_all = Map(cbind, unique_all, ensID = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique_all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]),
                 baseMean = lapply(unique_all, function(x) Ipres.down[match(x$nearestID, rownames(Ipres.down)),"baseMean"]), 
                 Prenatal.LFC = lapply(unique_all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "log2FoldChange"]), 
                 Prenatal.SE = lapply(unique_all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "lfcSE"]), 
                 Prenatal.padj = lapply(unique_all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "padj"]),
                 Adult.LFC = lapply(unique_all, function(x) Apres[match(x$nearestID, rownames(Apres)), "log2FoldChange"]), 
                 Adult.SE = lapply(unique_all, function(x) Apres[match(x$nearestID, rownames(Apres)), "lfcSE"]),
                 Adult.padj = lapply(unique_all, function(x) Apres[match(x$nearestID, rownames(Apres)), "padj"]),
                 Nucleus.LFC = lapply(unique_all, function(x) Npres[match(x$nearestID, rownames(Npres)), "log2FoldChange"]), 
                 Nucleus.SE = lapply(unique_all, function(x) Npres[match(x$nearestID, rownames(Npres)), "lfcSE"]), 
                 Nucleus.padj = lapply(unique_all, function(x) Npres[match(x$nearestID, rownames(Npres)), "padj"]),
                 Cytosol.LFC = lapply(unique_all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "log2FoldChange"]), 
                 Cytosol.SE = lapply(unique_all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "lfcSE"]),
                 Cytosol.padj = lapply(unique_all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "padj"]))
DEGnames = c(lapply(sig[1:8], function(x) as.character(x$geneID)), lapply(age.sig, function(x) as.character(x$geneID))) 
unique_all = Map(cbind, unique_all, both_retained = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_retained, "both_retained", "no")),
                 both_exported = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_exported, "both_exported", "no")),
                 Fet_retained = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
                 Ad_retained = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
                 Fet_exported = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
                 Ad_exported = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
                 ret_Ad_exp_Fet = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
                 ret_Fet_exp_Ad = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
                 interacting = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$interacting, "interacting", "no")),
                 both_decreasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
                 both_increasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$both_increasing, "both_increasing", "no")),
                 Cyt_decreasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
                 Nuc_decreasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
                 Cyt_increasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
                 Nuc_increasing = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
                 decr_Nuc_incr_Cyt = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
                 decr_Cyt_incr_Nuc = lapply(unique_all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))
elementNROWS(unique_all)
lapply(unique_all, head)

all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]),
          baseMean = lapply(all, function(x) Ipres.down[match(x$nearestID, rownames(Ipres.down)),"baseMean"]), 
          Prenatal.LFC = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "log2FoldChange"]), 
          Prenatal.SE = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "lfcSE"]), 
          Prenatal.padj = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "padj"]),
          Adult.LFC = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "log2FoldChange"]), 
          Adult.SE = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "lfcSE"]),
          Adult.padj = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "padj"]),
          Nucleus.LFC = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "log2FoldChange"]), 
          Nucleus.SE = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "lfcSE"]), 
          Nucleus.padj = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "padj"]),
          Cytosol.LFC = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "log2FoldChange"]), 
          Cytosol.SE = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "lfcSE"]),
          Cytosol.padj = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "padj"]),
          both_retained = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_retained, "both_retained", "no")),
          both_exported = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_exported, "both_exported", "no")),
          Fet_retained = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
          Ad_retained = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
          Fet_exported = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
          Ad_exported = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
          ret_Ad_exp_Fet = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
          ret_Fet_exp_Ad = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
          interacting = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$interacting, "interacting", "no")),
          both_decreasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
          both_increasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$both_increasing, "both_increasing", "no")),
          Cyt_decreasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
          Nuc_decreasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
          Cyt_increasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
          Nuc_increasing = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
          decr_Nuc_incr_Cyt = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
          decr_Cyt_incr_Nuc = lapply(all, function(x) ifelse(x$nearestID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))


### How many sites are there per gene, and is this different than the distribution of sites that aren't specific to a group and in all samples in a group?

## Number of sites per gene
numsites_bygene = lapply(unique_all, function(x) x[,length(unique(editingID)),by="nearestID"])
numsites_bygene_all = lapply(all, function(x) x[,length(unique(editingID)),by="nearestID"])
numsites_bygene = lapply(numsites_bygene,as.data.frame)
numsites_bygene_all = lapply(numsites_bygene_all, as.data.frame)


## Compare to the number of sites per gene in sites present in each compartment but not unique, and the total number of sites in the same genes as those of the unique sites
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
t.num.bygene = t.editedgenes.numsites = list()
for (i in 1:length(numsites_bygene)){
  t.num.bygene[[i]] = t.test(x = numsites_bygene[[group[i]]][,"V1"], y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"nearestID"] %in% numsites_bygene[[group[i]]][,"nearestID"]),"V1"])
  t.editedgenes.numsites[[i]] = t.test(x = numsites_bygene_all[[allgroup[i]]][which(numsites_bygene_all[[allgroup[i]]][,"nearestID"] %in% numsites_bygene[[group[i]]][,"nearestID"]),"V1"], 
                                       y = numsites_bygene_all[[allgroup[i]]][-which(numsites_bygene_all[[allgroup[i]]][,"nearestID"] %in% numsites_bygene[[group[i]]][,"nearestID"]),"V1"])
}
names(t.num.bygene) = names(t.editedgenes.numsites) = group
t.num.bygene.editing = rbind(Tstat = data.frame(lapply(t.num.bygene, function(x) round(x$statistic,3))), pval = data.frame(lapply(t.num.bygene, function(x) x$p.value)),
                             confInt = data.frame(lapply(t.num.bygene, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(t.num.bygene, function(x) round(x$estimate,3))),
                             max = data.frame(lapply(numsites_bygene, function(x) max(x$V1))), median = data.frame(lapply(numsites_bygene, function(x) median(x$V1))), 
                             sd = data.frame(lapply(numsites_bygene, function(x) round(sd(x$V1),3))))
#                     adultOnly  prenatalOnly       ACnotAN       ANnotAC       PCnotPN       PNnotPC       ACnotPC       PCnotAC       ANnotPN      PNnotAN
#Tstat              -1.29200e+01 -6.487000e+00 -1.498300e+01 -1.846700e+01 -1.581900e+01 -1.647300e+01 -6.995000e+00 -4.147000e+00 -9.159000e+00 -5.74400e+00
#pval                1.18511e-27  1.776006e-08  4.295059e-29  8.882124e-65  3.296526e-21  1.783768e-44  1.367934e-11  4.957483e-05  4.538716e-19  2.78473e-08
#confInt.1          -3.09000e+00 -2.317000e+00 -2.541000e+00 -2.694000e+00 -2.435000e+00 -2.214000e+00 -1.439000e+00 -1.113000e+00 -1.645000e+00 -1.26100e+00
#confInt.2          -2.27200e+00 -1.225000e+00 -1.948000e+00 -2.176000e+00 -1.887000e+00 -1.742000e+00 -8.070000e-01 -3.960000e-01 -1.064000e+00 -6.17000e-01
#estMeans.mean of x  1.68800e+00  1.935000e+00  1.200000e+00  1.359000e+00  1.091000e+00  1.275000e+00  1.845000e+00  1.962000e+00  1.951000e+00  1.90700e+00
#estMeans.mean of y  4.36900e+00  3.706000e+00  3.444000e+00  3.794000e+00  3.252000e+00  3.253000e+00  2.968000e+00  2.717000e+00  3.305000e+00  2.84600e+00
#max                 7.00000e+00  9.000000e+00  3.000000e+00  6.000000e+00  2.000000e+00  3.000000e+00  9.000000e+00  1.300000e+01  1.500000e+01  1.10000e+01
#median              1.00000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.00000e+00
#sd                  1.45300e+00  1.718000e+00  5.000000e-01  8.140000e-01  3.020000e-01  5.320000e-01  1.564000e+00  1.885000e+00  1.903000e+00  1.77400e+00
t.editedgenes.numsites.df = rbind(Tstat = data.frame(lapply(t.editedgenes.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedgenes.numsites, function(x) x$p.value)),
                                  confInt = data.frame(lapply(t.editedgenes.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedgenes.numsites, function(x) x$estimate)),
                                  max = mapply(function(x,y) max(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup),
                                  median = mapply(function(x,y) median(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup),
                                  sd = mapply(function(x,y) sd(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup))
#                     adultOnly prenatalOnly      ACnotAN      ANnotAC     PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#Tstat              6.833708e+00 7.925401e+00 4.714393e+00 9.918612e+00  2.67267044 8.504629e+00 1.005641e+01 1.054639e+01 1.365258e+01 1.206220e+01
#pval               1.705879e-09 4.222442e-10 8.480397e-05 3.086173e-17  0.02331645 2.415161e-11 1.448170e-18 2.334049e-19 2.996391e-33 5.930008e-24
#confInt.1          7.192684e+00 1.025659e+01 5.778373e+00 8.934508e+00  1.06451307 6.442718e+00 5.802645e+00 5.929871e+00 7.462534e+00 7.002663e+00
#confInt.2          1.310791e+01 1.724408e+01 1.477300e+01 1.339197e+01 11.70432431 1.042453e+01 8.640001e+00 8.667021e+00 9.976335e+00 9.745848e+00
#estMeans.mean of x 1.451948e+01 1.745652e+01 1.372000e+01 1.495726e+01  9.63636364 1.168627e+01 1.018919e+01 1.001504e+01 1.202465e+01 1.122000e+01
#estMeans.mean of y 4.369185e+00 3.706188e+00 3.444312e+00 3.794027e+00  3.25194494 3.252653e+00 2.967866e+00 2.716591e+00 3.305213e+00 2.845745e+00
#max                7.500000e+01 5.300000e+01 4.800000e+01 5.100000e+01 36.00000000 4.400000e+01 4.800000e+01 2.900000e+01 5.300000e+01 4.400000e+01
#median             2.000000e+00 2.000000e+00 1.000000e+00 2.000000e+00  2.00000000 2.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 2.000000e+00
#sd                 6.480045e+00 4.903550e+00 4.569718e+00 5.165724e+00  4.16809461 4.189500e+00 3.797138e+00 3.140198e+00 4.385060e+00 3.291195e+00


## What about at the exon level?

exonMap$exonID = rownames(exonMap)
exonMap_gr = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
ov = findOverlaps(editing_anno_gr, exonMap_gr)
tog = cbind(editing_anno[queryHits(ov),], exonMap[subjectHits(ov),])
unique_exons = lapply(unique_all, function(x) tog[(editingID %in% x$editingID),,])
all_exons = lapply(all, function(x) tog[(editingID %in% x$editingID),,])

numsites_byexon = lapply(unique_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon_all = lapply(all_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon = lapply(numsites_byexon,as.data.frame)
numsites_byexon_all = lapply(numsites_byexon_all, as.data.frame)

t.num.byexon = t.editedexons.numsites = list()
for (i in 1:length(numsites_byexon)){
  t.num.byexon[[i]] = t.test(x = numsites_byexon[[group[i]]][,"V1"], y = numsites_byexon_all[[allgroup[i]]][-which(numsites_byexon_all[[allgroup[i]]][,"exonID"] %in% numsites_byexon[[group[i]]][,"exonID"]),"V1"])
  t.editedexons.numsites[[i]] = t.test(x = numsites_byexon_all[[allgroup[i]]][which(numsites_byexon_all[[allgroup[i]]][,"exonID"] %in% numsites_byexon[[group[i]]][,"exonID"]),"V1"], 
                                       y = numsites_byexon_all[[allgroup[i]]][-which(numsites_byexon_all[[allgroup[i]]][,"exonID"] %in% numsites_byexon[[group[i]]][,"exonID"]),"V1"])
}
names(t.num.byexon) = names(t.editedexons.numsites) = group
t.num.byexon.editing = rbind(Tstat = data.frame(lapply(t.num.byexon, function(x) round(x$statistic,3))), pval = data.frame(lapply(t.num.byexon, function(x) x$p.value)),
                             confInt = data.frame(lapply(t.num.byexon, function(x) round(x$conf.int,3))), estMeans = data.frame(lapply(t.num.byexon, function(x) round(x$estimate,3))),
                             max = data.frame(lapply(numsites_byexon, function(x) max(x$V1))), median = data.frame(lapply(numsites_byexon, function(x) median(x$V1))), 
                             sd = data.frame(lapply(numsites_byexon, function(x) round(sd(x$V1),3))))
#                       adultOnly prenatalOnly       ACnotAN       ANnotAC      PCnotPN       PNnotPC       ACnotPC     PCnotAC      ANnotPN      PNnotAN
#Tstat              -1.450700e+01  -3.04800000 -1.395200e+01 -1.364100e+01 -1.14670e+01 -1.371900e+01 -6.264000e+00 -2.09400000 -6.41800e+00 -2.611000000
#pval                7.826887e-41   0.00383909  2.083555e-32  1.904047e-36  3.48763e-11  1.134369e-27  7.737729e-10  0.03800362  2.26629e-10  0.009898742
#confInt.1          -2.227000e+00  -1.48300000 -1.868000e+00 -1.879000e+00 -1.72700e+00 -1.588000e+00 -9.750000e-01 -0.72600000 -9.29000e-01 -0.807000000
#confInt.2          -1.696000e+00  -0.30300000 -1.406000e+00 -1.406000e+00 -1.20000e+00 -1.188000e+00 -5.090000e-01 -0.02100000 -4.94000e-01 -0.112000000
#estMeans.mean of x  1.336000e+00   1.87500000  1.167000e+00  1.333000e+00  1.10000e+00  1.118000e+00  1.656000e+00  1.84100000  1.77800e+00  1.818000000
#estMeans.mean of y  3.298000e+00   2.76800000  2.804000e+00  2.976000e+00  2.56300e+00  2.506000e+00  2.398000e+00  2.21400000  2.49000e+00  2.278000000
#max                 7.000000e+00   9.00000000  2.000000e+00  6.000000e+00  2.00000e+00  3.000000e+00  8.000000e+00 11.00000000  1.30000e+01 11.000000000
#median              1.000000e+00   1.00000000  1.000000e+00  1.000000e+00  1.00000e+00  1.000000e+00  1.000000e+00  1.00000000  1.00000e+00  1.000000000
#sd                  9.210000e-01   1.78600000  3.790000e-01  8.460000e-01  3.16000e-01  4.090000e-01  1.193000e+00  1.77600000  1.53600e+00  1.803000000
t.editedexons.numsites.df = rbind(Tstat = data.frame(lapply(t.editedexons.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedexons.numsites, function(x) x$p.value)),
                                  confInt = data.frame(lapply(t.editedexons.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedexons.numsites, function(x) x$estimate)),
                                  max = mapply(function(x,y) max(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                                  median = mapply(function(x,y) median(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                                  sd = mapply(function(x,y) sd(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup))
#                     adultOnly prenatalOnly      ACnotAN      ANnotAC    PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#Tstat              6.752579e+00 6.826467e+00 5.336839e+00 7.485853e+00  1.7336315 5.109071e+00 9.276707e+00 8.604770e+00 1.189163e+01 9.052067e+00
#pval               6.877949e-10 3.632571e-08 9.726992e-06 1.774729e-11  0.1168335 1.315697e-05 4.772373e-17 4.857557e-14 1.943060e-27 2.679405e-15
#confInt.1          4.716173e+00 8.080140e+00 5.076844e+00 6.019345e+00 -0.8009582 3.944150e+00 4.795662e+00 4.705898e+00 5.587090e+00 4.541880e+00
#confInt.2          8.633537e+00 1.488328e+01 1.138245e+01 1.035304e+01  6.0741083 9.161925e+00 7.386507e+00 7.520797e+00 7.801852e+00 7.084327e+00
#estMeans.mean of x 9.972727e+00 1.425000e+01 1.103333e+01 1.116216e+01  5.2000000 9.058824e+00 8.488889e+00 8.327434e+00 9.184615e+00 8.090909e+00
#estMeans.mean of y 3.297872e+00 2.768288e+00 2.803684e+00 2.975971e+00  2.5634249 2.505786e+00 2.397804e+00 2.214086e+00 2.490145e+00 2.277806e+00
#max                7.700000e+01 4.100000e+01 5.300000e+01 5.800000e+01 36.0000000 3.200000e+01 4.600000e+01 2.900000e+01 3.300000e+01 3.200000e+01
#median             1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00  1.0000000 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
#sd                 5.592686e+00 3.884978e+00 4.303558e+00 4.485847e+00  3.4486019 3.318764e+00 3.427085e+00 2.639297e+00 3.392297e+00 2.867206e+00



### Compare the proportion of editing site in each annotation in sites unique to a specific group and those found in all samples in a group or not
# Using whether the site overlaps an annotation group at all rather than the CDS -> UTR -> Intron -> Other heirarchy

anno.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
unique = lapply(unique_bySamp[3:12], function(x) editing_anno[which(editingID %in% x$EditingID),,])
for (i in 1:length(unique_all)){
  anno.site[[i]] = list(
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("CDS", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("CDS", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("CDS", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("Intron", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("Intron", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("Intron", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("UTR5", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("UTR5", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("UTR5", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][grep("UTR3", anno),list(unique(editingID)),]), nrow(unique_all[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][grep("UTR3", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][grep("UTR3", anno),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])-nrow(unique_all[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(inAll = c(nrow(unique_all[[group[i]]][(annotation=="Other"),list(unique(editingID)),]), nrow(unique_all[[group[i]]][(annotation!="Other"),list(unique(editingID)),])),
               notInAll = c(nrow(unique[[group[i]]][(annotation=="Other"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(annotation=="Other"),list(unique(editingID)),]),
                            nrow(unique[[group[i]]][(annotation!="Other"),list(unique(editingID)),])-nrow(unique_all[[group[i]]][(annotation!="Other"),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")))
  names(anno.site[[i]]) = c("CDS","Intron","UTR5","UTR3","Other")
}
names(anno.site) = group
fisher.anno.site = lapply(anno.site, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#         adultOnly prenatalOnly   ACnotAN    ANnotAC   PCnotPN   PNnotPC   ACnotPC    PCnotAC      ANnotPN      PNnotAN
# CDS    0.030786642  1.000000000 0.3500569 0.01411564 1.0000000 0.5776303 0.8352318 1.00000000 0.4127481405 0.0095894643
# Intron 0.032163580  0.511061230 0.8546793 0.01002666 0.3947897 0.7932542 1.0000000 0.13681736 0.4759762429 0.3787189264
# UTR5   1.000000000  1.000000000 1.0000000 1.00000000 1.0000000 0.5472894 1.0000000 1.00000000 0.0625514485 0.4371744270
# UTR3   0.009019958  0.000281643 0.1458043 0.31233307 0.3854719 0.2877337 0.2118689 0.03813118 0.0001286324 0.0008078246
# Other  0.572627450  0.039140302 0.8199467 0.53044538 1.0000000 0.4044935 0.4195773 0.13451049 0.0978798584 0.0021495928
anno.site.props = lapply(anno.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### Characterize the overlap with differentially expressed genes

## Are editing sites specific to a group and shared by all samples in a group more or less likely to fall in an exported or retained gene by fraction than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
frac.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
for (i in 1:length(unique_all)){
  frac.site[[i]] = list(
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),]),
                 nrow(all[[allgroup[i]]][(both_retained=="both_retained"),list(unique(editingID)),])-
                   nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),]),
                 nrow(all[[allgroup[i]]][(both_exported=="both_exported"),list(unique(editingID)),])-
                   nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),]),
                 nrow(all[[allgroup[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])-
                   nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),]),
                 nrow(all[[allgroup[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])-
                   nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),]),
                 nrow(all[[allgroup[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])-
                   nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),]),
                 nrow(all[[allgroup[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])-
                   nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(frac.site[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG")
}
names(frac.site) = group
fisher.frac.site = lapply(frac.site, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#            adultOnly prenatalOnly    ACnotAN     ANnotAC PCnotPN    PNnotPC    ACnotPC   PCnotAC   ANnotPN   PNnotAN
#bothagesDEG 1.0000000            1 1.00000000 1.000000000       1 1.00000000 1.00000000 1.0000000 1.0000000 1.0000000
#adultDEG    0.8614905            1 0.01806298 0.002091996       1 0.07770767 0.05007798 0.5867589 0.1335647 0.5643383
#prenatalDEG 1.0000000            1 1.00000000 1.000000000       1 1.00000000 1.00000000 1.0000000 1.0000000 1.0000000
frac.site.props = lapply(frac.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))


# number of significantly retained or exported genes that either contain an editing site specific to a group and in all samples, or don't
frac.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  frac.gene[[i]] = list(
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),]),
                            length(unique(sig$both_retained$geneID))-nrow(unique_all[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),]),
                            length(unique(sig$both_exported$geneID))-nrow(unique_all[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Ad_retained$geneID)))-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Ad_exported$geneID)))-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Fet_retained$geneID)))-
                              nrow(unique_all[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Fet_exported$geneID)))-
                              nrow(unique_all[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(frac.gene[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG")
}
names(frac.gene) = group
fisher.frac.gene = lapply(frac.gene, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
# adultOnly prenatalOnly    ACnotAN      ANnotAC PCnotPN     PNnotPC   ACnotPC   PCnotAC     ANnotPN     PNnotAN
# bothagesDEG 1.0000000     1.000000 1.00000000 0.5880699863       1 0.588069986 1.0000000 1.0000000 0.209325682 0.347362128
# adultDEG    0.5133018     0.267304 0.01968468 0.0000379302       1 0.001106946 0.6018692 0.2511116 0.002154625 0.001230031
# prenatalDEG 1.0000000     1.000000 1.00000000 0.5875735757       1 0.587573576 1.0000000 1.0000000 0.213019020 0.353090184
frac.gene.props = lapply(frac.gene, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



## Are editing sites specific to a group and shared by all samples in a group more or less likely to fall in an increasing or decreasing genes by age than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
age.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  age.site[[i]] = list(
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])-
                              nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(age.site[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(age.site) = group
fisher.age.site = lapply(age.site, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                    adultOnly prenatalOnly   ACnotAN   ANnotAC PCnotPN   PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#bothfractionsDEG 5.137825e-10 6.996984e-06 0.7623512 1.0000000       1 0.1032770 1.150492e-13 2.000945e-06 2.222938e-13 5.013174e-09
#CytosolDEG       5.831627e-12 2.655909e-06 0.4651877 0.9124600       1 0.0120951 8.649887e-16 1.599953e-07 5.994252e-15 2.532017e-12
#NucleusDEG       2.196999e-10 2.165893e-08 0.2525420 0.3934576       1 0.2781730 3.381863e-11 9.332538e-12 1.154839e-17 1.036203e-13
age.site.props = lapply(age.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))


# number of significantly increasing or decreasing genesby age that either contain an editing site specific to a group and in all samples, or don't
age.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  age.gene[[i]] = list(
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),]),
                            length(unique(age.sig$both_increasing$geneID))-nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),]),
                            length(unique(age.sig$both_decreasing$geneID))-nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),]),
                            sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Cyt_increasing$geneID)))-
                              nrow(unique_all[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),]),
                            sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Cyt_decreasing$geneID)))-
                              nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),]),
                            sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Nuc_increasing$geneID)))-
                              nrow(unique_all[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),])),
               decreasing = c(nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),]),
                            sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Nuc_decreasing$geneID)))-
                              nrow(unique_all[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(age.gene[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(age.gene) = group
fisher.age.gene = lapply(age.gene, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                    adultOnly prenatalOnly   ACnotAN      ANnotAC   PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#bothfractionsDEG 9.382517e-07 3.696026e-10 1.0000000 0.0077734607 0.2137893 1.736459e-05 1.677679e-07 8.024038e-19 5.908792e-10 7.191367e-19
#CytosolDEG       5.071823e-06 1.363399e-10 1.0000000 0.1042451089 0.1098819 3.470043e-07 8.787733e-07 1.978561e-19 5.629772e-06 2.136388e-21
#NucleusDEG       1.861462e-08 1.750548e-10 0.7667066 0.0003625822 0.2000911 2.170476e-04 1.525285e-07 7.439758e-21 3.063336e-14 1.452532e-20
age.gene.props = lapply(age.gene, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### Does editing rate correlate with gene expression in the group the editing site appears?
## correlate LFC with editing rate in editing sites unique to a group

corr.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
unique_all_df = lapply(unique_all, as.data.frame)
for (i in 1:length(unique_all)){
  corr.site[[i]] = list(Adult = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Adult.LFC"], use = "complete.obs"),3),
                             Prenatal = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Prenatal.LFC"], use = "complete.obs"),3),
                             Cytosol = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Cytosol.LFC"], use = "complete.obs"),3),
                             Nucleus = round(cor(x = unique_all_df[[group[i]]][,"rate"], y = unique_all_df[[group[i]]][,"Nucleus.LFC"], use = "complete.obs"),3))
}
names(corr.site) = group
data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))
#         adultOnly prenatalOnly ACnotAN ANnotAC PCnotPN PNnotPC ACnotPC PCnotAC ANnotPN PNnotAN
#Adult        0.141       -0.225  -0.192   0.034   0.138   0.071   0.168  -0.114   0.080  -0.204
#Prenatal     0.221       -0.022  -0.019   0.060   0.044   0.079   0.257   0.009   0.123  -0.022
#Cytosol     -0.416       -0.376  -0.038  -0.014   0.289  -0.032  -0.268  -0.239  -0.190  -0.237
#Nucleus     -0.417       -0.343   0.044  -0.007   0.346  -0.027  -0.276  -0.222  -0.196  -0.167
max(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))) # 0.346
min(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))) # -0.417
min(abs(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F))))) # 0.007


### Is gene expression greater in the compartment/age exhibiting the editing site than in the compared group?
## t test of LFC between unique sites in both groups

t.frac.site = t.age.site = list(list())
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.frac.site[[i]] = list(adult = t.test(x = unique_all_df[[comps[[i]][1]]][,"Adult.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Adult.LFC"]),
                          prenatal = t.test(x = unique_all_df[[comps[[i]][1]]][,"Prenatal.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Prenatal.LFC"]))
  t.age.site[[i]] = list(cytosol = t.test(x = unique_all_df[[comps[[i]][1]]][,"Cytosol.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Cytosol.LFC"]),
                         nucleus = t.test(x = unique_all_df[[comps[[i]][1]]][,"Nucleus.LFC"], y = unique_all_df[[comps[[i]][2]]][,"Nucleus.LFC"]))
}
names(t.frac.site) = names(t.age.site) = names(comps)
t.LFC.editing = rbind(frac.Tstat = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      age.Tstat = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      frac.pval = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      age.pval = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      frac.confInt = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                      age.confInt = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))),
                      frac.estMeans = data.frame(lapply(t.frac.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))),
                      age.estMeans = data.frame(lapply(t.age.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)



### In intronic editing sites, is the IR ratio greater in the compartment exhibiting the site than the compared group?
## t test of IR ratio of introns with editing sites by group

names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
  IRres[[i]][,"Chr"] = paste0("chr", IRres[[i]][,"Chr"])
}
names(IRres) = shortenedNames
IRres_gr = lapply(IRres, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr",start.field = "Start",end.field = "End",strand.field = "Direction",keep.extra.columns = T))
ov = lapply(IRres_gr, function(x,y) findOverlaps(x, editing_anno_gr))
togintron = list()
for (i in 1:length(IRres)){
  tmp = IRres[[i]]
  togintron[[i]] = cbind(editing_anno[subjectHits(ov[[i]]),], IRratio = tmp$IRratio[queryHits(ov[[i]])],
                   sampleIDintron = names(IRres)[i])
}
togintron = do.call(rbind, togintron)
IRratio_editing = lapply(unique_all_df, function(x) togintron[which(togintron$editingID %in% x$editingID),])

t.IRratio = list()
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.IRratio[[i]] = t.test(x = IRratio_editing[[comps[[i]][1]]][,"IRratio"], y = IRratio_editing[[comps[[i]][2]]][,"IRratio"])
}
names(t.IRratio) = names(comps)
t.IRratio.editing = rbind(Tstat = data.frame(lapply(t.IRratio, function(x) x$statistic)), pval = data.frame(lapply(t.IRratio, function(x) x$p.value)),
                          confInt = data.frame(lapply(t.IRratio, function(x) x$conf.int)), estMeans = data.frame(lapply(t.IRratio, function(x) x$estimate)))
write.csv(t.IRratio.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.IRratio.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)



### in 3'UTR editing sites, is the exon differentially expressed by group?

tog.3UTR = tog[which(tog$UTR3=="UTR3"),]
counts3UTR = exonCounts.down[which(rownames(exonCounts.down) %in% tog.3UTR$exonID), grep("polyA", colnames(exonCounts.down))]
match(rownames(pd[grep("polyA", rownames(pd)),]), colnames(counts3UTR))

## DESeq2 on exons, is the exon the last exon with the greatest count
exonsA = DESeqDataSetFromMatrix(countData = counts3UTR[,-grep("53", colnames(counts3UTR))], 
                                colData = pd[which(pd$Fetal == "Adult" & pd$Library =="polyA"),], design = ~ Zone)
exonsF = DESeqDataSetFromMatrix(countData = counts3UTR[,grep("53", colnames(counts3UTR))], 
                                colData = pd[which(pd$Fetal == "Prenatal" & pd$Library =="polyA"),], design = ~ Zone)
exonsC = DESeqDataSetFromMatrix(countData = counts3UTR[,grep("C", colnames(counts3UTR))], 
                                colData = pd[which(pd$Zone == "Cytosol" & pd$Library =="polyA"),], design = ~ Fetal)
exonsN = DESeqDataSetFromMatrix(countData = counts3UTR[,grep("N", colnames(counts3UTR))], 
                                colData = pd[which(pd$Zone == "Nucleus" & pd$Library =="polyA"),], design = ~ Fetal)
exonsA = DESeq(exonsA)
exonsF = DESeq(exonsF)
exonsC = DESeq(exonsC)
exonsN = DESeq(exonsN)
exonres = list(exonsA.res = results(exonsA), exonsF.res = results(exonsF), exonsC.res = results(exonsC), exonsN.res = results(exonsN))
elementNROWS(exonres)
elementNROWS(lapply(exonres, function(x) x[which(x$padj<=0.05),]))
tog.3UTR = cbind(tog.3UTR, A.LFC = exonsA.res[match(tog.3UTR$exonID, rownames(exonsA.res)),"log2FoldChange"],
                 A.padj = exonsA.res[match(tog.3UTR$exonID, rownames(exonsA.res)),"padj"],
                 P.LFC = exonsF.res[match(tog.3UTR$exonID, rownames(exonsF.res)),"log2FoldChange"],
                 P.padj = exonsF.res[match(tog.3UTR$exonID, rownames(exonsF.res)),"padj"],
                 C.LFC = exonsC.res[match(tog.3UTR$exonID, rownames(exonsC.res)),"log2FoldChange"],
                 C.padj = exonsC.res[match(tog.3UTR$exonID, rownames(exonsC.res)),"padj"],
                 N.LFC = exonsN.res[match(tog.3UTR$exonID, rownames(exonsN.res)),"log2FoldChange"],
                 N.padj = exonsN.res[match(tog.3UTR$exonID, rownames(exonsN.res)),"padj"])


## Compare LFC and significance of groups of editing sites in 3'UTR

uniqueAll_3UTR = lapply(unique_all, function(x) tog.3UTR[which(tog.3UTR$editingID %in% x$editingID),])
t.3UTR.site = list(list())
comps = list(byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.3UTR.site[[i]] = list(adult = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"A.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"P.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"C.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = uniqueAll_3UTR[[comps[[i]][1]]][,"N.LFC"], y = uniqueAll_3UTR[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.3UTR.site) = names(comps)
t.LFC.3UTR.editing = rbind(Tstat = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      pval = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      confInt = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                      estMeans = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.3UTR.editing, 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.3UTR.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.3UTR.site = list(list())
group = c("adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
for (i in 1:length(group)){
  corr.3UTR.site[[i]] = list(fraction = cor(x = uniqueAll_3UTR[[group[i]]][,"A.LFC"], y = uniqueAll_3UTR[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = uniqueAll_3UTR[[group[i]]][,"C.LFC"], y = uniqueAll_3UTR[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.3UTR.site) = group
data.frame(lapply(corr.3UTR.site, function(x) unlist(x, recursive=F)))
#         adultOnly prenatalOnly   ACnotAN   ANnotAC   PCnotPN   PNnotPC   ACnotPC   PCnotAC   ANnotPN   PNnotAN
#fraction 0.8730890    0.4996052 0.3140912 0.7291933 0.2640446 0.4316696 0.8124232 0.6049925 0.8140242 0.5064677
#age      0.9043124    0.9664685 0.9115609 0.9598261 0.9853123 0.8788988 0.8788835 0.9607036 0.9293007 0.9573088


## Compare the number of edited 3'UTRs that are significantly and non-significantly up- and down-expressed per group

fisher.3UTR = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique_all)){
  fisher.3UTR[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_3UTR[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_3UTR[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_3UTR[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_3UTR[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.3UTR[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.3UTR) = group
fisher.3UTR.editing = lapply(fisher.3UTR, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.3UTR.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                   adultOnly prenatalOnly   ACnotAN    ANnotAC PCnotPN   PNnotPC    ACnotPC    PCnotAC    ANnotPN      PNnotAN
#bothagesDEG      0.08259109    1.0000000 1.0000000 1.00000000       1 1.0000000 0.05520353 1.00000000 0.02656958 1.000000e+00
#AdultDEG         1.00000000    0.3106864 0.2854701 0.30098834       1 0.5764706 0.70051452 0.58521437 0.57841961 7.833541e-01
#PrenatalDEG      0.23972603    1.0000000 1.0000000 1.00000000       1 1.0000000 0.09005497 1.00000000 0.02875768 1.000000e+00
#bothfractionsDEG 0.43214195    1.0000000 1.0000000 1.00000000       1 0.2727273 0.57219722 1.00000000 0.68865682 1.156235e-02
#CytosolDEG       0.12935880    1.0000000 0.7063310 0.06733649       1 0.2450980 0.24832078 1.00000000 0.41151685 6.845838e-02
#NucleusDEG       0.06206983    1.0000000 0.6945843 0.14689203       1 0.2745098 0.02171784 0.02936896 0.02182513 4.126946e-05
fisher.3UTR.props = lapply(fisher.3UTR, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### in 3'UTR editing sites, is the exon differentially expressed by group?

tog.CDS = tog[which(tog$cds=="CDS"),]
countsCDS = exonCounts.down[which(rownames(exonCounts.down) %in% tog.CDS$exonID), grep("polyA", colnames(exonCounts.down))]
match(rownames(pd[grep("polyA", rownames(pd)),]), colnames(countsCDS))

## DESeq2 on exons, is the exon the last exon with the greatest count
exonsA = DESeqDataSetFromMatrix(countData = countsCDS[,-grep("53", colnames(countsCDS))], 
                                colData = pd[which(pd$Fetal == "Adult" & pd$Library =="polyA"),], design = ~ Zone)
exonsF = DESeqDataSetFromMatrix(countData = countsCDS[,grep("53", colnames(countsCDS))], 
                                colData = pd[which(pd$Fetal == "Prenatal" & pd$Library =="polyA"),], design = ~ Zone)
exonsC = DESeqDataSetFromMatrix(countData = countsCDS[,grep("C", colnames(countsCDS))], 
                                colData = pd[which(pd$Zone == "Cytosol" & pd$Library =="polyA"),], design = ~ Fetal)
exonsN = DESeqDataSetFromMatrix(countData = countsCDS[,grep("N", colnames(countsCDS))], 
                                colData = pd[which(pd$Zone == "Nucleus" & pd$Library =="polyA"),], design = ~ Fetal)
exonsA = DESeq(exonsA)
exonsF = DESeq(exonsF)
exonsC = DESeq(exonsC)
exonsN = DESeq(exonsN)
exonres = list(exonsA.res = results(exonsA), exonsF.res = results(exonsF), exonsC.res = results(exonsC), exonsN.res = results(exonsN))
elementNROWS(exonres)
elementNROWS(lapply(exonres, function(x) x[which(x$padj<=0.05),]))
tog.CDS = cbind(tog.CDS, A.LFC = exonsA.res[match(tog.CDS$exonID, rownames(exonsA.res)),"log2FoldChange"],
                A.padj = exonsA.res[match(tog.CDS$exonID, rownames(exonsA.res)),"padj"],
                P.LFC = exonsF.res[match(tog.CDS$exonID, rownames(exonsF.res)),"log2FoldChange"],
                P.padj = exonsF.res[match(tog.CDS$exonID, rownames(exonsF.res)),"padj"],
                C.LFC = exonsC.res[match(tog.CDS$exonID, rownames(exonsC.res)),"log2FoldChange"],
                C.padj = exonsC.res[match(tog.CDS$exonID, rownames(exonsC.res)),"padj"],
                N.LFC = exonsN.res[match(tog.CDS$exonID, rownames(exonsN.res)),"log2FoldChange"],
                N.padj = exonsN.res[match(tog.CDS$exonID, rownames(exonsN.res)),"padj"],
                A.basemean = exonsA.res[match(tog.CDS$exonID, rownames(exonsA.res)),"baseMean"],
                P.basemean = exonsF.res[match(tog.CDS$exonID, rownames(exonsF.res)),"baseMean"],
                C.basemean = exonsC.res[match(tog.CDS$exonID, rownames(exonsC.res)),"baseMean"],
                N.basemean = exonsN.res[match(tog.CDS$exonID, rownames(exonsN.res)),"baseMean"])


## Compare LFC and significance of groups of editing sites in 3'UTR

uniqueAll_CDS = lapply(unique_all, function(x) tog.CDS[which(tog.CDS$editingID %in% x$editingID),])
t.CDS.site = list(list())
comps = list(byFractionInAdult = c("ACnotAN","ANnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN")) # not enough observations for byAge = c("adultOnly","prenatalOnly"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), byAgeInCytosol = c("ACnotPC","PCnotAC")
for (i in 1:length(comps)){
  t.CDS.site[[i]] = list(adult = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"A.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"P.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"C.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = uniqueAll_CDS[[comps[[i]][1]]][,"N.LFC"], y = uniqueAll_CDS[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.CDS.site) = names(comps)
t.LFC.CDS.editing = rbind(Tstat = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                           pval = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                           confInt = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                           estMeans = data.frame(lapply(t.CDS.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.CDS.editing, 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.CDS.LFC.between.unique.editing.sites.shared.by.all.in.group.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.CDS.site = list(list())
for (i in 1:length(group)){
  corr.CDS.site[[i]] = list(fraction = cor(x = uniqueAll_CDS[[group[i]]][,"A.LFC"], y = uniqueAll_CDS[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = uniqueAll_CDS[[group[i]]][,"C.LFC"], y = uniqueAll_CDS[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.CDS.site) = group
data.frame(lapply(corr.CDS.site, function(x) unlist(x, recursive=F)))
#         adultOnly ACnotAN   ANnotAC PNnotPC   ACnotPC   ANnotPN PNnotAN
#fraction 0.9463874      NA 0.1542199      NA 0.9444664 0.8277404      NA
#age      0.8829940      NA 0.8139241      NA 0.8841078 0.6899047      NA


## Compare the number of edited CDSs that are significantly and non-significantly up- and down-expressed per group

fisher.CDS = list(list(),list(),list(),list(),list(),list(),list())
group = c("adultOnly","ACnotAN","ANnotAC","PNnotPC","ACnotPC","ANnotPN","PNnotAN") #not enough obs for "prenatalOnly","PCnotAC","PCnotPN"
for (i in 1:length(group)){
  fisher.CDS[[i]] = list(
    data.frame(exported = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(uniqueAll_CDS[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(uniqueAll_CDS[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(uniqueAll_CDS[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(uniqueAll_CDS[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(uniqueAll_CDS[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.CDS[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.CDS) = group
fisher.CDS.editing = lapply(fisher.CDS, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.CDS.editing, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                  adultOnly ACnotAN    ANnotAC PNnotPC    ACnotPC    ANnotPN PNnotAN
#bothagesDEG      0.10000000       1 1.00000000       1 0.10000000 0.06666667       1
#AdultDEG         0.04274402       1 1.00000000       1 0.01805986 0.06636156       1
#PrenatalDEG      1.00000000       1 1.00000000       1 1.00000000 0.49090909       1
#bothfractionsDEG 0.30769231       1 0.33333333       1 0.40000000 0.35294118       1
#CytosolDEG       0.47058824       1 0.04761905       1 1.00000000 0.10521739       1
#NucleusDEG       0.23529412       1 1.00000000       1 0.31578947 0.12000000       1
fisher.CDS.props = lapply(fisher.CDS, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



#showing:
# that genes w editing site in a group affect expression of that gene in that group
# # of editing sites per gene
# annotation of distribution of sites per group
# gene set enrichment of edited genes per group
# mechanisms behind results
#questions:
# is nuclear editing sequestering transcripts
#are nuclear and cytoplasmic editing serving different purposes, or reflecting different strategies
# is this affected by age?

## Are cytosolic/nuclear/adult/prenatal-specific editing sites enriched for DEG Fraction/Age overall, and broken down by group?
# retained or exported DEG and presence or absence of cytosolic-specific editing site
# sites within retained or exported DEG and cytosolic-specific or non-specific status

### Check gene enrichment by annotation,
# retained or exported DEG and presence or absence of cytosolic-specific editing site
# Is a nuclear-specific editing site higher expressed in nucleus, for both intronic and other annotations?
# Is there a relationship between having an editing site and expression by fraction and age?
# is this affected by where in the gene the editing site falls (annotation)?
