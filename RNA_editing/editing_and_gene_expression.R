library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/unique_editingSites_bySample.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


### Isolate the sites unique to each group

unique = lapply(unique_bySamp, function(x) editing_anno[which(editingID %in% x$EditingID),,])
elementNROWS(unique)
unique = Map(cbind, unique, ensID = lapply(unique, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
                 Type = lapply(unique, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]),
                 baseMean = lapply(unique, function(x) Ipres.down[match(x$nearestID, rownames(Ipres.down)),"baseMean"]), 
                 P.LFC = lapply(unique, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "log2FoldChange"]), 
                 P.SE = lapply(unique, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "lfcSE"]), 
                 P.padj = lapply(unique, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "padj"]),
                 A.LFC = lapply(unique, function(x) Apres[match(x$nearestID, rownames(Apres)), "log2FoldChange"]), 
                 A.SE = lapply(unique, function(x) Apres[match(x$nearestID, rownames(Apres)), "lfcSE"]),
                 A.padj = lapply(unique, function(x) Apres[match(x$nearestID, rownames(Apres)), "padj"]),
                 N.LFC = lapply(unique, function(x) Npres[match(x$nearestID, rownames(Npres)), "log2FoldChange"]), 
                 N.SE = lapply(unique, function(x) Npres[match(x$nearestID, rownames(Npres)), "lfcSE"]), 
                 N.padj = lapply(unique, function(x) Npres[match(x$nearestID, rownames(Npres)), "padj"]),
                 C.LFC = lapply(unique, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "log2FoldChange"]), 
                 C.SE = lapply(unique, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "lfcSE"]),
                 C.padj = lapply(unique, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "padj"]))
DEGnames = c(lapply(sig[1:8], function(x) as.character(x$geneID)), lapply(age.sig, function(x) as.character(x$geneID))) 
unique = Map(cbind, unique, both_retained = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$both_retained, "both_retained", "no")),
                 both_exported = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$both_exported, "both_exported", "no")),
                 Fet_retained = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Fet_retained, "Fet_retained", "no")),
                 Ad_retained = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Ad_retained, "Ad_retained", "no")),
                 Fet_exported = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Fet_exported, "Fet_exported", "no")),
                 Ad_exported = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Ad_exported, "Ad_exported", "no")),
                 ret_Ad_exp_Fet = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$ret_Ad_exp_Fet, "ret_Ad_exp_Fet", "no")),
                 ret_Fet_exp_Ad = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$ret_Fet_exp_Ad, "ret_Fet_exp_Ad", "no")),
                 interacting = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$interacting, "interacting", "no")),
                 both_decreasing = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$both_decreasing, "both_decreasing", "no")),
                 both_increasing = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$both_increasing, "both_increasing", "no")),
                 Cyt_decreasing = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_decreasing, "Cyt_decreasing", "no")),
                 Nuc_decreasing = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_decreasing, "Nuc_decreasing", "no")),
                 Cyt_increasing = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Cyt_increasing, "Cyt_increasing", "no")),
                 Nuc_increasing = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$Nuc_increasing, "Nuc_increasing", "no")),
                 decr_Nuc_incr_Cyt = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$decr_Nuc_incr_Cyt, "decr_Nuc_incr_Cyt", "no")),
                 decr_Cyt_incr_Nuc = lapply(unique, function(x) ifelse(x$nearestID %in% DEGnames$decr_Cyt_incr_Nuc, "decr_Cyt_incr_Nuc", "no")))
elementNROWS(unique)
lapply(unique, head)

all = Map(cbind, all, ensID = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"ensemblID"]),
          Type = lapply(all, function(x) geneMap[match(x$nearestID,geneMap$gencodeID),"gene_type"]),
          baseMean = lapply(all, function(x) Ipres.down[match(x$nearestID, rownames(Ipres.down)),"baseMean"]), 
          P.LFC = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "log2FoldChange"]), 
          P.SE = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "lfcSE"]), 
          P.padj = lapply(all, function(x) Fpres.down[match(x$nearestID, rownames(Fpres.down)), "padj"]),
          A.LFC = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "log2FoldChange"]), 
          A.SE = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "lfcSE"]),
          A.padj = lapply(all, function(x) Apres[match(x$nearestID, rownames(Apres)), "padj"]),
          N.LFC = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "log2FoldChange"]), 
          N.SE = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "lfcSE"]), 
          N.padj = lapply(all, function(x) Npres[match(x$nearestID, rownames(Npres)), "padj"]),
          C.LFC = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "log2FoldChange"]), 
          C.SE = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "lfcSE"]),
          C.padj = lapply(all, function(x) Cpres.down[match(x$nearestID, rownames(Cpres.down)), "padj"]),
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
numsites_bygene = lapply(unique, function(x) x[,length(unique(editingID)),by="nearestID"])
numsites_bygene_all = lapply(all, function(x) x[,length(unique(editingID)),by="nearestID"])
numsites_bygene = lapply(numsites_bygene,as.data.frame)
numsites_bygene_all = lapply(numsites_bygene_all, as.data.frame)
elementNROWS(numsites_bygene)
elementNROWS(numsites_bygene_all)

## Compare to the number of sites per gene in sites present in each compartment but not unique, and the total number of sites in the same genes as those of the unique sites
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
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
#                   cytosolOnly  nucleusOnly   adultOnly prenatalOnly      ACnotAN      ANnotAC      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#Tstat              1.082000e+01 2.037100e+01 1.55920e+01 4.953000e+00 6.613000e+00 1.781800e+01 7.597000e+00 1.182600e+01 1.019100e+01 4.822000e+00 1.571700e+01  2.968000000
#pval               1.037003e-26 4.262358e-87 1.09945e-49 8.796764e-07 5.324119e-11 1.357331e-66 5.046285e-14 3.891695e-31 2.650874e-22 1.852096e-06 5.555072e-49  0.003129269
#confInt.1          6.460000e-01 1.742000e+00 1.79100e+00 4.210000e-01 4.020000e-01 1.808000e+00 4.620000e-01 8.650000e-01 1.090000e+00 4.220000e-01 1.778000e+00  0.140000000
#confInt.2          9.330000e-01 2.114000e+00 2.30700e+00 9.740000e-01 7.420000e-01 2.256000e+00 7.840000e-01 1.209000e+00 1.610000e+00 1.002000e+00 2.286000e+00  0.686000000
#estMeans.mean of x 2.361000e+00 3.383000e+00 3.73500e+00 2.995000e+00 2.199000e+00 3.447000e+00 2.199000e+00 2.555000e+00 2.941000e+00 2.717000e+00 3.577000e+00  2.650000000
#estMeans.mean of y 1.572000e+00 1.455000e+00 1.68600e+00 2.297000e+00 1.627000e+00 1.415000e+00 1.576000e+00 1.519000e+00 1.591000e+00 2.005000e+00 1.545000e+00  2.238000000
#max                1.900000e+01 5.100000e+01 5.10000e+01 4.200000e+01 1.500000e+01 4.300000e+01 1.800000e+01 2.700000e+01 2.800000e+01 3.600000e+01 4.800000e+01 30.000000000
#median             1.000000e+00 2.000000e+00 2.00000e+00 1.000000e+00 1.000000e+00 2.000000e+00 1.000000e+00 2.000000e+00 1.000000e+00 1.000000e+00 2.000000e+00  1.000000000
#sd                 2.220000e+00 4.202000e+00 5.14500e+00 3.973000e+00 1.932000e+00 4.309000e+00 2.023000e+00 2.597000e+00 3.573000e+00 3.467000e+00 4.812000e+00  3.234000000
t.editedgenes.numsites.df = rbind(Tstat = data.frame(lapply(t.editedgenes.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedgenes.numsites, function(x) x$p.value)),
                                  confInt = data.frame(lapply(t.editedgenes.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedgenes.numsites, function(x) x$estimate)),
                                  max = mapply(function(x,y) max(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup),
                                  median = mapply(function(x,y) median(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup),
                                  sd = mapply(function(x,y) sd(numsites_bygene_all[[y]][-which(numsites_bygene_all[[y]][,"nearestID"] %in% numsites_bygene[[x]][,"nearestID"]),"V1"]), group, allgroup))
#                    cytosolOnly   nucleusOnly    adultOnly prenatalOnly      ACnotAN      ANnotAC      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#Tstat              2.227659e+01  2.597561e+01 1.953522e+01 1.156446e+01 1.657753e+01 2.217149e+01 1.686707e+01 1.877149e+01 1.371927e+01 8.996844e+00 1.909397e+01 8.733816e+00
#pval               7.439237e-98 1.210852e-134 1.518419e-77 1.199959e-29 1.504480e-56 4.564877e-99 6.698610e-58 2.468235e-72 5.014523e-39 1.823673e-18 2.711112e-72 1.256083e-17
#confInt.1          3.387897e+00  3.606980e+00 2.887856e+00 1.625644e+00 2.633570e+00 3.267953e+00 2.404351e+00 2.361264e+00 1.888094e+00 1.159525e+00 2.669316e+00 1.091598e+00
#confInt.2          4.042022e+00  4.195976e+00 3.532414e+00 2.289795e+00 3.340488e+00 3.902115e+00 3.037256e+00 2.912224e+00 2.518475e+00 1.806760e+00 3.280602e+00 1.724431e+00
#estMeans.mean of x 5.286683e+00  5.356375e+00 4.895850e+00 4.254455e+00 4.613535e+00 5.000000e+00 4.296610e+00 4.155541e+00 3.794194e+00 3.487688e+00 4.519787e+00 3.645562e+00
#estMeans.mean of y 1.571723e+00  1.454897e+00 1.685714e+00 2.296736e+00 1.626506e+00 1.414966e+00 1.575806e+00 1.518797e+00 1.590909e+00 2.004545e+00 1.544828e+00 2.237548e+00
#max                1.500000e+01  1.700000e+01 6.000000e+00 1.300000e+01 1.600000e+01 1.600000e+01 1.200000e+01 1.500000e+01 8.000000e+00 1.000000e+01 6.000000e+00 1.300000e+01
#median             1.000000e+00  1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
#sd                 1.487991e+00  1.299340e+00 1.204571e+00 2.013590e+00 1.554721e+00 1.242305e+00 1.333935e+00 1.301030e+00 1.197020e+00 1.729406e+00 9.645041e-01 1.870080e+00


## What about at the exon level?

exonMap$exonID = rownames(exonMap)
exonMap_gr = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
editing_anno_gr = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
ov = findOverlaps(editing_anno_gr, exonMap_gr)
tog = cbind(editing_anno[queryHits(ov),], exonMap[subjectHits(ov),])
unique_exons = lapply(unique, function(x) tog[(editingID %in% x$editingID),,])
all_exons = lapply(all, function(x) tog[(editingID %in% x$editingID),,])

numsites_byexon = lapply(unique_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon_all = lapply(all_exons, function(x) x[,length(unique(editingID)),by="exonID"])
numsites_byexon = lapply(numsites_byexon,as.data.frame)
numsites_byexon_all = lapply(numsites_byexon_all, as.data.frame)
elementNROWS(numsites_byexon)
elementNROWS(numsites_byexon_all)

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
#                    cytosolOnly  nucleusOnly    adultOnly prenatalOnly      ACnotAN      ANnotAC      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN     PNnotAN
#Tstat              1.129400e+01 1.800700e+01 1.424000e+01    0.3970000 9.041000e+00 1.673000e+01 6.770000e+00 1.046800e+01 9.638000e+00  2.773000000 1.286600e+01 -1.68100000
#pval               7.021997e-29 9.704584e-69 4.599877e-44    0.6917896 3.745307e-19 2.216909e-59 1.749433e-11 5.726088e-25 5.569381e-21  0.005686943 9.356659e-36  0.09327819
#confInt.1          5.490000e-01 1.108000e+00 1.150000e+00   -0.1680000 4.700000e-01 1.266000e+00 3.370000e-01 5.570000e-01 7.940000e-01  0.086000000 1.060000e+00 -0.39000000
#confInt.2          7.800000e-01 1.378000e+00 1.517000e+00    0.2520000 7.300000e-01 1.602000e+00 6.110000e-01 8.140000e-01 1.199000e+00  0.501000000 1.442000e+00  0.03000000
#estMeans.mean of x 2.040000e+00 2.504000e+00 2.788000e+00    2.0950000 1.982000e+00 2.694000e+00 1.908000e+00 2.028000e+00 2.398000e+00  2.039000000 2.675000e+00  1.86800000
#estMeans.mean of y 1.375000e+00 1.261000e+00 1.455000e+00    2.0520000 1.382000e+00 1.260000e+00 1.434000e+00 1.343000e+00 1.402000e+00  1.746000000 1.424000e+00  2.04800000
#max                1.900000e+01 3.200000e+01 4.800000e+01   42.0000000 1.900000e+01 3.700000e+01 1.600000e+01 1.900000e+01 3.400000e+01 36.000000000 4.500000e+01 28.00000000
#median             1.000000e+00 1.000000e+00 1.000000e+00    1.0000000 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00  1.000000000 1.000000e+00  1.00000000
#sd                 1.966000e+00 3.067000e+00 4.058000e+00    2.8080000 1.837000e+00 3.310000e+00 1.763000e+00 1.949000e+00 3.137000e+00  2.607000000 3.759000e+00  2.22200000
t.editedexons.numsites.df = rbind(Tstat = data.frame(lapply(t.editedexons.numsites, function(x) x$statistic)), pval = data.frame(lapply(t.editedexons.numsites, function(x) x$p.value)),
                                  confInt = data.frame(lapply(t.editedexons.numsites, function(x) x$conf.int)), estMeans = data.frame(lapply(t.editedexons.numsites, function(x) x$estimate)),
                                  max = mapply(function(x,y) max(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                                  median = mapply(function(x,y) median(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup),
                                  sd = mapply(function(x,y) sd(numsites_byexon_all[[y]][-which(numsites_byexon_all[[y]][,"exonID"] %in% numsites_byexon[[x]][,"exonID"]),"V1"]), group, allgroup))
#                    cytosolOnly   nucleusOnly    adultOnly prenatalOnly      ACnotAN      ANnotAC      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#Tstat              2.040378e+01  2.292049e+01 1.792096e+01 8.340680e+00 1.648868e+01 2.045878e+01 1.520728e+01 1.715550e+01 1.306292e+01 8.051233e+00 1.646333e+01 5.332036e+00
#pval               1.495286e-82 2.820381e-105 2.623490e-68 1.320516e-16 3.116779e-55 2.107017e-84 1.278989e-47 6.325068e-60 4.068688e-37 1.848075e-15 1.895994e-57 1.147852e-07
#confInt.1          2.973827e+00  2.889545e+00 2.087489e+00 8.776632e-01 2.567003e+00 2.762602e+00 1.964866e+00 1.967063e+00 1.470248e+00 7.631079e-01 1.836466e+00 4.392564e-01
#confInt.2          3.606387e+00  3.430232e+00 2.600397e+00 1.417254e+00 3.260386e+00 3.348405e+00 2.546999e+00 2.475003e+00 1.989789e+00 1.254803e+00 2.333139e+00 9.506463e-01
#estMeans.mean of x 4.665565e+00  4.420959e+00 3.798754e+00 3.199795e+00 4.295660e+00 4.315551e+00 3.689834e+00 3.563953e+00 3.131593e+00 2.754464e+00 3.509045e+00 2.742690e+00
#estMeans.mean of y 1.375458e+00  1.261071e+00 1.454810e+00 2.052336e+00 1.381966e+00 1.260047e+00 1.433902e+00 1.342920e+00 1.401575e+00 1.745509e+00 1.424242e+00 2.047739e+00
#max                1.700000e+01  1.700000e+01 6.000000e+00 1.500000e+01 1.600000e+01 1.600000e+01 1.200000e+01 1.500000e+01 1.000000e+01 8.000000e+00 6.000000e+00 1.300000e+01
#median             1.000000e+00  1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
#sd                 1.220897e+00  9.434672e-01 9.688867e-01 1.991807e+00 1.155724e+00 9.771287e-01 1.251851e+00 1.011394e+00 1.164426e+00 1.514066e+00 9.003001e-01 1.846839e+00



### Compare the proportion of editing site in each annotation in sites unique to a specific group and those that aren't
# Using whether the site overlaps an annotation group at all rather than the CDS -> UTR -> Intron -> Other heirarchy

anno.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
for (i in 1:length(unique)){
  anno.site[[i]] = list(
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("CDS", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("CDS", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("CDS", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("CDS", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("CDS", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("Intron", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("Intron", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("Intron", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("Intron", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("Intron", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("UTR5", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("UTR5", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("UTR5", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("UTR5", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("UTR5", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][grep("UTR3", anno),list(unique(editingID)),]), nrow(unique[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][grep("UTR3", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][grep("UTR3", anno),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][-grep("UTR3", anno),list(unique(editingID)),])-nrow(unique[[group[i]]][-grep("UTR3", anno),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")),
    data.frame(Unique = c(nrow(unique[[group[i]]][(annotation=="Other"),list(unique(editingID)),]), nrow(unique[[group[i]]][(annotation!="Other"),list(unique(editingID)),])),
               notUnique = c(nrow(all[[allgroup[i]]][(annotation=="Other"),list(unique(editingID)),])-nrow(unique[[group[i]]][(annotation=="Other"),list(unique(editingID)),]),
                             nrow(all[[allgroup[i]]][(annotation!="Other"),list(unique(editingID)),])-nrow(unique[[group[i]]][(annotation!="Other"),list(unique(editingID)),])), row.names = c("inAnno","notInAnno")))
  names(anno.site[[i]]) = c("CDS","Intron","UTR5","UTR3","Other")
}
names(anno.site) = group
fisher.anno.site = lapply(anno.site, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.anno.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#        cytosolOnly   nucleusOnly    adultOnly prenatalOnly      ACnotAN      ANnotAC      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#CDS    7.479212e-09  9.405109e-15 7.640423e-01 4.052958e-01 7.396955e-06 2.295851e-11 4.904124e-08 3.257281e-10 2.589635e-01 4.558662e-01 3.486063e-01 5.587548e-01
#Intron 3.448633e-01  9.082804e-49 5.026308e-01 2.981328e-07 3.900058e-01 2.494848e-40 9.781940e-01 2.658615e-20 7.697652e-01 6.756008e-08 5.847943e-01 2.231538e-03
#UTR5   5.218612e-01  1.206456e-01 6.383796e-02 4.035939e-04 2.957170e-01 1.183810e-01 7.895968e-01 7.367531e-01 6.147130e-01 7.729938e-02 1.279998e-01 1.855744e-03
#UTR3   1.751213e-10 6.379695e-135 1.415903e-30 5.001359e-61 5.254058e-07 1.221598e-98 6.873078e-05 1.496410e-44 2.251906e-06 2.422827e-29 2.486909e-25 9.635696e-45
#Other  1.991734e-08  7.062937e-05 3.416870e-12 9.000582e-11 1.130273e-08 1.184954e-04 2.047083e-03 1.129417e-01 4.471699e-05 9.645687e-05 1.437100e-12 3.490934e-12
anno.site.props = lapply(anno.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### Characterize the overlap with differentially expressed genes

## Are editing sites specific to a group more or less likely to fall in an exported or retained gene by fraction than chance?

# number of editing sites falling within significantly retained or imported genes that are group-specific, or those found in the same group that aren't specific
frac.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
for (i in 1:length(unique)){
  frac.site[[i]] = list(
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(editingID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),]),
                            nrow(all[[allgroup[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])-
                              nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(frac.site[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG")
}
names(frac.site) = group
fisher.frac.site = lapply(frac.site, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.frac.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#            cytosolOnly  nucleusOnly adultOnly prenatalOnly      ACnotAN      ANnotAC    PCnotPN      PNnotPC   ACnotPC      PCnotAC   ANnotPN   PNnotAN
#bothagesDEG 0.276452942 1.955890e-02 0.6962943    1.0000000 0.4254588350 6.685020e-03 0.54923469 4.330159e-01 0.4787456 1.727041e-01 0.5969333 1.0000000
#adultDEG    0.000012599 2.497990e-15 0.7880560    0.8303155 0.0000342682 6.056706e-20 0.08030351 1.783635e-06 0.2637201 3.623602e-07 0.4599105 0.3268008
#prenatalDEG 0.285831785 4.355389e-02 1.0000000    1.0000000 0.4483603415 1.492129e-02 0.28787879 5.528326e-01 0.3393518 1.549062e-01 0.5979982 1.0000000
frac.site.props = lapply(frac.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))


# number of significantly retained or exported genes that either contain an editing site specific to a group, or don't
frac.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  frac.gene[[i]] = list(
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),]),
                            length(unique(sig$both_retained$geneID))-nrow(unique[[group[i]]][(both_retained=="both_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),]),
                            length(unique(sig$both_exported$geneID))-nrow(unique[[group[i]]][(both_exported=="both_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Ad_retained$geneID)))-
                              nrow(unique[[group[i]]][(both_retained=="both_retained" | Ad_retained=="Ad_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Ad_exported$geneID)))-
                              nrow(unique[[group[i]]][(both_exported=="both_exported" | Ad_exported=="Ad_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(retained = c(nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_retained$geneID)),length(unique(sig$Fet_retained$geneID)))-
                              nrow(unique[[group[i]]][(both_retained=="both_retained"| Fet_retained=="Fet_retained"),list(unique(nearestID)),])),
               exported = c(nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),]),
                            sum(length(unique(sig$both_exported$geneID)),length(unique(sig$Fet_exported$geneID)))-
                              nrow(unique[[group[i]]][(both_exported=="both_exported"| Fet_exported=="Fet_exported"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(frac.gene[[i]]) = c("bothagesDEG","adultDEG","prenatalDEG")
}
names(frac.gene) = group
fisher.frac.gene = lapply(frac.gene, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.frac.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#            cytosolOnly  nucleusOnly    adultOnly prenatalOnly   ACnotAN      ANnotAC   PCnotPN      PNnotPC    ACnotPC     PCnotAC      ANnotPN      PNnotAN
#bothagesDEG   1.0000000 1.587033e-04 1.649457e-02 0.0364227926 0.7233423 4.783678e-04 1.0000000 1.361014e-02 1.00000000 0.264831900 2.836746e-03 3.584345e-02
#adultDEG      0.5870878 4.332115e-16 6.047586e-12 0.0002520483 0.5136990 1.312404e-16 0.2282471 1.373863e-06 0.01142585 0.002617012 8.948925e-16 5.403304e-06
#prenatalDEG   1.0000000 1.229956e-04 3.833183e-03 0.0234104005 1.0000000 8.478201e-05 1.0000000 1.596584e-02 0.64525008 0.187657965 5.464908e-04 3.528505e-02
frac.gene.props = lapply(frac.gene, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))


## Are editing sites specific to a group more or less likely to fall in an increasing or decreasing genes by age than chance?

# number of editing sites falling within significantly retained or imported genes that are specific and in all samples, or those found in the same group that aren't specific
age.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  age.site[[i]] = list(
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(editingID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),]),
                              nrow(all[[allgroup[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(editingID)),])), row.names = c("unique","notUnique")))
  names(age.site[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(age.site) = group
fisher.age.site = lapply(age.site, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.age.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                  cytosolOnly  nucleusOnly     adultOnly prenatalOnly      ACnotAN   ANnotAC   PCnotPN    PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#bothfractionsDEG 0.0053083187 2.376585e-14 1.181044e-115 2.423966e-56 4.362771e-03 0.5434692 0.5755423 0.06343512 3.307034e-52 9.115417e-29 1.528439e-89 4.108609e-34
#CytosolDEG       0.0253018141 1.226418e-11 1.076269e-105 1.252251e-55 1.305056e-01 0.7960983 0.5482321 0.01624782 2.116403e-46 3.822950e-35 4.220221e-84 2.453798e-29
#NucleusDEG       0.0004868032 5.677654e-15 4.316760e-106 3.415726e-79 3.446854e-05 0.2585892 0.2740360 0.01718529 2.461866e-47 9.046786e-37 1.356967e-84 2.346243e-54
age.site.props = lapply(age.site, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))


# number of significantly increasing or decreasing genesby age that either contain an editing site specific to a group and in all samples, or don't
age.gene = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  age.gene[[i]] = list(
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),]),
                              length(unique(age.sig$both_increasing$geneID))-nrow(unique[[group[i]]][(both_increasing=="both_increasing"),list(unique(nearestID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),]),
                              length(unique(age.sig$both_decreasing$geneID))-nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Cyt_increasing$geneID)))-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing" | Cyt_increasing=="Cyt_increasing"),list(unique(nearestID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Cyt_decreasing$geneID)))-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing" | Cyt_decreasing=="Cyt_decreasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")),
    data.frame(increasing = c(nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_increasing$geneID)),length(unique(age.sig$Nuc_increasing$geneID)))-
                                nrow(unique[[group[i]]][(both_increasing=="both_increasing"| Nuc_increasing=="Nuc_increasing"),list(unique(nearestID)),])),
               decreasing = c(nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),]),
                              sum(length(unique(age.sig$both_decreasing$geneID)),length(unique(age.sig$Nuc_decreasing$geneID)))-
                                nrow(unique[[group[i]]][(both_decreasing=="both_decreasing"| Nuc_decreasing=="Nuc_decreasing"),list(unique(nearestID)),])), row.names = c("ContainsUnique","NoUniquePresent")))
  names(age.gene[[i]]) = c("bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(age.gene) = group
fisher.age.gene = lapply(age.gene, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.age.gene, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F)))
#                  cytosolOnly  nucleusOnly    adultOnly prenatalOnly   ACnotAN      ANnotAC      PCnotPN      PNnotPC      ACnotPC      PCnotAC      ANnotPN      PNnotAN
#bothfractionsDEG 1.071982e-12 1.416789e-04 3.179605e-10 2.219765e-79 0.5113259 1.578216e-06 3.116086e-40 3.687691e-47 2.697564e-06 7.560092e-63 1.156650e-11 1.960528e-66
#CytosolDEG       1.930731e-16 4.453167e-10 2.452306e-04 5.576762e-90 0.6641030 1.629152e-01 2.696275e-49 1.728221e-52 4.118452e-04 2.444984e-74 1.809030e-04 1.431815e-72
#NucleusDEG       4.747446e-10 1.437169e-02 6.882808e-15 4.294430e-75 0.2401903 3.649475e-10 7.579455e-35 2.023333e-42 1.666460e-08 2.548463e-54 3.034939e-16 3.952854e-64
age.gene.props = lapply(age.gene, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))



### Does editing rate correlate with gene expression in the group the editing site appears?
## correlate LFC with editing rate in editing rates shared between fraction

corr.shared = list(list(),list(),list(),list(),list(),list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
allgroup = c("cytosolAll","nucleusAll","adultAll","prenatalAll","allAC","allAN","allPC","allPN","allAC","allPC","allAN","allPN")
unique_df = lapply(unique, as.data.frame)
all_df = lapply(all, as.data.frame)
for (i in 1:length(unique)){
  corr.shared[[i]] = list(Adult = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                            y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"A.LFC"], use = "complete.obs"),3),
                          Prenatal = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                               y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"P.LFC"], use = "complete.obs"),3),
                          Cytosol = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                              y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"C.LFC"], use = "complete.obs"),3),
                          Nucleus = round(cor(x = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"rate"],
                                              y = all_df[[allgroup[i]]][-which(all_df[[allgroup[i]]][,"editingID"] %in% unique_df[[group[i]]][,"editingID"]),"N.LFC"], use = "complete.obs"),3))
}
names(corr.shared) = c("sharedbyFrac:inCyt","sharedbyFrac:inNuc","sharedbyAge:inAd","sharedbyAge:inPren","sharedbyFrac:inAC","sharedbyFrac:inAN","sharedbyFrac:inPC","sharedbyFrac:inPN",
                       "sharedbyAge:inAC","sharedbyAge:inPC","sharedbyAge:inAN","sharedbyAge:inPN")
data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F)))
#         sharedbyFrac.inCyt sharedbyFrac.inNuc sharedbyAge.inAd sharedbyAge.inPren sharedbyFrac.inAC sharedbyFrac.inAN sharedbyFrac.inPC sharedbyFrac.inPN sharedbyAge.inAC sharedbyAge.inPC
#Adult                 0.009              0.004            0.071             -0.006             0.065             0.020            -0.036            -0.056            0.067           -0.018
#Prenatal              0.070              0.039            0.094              0.063             0.074             0.023             0.054             0.020            0.087            0.079
#Cytosol              -0.009              0.015            0.046             -0.036             0.056             0.066            -0.079            -0.074            0.069           -0.017
#Nucleus               0.002              0.019            0.024             -0.008             0.034             0.054            -0.049            -0.041            0.047            0.020
#         sharedbyAge.inAN sharedbyAge.inPN
#Adult               0.067            0.031
#Prenatal            0.069            0.061
#Cytosol             0.049           -0.052
#Nucleus             0.019           -0.059
max(data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F)))) # 0.094
min(data.frame(lapply(corr.shared, function(x) unlist(x, recursive=F)))) # -0.079

## correlate LFC with editing rate in editing sites unique to a group

corr.site = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  corr.site[[i]] = list(Adult = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"A.LFC"], use = "complete.obs"),3),
                        Prenatal = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"P.LFC"], use = "complete.obs"),3),
                        Cytosol = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"C.LFC"], use = "complete.obs"),3),
                        Nucleus = round(cor(x = unique_df[[group[i]]][,"rate"], y = unique_df[[group[i]]][,"N.LFC"], use = "complete.obs"),3))
}
names(corr.site) = group
data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))
#         cytosolOnly nucleusOnly adultOnly prenatalOnly ACnotAN ANnotAC PCnotPN PNnotPC ACnotPC PCnotAC ANnotPN PNnotAN
#Adult         -0.029      -0.005     0.036       -0.107  -0.117   0.052  -0.027  -0.015   0.039  -0.034   0.035  -0.153
#Prenatal       0.065       0.038     0.068        0.001   0.034   0.072   0.086   0.061   0.091   0.059   0.070  -0.028
#Cytosol        0.009      -0.044    -0.058       -0.131   0.079   0.005  -0.078  -0.130  -0.058  -0.112  -0.072  -0.145
#Nucleus        0.034      -0.031    -0.062       -0.089   0.165  -0.002  -0.050  -0.111  -0.060  -0.086  -0.077  -0.079
max(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))) # 0.165
min(data.frame(lapply(corr.site, function(x) unlist(x, recursive=F)))) # -0.153



### Is gene expression greater in the compartment/age exhibiting the editing site than in the compared group?
## t test of LFC between unique sites in both groups

t.site = list(list(),list(),list(),list(),list(),list())
comps = list(byFraction = c("cytosolOnly","nucleusOnly"), byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.site[[i]] = list(adult = t.test(x = unique_df[[comps[[i]][1]]][,"A.LFC"], y = unique_df[[comps[[i]][2]]][,"A.LFC"]),
                     prenatal = t.test(x = unique_df[[comps[[i]][1]]][,"P.LFC"], y = unique_df[[comps[[i]][2]]][,"P.LFC"]),
                     cytosol = t.test(x = unique_df[[comps[[i]][1]]][,"C.LFC"], y = unique_df[[comps[[i]][2]]][,"C.LFC"]),
                     nucleus = t.test(x = unique_df[[comps[[i]][1]]][,"N.LFC"], y = unique_df[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.site) = names(comps)
t.LFC.editing = rbind(Tstat = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                      pval = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                      confInt = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                      estMeans = data.frame(lapply(t.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.LFC.between.editing.sites.unique.to.group.csv", quote=F)



### In intronic editing sites, is the IR ratio greater in the compartment exhibiting the site than the compared group?

names = scan("./Dropbox/BrainRNACompartments/SampleIDs.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0("./Dropbox/sorted_figures/IRfinder/PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
  IRres[[i]][,"Chr"] = paste0("chr", IRres[[i]][,"Chr"])
}
names(IRres) = shortenedNames
IRres_gr = lapply(IRres, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr", start.field = "Start", end.field = "End", strand.field = "Direction", keep.extra.columns = T))
ov = lapply(IRres_gr, function(x) findOverlaps(editing_anno_gr, x))
togintron = list()
for (i in 1:length(IRres)){
  tmp = IRres[[i]]
  togintron[[i]] = cbind(editing_anno[queryHits(ov[[i]]),], IRratio = tmp$IRratio[subjectHits(ov[[i]])],
                   sampleIDintron = names(IRres)[i])
}
togintron = do.call(rbind, togintron)
IRratio_editing = lapply(unique_df, function(x) togintron[which(togintron$editingID %in% x$editingID),])

t.IRratio = list()
for (i in 1:length(comps)){
  t.IRratio[[i]] = t.test(x = IRratio_editing[[comps[[i]][1]]][,"IRratio"], y = IRratio_editing[[comps[[i]][2]]][,"IRratio"])
}
names(t.IRratio) = names(comps)
t.IRratio.editing = rbind(Tstat = data.frame(lapply(t.IRratio, function(x) x$statistic)), pval = data.frame(lapply(t.IRratio, function(x) x$p.value)),
                          confInt = data.frame(lapply(t.IRratio, function(x) x$conf.int)), estMeans = data.frame(lapply(t.IRratio, function(x) x$estimate)))
write.csv(t.IRratio.editing, file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.IRratio.between.unique.editing.sites.to.group.csv", quote=F)



### in 3'UTR editing sites, is the exon differentially expressed by fraction or group?

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
exonsA.res = results(exonsA)
exonsF.res = results(exonsF)
exonsC.res = results(exonsC)
exonsN.res = results(exonsN)
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

unique_3UTR = lapply(unique, function(x) tog.3UTR[which(tog.3UTR$editingID %in% x$editingID),])
t.3UTR.site = list(list())
comps = list(byFraction = c("cytosolOnly","nucleusOnly"), byAge = c("adultOnly","prenatalOnly"), byFractionInAdult = c("ACnotAN","ANnotAC"), byFractionInPrenatal = c("PCnotPN","PNnotPC"), 
             byAgeInCytosol = c("ACnotPC","PCnotAC"), byAgeInNucleus = c("ANnotPN","PNnotAN"))
for (i in 1:length(comps)){
  t.3UTR.site[[i]] = list(adult = t.test(x = unique_3UTR[[comps[[i]][1]]][,"A.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"A.LFC"]),
                          prenatal = t.test(x = unique_3UTR[[comps[[i]][1]]][,"P.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"P.LFC"]),
                          cytosol = t.test(x = unique_3UTR[[comps[[i]][1]]][,"C.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"C.LFC"]),
                          nucleus = t.test(x = unique_3UTR[[comps[[i]][1]]][,"N.LFC"], y = unique_3UTR[[comps[[i]][2]]][,"N.LFC"]))
}
names(t.3UTR.site) = names(comps)
t.LFC.3UTR.editing = rbind(Tstat = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$statistic), recursive=F))),
                           pval = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$p.value), recursive=F))),
                           confInt = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$conf.int), recursive=F))), 
                           estMeans = data.frame(lapply(t.3UTR.site, function(x) unlist(lapply(x, function(y) y$estimate), recursive=F))))
write.csv(t.LFC.3UTR.editing, 
          file = "./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/t.test.of.3UTR.LFC.between.unique.editing.sites.to.group.csv", quote=F)


## Correlate the LFC between age and fraction comparisons in different groups
corr.3UTR.site = list(list())
group = c("cytosolOnly","nucleusOnly","adultOnly","prenatalOnly","ACnotAN","ANnotAC","PCnotPN","PNnotPC","ACnotPC","PCnotAC","ANnotPN","PNnotAN")
for (i in 1:length(group)){
  corr.3UTR.site[[i]] = list(fraction = cor(x = unique_3UTR[[group[i]]][,"A.LFC"], y = unique_3UTR[[group[i]]][,"P.LFC"], use = "complete.obs"),
                             age = cor(x = unique_3UTR[[group[i]]][,"C.LFC"], y = unique_3UTR[[group[i]]][,"N.LFC"], use = "complete.obs"))
}
names(corr.3UTR.site) = group
data.frame(lapply(corr.3UTR.site, function(x) unlist(x, recursive=F)))
#         cytosolOnly nucleusOnly adultOnly prenatalOnly   ACnotAN   ANnotAC   PCnotPN   PNnotPC   ACnotPC   PCnotAC   ANnotPN   PNnotAN
#fraction   0.5231701   0.6678311 0.7007121    0.5023542 0.5429321 0.7208554 0.6083703 0.6499986 0.6971967 0.6105784 0.7143684 0.5119560
#age        0.9358370   0.9437012 0.9376431    0.9601124 0.9132981 0.9278078 0.9308891 0.9425266 0.9397684 0.9501729 0.9393015 0.9601277


## Compare the number of edited 3'UTRs that are significantly and non-significantly up- and down-expressed per group

fisher.3UTR = list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list())
for (i in 1:length(unique)){
  fisher.3UTR[[i]] = list(
    data.frame(exported = c(nrow(unique_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(unique_3UTR[[group[i]]][(A.LFC<0 & P.LFC<0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj<=0.05 & P.padj<=0.05),list(unique(exonID)),]),
                            nrow(unique_3UTR[[group[i]]][(A.LFC>0 & P.LFC>0 & A.padj>0.05 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(unique_3UTR[[group[i]]][(A.LFC<0 & A.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(A.LFC<0 & A.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_3UTR[[group[i]]][(A.LFC>0 & A.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(A.LFC>0 & A.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(exported = c(nrow(unique_3UTR[[group[i]]][(P.LFC<0 & P.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(P.LFC<0 & P.padj>0.05),list(unique(exonID)),])),
               retained = c(nrow(unique_3UTR[[group[i]]][(P.LFC>0 & P.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(P.LFC>0 & P.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(unique_3UTR[[group[i]]][(C.LFC<0 & N.LFC<0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj<=0.05 & N.padj<=0.05),list(unique(exonID)),]),
                              nrow(unique_3UTR[[group[i]]][(C.LFC>0 & N.LFC>0 & C.padj>0.05 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC<0 & C.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(C.LFC<0 & C.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_3UTR[[group[i]]][(C.LFC>0 & C.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(C.LFC>0 & C.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")),
    data.frame(increasing = c(nrow(unique_3UTR[[group[i]]][(N.LFC<0 & N.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(N.LFC<0 & N.padj>0.05),list(unique(exonID)),])),
               decreasing = c(nrow(unique_3UTR[[group[i]]][(N.LFC>0 & N.padj<=0.05),list(unique(exonID)),]),nrow(unique_3UTR[[group[i]]][(N.LFC>0 & N.padj>0.05),list(unique(exonID)),])), row.names = c("sig","non-sig")))
  names(fisher.3UTR[[i]]) = c("bothagesDEG","AdultDEG","PrenatalDEG","bothfractionsDEG","CytosolDEG","NucleusDEG")
}
names(fisher.3UTR) = group
fisher.3UTR.editing = lapply(fisher.3UTR, function(x) lapply(x, fisher.test))
data.frame(lapply(fisher.3UTR.editing, function(x) unlist(lapply(x, function(y) round(y$p.value,3)), recursive=F)))
#                 cytosolOnly nucleusOnly adultOnly prenatalOnly ACnotAN ANnotAC PCnotPN PNnotPC ACnotPC PCnotAC ANnotPN PNnotAN
#bothagesDEG            0.520       0.000     0.000        0.004   0.551   0.023   1.000   0.010   0.002   0.146   0.000   0.004
#AdultDEG               0.005       0.947     0.331        0.097   0.007   0.887   0.010   0.636   0.113   0.321   0.688   0.755
#PrenatalDEG            0.626       0.000     0.000        0.000   0.610   0.017   1.000   0.001   0.000   0.035   0.000   0.000
#bothfractionsDEG       0.078       0.531     0.007        0.000   0.857   0.197   0.007   0.000   0.030   0.000   0.007   0.000
#CytosolDEG             0.001       0.068     0.099        0.000   0.260   0.637   0.000   0.000   0.240   0.000   0.055   0.000
#NucleusDEG             0.141       0.670     0.000        0.000   0.442   0.068   0.003   0.001   0.000   0.000   0.000   0.000
fisher.3UTR.props = lapply(fisher.3UTR, function(y) lapply(y, function(x) 
  c(row1prop = x[1,1]/rowSums(x[1,]), row2prop = x[2,1]/rowSums(x[2,]), col1prop = x[1,1]/sum(x[,1]), col2prop = x[1,2]/sum(x[,2]))))

