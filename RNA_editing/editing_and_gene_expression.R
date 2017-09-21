library(GenomicRanges)
library(data.table)
library(GenomicFeatures)
library(reshape2)
library(ggplot2)
library(plyr)

load("./Dropbox/sorted_figures/new/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
load("./Dropbox/sorted_figures/new/github_controlled/rna_editing/data/rna_editing_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/characterize_fractioned_transcriptome/data/DESeq2_results.rda")
load("./Dropbox/sorted_figures/new/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")

# Explore results
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

# Annotate editing sites to features in the genome
txdb = loadDb("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/SGSeq_out/gencode.v25lift37.annotation.sqlite")
features = list(CDS = cdsBy(txdb, by="tx", use.names=T), Introns = intronsByTranscript(txdb, use.names=T), 
                UTR5 = fiveUTRsByTranscript(txdb, use.names=T), UTR3 = threeUTRsByTranscript(txdb, use.names=T))
features = lapply(features, function(x) unlist(x, recursive = TRUE, use.names = TRUE))
lapply(features, head)

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

## Editing sites present in one group but not another
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
all = list(cytosolAll = cyt, nucleusAll = nuc, adultAll = ad, prenatalAll = pren, allAC = AC, allAN = AN, allPC = PC, allPN = PN)



### Characterize the overlap with differentially expressed genes

editing_ranges = makeGRangesFromDataFrame(editing_anno, keep.extra.columns = T)
DEG_hits = findOverlaps(geneMapGR, reduce(editing_ranges))
editing_anno = as.data.frame(editing_anno)
DEG_editing = lapply(sig, function(x) editing_anno[which(editing_anno$collapsedconversion=="A:G / T:C" & 
                                                           editing_anno$nearestID %in% as.character(x$geneID)),])

## Are cytosolic-specific editing sites enriched for DEG Fraction?
names(unique)
cyt = unique[["cytosolOnly"]]
dim(cyt) #4285
cytonly.deg = lapply(sig, function(x) cyt[which(cyt$collapsedconversion=="A:G / T:C" & cyt$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cytonly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.056533 2.167616
#sample estimates:
#  odds ratio 
#1.508694 
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.206004 1.761161
#sample estimates:
#  odds ratio 
#1.455631
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0429
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5209447 0.9903415
#sample estimates:
#  odds ratio 
#0.7196827 

names(all)
cytosolAll = all[["cytosolAll"]]
dim(cytosolAll) #16513
cytosolAll.deg = lapply(sig, function(x) cytosolAll[which(cytosolAll$collapsedconversion=="A:G / T:C" & cytosolAll$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(cytosolAll.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#425           653            82          1851            76           588 
elementNROWS(lapply(cytonly.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#128           219            34           645            30           212 

fisher.test(data.frame(c(128,425-128), c(219,653-219))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 0.257
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6500174 1.1202317
#sample estimates:
#  odds ratio 
#0.8542065 
fisher.test(data.frame(c(128+645,425+1851-(128+645)), c(219+212,653+588-(219+212)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 0.6555
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8335815 1.1213151
#sample estimates:
#  odds ratio 
#0.966544 
fisher.test(data.frame(c(128+34,425+82-(128+34)), c(219+30,653+76-(219+30)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.4256
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7054002 1.1603049
#sample estimates:
#  odds ratio 
#0.9052591 


## Are nuclear-specific editing sites enriched for DEG Fraction?
nuc = unique[["nucleusOnly"]]
nuconly.deg = lapply(sig, function(x) nuc[which(nuc$collapsedconversion=="A:G / T:C" & nuc$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuconly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#167           112            47           546            19           136

fisher.test(data.frame(c(167,975-167), c(112,1010-112))) # both ages: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value = 0.0001323
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.271346 2.164859
#sample estimates:
#  odds ratio 
#1.656712
fisher.test(data.frame(c(167+546,975+3427-(167+546)), c(112+136,1010+2442-(112+136)))) # adult: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.140855 2.919265
#sample estimates:
#  odds ratio 
#2.496752
fisher.test(data.frame(c(167+47,975+354-(167+47)), c(112+19,1010+350-(112+19)))) # prenatal: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value = 6.311e-07
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.419394 2.289464
#sample estimates:
#  odds ratio 
#1.800193

names(all)
nucleusAll = all[["nucleusAll"]]
dim(nucleusAll) #24407
nucleusAll.deg = lapply(sig, function(x) nucleusAll[which(nucleusAll$collapsedconversion=="A:G / T:C" & nucleusAll$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(nucleusAll.deg[1:6], function(x) unique(x$editingID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported
#1013           775           200          3376           116           812
elementNROWS(lapply(nuconly.deg[1:6], function(x) unique(x$editingID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported
#716           341           152          2170            70           436

fisher.test(data.frame(c(716,1013-716), c(341,775-341))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.509839 3.751007
#sample estimates:
#  odds ratio 
#3.06619
fisher.test(data.frame(c(716+2170,1013+3376-(716+2170)), c(341+436,775+812-(341+436)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
# p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.778337 2.252998
#sample estimates:
#  odds ratio 
#2.001491
fisher.test(data.frame(c(716+152,1013+200-(716+152)), c(341+70,775+116-(341+70)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.440620 3.537815
#sample estimates:
#  odds ratio 
#2.936724 

## Are adult-specific editing sites enriched for DEG Fraction?
ad = unique[["adultOnly"]]
adonly.deg = lapply(sig, function(x) ad[which(ad$collapsedconversion=="A:G / T:C" & ad$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(adonly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#154           127            44           509            18           157 

fisher.test(data.frame(c(154,975-154), c(127,1010-127))) # both ages: retained or exported DEG and presence or absence of adult-specific editing site
#p-value = 0.04576
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.004371 1.695185
#sample estimates:
#  odds ratio 
#1.303989
fisher.test(data.frame(c(154+509,975+3427-(154+509)), c(127+157,1010+2442-(127+157)))) # adult: retained or exported DEG and presence or absence of adult-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.704538 2.299308
#sample estimates:
#  odds ratio 
#1.977819 
fisher.test(data.frame(c(154+44,975+354-(154+44)), c(127+18,1010+350-(127+18)))) # prenatal: retained or exported DEG and presence or absence of adult-specific editing site
#p-value = 0.001182
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.159608 1.858403
#sample estimates:
#  odds ratio 
#1.466712

names(all)
adultAll = all[["adultAll"]]
dim(adultAll) # 23144
adultAll.deg = lapply(sig, function(x) adultAll[which(adultAll$collapsedconversion=="A:G / T:C" & adultAll$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(adultAll.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#770           758           160          2627            98           699 
elementNROWS(lapply(adonly.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#582           545           138          1939            68           530 

fisher.test(data.frame(c(582,770-582), c(545,758-545))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 0.1038
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9566833 1.5305059
#sample estimates:
#  odds ratio 
#1.209741
fisher.test(data.frame(c(582+1939,770+2627-(582+1939)), c(545+530,758+699-(545+530)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
# p-value = 0.775
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.886572 1.178452
#sample estimates:
#  odds ratio 
#1.022653
fisher.test(data.frame(c(582+138,770+160-(582+138)), c(545+68,758+98-(545+68)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.005474
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.091273 1.693146
#sample estimates:
#  odds ratio 
#1.358873


## Are prenatal-specific editing sites enriched for DEG Fraction?
pren = unique[["prenatalOnly"]]
prenonly.deg = lapply(sig, function(x) pren[which(pren$collapsedconversion=="A:G / T:C" & pren$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(prenonly.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#104           105            24           380            15           106 

fisher.test(data.frame(c(104,975-104), c(105,1010-105))) # both ages: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 0.8838
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7643917 1.3854113
#sample estimates:
#  odds ratio 
#1.029133
fisher.test(data.frame(c(104+380,975+3427-(104+380)), c(105+106,1010+2442-(105+106)))) # adult: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 2.018e-14
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.599683 2.256263
#sample estimates:
#  odds ratio 
#1.897326
fisher.test(data.frame(c(104+24,975+354-(104+24)), c(105+15,1010+350-(105+15)))) # prenatal: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 0.5052
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8406777 1.4434121
#sample estimates:
#  odds ratio 
#1.101264

names(all)
prenatalAll = all[["prenatalAll"]]
dim(prenatalAll) # 23144
prenatalAll.deg = lapply(sig, function(x) prenatalAll[which(prenatalAll$collapsedconversion=="A:G / T:C" & prenatalAll$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(prenatalAll.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          559           449            96          2082            78           494 
elementNROWS(lapply(prenonly.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          371           236            74          1394            48           325

fisher.test(data.frame(c(371,559-371), c(236,449-236))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 1.037e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.368608 2.317642
#sample estimates:
#  odds ratio 
#1.780072
fisher.test(data.frame(c(371+1394,559+2082-(371+1394)), c(236+325,449+494-(236+325)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 5.91e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.173150 1.603399
#sample estimates:
#  odds ratio 
#1.371788
fisher.test(data.frame(c(371+74,559+96-(371+74)), c(236+48,449+78-(236+48)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 1.036e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.420105 2.314834
#sample estimates:
#  odds ratio 
#1.812226

## In Adult: Are cytosolic-specific editing sites enriched for DEG Fraction?
names(unique)
ACnotAN = unique[["ACnotAN"]]
dim(ACnotAN) #2725
ACnotAN.deg = lapply(sig, function(x) ACnotAN[which(ACnotAN$collapsedconversion=="A:G / T:C" & ACnotAN$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotAN.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#32            64            15           170             9            70            

fisher.test(data.frame(c(32+170,975+3427-(32+170)), c(64+70,1010+2442-(64+70)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
# p-value = 0.1295
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9477982 1.5000514
#sample estimates:
#  odds ratio 
#1.190869
fisher.test(data.frame(c(32+15,975+354-(32+15)), c(64+9,1010+350-(64+9)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
# p-value = 0.0247
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4345620 0.9538287
#sample estimates:
#  odds ratio 
#0.6464544

names(all)
allAC = all[["allAC"]]
dim(allAC) #16513
allAC.deg = lapply(sig, function(x) allAC[which(allAC$collapsedconversion=="A:G / T:C" & allAC$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allAC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          204           462            61           929            41           369
elementNROWS(lapply(ACnotAN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#           58           157            28           314            17           148

fisher.test(data.frame(c(58+314,204+929-(58+314)), c(157+148,462+369-(157+148)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
# p-value = 0.07567
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.69569 1.02199
#sample estimates:
#  odds ratio 
#0.8431107
fisher.test(data.frame(c(58+28,204+61-(58+28)), c(157+17,462+41-(157+17)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
# p-value = 0.5751
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6528906 1.2599697
#sample estimates:
#  odds ratio 
#0.9085604 


## Are nuclear-specific editing sites enriched for DEG Fraction?
ANnotAC = unique[["ANnotAC"]]
ANnotAC.deg = lapply(sig, function(x) ANnotAC[which(ANnotAC$collapsedconversion=="A:G / T:C" & ANnotAC$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotAC.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          132            88            34           429            17           109 

fisher.test(data.frame(c(132+429,975+3427-(132+429)), c(88+109,1010+2442-(88+109)))) # adult: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.033639 2.872376
#sample estimates:
#  odds ratio 
#2.412831
fisher.test(data.frame(c(132+34,975+354-(132+34)), c(88+17,1010+350-(88+17)))) # prenatal: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value = 3.968e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.309963 2.228546
#sample estimates:
#  odds ratio 
#1.705658

allAN = all[["allAN"]]
dim(allAN) #24407
allAN.deg = lapply(sig, function(x) allAN[which(allAN$collapsedconversion=="A:G / T:C" & allAN$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allAN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#712           601           132          2313            81           551
elementNROWS(lapply(ANnotAC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#566           296            99          1698            57           330 

fisher.test(data.frame(c(566+1698,712+2313-(566+1698)), c(296+330,601+551-(296+330)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.162781 2.888352
#sample estimates:
#  odds ratio 
#2.499002
fisher.test(data.frame(c(566+99,712+132-(566+99)), c(296+57,601+81-(296+57)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.753205 4.356842
#sample estimates:
#  odds ratio 
#3.459475 

## In Prenatal: Are cytosolic-specific editing sites enriched for DEG Fraction?
names(unique)
PCnotPN = unique[["PCnotPN"]]
dim(PCnotPN) #2725
PCnotPN.deg = lapply(sig, function(x) PCnotPN[which(PCnotPN$collapsedconversion=="A:G / T:C" & PCnotPN$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotPN.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#           45            57             3           192            12            61            

fisher.test(data.frame(c(45+192,975+3427-(45+192)), c(57+61,1010+2442-(57+61)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 3.013e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.277397 2.032065
#sample estimates:
#  odds ratio 
#1.607646
fisher.test(data.frame(c(45+3,975+354-(45+3)), c(57+12,1010+350-(57+12)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.07221
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4707109 1.0372880
#sample estimates:
#  odds ratio 
#0.701178

names(all)
allPC = all[["allPC"]]
allPC.deg = lapply(sig, function(x) allPC[which(allPC$collapsedconversion=="A:G / T:C" & allPC$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allPC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          262           324            26          1155            52           305
elementNROWS(lapply(PCnotPN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          100           129             8           471            22           117

fisher.test(data.frame(c(100+471,262+1155-(100+471)), c(129+117,324+305-(129+117)))) # adult: # sites within retained or exported DEG PNd cytosolic-specific or non-specific status
#p-value = 0.625
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8634059 1.2799323
#sample estimates:
#  odds ratio 
#1.050815
fisher.test(data.frame(c(100+8,262+26-(100+8)), c(129+22,324+52-(129+22)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.521
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6438427 1.2402827
#sample estimates:
#  odds ratio 
#0.8941912


## In prenatal: Are nuclear-specific editing sites enriched for DEG Fraction?
PNnotPC = unique[["PNnotPC"]]
PNnotPC.deg = lapply(sig, function(x) PNnotPC[which(PNnotPC$collapsedconversion=="A:G / T:C" & PNnotPC$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotPC.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#           84            67            22           303             8            68 

fisher.test(data.frame(c(84+303,975+3427-(84+303)), c(67+68,1010+2442-(67+68)))) # adult: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.930785 2.918211
#sample estimates:
#  odds ratio 
#2.367956
fisher.test(data.frame(c(84+22,975+354-(84+22)), c(67+8,1010+350-(67+8)))) # prenatal: retained or exported DEG and presence or absence of nuclear-specific editing site
#p-value = 0.01111
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.082317 2.045061
#sample estimates:
#  odds ratio 
#1.484753

allPN = all[["allPN"]]
dim(allPN) #24407
allPN.deg = lapply(sig, function(x) allPN[which(allPN$collapsedconversion=="A:G / T:C" & allPN$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allPN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          459           320            88          1611            56           377
elementNROWS(lapply(PNnotPC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          297           125            70           927            26           189 

fisher.test(data.frame(c(297+927,459+1611-(297+927)), c(125+189,320+377-(125+189)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 1.134e-10
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.478673 2.106195
#sample estimates:
#  odds ratio 
#1.764408 
fisher.test(data.frame(c(297+70,459+88-(297+70)), c(125+26,320+56-(125+26)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 6.179e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.291805 4.027934
#sample estimates:
#  odds ratio 
#3.03419 

## In cytosol: Are adult-specific editing sites enriched for DEG Fraction?
ACnotPC = unique[["ACnotPC"]]
ACnotPC.deg = lapply(sig, function(x) ACnotPC[which(ACnotPC$collapsedconversion=="A:G / T:C" & ACnotPC$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442 
elementNROWS(lapply(ACnotPC.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#           69            88            24           269            10            94

fisher.test(data.frame(c(69+269,975+3427-(69+269)), c(88+94,1010+2442-(88+94)))) # adult: retained or exported DEG and presence or absence of adult-specific editing site
#p-value = 2.031e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.236623 1.810051
#sample estimates:
#  odds ratio 
#1.494227 
fisher.test(data.frame(c(69+24,975+354-(69+24)), c(88+10,1010+350-(88+10)))) # prenatal: retained or exported DEG and presence or absence of adult-specific editing site
#p-value = 0.8807
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7135514 1.3150579
#sample estimates:
#  odds ratio 
#0.9689459

names(all)
allAC = all[["allAC"]]
allAC.deg = lapply(sig, function(x) allAC[which(allAC$collapsedconversion=="A:G / T:C" & allAC$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allAC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          204           462            61           929            41           369  
elementNROWS(lapply(ACnotPC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          163           329            56           696            24           283 

fisher.test(data.frame(c(163+696,204+929-(163+696)), c(329+283,462+369-(329+283)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 0.2922
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.908016 1.385111
#sample estimates:
#  odds ratio 
#1.121783 
fisher.test(data.frame(c(163+56,204+61-(163+56)), c(329+24,462+41-(329+24)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0001712
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.379538 3.000998
#sample estimates:
#  odds ratio 
#2.021264


## In cytosol: Are prenatal-specific editing sites enriched for DEG Fraction?
PCnotAC = unique[["PCnotAC"]]
PCnotAC.deg = lapply(sig, function(x) PCnotAC[which(PCnotAC$collapsedconversion=="A:G / T:C" & PCnotAC$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotAC.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#           69            85            10           283            12            82 

fisher.test(data.frame(c(69+283,975+3427-(69+283)), c(85+82,1010+2442-(85+82)))) # adult: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 1.606e-08
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.409500 2.079881
#sample estimates:
#  odds ratio 
#1.709529
fisher.test(data.frame(c(69+10,975+354-(69+10)), c(85+12,1010+350-(85+12)))) # prenatal: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 0.2421
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5971791 1.1313507
#sample estimates:
#  odds ratio 
#0.8229641

names(all)
allPC = all[["allPC"]]
allPC.deg = lapply(sig, function(x) allPC[which(allPC$collapsedconversion=="A:G / T:C" & allPC$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allPC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          262           324            26          1155            52           305  
elementNROWS(lapply(PCnotAC.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          221           191            21           922            35           219 

fisher.test(data.frame(c(221+922,262+1155-(221+922)), c(191+219,324+305-(191+219)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 1.486e-13
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.794137 2.764096
#sample estimates:
#  odds ratio 
#2.227265 
fisher.test(data.frame(c(221+21,262+26-(221+21)), c(191+35,324+52-(191+35)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 1.298e-11
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  2.363788 5.207554
#sample estimates:
#  odds ratio 
#3.485251


## In nucleus: Are adult-specific editing sites enriched for DEG Fraction?
ANnotPN = unique[["ANnotPN"]]
ANnotPN.deg = lapply(sig, function(x) ANnotPN[which(ANnotPN$collapsedconversion=="A:G / T:C" & ANnotPN$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442 
elementNROWS(lapply(ANnotPN.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          148           111            39           472            16           135 

fisher.test(data.frame(c(148+472,975+3427-(148+472)), c(111+135,1010+2442-(111+135)))) # adult: retained or exported DEG and presence or absence of adult-specific editing site
# p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.825913 2.505054
#sample estimates:
#  odds ratio 
#2.136273
fisher.test(data.frame(c(148+39,975+354-(148+39)), c(111+16,1010+350-(111+16)))) # prenatal: retained or exported DEG and presence or absence of adult-specific editing site
# p-value = 0.000149
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.243399 2.036752
#sample estimates:
#  odds ratio 
#1.58948

allAN = all[["allAN"]]
allAN.deg = lapply(sig, function(x) allAN[which(allAN$collapsedconversion=="A:G / T:C" & allAN$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allAN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          712           601           132          2313            81           551  
elementNROWS(lapply(ANnotPN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          554           455           112          1765            60           435 

fisher.test(data.frame(c(554+1765,712+2313-(554+1765)), c(455+435,601+551-(455+435)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 0.712
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8193136 1.1392673
#sample estimates:
#  odds ratio 
#0.9669628
fisher.test(data.frame(c(554+112,712+132-(554+112)), c(455+60,601+81-(455+60)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.1239
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9467356 1.5542015
#sample estimates:
#  odds ratio 
#1.213125


## In nucleus: Are prenatal-specific editing sites enriched for DEG Fraction?
PNnotAN = unique[["PNnotAN"]]
PNnotAN.deg = lapply(sig, function(x) PNnotAN[which(PNnotAN$collapsedconversion=="A:G / T:C" & PNnotAN$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotAN.deg[1:6], function(x) unique(x$nearestID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#           95            86            23           346             9            88 

fisher.test(data.frame(c(95+346,975+3427-(95+346)), c(86+88,1010+2442-(86+88)))) # adult: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.744678 2.529655
#sample estimates:
#  odds ratio 
#2.097265
fisher.test(data.frame(c(95+23,975+354-(95+23)), c(86+9,1010+350-(86+9)))) # prenatal: retained or exported DEG and presence or absence of prenatal-specific editing site
#p-value = 0.07423
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.970067 1.738612
#sample estimates:
#  odds ratio 
#1.297361

allPN = all[["allPN"]]
allPN.deg = lapply(sig, function(x) allPN[which(allPN$collapsedconversion=="A:G / T:C" & allPN$nearestID %in% as.character(x$geneID)),])
elementNROWS(lapply(allPN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          459           320            88          1611            35           377   
elementNROWS(lapply(PNnotAN.deg[1:6], function(x) unique(x$editingID)))
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#          301           174            68          1063            35           261 

fisher.test(data.frame(c(301+1063,459+1611-(301+1063)), c(174+261,320+377-(174+261)))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value = 0.09846
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9692207 1.3955462
#sample estimates:
#  odds ratio 
#1.163575
fisher.test(data.frame(c(301+68,459+88-(301+68)), c(174+35,320+35-(174+35)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0105
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.086809 1.927998
#sample estimates:
#  odds ratio 
#1.4475

## Looking into gene expression by age and editing sites
AgeList = list(Cres = as.data.frame(Cres.down), Nres = as.data.frame(Nres))
AgeList = Map(cbind, AgeList, lapply(AgeList, function(x) geneMap[match(rownames(x),rownames(geneMap)),]))
SigAgeList = lapply(AgeList, function(x) x[which(x$padj<=0.05 & abs(x$log2FoldChange) >=1),])
elementNROWS(SigAgeList)
SigAgeList = Map(cbind, SigAgeList, Sign = lapply(SigAgeList, function(x) ifelse(x$log2FoldChange > 0,"UpPrenatal", "DownPrenatal")))
SigAgeList = lapply(SigAgeList, function(x) split(x, x$Sign))
SigList = unlist(SigAgeList, recursive = F) 
lapply(SigList, head)
DEG_editing = lapply(SigList, function(x) editing_anno[which(editing_anno$collapsedconversion=="A:G / T:C" & 
                                                               editing_anno$nearestID %in% rownames(x)),])
elementNROWS(DEG_editing)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5642              8525              5597              6189 
elementNROWS(lapply(DEG_editing, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#2815              3764              2854              2693 
elementNROWS(lapply(DEG_editing, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#641               760               661               607 


## Are cytosolic-specific editing sites enriched for DEG Age?
cytonly.deg = lapply(SigList, function(x) cyt[which(cyt$collapsedconversion=="A:G / T:C" & cyt$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(cytonly.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.3622105 0.5144960
#sample estimates:
#  odds ratio 
#0.432008
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value = 4.254e-13
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.4221270 0.6132562
#sample estimates:
#  odds ratio 
#0.5089893

cytosolAll.deg = lapply(SigList, function(x) cytosolAll[which(cytosolAll$collapsedconversion=="A:G / T:C" & cytosolAll$nearestID %in% rownames(x)),])
elementNROWS(lapply(cytosolAll.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#1360              2192              1293              1648 
elementNROWS(lapply(cytonly.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#481               854               437               649 

fisher.test(data.frame(c(481,1360-481), c(854,2192-854))) # cytosol: # sites within increasing or decreasing age DEG and cytosolic-specific or non-specific status
# p-value = 0.03254
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7429849 0.9889736
#sample estimates:
#  odds ratio 
#0.857381
fisher.test(data.frame(c(437,1293-437), c(649,1648-649))) # nucleus: # sites within increasing or decreasing age DEG and cytosolic-specific or non-specific status
#p-value = 0.002069
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6730767 0.9172312
#sample estimates:
#  odds ratio 
#0.7858987


## Are nuclear-specific editing sites enriched for DEG Age?
nuconly.deg = lapply(SigList, function(x) nuc[which(nuc$collapsedconversion=="A:G / T:C" & nuc$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(nuconly.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#430               513               459               383

fisher.test(data.frame(c(430,5634-430), c(513,4019-513))) # cytosol: increasing or decreasing age DEG and presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.492096 0.647859
#sample estimates:
#  odds ratio 
#0.5647625
fisher.test(data.frame(c(459,5137-459), c(383,3347-383))) # nucleus: increasing or decreasing age DEG and  presence or absence of nuclear-specific editing site
#p-value = 0.0001743
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6563835 0.8788090
#sample estimates:
#  odds ratio 
#0.7593713

nucleusAll.deg = lapply(SigList, function(x) nucleusAll[which(nucleusAll$collapsedconversion=="A:G / T:C" & nucleusAll$nearestID %in% rownames(x)),])
elementNROWS(lapply(nucleusAll.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#2334              2910              2417              2044
elementNROWS(lapply(nuconly.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#1455              1572              1561              1045

fisher.test(data.frame(c(1455,2334-1455), c(1572,2910-1572))) # cytosol: # sites within increasing or decreasing age DEG and nuclear-specific or non-specific status
#p-value = 1.445e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.258900 1.576871
#sample estimates:
#  odds ratio 
#1.408778
fisher.test(data.frame(c(1561,2417-1561), c(1045,2044-1045))) # nucleus: # sites within increasing or decreasing age DEG and nuclear-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.542882 1.969795
#sample estimates:
#  odds ratio 
#1.743135

## Are adult-specific editing sites enriched for DEG Age?
adonly.deg = lapply(SigList, function(x) ad[which(ad$collapsedconversion=="A:G / T:C" & ad$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347 
elementNROWS(lapply(adonly.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#562               323               584               229 

fisher.test(data.frame(c(562,5634-562), c(323,4019-323))) # cytosol: increasing or decreasing age DEG and presence or absence of adult-specific editing site
#p-value = 0.001126
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.096405 1.467928
#sample estimates:
#  odds ratio 
#1.267873
fisher.test(data.frame(c(584,5137-584), c(229,3347-229))) # nucleus: increasing or decreasing age DEG and presence or absence of adult-specific editing site
# p-value = 1.853e-12
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.485841 2.057601
#sample estimates:
#  odds ratio 
#1.746337

adultAll.deg = lapply(SigList, function(x) adultAll[which(adultAll$collapsedconversion=="A:G / T:C" & adultAll$nearestID %in% rownames(x)),])
elementNROWS(lapply(adultAll.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#2592              1514              2624               987 
elementNROWS(lapply(adonly.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#2320               879              2376               583 

fisher.test(data.frame(c(2320,2592-2320), c(879,1514-879))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.225823 7.272376
#sample estimates:
#  odds ratio 
#6.158586
fisher.test(data.frame(c(2376,2624-2376), c(583,987-583))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.509525 7.999610
#sample estimates:
#  odds ratio 
#6.634573


## Are prenatal-specific editing sites enriched for DEG Age?
prenonly.deg = lapply(SigList, function(x) pren[which(pren$collapsedconversion=="A:G / T:C" & pren$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(prenonly.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#130               582               131               465 

fisher.test(data.frame(c(130,5634-130), c(582,4019-582))) # cytosol: increasing/decreasing age DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1138996 0.1698730
#sample estimates:
#  odds ratio 
#0.1395105
fisher.test(data.frame(c(131,5137-131), c(465,3347-465))) # nucleus: increasing/decreasing age DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1318666 0.1984455
#sample estimates:
#  odds ratio 
#0.1622252

prenatalAll.deg = lapply(SigList, function(x) prenatalAll[which(prenatalAll$collapsedconversion=="A:G / T:C" & prenatalAll$nearestID %in% rownames(x)),])
elementNROWS(lapply(prenatalAll.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#495              2885               478              2110 
elementNROWS(lapply(prenonly.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#223              2250               230              1706

fisher.test(data.frame(c(223,495-223), c(2250,2885-2250))) # cytosol: # sites within increasing/decreasing age DEG and prenatal-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1889044 0.2834098
#sample estimates:
#  odds ratio 
#0.2315041 
fisher.test(data.frame(c(230,478-230), c(1706,2110-1706))) # nucleus: # sites within increasing/decreasing age DEG and prenatal-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1771500 0.2724677
#sample estimates:
#  odds ratio 
#0.2197838


## In adult: Are cytosolic-specific editing sites enriched for DEG Age?
ACnotAN.deg = lapply(SigList, function(x) ACnotAN[which(ACnotAN$collapsedconversion=="A:G / T:C" & ACnotAN$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             5634              4019              5137              3347 
elementNROWS(lapply(ACnotAN.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              208               156               202               126 

fisher.test(data.frame(c(208,5634-208), c(156,4019-156))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.6262
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7643087 1.1809113
#sample estimates:
#  odds ratio 
#0.9492834 
fisher.test(data.frame(c(202,5137-202), c(126,3347-126))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.7297
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8296287 1.3236924
#sample estimates:
#  odds ratio 
#1.046368

allAC.deg = lapply(SigList, function(x) allAC[which(allAC$collapsedconversion=="A:G / T:C" & allAC$nearestID %in% rownames(x)),])
elementNROWS(lapply(allAC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             1214               665              1144               508
elementNROWS(lapply(ACnotAN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              452               281               398               227

fisher.test(data.frame(c(452,1214-452), c(281,665-281))) # cytosol: # sites within increasing or decreasing age DEG and cytosolic-specific or non-specific status
#p-value = 0.03358
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6653819 0.9880163
#sample estimates:
#  odds ratio 
#0.8107035
fisher.test(data.frame(c(398,1144-398), c(227,508-227))) # nucleus: # sites within increasing or decreasing age DEG and cytosolic-specific or non-specific status
#p-value = 0.0001478
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5308015 0.8222648
#sample estimates:
#  odds ratio 
#0.6606228


## In adult: Are nuclear-specific editing sites enriched for DEG Age?
ANnotAC.deg = lapply(SigList, function(x) ANnotAC[which(ANnotAC$collapsedconversion=="A:G / T:C" & ANnotAC$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotAC.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              399               266               429               172

fisher.test(data.frame(c(399,5634-399), c(266,4019-266))) # cytosol: increasing or decreasing age DEG and presence or absence of nuclear-specific editing site
#p-value = 0.3921
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9132039 1.2678446
#sample estimates:
#  odds ratio 
#1.075369
fisher.test(data.frame(c(429,5137-429), c(172,3347-172))) # nucleus: increasing or decreasing age DEG and  presence or absence of nuclear-specific editing site
#p-value = 9.9e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.397662 2.031085
#sample estimates:
#  odds ratio 
#1.681939

allAN.deg = lapply(SigList, function(x) allAN[which(allAN$collapsedconversion=="A:G / T:C" & allAN$nearestID %in% rownames(x)),])
elementNROWS(lapply(allAN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             2140              1233              2226               760
elementNROWS(lapply(ANnotAC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#            1378               849              1480               479

fisher.test(data.frame(c(1378,2140-1378), c(849,1233-849))) # cytosol: # sites within increasing or decreasing age DEG and nuclear-specific or non-specific status
#p-value = 0.009186
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.7020910 0.9523032
#sample estimates:
#  odds ratio 
#0.8179831 
fisher.test(data.frame(c(1480,2226-1480), c(479,760-479))) # nucleus: # sites within increasing or decreasing age DEG and nuclear-specific or non-specific status
#p-value = 0.08473
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9761705 1.3862375
#sample estimates:
#  odds ratio 
#1.163773 


## In prenatal: Are cytosolic-specific editing sites enriched for DEG Age?
PCnotPN.deg = lapply(SigList, function(x) PCnotPN[which(PCnotPN$collapsedconversion=="A:G / T:C" & PCnotPN$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             5634              4019              5137              3347 
elementNROWS(lapply(PCnotPN.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#               67               304                64               232

fisher.test(data.frame(c(67,5634-67), c(304,4019-304))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1108179 0.1928972
#sample estimates:
#  odds ratio 
#0.1470929 
fisher.test(data.frame(c(64,5137-64), c(232,3347-232))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1259313 0.2252244
#sample estimates:
#  odds ratio 
#0.1694239

allPC.deg = lapply(SigList, function(x) allPC[which(allPC$collapsedconversion=="A:G / T:C" & allPC$nearestID %in% rownames(x)),])
elementNROWS(lapply(allPC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              278              1817               262              1364 
elementNROWS(lapply(PCnotPN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              103               749                98               557 

fisher.test(data.frame(c(103,278-103), c(749,1817-749))) # cytosol: # sites within increasing or decreasing age DEG and cytosolic-specific or non-specific status
#p-value = 0.1906
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6398391 1.0969241
#sample estimates:
#  odds ratio 
#0.8393145
fisher.test(data.frame(c(98,262-98), c(557,1364-557))) # nucleus: # sites within increasing or decreasing age DEG and cytosolic-specific or non-specific status
#p-value = 0.3357
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6518035 1.1459588
#sample estimates:
#  odds ratio 
#0.8658424 


## In prenatal: Are nuclear-specific editing sites enriched for DEG Age?
PNnotPC.deg = lapply(SigList, function(x) PNnotPC[which(PNnotPC$collapsedconversion=="A:G / T:C" & PNnotPC$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotPC.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              107               398               102               311

fisher.test(data.frame(c(107,5634-107), c(398,4019-398))) # cytosol: increasing or decreasing age DEG and presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.140369 0.219578
#sample estimates:
#  odds ratio 
#0.1761634
fisher.test(data.frame(c(102,5137-102), c(311,3347-311))) # nucleus: increasing or decreasing age DEG and  presence or absence of nuclear-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1558625 0.2493482
#sample estimates:
#  odds ratio 
#0.197801

allPN.deg = lapply(SigList, function(x) allPN[which(allPN$collapsedconversion=="A:G / T:C" & allPN$nearestID %in% rownames(x)),])
elementNROWS(lapply(allPN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              392              2136               380              1553 
elementNROWS(lapply(PNnotPC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              217              1068               216               746 

fisher.test(data.frame(c(217,392-217), c(1068,2136-1068))) # cytosol: # sites within increasing or decreasing age DEG and nuclear-specific or non-specific status
#p-value = 0.05442
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9930726 1.5496486
#sample estimates:
#  odds ratio 
#1.239891
fisher.test(data.frame(c(216,380-216), c(746,1553-746))) # nucleus: # sites within increasing or decreasing age DEG and nuclear-specific or non-specific status
#p-value = 0.002392
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.129529 1.798900
#sample estimates:
#  odds ratio 
#1.424493 


## in Cytosol: Are adult-specific editing sites enriched for DEG Age?
ACnotPC.deg = lapply(SigList, function(x) ACnotPC[which(ACnotPC$collapsedconversion=="A:G / T:C" & ACnotPC$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             5634              4019              5137              3347
elementNROWS(lapply(ACnotPC.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              344               186               344               144 

fisher.test(data.frame(c(344,5634-344), c(186,4019-186))) # cytosol: increasing or decreasing age DEG and presence or absence of adult-specific editing site
#p-value = 0.001744
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.112245 1.618421
#sample estimates:
#  odds ratio 
#1.340033
fisher.test(data.frame(c(344,5137-344), c(144,3347-144))) # nucleus: increasing or decreasing age DEG and presence or absence of adult-specific editing site
#p-value = 2.767e-06
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.303285 1.963161
#sample estimates:
#  odds ratio 
#1.596328

allAC.deg = lapply(SigList, function(x) allAC[which(allAC$collapsedconversion=="A:G / T:C" & allAC$nearestID %in% rownames(x)),])
elementNROWS(lapply(allAC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             1214               665              1144               508 
elementNROWS(lapply(ACnotPC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             1082               375              1031               284

fisher.test(data.frame(c(1082,1214-1082), c(375,665-375))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  4.969351 8.096762
#sample estimates:
#  odds ratio 
#6.331615
fisher.test(data.frame(c(1031,1144-1031), c(284,508-284))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.495315 9.433971
#sample estimates:
#  odds ratio 
#7.184923


## in Cytosol: Are prenatal-specific editing sites enriched for DEG Age?
PCnotAC.deg = lapply(SigList, function(x) PCnotAC[which(PCnotAC$collapsedconversion=="A:G / T:C" & PCnotAC$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotAC.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#               87               430                90               338

fisher.test(data.frame(c(87,5634-87), c(430,4019-430))) # cytosol: increasing/decreasing age DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1023394 0.1658913
#sample estimates:
#  odds ratio 
#0.1309208
fisher.test(data.frame(c(90,5137-90), c(338,3347-338))) # nucleus: increasing/decreasing age DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1238817 0.2018565
#sample estimates:
#  odds ratio 
#0.1587847

allPC.deg = lapply(SigList, function(x) allPC[which(allPC$collapsedconversion=="A:G / T:C" & allPC$nearestID %in% rownames(x)),])
elementNROWS(lapply(allPC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              278              1817               262              1364 
elementNROWS(lapply(PCnotAC.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             146              1527               149              1140 

fisher.test(data.frame(c(146,278-146), c(1527,1817-1527))) # cytosol: # sites within increasing/decreasing age DEG and prenatal-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1595685 0.2770385
#sample estimates:
#  odds ratio 
#0.2102651
fisher.test(data.frame(c(149,262-149), c(1140,1364-1140))) # nucleus: # sites within increasing/decreasing age DEG and prenatal-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1934672 0.3480022
#sample estimates:
#  odds ratio 
#0.2593643

## in Nucleus: Are adult-specific editing sites enriched for DEG Age?
ANnotPN.deg = lapply(SigList, function(x) ANnotPN[which(ANnotPN$collapsedconversion=="A:G / T:C" & ANnotPN$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347 
elementNROWS(lapply(ANnotPN.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              500               280               526               189 

fisher.test(data.frame(c(500,5634-500), c(280,4019-280))) # cytosol: increasing or decreasing age DEG and presence or absence of adult-specific editing site
#p-value = 0.0007392
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.114289 1.520040
#sample estimates:
#  odds ratio 
#1.300468
fisher.test(data.frame(c(526,5137-526), c(189,3347-189))) # nucleus: increasing or decreasing age DEG and presence or absence of adult-specific editing site
#p-value = 2.988e-14
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.600780 2.276763
#sample estimates:
#  odds ratio 
#1.905928

allAN.deg = lapply(SigList, function(x) allAN[which(allAN$collapsedconversion=="A:G / T:C" & allAN$nearestID %in% rownames(x)),])
elementNROWS(lapply(allAN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             2140              1233              2226               760 
elementNROWS(lapply(ANnotPN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#             1942               774              2037               491 

fisher.test(data.frame(c(1942,2140-1942), c(774,1233-774))) # both ages: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  4.809596 7.044141
#sample estimates:
#  odds ratio 
#5.81295
fisher.test(data.frame(c(2037,2226-2037), c(491,760-491))) # adult: # sites within retained or exported DEG and cytosolic-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  4.756517 7.329345
#sample estimates:
#  odds ratio 
#5.90018


## in Nucleus: Are prenatal-specific editing sites enriched for DEG Age?
PNnotAN.deg = lapply(SigList, function(x) PNnotAN[which(PNnotAN$collapsedconversion=="A:G / T:C" & PNnotAN$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotAN.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              117               500               119               407 

fisher.test(data.frame(c(117,5634-117), c(500,4019-500))) # cytosol: increasing/decreasing age DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1204461 0.1837764
#sample estimates:
#  odds ratio 
#0.1492844 
fisher.test(data.frame(c(119,5137-119), c(407,3347-407))) # nucleus: increasing/decreasing age DEG and presence or absence of prenatal-specific editing site
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1377963 0.2118020
#sample estimates:
#  odds ratio 
#0.1713415

allPN.deg = lapply(SigList, function(x) allPN[which(allPN$collapsedconversion=="A:G / T:C" & allPN$nearestID %in% rownames(x)),])
elementNROWS(lapply(allPN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              392              2136               380              1553
elementNROWS(lapply(PNnotAN.deg, function(x) unique(x$editingID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#              194              1677               191              1284 

fisher.test(data.frame(c(194,392-194), c(1677,2136-1677))) # cytosol: # sites within increasing/decreasing age DEG and prenatal-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.2132424 0.3375627
#sample estimates:
#  odds ratio 
#0.2683435
fisher.test(data.frame(c(191,380-191), c(1284,1553-1284))) # nucleus: # sites within increasing/decreasing age DEG and prenatal-specific or non-specific status
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.1653845 0.2713496
#sample estimates:
#  odds ratio 
#0.2119319


## Assess DEG patterns by annotation of editing site and whether the site is unique to a specific group
split.anno = lapply(unique, function(x) split(x, f = x$annotation))
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


## Are cytosolic-specific editing sites in 3'UTR enriched for DEG Fraction?
cyt.UTR3 = cytosolOnly[["UTR3"]]
cyt.UTR3.deg = lapply(sig, function(x) cyt.UTR3[which(cyt.UTR3$collapsedconversion=="A:G / T:C" & cyt.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cyt.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0429

cyt.UTR5 = cytosolOnly[["UTR5"]]
cyt.UTR5.deg = lapply(sig, function(x) cyt.UTR5[which(cyt.UTR5$collapsedconversion=="A:G / T:C" & cyt.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cyt.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0429

cyt.Intron = cytosolOnly[["Intron"]]
cyt.Intron.deg = lapply(sig, function(x) cyt.Intron[which(cyt.Intron$collapsedconversion=="A:G / T:C" & cyt.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cyt.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0429

cyt.Other = cytosolOnly[["Other"]]
cyt.Other.deg = lapply(sig, function(x) cyt.Other[which(cyt.Other$collapsedconversion=="A:G / T:C" & cyt.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cyt.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0429

cyt.CDS = cytosolOnly[["CDS"]]
cyt.CDS.deg = lapply(sig, function(x) cyt.CDS[which(cyt.CDS$collapsedconversion=="A:G / T:C" & cyt.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(cyt.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of cytosolic-specific editing site
#p-value = 0.0429


## Are nuclear-specific editing sites enriched for DEG Fraction?
nuc.UTR3 = nucleusOnly[["UTR3"]]
nuc.UTR3.deg = lapply(sig, function(x) nuc.UTR3[which(nuc.UTR3$collapsedconversion=="A:G / T:C" & nuc.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuc.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.0429

nuc.UTR5 = nucleusOnly[["UTR5"]]
nuc.UTR5.deg = lapply(sig, function(x) nuc.UTR5[which(nuc.UTR5$collapsedconversion=="A:G / T:C" & nuc.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuc.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.0429

nuc.Intron = nucleusOnly[["Intron"]]
nuc.Intron.deg = lapply(sig, function(x) nuc.Intron[which(nuc.Intron$collapsedconversion=="A:G / T:C" & nuc.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuc.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.0429

nuc.Other = nucleusOnly[["Other"]]
nuc.Other.deg = lapply(sig, function(x) nuc.Other[which(nuc.Other$collapsedconversion=="A:G / T:C" & nuc.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuc.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.0429

nuc.CDS = nucleusOnly[["CDS"]]
nuc.CDS.deg = lapply(sig, function(x) nuc.CDS[which(nuc.CDS$collapsedconversion=="A:G / T:C" & nuc.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(nuc.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of nucosolic-specific editing site
#p-value = 0.0429


## Are adult-specific editing sites enriched for DEG Fraction?
ad.UTR3 = adultOnly[["UTR3"]]
ad.UTR3.deg = lapply(sig, function(x) ad.UTR3[which(ad.UTR3$collapsedconversion=="A:G / T:C" & ad.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ad.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ad.UTR5 = adultOnly[["UTR5"]]
ad.UTR5.deg = lapply(sig, function(x) ad.UTR5[which(ad.UTR5$collapsedconversion=="A:G / T:C" & ad.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ad.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ad.Intron = adultOnly[["Intron"]]
ad.Intron.deg = lapply(sig, function(x) ad.Intron[which(ad.Intron$collapsedconversion=="A:G / T:C" & ad.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ad.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ad.Other = adultOnly[["Other"]]
ad.Other.deg = lapply(sig, function(x) ad.Other[which(ad.Other$collapsedconversion=="A:G / T:C" & ad.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ad.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ad.CDS = adultOnly[["CDS"]]
ad.CDS.deg = lapply(sig, function(x) ad.CDS[which(ad.CDS$collapsedconversion=="A:G / T:C" & ad.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ad.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429


## Are prenatal-specific editing sites enriched for DEG Fraction?
pren = unique[["prenatalOnly"]]
prenonly.deg = lapply(sig, function(x) pren[which(pren$collapsedconversion=="A:G / T:C" & pren$nearestID %in% as.character(x$geneID)),])
pren.UTR3 = prenatalOnly[["UTR3"]]
pren.UTR3.deg = lapply(sig, function(x) pren.UTR3[which(pren.UTR3$collapsedconversion=="A:G / T:C" & pren.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(pren.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

pren.UTR5 = prenatalOnly[["UTR5"]]
pren.UTR5.deg = lapply(sig, function(x) pren.UTR5[which(pren.UTR5$collapsedconversion=="A:G / T:C" & pren.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(pren.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

pren.Intron = prenatalOnly[["Intron"]]
pren.Intron.deg = lapply(sig, function(x) pren.Intron[which(pren.Intron$collapsedconversion=="A:G / T:C" & pren.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(pren.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

pren.Other = prenatalOnly[["Other"]]
pren.Other.deg = lapply(sig, function(x) pren.Other[which(pren.Other$collapsedconversion=="A:G / T:C" & pren.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(pren.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

pren.CDS = prenatalOnly[["CDS"]]
pren.CDS.deg = lapply(sig, function(x) pren.CDS[which(pren.CDS$collapsedconversion=="A:G / T:C" & pren.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(pren.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

## In Adult: Are cytosolic-specific editing sites enriched for DEG Fraction?
names(unique)
ACnotAN = unique[["ACnotAN"]]
ACnotAN.UTR3 = ACnotAN[["UTR3"]]
ACnotAN.UTR3.deg = lapply(sig, function(x) ACnotAN.UTR3[which(ACnotAN.UTR3$collapsedconversion=="A:G / T:C" & ACnotAN.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotAN.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotAN.UTR5 = ACnotAN[["UTR5"]]
ACnotAN.UTR5.deg = lapply(sig, function(x) ACnotAN.UTR5[which(ACnotAN.UTR5$collapsedconversion=="A:G / T:C" & ACnotAN.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotAN.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotAN.Intron = ACnotAN[["Intron"]]
ACnotAN.Intron.deg = lapply(sig, function(x) ACnotAN.Intron[which(ACnotAN.Intron$collapsedconversion=="A:G / T:C" & ACnotAN.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotAN.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotAN.Other = ACnotAN[["Other"]]
ACnotAN.Other.deg = lapply(sig, function(x) ACnotAN.Other[which(ACnotAN.Other$collapsedconversion=="A:G / T:C" & ACnotAN.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotAN.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotAN.CDS = ACnotAN[["CDS"]]
ACnotAN.CDS.deg = lapply(sig, function(x) ACnotAN.CDS[which(ACnotAN.CDS$collapsedconversion=="A:G / T:C" & ACnotAN.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotAN.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429


## In Adult: Are nuclear-specific editing sites enriched for DEG Fraction?
ANnotAC = unique[["ANnotAC"]]
ANnotAC.deg = lapply(sig, function(x) ANnotAC[which(ANnotAC$collapsedconversion=="A:G / T:C" & ANnotAC$nearestID %in% as.character(x$geneID)),])
ANnotAC.UTR3 = ANnotAC[["UTR3"]]
ANnotAC.UTR3.deg = lapply(sig, function(x) ANnotAC.UTR3[which(ANnotAC.UTR3$collapsedconversion=="A:G / T:C" & ANnotAC.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotAC.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotAC.UTR5 = ANnotAC[["UTR5"]]
ANnotAC.UTR5.deg = lapply(sig, function(x) ANnotAC.UTR5[which(ANnotAC.UTR5$collapsedconversion=="A:G / T:C" & ANnotAC.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotAC.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotAC.Intron = ANnotAC[["Intron"]]
ANnotAC.Intron.deg = lapply(sig, function(x) ANnotAC.Intron[which(ANnotAC.Intron$collapsedconversion=="A:G / T:C" & ANnotAC.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotAC.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotAC.Other = ANnotAC[["Other"]]
ANnotAC.Other.deg = lapply(sig, function(x) ANnotAC.Other[which(ANnotAC.Other$collapsedconversion=="A:G / T:C" & ANnotAC.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotAC.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotAC.CDS = ANnotAC[["CDS"]]
ANnotAC.CDS.deg = lapply(sig, function(x) ANnotAC.CDS[which(ANnotAC.CDS$collapsedconversion=="A:G / T:C" & ANnotAC.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotAC.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429


## In Prenatal: Are cytosolic-specific editing sites enriched for DEG Fraction?
PCnotPN = unique[["PCnotPN"]]
PCnotPN.deg = lapply(sig, function(x) PCnotPN[which(PCnotPN$collapsedconversion=="A:G / T:C" & PCnotPN$nearestID %in% as.character(x$geneID)),])
PCnotPN.UTR3 = PCnotPN[["UTR3"]]
PCnotPN.UTR3.deg = lapply(sig, function(x) PCnotPN.UTR3[which(PCnotPN.UTR3$collapsedconversion=="A:G / T:C" & PCnotPN.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotPN.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotPN.UTR5 = PCnotPN[["UTR5"]]
PCnotPN.UTR5.deg = lapply(sig, function(x) PCnotPN.UTR5[which(PCnotPN.UTR5$collapsedconversion=="A:G / T:C" & PCnotPN.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotPN.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotPN.Intron = PCnotPN[["Intron"]]
PCnotPN.Intron.deg = lapply(sig, function(x) PCnotPN.Intron[which(PCnotPN.Intron$collapsedconversion=="A:G / T:C" & PCnotPN.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotPN.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotPN.Other = PCnotPN[["Other"]]
PCnotPN.Other.deg = lapply(sig, function(x) PCnotPN.Other[which(PCnotPN.Other$collapsedconversion=="A:G / T:C" & PCnotPN.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotPN.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotPN.CDS = PCnotPN[["CDS"]]
PCnotPN.CDS.deg = lapply(sig, function(x) PCnotPN.CDS[which(PCnotPN.CDS$collapsedconversion=="A:G / T:C" & PCnotPN.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotPN.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

## In prenatal: Are nuclear-specific editing sites enriched for DEG Fraction?
PNnotPC = unique[["PNnotPC"]]
PNnotPC.deg = lapply(sig, function(x) PNnotPC[which(PNnotPC$collapsedconversion=="A:G / T:C" & PNnotPC$nearestID %in% as.character(x$geneID)),])
PNnotPC.UTR3 = PNnotPC.[["UTR3"]]
PNnotPC.UTR3.deg = lapply(sig, function(x) PNnotPC.UTR3[which(PNnotPC.UTR3$collapsedconversion=="A:G / T:C" & PNnotPC.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotPC.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotPC.UTR5 = PNnotPC.[["UTR5"]]
PNnotPC.UTR5.deg = lapply(sig, function(x) PNnotPC.UTR5[which(PNnotPC.UTR5$collapsedconversion=="A:G / T:C" & PNnotPC.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotPC.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotPC.Intron = PNnotPC.[["Intron"]]
PNnotPC.Intron.deg = lapply(sig, function(x) PNnotPC.Intron[which(PNnotPC.Intron$collapsedconversion=="A:G / T:C" & PNnotPC.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotPC.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotPC.Other = PNnotPC.[["Other"]]
PNnotPC.Other.deg = lapply(sig, function(x) PNnotPC.Other[which(PNnotPC.Other$collapsedconversion=="A:G / T:C" & PNnotPC.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotPC.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotPC.CDS = PNnotPC.[["CDS"]]
PNnotPC.CDS.deg = lapply(sig, function(x) PNnotPC.CDS[which(PNnotPC.CDS$collapsedconversion=="A:G / T:C" & PNnotPC.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotPC.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

## In cytosol: Are adult-specific editing sites enriched for DEG Fraction?
ACnotPC = unique[["ACnotPC"]]
ACnotPC.deg = lapply(sig, function(x) ACnotPC[which(ACnotPC$collapsedconversion=="A:G / T:C" & ACnotPC$nearestID %in% as.character(x$geneID)),])
ACnotPC.UTR3 = ACnotPC[["UTR3"]]
ACnotPC.UTR3.deg = lapply(sig, function(x) ACnotPC.UTR3[which(ACnotPC.UTR3$collapsedconversion=="A:G / T:C" & ACnotPC.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotPC.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotPC.UTR5 = ACnotPC[["UTR5"]]
ACnotPC.UTR5.deg = lapply(sig, function(x) ACnotPC.UTR5[which(ACnotPC.UTR5$collapsedconversion=="A:G / T:C" & ACnotPC.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotPC.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotPC.Intron = ACnotPC[["Intron"]]
ACnotPC.Intron.deg = lapply(sig, function(x) ACnotPC.Intron[which(ACnotPC.Intron$collapsedconversion=="A:G / T:C" & ACnotPC.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotPC.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotPC.Other = ACnotPC[["Other"]]
ACnotPC.Other.deg = lapply(sig, function(x) ACnotPC.Other[which(ACnotPC.Other$collapsedconversion=="A:G / T:C" & ACnotPC.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotPC.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ACnotPC.CDS = ACnotPC[["CDS"]]
ACnotPC.CDS.deg = lapply(sig, function(x) ACnotPC.CDS[which(ACnotPC.CDS$collapsedconversion=="A:G / T:C" & ACnotPC.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ACnotPC.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

## In cytosol: Are prenatal-specific editing sites enriched for DEG Fraction?
PCnotAC = unique[["PCnotAC"]]
PCnotAC.deg = lapply(sig, function(x) PCnotAC[which(PCnotAC$collapsedconversion=="A:G / T:C" & PCnotAC$nearestID %in% as.character(x$geneID)),])
PCnotAC.UTR3 = PCnotAC[["UTR3"]]
PCnotAC.UTR3.deg = lapply(sig, function(x) PCnotAC.UTR3[which(PCnotAC.UTR3$collapsedconversion=="A:G / T:C" & PCnotAC.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotAC.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotAC.UTR5 = PCnotAC[["UTR5"]]
PCnotAC.UTR5.deg = lapply(sig, function(x) PCnotAC.UTR5[which(PCnotAC.UTR5$collapsedconversion=="A:G / T:C" & PCnotAC.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotAC.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotAC.Intron = PCnotAC[["Intron"]]
PCnotAC.Intron.deg = lapply(sig, function(x) PCnotAC.Intron[which(PCnotAC.Intron$collapsedconversion=="A:G / T:C" & PCnotAC.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotAC.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotAC.Other = PCnotAC[["Other"]]
PCnotAC.Other.deg = lapply(sig, function(x) PCnotAC.Other[which(PCnotAC.Other$collapsedconversion=="A:G / T:C" & PCnotAC.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotAC.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PCnotAC.CDS = PCnotAC[["CDS"]]
PCnotAC.CDS.deg = lapply(sig, function(x) PCnotAC.CDS[which(PCnotAC.CDS$collapsedconversion=="A:G / T:C" & PCnotAC.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PCnotAC.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

## In nucleus: Are adult-specific editing sites enriched for DEG Fraction?
ANnotPN = unique[["ANnotPN"]]
ANnotPN.deg = lapply(sig, function(x) ANnotPN[which(ANnotPN$collapsedconversion=="A:G / T:C" & ANnotPN$nearestID %in% as.character(x$geneID)),])
ANnotPN.UTR3 = ANnotPN[["UTR3"]]
ANnotPN.UTR3.deg = lapply(sig, function(x) ANnotPN.UTR3[which(ANnotPN.UTR3$collapsedconversion=="A:G / T:C" & ANnotPN.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotPN.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotPN.UTR5 = ANnotPN[["UTR5"]]
ANnotPN.UTR5.deg = lapply(sig, function(x) ANnotPN.UTR5[which(ANnotPN.UTR5$collapsedconversion=="A:G / T:C" & ANnotPN.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotPN.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotPN.Intron = ANnotPN[["Intron"]]
ANnotPN.Intron.deg = lapply(sig, function(x) ANnotPN.Intron[which(ANnotPN.Intron$collapsedconversion=="A:G / T:C" & ANnotPN.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotPN.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotPN.Other = ANnotPN[["Other"]]
ANnotPN.Other.deg = lapply(sig, function(x) ANnotPN.Other[which(ANnotPN.Other$collapsedconversion=="A:G / T:C" & ANnotPN.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotPN.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

ANnotPN.CDS = ANnotPN[["CDS"]]
ANnotPN.CDS.deg = lapply(sig, function(x) ANnotPN.CDS[which(ANnotPN.CDS$collapsedconversion=="A:G / T:C" & ANnotPN.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(ANnotPN.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429


## In nucleus: Are prenatal-specific editing sites enriched for DEG Fraction?
PNnotAN = unique[["PNnotAN"]]
PNnotAN.deg = lapply(sig, function(x) PNnotAN[which(PNnotAN$collapsedconversion=="A:G / T:C" & PNnotAN$nearestID %in% as.character(x$geneID)),])
PNnotAN.UTR3 = PNnotAN[["UTR3"]]
PNnotAN.UTR3.deg = lapply(sig, function(x) PNnotAN.UTR3[which(PNnotAN.UTR3$collapsedconversion=="A:G / T:C" & PNnotAN.UTR3$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotAN.UTR3.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotAN.UTR5 = PNnotAN[["UTR5"]]
PNnotAN.UTR5.deg = lapply(sig, function(x) PNnotAN.UTR5[which(PNnotAN.UTR5$collapsedconversion=="A:G / T:C" & PNnotAN.UTR5$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotAN.UTR5.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotAN.Intron = PNnotAN[["Intron"]]
PNnotAN.Intron.deg = lapply(sig, function(x) PNnotAN.Intron[which(PNnotAN.Intron$collapsedconversion=="A:G / T:C" & PNnotAN.Intron$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotAN.Intron.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotAN.Other = PNnotAN[["Other"]]
PNnotAN.Other.deg = lapply(sig, function(x) PNnotAN.Other[which(PNnotAN.Other$collapsedconversion=="A:G / T:C" & PNnotAN.Other$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotAN.Other.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

PNnotAN.CDS = PNnotAN[["CDS"]]
PNnotAN.CDS.deg = lapply(sig, function(x) PNnotAN.CDS[which(PNnotAN.CDS$collapsedconversion=="A:G / T:C" & PNnotAN.CDS$nearestID %in% as.character(x$geneID)),])
elementNROWS(sig[1:6])
#both_retained both_exported  Fet_retained   Ad_retained  Fet_exported   Ad_exported 
#975          1010           354          3427           350          2442
elementNROWS(lapply(PNnotAN.CDS.deg[1:6], function(x) unique(x$nearestID)))
#both_retained  both_exported   Fet_retained    Ad_retained   Fet_exported    Ad_exported 
#58             88             16            279             15             98            

fisher.test(data.frame(c(88,1010-88), c(58,975-58))) # both ages: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.02011
fisher.test(data.frame(c(58+279,975+3427-(58+279)), c(88+98,1010+2442-(88+98)))) # adult: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 5.862e-05
fisher.test(data.frame(c(58+16,975+354-(58+16)), c(88+15,1010+350-(88+15)))) # prenatal: retained or exported DEG and presence or absence of adosolic-specific editing site
#p-value = 0.0429

## Are cytosolic-specific editing sites enriched for DEG Age?
cyt.UTR3.deg = lapply(SigList, function(x) cyt.UTR3[which(cyt.UTR3$collapsedconversion=="A:G / T:C" & cyt.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(cyt.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value = 4.254e-13

cyt.UTR5.deg = lapply(SigList, function(x) cyt.UTR5[which(cyt.UTR5$collapsedconversion=="A:G / T:C" & cyt.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(cyt.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value = 4.254e-13

cyt.Intron.deg = lapply(SigList, function(x) cyt.Intron[which(cyt.Intron$collapsedconversion=="A:G / T:C" & cyt.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(cyt.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value = 4.254e-13

cyt.CDS.deg = lapply(SigList, function(x) cyt.CDS[which(cyt.CDS$collapsedconversion=="A:G / T:C" & cyt.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(cyt.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value = 4.254e-13

cyt.Other.deg = lapply(SigList, function(x) cyt.Other[which(cyt.Other$collapsedconversion=="A:G / T:C" & cyt.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(cyt.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # cytosol: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of cytosolic-specific editing site
# p-value = 4.254e-13


## Are nuclear-specific editing sites enriched for DEG Age?
nuc.UTR3.deg = lapply(SigList, function(x) nuc.UTR3[which(nuc.UTR3$collapsedconversion=="A:G / T:C" & nuc.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(nuc.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # nucosol: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value = 4.254e-13

nuc.UTR5.deg = lapply(SigList, function(x) nuc.UTR5[which(nuc.UTR5$collapsedconversion=="A:G / T:C" & nuc.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(nuc.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # nucosol: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value = 4.254e-13

nuc.Intron.deg = lapply(SigList, function(x) nuc.Intron[which(nuc.Intron$collapsedconversion=="A:G / T:C" & nuc.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(nuc.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # nucosol: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value = 4.254e-13

nuc.CDS.deg = lapply(SigList, function(x) nuc.CDS[which(nuc.CDS$collapsedconversion=="A:G / T:C" & nuc.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(nuc.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # nucosol: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value = 4.254e-13

nuc.Other.deg = lapply(SigList, function(x) nuc.Other[which(nuc.Other$collapsedconversion=="A:G / T:C" & nuc.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(nuc.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # nucosol: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of nucosolic-specific editing site
# p-value = 4.254e-13


## Are adult-specific editing sites enriched for DEG Age?
ad.UTR3.deg = lapply(SigList, function(x) ad.UTR3[which(ad.UTR3$collapsedconversion=="A:G / T:C" & ad.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ad.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # adosol: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value = 4.254e-13

ad.UTR5.deg = lapply(SigList, function(x) ad.UTR5[which(ad.UTR5$collapsedconversion=="A:G / T:C" & ad.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ad.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # adosol: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value = 4.254e-13

ad.Intron.deg = lapply(SigList, function(x) ad.Intron[which(ad.Intron$collapsedconversion=="A:G / T:C" & ad.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ad.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # adosol: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value = 4.254e-13

ad.CDS.deg = lapply(SigList, function(x) ad.CDS[which(ad.CDS$collapsedconversion=="A:G / T:C" & ad.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ad.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # adosol: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value = 4.254e-13

ad.Other.deg = lapply(SigList, function(x) ad.Other[which(ad.Other$collapsedconversion=="A:G / T:C" & ad.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ad.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # adosol: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of adosolic-specific editing site
# p-value = 4.254e-13


## Are prenatal-specific editing sites enriched for DEG Age?
pren.UTR3.deg = lapply(SigList, function(x) pren.UTR3[which(pren.UTR3$collapsedconversion=="A:G / T:C" & pren.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(pren.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

pren.UTR5.deg = lapply(SigList, function(x) pren.UTR5[which(pren.UTR5$collapsedconversion=="A:G / T:C" & pren.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(pren.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

pren.Intron.deg = lapply(SigList, function(x) pren.Intron[which(pren.Intron$collapsedconversion=="A:G / T:C" & pren.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(pren.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

pren.CDS.deg = lapply(SigList, function(x) pren.CDS[which(pren.CDS$collapsedconversion=="A:G / T:C" & pren.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(pren.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

pren.Other.deg = lapply(SigList, function(x) pren.Other[which(pren.Other$collapsedconversion=="A:G / T:C" & pren.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(pren.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

## In adult: Are cytosolic-specific editing sites enriched for DEG Age?
ACnotAN.UTR3.deg = lapply(SigList, function(x) ACnotAN.UTR3[which(ACnotAN.UTR3$collapsedconversion=="A:G / T:C" & ACnotAN.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotAN.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotAN.UTR5.deg = lapply(SigList, function(x) ACnotAN.UTR5[which(ACnotAN.UTR5$collapsedconversion=="A:G / T:C" & ACnotAN.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotAN.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotAN.Intron.deg = lapply(SigList, function(x) ACnotAN.Intron[which(ACnotAN.Intron$collapsedconversion=="A:G / T:C" & ACnotAN.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotAN.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotAN.CDS.deg = lapply(SigList, function(x) ACnotAN.CDS[which(ACnotAN.CDS$collapsedconversion=="A:G / T:C" & ACnotAN.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotAN.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotAN.Other.deg = lapply(SigList, function(x) ACnotAN.Other[which(ACnotAN.Other$collapsedconversion=="A:G / T:C" & ACnotAN.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotAN.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

## In adult: Are nuclear-specific editing sites enriched for DEG Age?
ANnotAC.UTR3.deg = lapply(SigList, function(x) ANnotAC.UTR3[which(ANnotAC.UTR3$collapsedconversion=="A:G / T:C" & ANnotAC.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotAC.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotAC.UTR5.deg = lapply(SigList, function(x) ANnotAC.UTR5[which(ANnotAC.UTR5$collapsedconversion=="A:G / T:C" & ANnotAC.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotAC.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotAC.Intron.deg = lapply(SigList, function(x) ANnotAC.Intron[which(ANnotAC.Intron$collapsedconversion=="A:G / T:C" & ANnotAC.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotAC.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotAC.CDS.deg = lapply(SigList, function(x) ANnotAC.CDS[which(ANnotAC.CDS$collapsedconversion=="A:G / T:C" & ANnotAC.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotAC.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotAC.Other.deg = lapply(SigList, function(x) ANnotAC.Other[which(ANnotAC.Other$collapsedconversion=="A:G / T:C" & ANnotAC.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotAC.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

## In prenatal: Are cytosolic-specific editing sites enriched for DEG Age?
PCnotPN.UTR3.deg = lapply(SigList, function(x) PCnotPN.UTR3[which(PCnotPN.UTR3$collapsedconversion=="A:G / T:C" & PCnotPN.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotPN.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotPN.UTR5.deg = lapply(SigList, function(x) PCnotPN.UTR5[which(PCnotPN.UTR5$collapsedconversion=="A:G / T:C" & PCnotPN.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotPN.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotPN.Intron.deg = lapply(SigList, function(x) PCnotPN.Intron[which(PCnotPN.Intron$collapsedconversion=="A:G / T:C" & PCnotPN.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotPN.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotPN.CDS.deg = lapply(SigList, function(x) PCnotPN.CDS[which(PCnotPN.CDS$collapsedconversion=="A:G / T:C" & PCnotPN.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotPN.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotPN.Other.deg = lapply(SigList, function(x) PCnotPN.Other[which(PCnotPN.Other$collapsedconversion=="A:G / T:C" & PCnotPN.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotPN.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

## In prenatal: Are nuclear-specific editing sites enriched for DEG Age?
PNnotPC.UTR3.deg = lapply(SigList, function(x) PNnotPC.UTR3[which(PNnotPC.UTR3$collapsedconversion=="A:G / T:C" & PNnotPC.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotPC.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotPC.UTR5.deg = lapply(SigList, function(x) PNnotPC.UTR5[which(PNnotPC.UTR5$collapsedconversion=="A:G / T:C" & PNnotPC.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotPC.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotPC.Intron.deg = lapply(SigList, function(x) PNnotPC.Intron[which(PNnotPC.Intron$collapsedconversion=="A:G / T:C" & PNnotPC.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotPC.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotPC.CDS.deg = lapply(SigList, function(x) PNnotPC.CDS[which(PNnotPC.CDS$collapsedconversion=="A:G / T:C" & PNnotPC.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotPC.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotPC.Other.deg = lapply(SigList, function(x) PNnotPC.Other[which(PNnotPC.Other$collapsedconversion=="A:G / T:C" & PNnotPC.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotPC.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

## in Cytosol: Are adult-specific editing sites enriched for DEG Age?
ACnotPC.UTR3.deg = lapply(SigList, function(x) ACnotPC.UTR3[which(ACnotPC.UTR3$collapsedconversion=="A:G / T:C" & ACnotPC.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotPC.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotPC.UTR5.deg = lapply(SigList, function(x) ACnotPC.UTR5[which(ACnotPC.UTR5$collapsedconversion=="A:G / T:C" & ACnotPC.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotPC.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotPC.Intron.deg = lapply(SigList, function(x) ACnotPC.Intron[which(ACnotPC.Intron$collapsedconversion=="A:G / T:C" & ACnotPC.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotPC.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotPC.CDS.deg = lapply(SigList, function(x) ACnotPC.CDS[which(ACnotPC.CDS$collapsedconversion=="A:G / T:C" & ACnotPC.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotPC.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ACnotPC.Other.deg = lapply(SigList, function(x) ACnotPC.Other[which(ACnotPC.Other$collapsedconversion=="A:G / T:C" & ACnotPC.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ACnotPC.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13


## in Cytosol: Are prenatal-specific editing sites enriched for DEG Age?
PCnotAC.UTR3.deg = lapply(SigList, function(x) PCnotAC.UTR3[which(PCnotAC.UTR3$collapsedconversion=="A:G / T:C" & PCnotAC.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotAC.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotAC.UTR5.deg = lapply(SigList, function(x) PCnotAC.UTR5[which(PCnotAC.UTR5$collapsedconversion=="A:G / T:C" & PCnotAC.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotAC.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotAC.Intron.deg = lapply(SigList, function(x) PCnotAC.Intron[which(PCnotAC.Intron$collapsedconversion=="A:G / T:C" & PCnotAC.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotAC.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotAC.CDS.deg = lapply(SigList, function(x) PCnotAC.CDS[which(PCnotAC.CDS$collapsedconversion=="A:G / T:C" & PCnotAC.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotAC.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PCnotAC.Other.deg = lapply(SigList, function(x) PCnotAC.Other[which(PCnotAC.Other$collapsedconversion=="A:G / T:C" & PCnotAC.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PCnotAC.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

## in Nucleus: Are adult-specific editing sites enriched for DEG Age?
ANnotPN.UTR3.deg = lapply(SigList, function(x) ANnotPN.UTR3[which(ANnotPN.UTR3$collapsedconversion=="A:G / T:C" & ANnotPN.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotPN.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotPN.UTR5.deg = lapply(SigList, function(x) ANnotPN.UTR5[which(ANnotPN.UTR5$collapsedconversion=="A:G / T:C" & ANnotPN.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotPN.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotPN.Intron.deg = lapply(SigList, function(x) ANnotPN.Intron[which(ANnotPN.Intron$collapsedconversion=="A:G / T:C" & ANnotPN.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotPN.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotPN.CDS.deg = lapply(SigList, function(x) ANnotPN.CDS[which(ANnotPN.CDS$collapsedconversion=="A:G / T:C" & ANnotPN.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotPN.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

ANnotPN.Other.deg = lapply(SigList, function(x) ANnotPN.Other[which(ANnotPN.Other$collapsedconversion=="A:G / T:C" & ANnotPN.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(ANnotPN.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13


## in Nucleus: Are prenatal-specific editing sites enriched for DEG Age?
PNnotAN.UTR3.deg = lapply(SigList, function(x) PNnotAN.UTR3[which(PNnotAN.UTR3$collapsedconversion=="A:G / T:C" & PNnotAN.UTR3$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotAN.UTR3.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotAN.UTR5.deg = lapply(SigList, function(x) PNnotAN.UTR5[which(PNnotAN.UTR5$collapsedconversion=="A:G / T:C" & PNnotAN.UTR5$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotAN.UTR5.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotAN.Intron.deg = lapply(SigList, function(x) PNnotAN.Intron[which(PNnotAN.Intron$collapsedconversion=="A:G / T:C" & PNnotAN.Intron$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotAN.Intron.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotAN.CDS.deg = lapply(SigList, function(x) PNnotAN.CDS[which(PNnotAN.CDS$collapsedconversion=="A:G / T:C" & PNnotAN.CDS$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotAN.CDS.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13

PNnotAN.Other.deg = lapply(SigList, function(x) PNnotAN.Other[which(PNnotAN.Other$collapsedconversion=="A:G / T:C" & PNnotAN.Other$nearestID %in% rownames(x)),])
elementNROWS(SigList)
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#5634              4019              5137              3347
elementNROWS(lapply(PNnotAN.Other.deg, function(x) unique(x$nearestID)))
#Cres.DownPrenatal   Cres.UpPrenatal Nres.DownPrenatal   Nres.UpPrenatal 
#227               356               223               274 

fisher.test(data.frame(c(227,5634-227), c(356,4019-356))) # prenosol: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value < 2.2e-16
fisher.test(data.frame(c(223,5137-223), c(274,3347-274))) # nucleus: increasing or decreasing age DEG and presence or absence of prenosolic-specific editing site
# p-value = 4.254e-13


### does editing rate correlate with gene expression in the group the editing site appears?

AGonly.DEG = as.data.frame(cbind(AGonly, LFC.Adult = Apres[match(AGonly$nearestID, rownames(Apres)), "log2FoldChange"],
                                 padj.Adult = Apres[match(AGonly$nearestID, rownames(Apres)), "padj"],
                                 LFC.Prenatal = Fpres[match(AGonly$nearestID, rownames(Fpres)), "log2FoldChange"],
                                 padj.Prenatal = Fpres[match(AGonly$nearestID, rownames(Fpres)), "padj"]))
cyt = AGonly[Fraction=="Cytosol",,]
nuc = AGonly[Fraction=="Nucleus",,]
ad = AGonly[Age=="Adult",,]
pren = AGonly[Age=="Prenatal",,]
AC = AGonly[Group=="Adult:Cytosol",,]
AN = AGonly[Group=="Adult:Nucleus",,]
PC = AGonly[Group=="Prenatal:Cytosol",,]
PN = AGonly[Group=="Prenatal:Nucleus",,]


## correlate LFC with editing rate in editing rates shared between fraction

# Shared between nucleus and cytosol
cor(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% as.character(nuc$editingID)),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% as.character(nuc$editingID)),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% as.character(nuc$editingID)),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% as.character(nuc$editingID)),"LFC.Prenatal"], use = "complete.obs")

# Shared between nucleus and cytosol in adult
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% as.character(AN$editingID)),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% as.character(AN$editingID)),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% as.character(AN$editingID)),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% as.character(AN$editingID)),"LFC.Prenatal"], use = "complete.obs")

# Shared between nucleus and cytosol in prenatal
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% as.character(PN$editingID)),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% as.character(PN$editingID)),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% as.character(PN$editingID)),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% as.character(PN$editingID)),"LFC.Prenatal"], use = "complete.obs")


## correlate LFC with editing rate in editing sites unique to a group

ids = lapply(unique, function(x) as.character(x$editingID))
# cytosolOnly
cor(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% ids[["cytosolOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% ids[["cytosolOnly"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% ids[["cytosolOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% ids[["cytosolOnly"]]),"LFC.Prenatal"], use = "complete.obs")
# nucleusOnly
cor(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Nucleus" & AGonly.DEG$editingID %in% ids[["nucleusOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Fraction=="Nucleus" & AGonly.DEG$editingID %in% ids[["nucleusOnly"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Nucleus" & AGonly.DEG$editingID %in% ids[["nucleusOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Fraction=="Nucleus" & AGonly.DEG$editingID %in% ids[["nucleusOnly"]]),"LFC.Prenatal"], use = "complete.obs")
# adultOnly
cor(x = AGonly.DEG[which(AGonly.DEG$Age=="Adult" & AGonly.DEG$editingID %in% ids[["adultOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Age=="Adult" & AGonly.DEG$editingID %in% ids[["adultOnly"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Age=="Adult" & AGonly.DEG$editingID %in% ids[["adultOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Age=="Adult" & AGonly.DEG$editingID %in% ids[["adultOnly"]]),"LFC.Prenatal"], use = "complete.obs")
# prenatalOnly
cor(x = AGonly.DEG[which(AGonly.DEG$Age=="Prenatal" & AGonly.DEG$editingID %in% ids[["prenatalOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Age=="Prenatal" & AGonly.DEG$editingID %in% ids[["prenatalOnly"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Age=="Prenatal" & AGonly.DEG$editingID %in% ids[["prenatalOnly"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Age=="Prenatal" & AGonly.DEG$editingID %in% ids[["prenatalOnly"]]),"LFC.Prenatal"], use = "complete.obs")
# ANnotAC 
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotAC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotAC"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotAC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotAC"]]),"LFC.Prenatal"], use = "complete.obs")
# ACnotAN
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotAN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotAN"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotAN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotAN"]]),"LFC.Prenatal"], use = "complete.obs")
# ANnotPN 
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotPN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotPN"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotPN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotPN"]]),"LFC.Prenatal"], use = "complete.obs")
# PNnotAN
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotAN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotAN"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotAN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotAN"]]),"LFC.Prenatal"], use = "complete.obs")
# ACnotPC 
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotPC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotPC"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotPC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotPC"]]),"LFC.Prenatal"], use = "complete.obs")
# PCnotAC
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotAC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotAC"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotAC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotAC"]]),"LFC.Prenatal"], use = "complete.obs")
# PCnotPN
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotPN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotPN"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotPN"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotPN"]]),"LFC.Prenatal"], use = "complete.obs")
# PNnotPC
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotPC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotPC"]]),"LFC.Adult"], use = "complete.obs")
cor(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotPC"]]),"rate"],
    y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotPC"]]),"LFC.Prenatal"], use = "complete.obs")


### Is gene expression greater in the compartment/age exhibiting the editing site than in the compared group?

## t test of LFC between unique sites in both groups

# by fraction
t.test(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% ids[["cytosolOnly"]]),"LFC.Adult"],
       y = AGonly.DEG[which(AGonly.DEG$Fraction=="Nucleus" & AGonly.DEG$editingID %in% ids[["nucleusOnly"]]),"LFC.Adult"])
#t = -21.598, df = 9029, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2518054 -0.2099014
#sample estimates:
#  mean of x    mean of y 
#-0.001810615  0.229042752
t.test(x = AGonly.DEG[which(AGonly.DEG$Fraction=="Cytosol" & AGonly.DEG$editingID %in% ids[["cytosolOnly"]]),"LFC.Prenatal"],
       y = AGonly.DEG[which(AGonly.DEG$Fraction=="Nucleus" & AGonly.DEG$editingID %in% ids[["nucleusOnly"]]),"LFC.Prenatal"])
#t = -29.129, df = 9823.2, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1350161 -0.1179904
#sample estimates:
#  mean of x   mean of y 
#-0.04433652  0.08216672

# by fraction in adult
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotAN"]]),"LFC.Adult"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotAC"]]),"LFC.Adult"])
#t = -25.049, df = 5594.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3507228 -0.2998108
#sample estimates:
#  mean of x   mean of y 
#-0.08080479  0.24446203
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotAN"]]),"LFC.Prenatal"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotAC"]]),"LFC.Prenatal"])
#t = -24.673, df = 6010.6, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1387849 -0.1183539
#sample estimates:
#  mean of x   mean of y 
#-0.04841949  0.08014994

# by fraction in prenatal
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotPN"]]),"LFC.Adult"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotPC"]]),"LFC.Adult"])
#t = -8.3382, df = 5475.3, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.15082831 -0.09340605
#sample estimates:
#  mean of x  mean of y 
#0.06780659 0.18992377
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotPN"]]),"LFC.Prenatal"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotPC"]]),"LFC.Prenatal"])
#t = -17.865, df = 5949.7, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.12220644 -0.09803906
#sample estimates:
#  mean of x   mean of y 
#-0.04134723  0.06877552

# by age
t.test(x = AGonly.DEG[which(AGonly.DEG$Age=="Adult" & AGonly.DEG$editingID %in% ids[["adultOnly"]]),"LFC.Adult"],
       y = AGonly.DEG[which(AGonly.DEG$Age=="Prenatal" & AGonly.DEG$editingID %in% ids[["prenatalOnly"]]),"LFC.Adult"])
#t = -4.4131, df = 20237, p-value = 1.024e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.05149456 -0.01982028
#sample estimates:
#  mean of x mean of y 
#0.1102599 0.1459173
t.test(x = AGonly.DEG[which(AGonly.DEG$Age=="Adult" & AGonly.DEG$editingID %in% ids[["adultOnly"]]),"LFC.Prenatal"],
       y = AGonly.DEG[which(AGonly.DEG$Age=="Prenatal" & AGonly.DEG$editingID %in% ids[["prenatalOnly"]]),"LFC.Prenatal"])
#t = 3.9966, df = 20245, p-value = 6.448e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.006993799 0.020456280
#sample estimates:
#  mean of x  mean of y 
#0.02441815 0.01069311 

# by age in cytosol
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotAC"]]),"LFC.Adult"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotPC"]]),"LFC.Adult"])
#t = 12.944, df = 10242, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.1220085 0.1655559
#sample estimates:
#  mean of x  mean of y 
#0.1324222 -0.0113600
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Cytosol" & AGonly.DEG$editingID %in% ids[["PCnotAC"]]),"LFC.Prenatal"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Cytosol" & AGonly.DEG$editingID %in% ids[["ACnotPC"]]),"LFC.Prenatal"])
#t = 2.9228, df = 10366, p-value = 0.003477
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.004401435 0.022327306
#sample estimates:
#  mean of x   mean of y 
#-0.01388109 -0.02724546

# by age in nucleus
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotAN"]]),"LFC.Adult"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotPN"]]),"LFC.Adult"])
#t = -2.0137, df = 12699, p-value = 0.04406
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0404629667 -0.0005450153
#sample estimates:
#  mean of x mean of y 
#0.144934  0.165438
t.test(x = AGonly.DEG[which(AGonly.DEG$Group=="Prenatal:Nucleus" & AGonly.DEG$editingID %in% ids[["PNnotAN"]]),"LFC.Prenatal"],
       y = AGonly.DEG[which(AGonly.DEG$Group=="Adult:Nucleus" & AGonly.DEG$editingID %in% ids[["ANnotPN"]]),"LFC.Prenatal"])
#t = -3.7779, df = 12588, p-value = 0.0001589
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.025051555 -0.007935934
#sample estimates:
#  mean of x  mean of y 
#0.02931815 0.04581189


### In intronic editing sites, is the IR ratio greater in the compartment exhibiting the site than the compared group?

names = scan("./Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "./Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)
}
names(IRres) = shortenedNames
for (i in 1:length(IRres)){ 
  tmp = IRres[[i]]
  tmp$Chr = paste0("chr", tmp$Chr)
  IRres[[i]] = tmp
  }
IRres_ranges = lapply(IRres, function(x) makeGRangesFromDataFrame(x, seqnames.field = "Chr", 
                                                                  start.field = "Start",
                                                                  end.field = "End", 
                                                                  strand.field = "Direction",
                                                                  keep.extra.columns = T))
ov = lapply(IRres_ranges, function(x) findOverlaps(editing_ranges, x))
tog = list()
for (i in 1:length(IRres)){
  tmp = IRres[[i]]
  tog[[i]] = cbind(editing_anno[queryHits(ov[[i]]),], IRratio = tmp$IRratio[subjectHits(ov[[i]])],
                   sampleIDintron = names(IRres)[i])
}
tog = do.call(rbind, tog)

## t test of IR ratio of introns with editing sites by group    

# cytosol only vs nucleus only
t.test(x = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Fraction=="Cytosol" & tog$annotation=="Intron" & tog$editingID %in% ids[["cytosolOnly"]]),"IRratio"],
       y = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Fraction=="Nucleus" & tog$annotation=="Intron" & tog$editingID %in% ids[["nucleusOnly"]]),"IRratio"])
#t = 10.412, df = 27155, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.02451516 0.03588493
#sample estimates:
#  mean of x mean of y 
#0.245736  0.215536

# cytosol only vs nucleus only in adult
t.test(x = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Adult:Cytosol" & tog$annotation=="Intron" & tog$editingID %in% ids[["ACnotAN"]]),"IRratio"],
       y = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Adult:Nucleus" & tog$annotation=="Intron" & tog$editingID %in% ids[["ANnotAC"]]),"IRratio"])
#t = -7.7611, df = 13659, p-value = 9.015e-15
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.03575389 -0.02133546
#sample estimates:
#  mean of x mean of y 
#0.1898456 0.2183903

# cytosol only vs nucleus only in prenatal
t.test(x = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Prenatal:Cytosol" & tog$annotation=="Intron" & tog$editingID %in% ids[["PCnotPN"]]),"IRratio"],
       y = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Prenatal:Nucleus" & tog$annotation=="Intron" & tog$editingID %in% ids[["PNnotPC"]]),"IRratio"])
#t = 10.064, df = 22215, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.03168217 0.04700698
#sample estimates:
#  mean of x mean of y 
#0.2904688 0.2511242


# adult only vs prenatal only
t.test(x = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Age=="Adult" & tog$annotation=="Intron" & tog$editingID %in% ids[["adultOnly"]]),"IRratio"],
       y = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Age=="Prenatal" & tog$annotation=="Intron" & tog$editingID %in% ids[["prenatalOnly"]]),"IRratio"])
#t = -37.585, df = 116730, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.07793774 -0.07021209
#sample estimates:
#  mean of x mean of y 
#0.1892185 0.2632934 

# adult only vs prenatal only in cytosol
t.test(x = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Adult:Cytosol" & tog$annotation=="Intron" & tog$editingID %in% ids[["ACnotPC"]]),"IRratio"],
       y = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Prenatal:Cytosol" & tog$annotation=="Intron" & tog$editingID %in% ids[["PCnotAC"]]),"IRratio"])
#t = -24.442, df = 35809, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.09489925 -0.08080893
#sample estimates:
#  mean of x mean of y 
#0.2168146 0.3046687

# adult only vs prenatal only in nucleus
t.test(x = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Adult:Nucleus" & tog$annotation=="Intron" & tog$editingID %in% ids[["ANnotPN"]]),"IRratio"],
       y = tog[which(tog$collapsedconversion=="A:G / T:C" & tog$Group=="Prenatal:Nucleus" & tog$annotation=="Intron" & tog$editingID %in% ids[["PNnotAN"]]),"IRratio"])
#t = -22.329, df = 81832, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.05525209 -0.04633493
#sample estimates:
#  mean of x mean of y 
#0.1932618 0.2440553


### in 3'UTR editing sites, is the exon differentially expressed by fraction or group?

exonMap$exonID = rownames(exonMap)
exonMap_ranges = makeGRangesFromDataFrame(exonMap, keep.extra.columns = T)
ov = findOverlaps(editing_ranges, exonMap_ranges)
tog = cbind(editing_anno[queryHits(ov),], exonMap[subjectHits(ov),])
tog.3UTR = tog[which(tog$annotation=="UTR3"),]
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
lapply(exonres, function(x) dim(x[which(x$padj<=0.05),]))

tog.3UTR = cbind(tog.3UTR, A.LFC = exonsA.res[match(tog.3UTR$exonID, rownames(exonsA.res)),"log2FoldChange"],
                 A.padj = exonsA.res[match(tog.3UTR$exonID, rownames(exonsA.res)),"padj"],
                 P.LFC = exonsF.res[match(tog.3UTR$exonID, rownames(exonsF.res)),"log2FoldChange"],
                 P.padj = exonsF.res[match(tog.3UTR$exonID, rownames(exonsF.res)),"padj"],
                 C.LFC = exonsC.res[match(tog.3UTR$exonID, rownames(exonsC.res)),"log2FoldChange"],
                 C.padj = exonsC.res[match(tog.3UTR$exonID, rownames(exonsC.res)),"padj"],
                 N.LFC = exonsN.res[match(tog.3UTR$exonID, rownames(exonsN.res)),"log2FoldChange"],
                 N.padj = exonsN.res[match(tog.3UTR$exonID, rownames(exonsN.res)),"padj"])

## Compare LFC and significance of groups of editing sites in 3'UTR

# LFC between cytosol only and nucleus only 
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]]),"A.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]]),"A.LFC"])
#t = -20.725, df = 5993.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2676177 -0.2213662
#sample estimates:
#  mean of x   mean of y 
#-0.04224752  0.20224439
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]]),"P.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]]),"P.LFC"])
#t = -19.528, df = 6645.7, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.10307914 -0.08427196
#sample estimates:
#  mean of x   mean of y 
#-0.05419901  0.03947654
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]]),"C.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]]),"C.LFC"])
#t = -0.02454, df = 5758.9, p-value = 0.9804
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.06697281  0.06531684
#sample estimates:
#  mean of x mean of y 
#0.1549564 0.1557844
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]]),"N.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]]),"N.LFC"])
#t = 4.3936, df = 5861.1, p-value = 1.134e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.07409559 0.19348591
#sample estimates:
#  mean of x  mean of y 
#0.1126460 -0.0211448

# LFC between cytosol only and nucleus only in adult
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]]),"A.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]]),"A.LFC"])
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3638747 -0.3107056
#sample estimates:
#  mean of x  mean of y 
#-0.1270911  0.2101990
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]]),"P.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]]),"P.LFC"])
#t = -17.42, df = 4677, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1065198 -0.0849698
#sample estimates:
#  mean of x   mean of y 
#-0.05411592  0.04162890
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]]),"C.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]]),"C.LFC"])
#t = -8.8735, df = 4281.8, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3565655 -0.2275174
#sample estimates:
#  mean of x   mean of y 
#-0.30261834 -0.01057685 
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]]),"N.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]]),"N.LFC"])
#t = -1.5209, df = 4051.4, p-value = 0.1284
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.10773703  0.01360584
#sample estimates:
#  mean of x  mean of y 
#-0.2418584 -0.1947928

# LFC between cytosol only and nucleus only in prenatal

t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]]),"A.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]]),"A.LFC"])
#t = -7.0904, df = 3371.9, p-value = 1.621e-12
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.14440399 -0.08184188
#sample estimates:
#  mean of x  mean of y 
#0.06419681 0.17731974
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]]),"P.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]]),"P.LFC"])
#t = -13.369, df = 3481.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.10415322 -0.07751018
#sample estimates:
#  mean of x   mean of y 
#-0.05733782  0.03349388
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]]),"C.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]]),"C.LFC"])
#t = 1.9146, df = 3405.1, p-value = 0.05563
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.002072212  0.174301025
#sample estimates:
#  mean of x mean of y 
#0.7066692 0.6205548
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]]),"N.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]]),"N.LFC"])
#t = 2.0179, df = 3440.8, p-value = 0.04368
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.002274204 0.158174215
#sample estimates:
#  mean of x mean of y 
#0.5345192 0.4542950

# LFC between adult and prenatal only
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]]),"A.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]]),"A.LFC"])
#t = -19.584, df = 10008, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1871397 -0.1530851
#sample estimates:
#  mean of x   mean of y 
#0.003193738 0.173306155 
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]]),"P.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]]),"P.LFC"])
#t = 3.1428, df = 9761.9, p-value = 0.001679
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.004563779 0.019693498
#sample estimates:
#  mean of x    mean of y 
#-0.005400458 -0.017529096
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]]),"C.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]]),"C.LFC"])
#t = -67.267, df = 7601.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.920432 -1.811672
#sample estimates:
#  mean of x  mean of y 
#-0.4877324  1.3783193
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]]),"N.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]]),"N.LFC"])
#t = -66.34, df = 8017.8, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.691864 -1.594749
#sample estimates:
#  mean of x  mean of y 
#-0.5006527  1.1426537

# LFC between adult and prenatal only in cytosol
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]]),"A.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]]),"A.LFC"])
#t = -23.411, df = 6304.6, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.2918628 -0.2467609
#sample estimates:
#  mean of x   mean of y 
#-0.08777889  0.18153293
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]]),"P.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]]),"P.LFC"])
#t = -2.4179, df = 6341.6, p-value = 0.01564
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.021885067 -0.002287071
#sample estimates:
#  mean of x   mean of y 
#-0.02976957 -0.01768350 
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]]),"C.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]]),"C.LFC"])
#t = -53.593, df = 5396.2, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.847335 -1.716955
#sample estimates:
#  mean of x  mean of y 
#-0.5237457  1.2583993
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]]),"N.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]]),"N.LFC"])
#t = -49.432, df = 5759.5, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.538207 -1.420856
#sample estimates:
#  mean of x  mean of y 
#-0.4708156  1.0087154

# LFC between adult and prenatal only in nucleus
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]]),"A.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]]),"A.LFC"])
#t = -8.7875, df = 5748.6, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.12012272 -0.07630259
#sample estimates:
#  mean of x  mean of y 
#0.06278712 0.16099977
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]]),"P.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]]),"P.LFC"])
#t = 2.4273, df = 5687.2, p-value = 0.01524
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.002316646 0.021770384
#sample estimates:
#  mean of x    mean of y 
#0.006403311 -0.005640204 
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]]),"C.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]]),"C.LFC"])
#t = -46.789, df = 4409.4, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.706939 -1.569648
#sample estimates:
#  mean of x  mean of y 
#-0.3954073  1.2428860
t.test(x=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]]),"N.LFC"],
       y=tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]]),"N.LFC"])
#t = -48.008, df = 4611.4, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.56127 -1.43876
#sample estimates:
#  mean of x  mean of y 
#-0.4587846  1.0412301

# Compare the number of 3'UTRs with and without an editing site and whether or not they are significantly DE

length(unique(tog.3UTR$editingID)) # 9837
length(unique(tog.3UTR$exonID)) # 5133
head(tog.3UTR)

# between cytosol-only and nucleus-only, compare proportion of differentially expressed 3'UTRs by fraction in adult
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]] & tog.3UTR$A.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]] & tog.3UTR$A.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]] & tog.3UTR$A.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]] & tog.3UTR$A.padj<=0.05),"exonID"])))))
#p-value = 0.6502
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8577487 1.2873497
#sample estimates:
#  odds ratio 
#1.050377

# between cytosol-only and nucleus-only, compare proportion of differentially expressed 3'UTRs by fraction in prenatal
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]] & tog.3UTR$P.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Cytosol" & tog.3UTR$editingID %in% ids[["cytosolOnly"]] & tog.3UTR$P.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]] & tog.3UTR$P.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Fraction=="Nucleus" & tog.3UTR$editingID %in% ids[["nucleusOnly"]] & tog.3UTR$P.padj<=0.05),"exonID"])))))
#p-value = 0.01625
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.248422 371.081454
#sample estimates:
#  odds ratio 
#8.621243

# in adult: between cytosol-only and nucleus-only, compare proportion of differentially expressed 3'UTRs by fraction in adult
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]] & tog.3UTR$A.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]] & tog.3UTR$A.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]] & tog.3UTR$A.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]] & tog.3UTR$A.padj<=0.05),"exonID"])))))
#p-value = 0.3465
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8900516 1.3945858
#sample estimates:
#  odds ratio 
#1.113256

# in adult: between cytosol-only and nucleus-only, compare proportion of differentially expressed 3'UTRs by fraction in prenatal
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]] & tog.3UTR$P.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotAN"]] & tog.3UTR$P.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]] & tog.3UTR$P.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotAC"]] & tog.3UTR$P.padj<=0.05),"exonID"])))))
#p-value = 0.6459
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.170731 112.063929
#sample estimates:
#  odds ratio 
#2.131672

# in prenatal: between cytosol-only and nucleus-only, compare proportion of differentially expressed 3'UTRs by fraction in adult
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]] & tog.3UTR$A.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]] & tog.3UTR$A.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]] & tog.3UTR$A.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]] & tog.3UTR$A.padj<=0.05),"exonID"])))))
#p-value = 0.2559
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.6541991 1.1241311
#sample estimates:
#  odds ratio 
#0.8575531

# in prenatal: between cytosol-only and nucleus-only, compare proportion of differentially expressed 3'UTRs by fraction in prenatal
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]] & tog.3UTR$P.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotPN"]] & tog.3UTR$P.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]] & tog.3UTR$P.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotPC"]] & tog.3UTR$P.padj<=0.05),"exonID"])))))
#p-value = 0.04178
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9488636 316.0353396
#sample estimates:
#  odds ratio 
#7.112278


# between adult-only and prenatal-only, compare proportion of differentially expressed 3'UTRs by age in cytosol
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]] & tog.3UTR$C.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]] & tog.3UTR$C.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]] & tog.3UTR$C.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]] & tog.3UTR$C.padj<=0.05),"exonID"])))))
#p-value = 0.01603
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.034176 1.413109
#sample estimates:
#  odds ratio 
#1.208582

# between adult-only and prenatal-only, compare proportion of differentially expressed 3'UTRs by age in nucleus
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]] & tog.3UTR$N.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Adult" & tog.3UTR$editingID %in% ids[["adultOnly"]] & tog.3UTR$N.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]] & tog.3UTR$N.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Age=="Prenatal" & tog.3UTR$editingID %in% ids[["prenatalOnly"]] & tog.3UTR$N.padj<=0.05),"exonID"])))))
#p-value = 0.1903
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.9503293 1.2937844
#sample estimates:
#  odds ratio 
#1.108671

# in cytosol: between adult-only and prenatal-only, compare proportion of differentially expressed 3'UTRs by age in cytosol
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]] & tog.3UTR$C.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]] & tog.3UTR$C.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]] & tog.3UTR$C.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]] & tog.3UTR$C.padj<=0.05),"exonID"])))))
#p-value = 0.01327
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.045332 1.483054
#sample estimates:
#  odds ratio 
#1.24477

# in cytosol: between adult-only and prenatal-only, compare proportion of differentially expressed 3'UTRs by age in nucleus
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]] & tog.3UTR$N.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Cytosol" & tog.3UTR$editingID %in% ids[["ACnotPC"]] & tog.3UTR$N.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]] & tog.3UTR$N.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Cytosol" & tog.3UTR$editingID %in% ids[["PCnotAC"]] & tog.3UTR$N.padj<=0.05),"exonID"])))))
#p-value = 0.4903
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.8959501 1.2659675
#sample estimates:
#  odds ratio 
#1.064832

# in nucleus: between adult-only and prenatal-only, compare proportion of differentially expressed 3'UTRs by age in cytosol
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]] & tog.3UTR$C.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]] & tog.3UTR$C.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]] & tog.3UTR$C.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]] & tog.3UTR$C.padj<=0.05),"exonID"])))))
#p-value = 0.04166
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.003301 1.401989
#sample estimates:
#  odds ratio 
#1.185712

# in nucleus: between adult-only and prenatal-only, compare proportion of differentially expressed 3'UTRs by age in nucleus
fisher.test(data.frame(c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]] & tog.3UTR$N.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Adult:Nucleus" & tog.3UTR$editingID %in% ids[["ANnotPN"]] & tog.3UTR$N.padj<=0.05),"exonID"]))),
                       c(length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]] & tog.3UTR$N.padj>0.05),"exonID"])),
                         length(unique(tog.3UTR[which(tog.3UTR$collapsedconversion=="A:G / T:C" & tog.3UTR$Group=="Prenatal:Nucleus" & tog.3UTR$editingID %in% ids[["PNnotAN"]] & tog.3UTR$N.padj<=0.05),"exonID"])))))
#p-value = 0.7105
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.875978 1.218833
#sample estimates:
#  odds ratio 
#1.033162



## Are cytosolic/nuclear/adult/prenatal-specific editing sites enriched for DEG Fraction/Age overall, and broken down by group?
# retained or exported DEG and presence or absence of cytosolic-specific editing site
# sites within retained or exported DEG and cytosolic-specific or non-specific status

### Check gene enrichment by annotation,
# retained or exported DEG and presence or absence of cytosolic-specific editing site
# Is a nuclear-specific editing site higher expressed in nucleus, for both intronic and other annotations?
# Is there a relationship between having an editing site and expression by fraction and age?
# is this affected by where in the gene the editing site falls (annotation)?
