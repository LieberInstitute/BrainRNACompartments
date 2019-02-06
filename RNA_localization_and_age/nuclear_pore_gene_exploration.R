library(readxl)
library(DESeq2)
library(GenomicRanges)
library("ggplot2")


load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


## which NPC genes are differentially expressed by fraction or age?
# (Nucleoporin genes from http://amigo.geneontology.org/grebe)

nups = read.table("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/nuclear_pore_genes2.txt", sep = "\t", header = F)
unique(nups$V5)
nupGenes = as.character(unique(nups$V3))
geneMap[which(geneMap$ensemblID=="ENSG00000272391"),"Symbol"] = "POM121C"
geneMap[which(geneMap$ensemblID=="ENSG00000139496"),"Symbol"] = "NUP58"
geneMap[which(geneMap$ensemblID=="ENSG00000273398"),"Symbol"] = "H7BYZ3"
geneMap[which(geneMap$ensemblID=="ENSG00000251184"),"Symbol"] = "H7C4K7"
nups = geneMap[which(geneMap$Symbol %in% nupGenes),]

DEnups.Frac = lapply(sig, function(x) x[which(x$geneID %in% nups$gencodeID),])[which(elementNROWS(lapply(sig, function(x) x[which(x$geneID %in% nups$gencodeID),]))>0)]
elementNROWS(DEnups.Frac)
#Ad_retained Ad_exported interacting 
#          8          11           2 

DEnups.Age = lapply(age.sig, function(x) x[which(x$geneID %in% nups$gencodeID),])[which(elementNROWS(lapply(age.sig, function(x) x[which(x$geneID %in% nups$gencodeID),]))>0)]
elementNROWS(DEnups.Age)
#both_decreasing both_increasing  Cyt_decreasing  Nuc_decreasing  Cyt_increasing  Nuc_increasing     interacting 
#             31               7               7               4               1               1               2 


write.csv(do.call(rbind, Map(cbind, Group = as.list(names(DEnups.Frac)), DEnups.Frac)), quote = F,
          file = "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/nuclear_pore_genes_fraction_expr.csv")

write.csv(do.call(rbind, Map(cbind, Group = as.list(names(DEnups.Age)), DEnups.Age)), quote = F,
          file = "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/nuclear_pore_genes_age_expr.csv")


## are the NPC genes retained or exported in adult increasing or decreasing over aging?

lapply(DEnups.Frac, function(x) as.character(x$geneID))
elementNROWS(DEnups.Age)
elementNROWS(lapply(DEnups.Age, function(a) a[which(as.character(a$geneID) %in% as.character(DEnups.Frac$Ad_retained$geneID)),]))
# both_decreasing both_increasing  Cyt_decreasing  Nuc_decreasing  Cyt_increasing  Nuc_increasing     interacting 
#               4               0               2               0               0               1               2 

elementNROWS(lapply(DEnups.Age, function(a) a[which(as.character(a$geneID) %in% as.character(DEnups.Frac$Ad_exported$geneID)),]))
# both_decreasing both_increasing  Cyt_decreasing  Nuc_decreasing  Cyt_increasing  Nuc_increasing     interacting 
#               1               1               0               3               0               0               0 


# Which are increasing?

c(as.character(DEnups.Age$both_increasing$Symbol), as.character(DEnups.Age$Cyt_increasing$Symbol), as.character(DEnups.Age$Nuc_increasing$Symbol))
# "MVP"     "RANGAP1" "TMEM33"  "EIF5A2"  "SENP2"   "RANBP3L" "MX2"     "SEH1L"   "NPIPA1"

do.call(rbind, DEnups.Age[grep("increasing", names(DEnups.Age))])
#               geneID   baseMean Cytosol.LFC Cytosol.SE Cytosol.padj Nucleus.LFC Nucleus.SE Nucleus.padj           ensID
# ENSG00000013364.18_2  116.73180  -2.0543641  0.3038078 1.187125e-10  -2.1284984  0.3007880 1.619316e-11 ENSG00000013364
# ENSG00000100401.19_2 2077.55150  -1.6310815  0.1603928 6.351213e-23  -1.2742009  0.1704878 9.581809e-13 ENSG00000100401
# ENSG00000109133.12_1  537.90908  -0.7652196  0.1896346 2.034263e-04  -0.6755769  0.2074309 3.679493e-03 ENSG00000109133
#  ENSG00000163577.7_1  185.22555  -2.4825703  0.2515881 1.217457e-21  -2.4544543  0.2905638 4.858780e-16 ENSG00000163577
# ENSG00000163904.12_2  548.88029  -1.0271569  0.1789380 6.015841e-08  -0.9806770  0.1909823 1.697680e-06 ENSG00000163904
#  ENSG00000164188.8_2   53.94444  -2.9254009  0.6925176 9.474902e-05  -2.1915462  0.5563918 3.367397e-04 ENSG00000164188
# ENSG00000183486.12_2   10.38979  -3.7311926  0.8473127 4.432310e-05  -2.9133135  0.8627760 2.493087e-03 ENSG00000183486
# ENSG00000085415.15_1  593.21746  -0.5253555  0.1912986 1.485254e-02  -0.5652341  0.2554611 6.179550e-02 ENSG00000085415
# ENSG00000183426.16_2  205.28277  -0.1034150  0.4618783 8.765670e-01  -1.0474395  0.3221878 3.750040e-03 ENSG00000183426

#  Symbol  EntrezID           Type
#     MVP      9961 protein_coding
# RANGAP1      5905 protein_coding
#  TMEM33     55161 protein_coding
#  EIF5A2     56648 protein_coding
#   SENP2     59343 protein_coding
# RANBP3L    202151 protein_coding
#     MX2      4600 protein_coding
#   SEH1L     81929 protein_coding
#  NPIPA1 101059953 protein_coding

## Add this info to Table S3

age = read.csv("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/nuclear_pore_genes_age_expr.csv")
df = do.call(rbind, Map(cbind, age.sig, Group = as.list(names(age.sig))))
df$Nup = ifelse(df$geneID %in% age$geneID, "NPC Gene", "")
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/updated_table_S3.csv")


# Which are decreasing?

unique(as.character(age[which(age$Group %in% c("both_decreasing", "Cyt_decreasing", "Nuc_decreasing")),"Symbol"]))
length(unique(as.character(age[which(age$Group %in% c("both_decreasing", "Cyt_decreasing", "Nuc_decreasing")),"Symbol"])))
# 42




