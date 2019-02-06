library(readxl)
library(DESeq2)
library("ggplot2")


load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")


## which RBP genes are differentially expressed by fraction or age?
# (RBP genes from https://www.biorxiv.org/content/early/2018/02/22/269043)

rbps = scan(file="./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/RBP_list.txt",
            sep=",", what="character")
rbps = gsub(" ", "", rbps)

geneMap[which(geneMap$ensemblID=="ENSG00000278053"),"Symbol"] = "DDX52"
geneMap[which(geneMap$ensemblID=="ENSG00000278619"),"Symbol"] = "MRM1"
geneMap[which(geneMap$ensemblID=="ENSG00000278845"),"Symbol"] = "MRPL45"
geneMap[which(geneMap$ensemblID=="ENSG00000247315"),"Symbol"] = "ZCCHC3"
geneMap[which(geneMap$ensemblID=="ENSG00000275700"),"Symbol"] = "AATF"
geneMap[which(geneMap$ensemblID=="ENSG00000091436"),"Symbol"] = "MLTK"
geneMap[which(geneMap$ensemblID=="ENSG00000115128"),"Symbol"] = "SF3B14"
geneMap[which(geneMap$ensemblID=="ENSG00000274523"),"Symbol"] = "WBSCR16"
geneMap[which(geneMap$ensemblID=="ENSG00000180574"),"Symbol"] = "EIF2S3L"
geneMap[which(geneMap$ensemblID=="ENSG00000275183"),"Symbol"] = "LENG9"
geneMap[which(geneMap$ensemblID=="ENSG00000266472"),"Symbol"] = "MRPS21"
geneMap[which(geneMap$ensemblID=="ENSG00000265241"),"Symbol"] = "RBM8A"

not = rbps[-which(rbps %in% geneMap$Symbol)]
not
"CHD2"

rbps = geneMap[which(geneMap$Symbol %in% rbps),]
dim(rbps)

DErbps.Frac = lapply(sig, function(x) x[which(x$geneID %in% rbps$gencodeID),])[which(elementNROWS(lapply(sig, function(x) x[which(x$geneID %in% rbps$gencodeID),]))>0)]
elementNROWS(DErbps.Frac)
#both_retained both_exported  Fet_retained   Ad_retained   Ad_exported   interacting 
#            6             3             2            50           311            29 

DErbps.Age = lapply(age.sig, function(x) x[which(x$geneID %in% rbps$gencodeID),])[which(elementNROWS(lapply(age.sig, function(x) x[which(x$geneID %in% rbps$gencodeID),]))>0)]
elementNROWS(DErbps.Age)
#  both_decreasing   both_increasing    Cyt_decreasing    Nuc_decreasing    Cyt_increasing 
#              340                61                62               100                59 
#Nuc_increasing decr_Nuc_incr_Cyt       interacting 
#            16                 1                29 


write.csv(do.call(rbind, Map(cbind, Group = as.list(names(DErbps.Frac)), DErbps.Frac)), quote = F,
          file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/RBP_DEGs_fraction_expr.csv")

write.csv(do.call(rbind, Map(cbind, Group = as.list(names(DErbps.Age)), DErbps.Age)), quote = F,
          file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/RBP_DEGs_age_expr.csv")


# Which are increasing?

c(as.character(DErbps.Age$both_increasing$Symbol), as.character(DErbps.Age$Cyt_increasing$Symbol), as.character(DErbps.Age$Nuc_increasing$Symbol))

do.call(rbind, DErbps.Age[grep("increasing", names(DErbps.Age))])