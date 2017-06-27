library(GenomicRanges)
library(ggplot2)
library(plyr)
library(data.table)

### Prepare Intron Lists

# Load IRFinder Results
names = scan("/Users/amanda/Dropbox/NucVsCytosol/names.txt", what = "character")
shortenedNames = unique(gsub( "_.*$", "", names))
path = "/Users/amanda/Dropbox/sorted_figures/IRfinder/"
IRres = list()
for (i in 1:length(shortenedNames)){
  IRres[[i]] = read.table(paste0(path,"PolyA/",shortenedNames[i],"/IRFinder-IR-nondir.txt"), header = TRUE)}
names(IRres) = shortenedNames
lapply(IRres, head)
# Filter introns
allIntrons = IRres[[1]]
allIntrons = allIntrons[grep("clean", allIntrons$GeneIntronDetails, fixed=T),]
string = unlist(strsplit(as.character(allIntrons$GeneIntronDetails), "/", fixed = TRUE), recursive = FALSE)
allIntrons = data.frame(allIntrons[,c(1:4,6)], genes = string[grep("ENSG", string)],
                        intronID = paste0("chr",allIntrons$Chr,":",allIntrons$Start,"-",allIntrons$End))
# read in Differential IR results files
path = "./Dropbox/sorted_figures/IRfinder/"
comps = c("Adult_PolyA_Zone","Fetal_PolyA_Zone","Cytosol_PolyA_Age","Nuclear_PolyA_Age","PolyA_Zone","PolyA_Age")
nonconst = list()
for (i in 1:length(comps)){
  nonconst[[i]] = read.table(paste0(path, "PolyA/",comps[i],"_nonconst.tab"), header = TRUE, comment.char="#")
}
names(nonconst) = comps
elementNROWS(nonconst)
string = lapply(nonconst, function(y) unlist(strsplit(as.character(y$Intron.GeneName.GeneID),"/", fixed = TRUE),recursive = FALSE))
IR.diff = lapply(nonconst, function(y) y$A.IRratio - y$B.IRratio)
dIR = Map(cbind, nonconst, ensID = lapply(string, function(y) y[grep("ENSG", y)]), 
          comments = lapply(string, function(y) as.character(y[seq.int(from = 3, to=length(y), by=3)])), 
          IR.diff = IR.diff, Sign = lapply(IR.diff, function(y) ifelse(y < 0,"MoreIRInNuc.Fetal", "MoreIRInCyt.Adult")), 
          intronID = lapply(nonconst, function(x) paste0("chr",x$Chr,":",x$Start,"-",x$End)))
dIRclean = lapply(dIR, function(y) y[which(y$A.warnings!="LowCover" & y$A.warnings!="LowSplicing" & y$A.warnings!="NonUniformIntronCover" &
                                           y$B.warnings!="LowCover" & y$B.warnings!="LowSplicing" & y$B.warnings!="NonUniformIntronCover" &
                                           y$comments=="clean"),])
sigdIR = lapply(dIRclean, function(x) x[which(x$p.diff<=0.05),])
sigdIR = lapply(sigdIR, function(x) split(x, x$Sign))
sigdIR = unlist(sigdIR, recursive=F)
names(sigdIR) = c("Adult:Cytosol-Increased","Adult:Nucleus-Increased","Prenatal:Cytosol-Increased","Prenatal:Nucleus-Increased",
                  "Cytosol:Adult-Increased","Cytosol:Prenatal-Increased","Nucleus:Adult-Increased","Nucleus:Prenatal-Increased",
                  "Cytosol-Increased", "Nucleus-Increased", "Adult-Increased", "Prenatal-Increased")
sigdIR = c(lapply(sigdIR[1:8], function(x) x[,c(1:4,6,8,30,32,33,34)]),lapply(sigdIR[9:12], function(x) x[,c(1:4,6,8,36,38,39,40)]))
introns = c("All Introns" = list(allIntrons), "Introns (Fraction)" = list(do.call(rbind, lapply(dIRclean[1:2], function(x) x[,c(1:4,6,8,30,32,33,34)]))), 
            "Introns (Age)" = list(do.call(rbind, lapply(dIRclean[3:4], function(x) x[,c(1:4,6,8,30,32,33,34)]))), sigdIR)
introns[["Introns (Fraction)"]] = introns[["Introns (Fraction)"]][unique(introns[["Introns (Fraction)"]][,"intronID"]),]
introns[["Introns (Age)"]] = introns[["Introns (Age)"]][unique(introns[["Introns (Age)"]][,"intronID"]),]

# For RNA Binding protein analysis, write coordinates in a way recognized by RBPMap
rnacoord = do.call(rbind, lapply(dIRclean, function(x) x[,c(1:3,6)]))
rnacoord = rnacoord[!duplicated(rnacoord), ]
rnacoord$Chr = paste0("chr",rnacoord$Chr)
rnacoord$length = abs(rnacoord$End - rnacoord$Start)
splitlonger1 = splitlonger2 = rnacoord[which(rnacoord$length>10000 & rnacoord$length<20000),]
splitlonger1$End = splitlonger1$End - ceiling(splitlonger1$length/2) 
splitlonger2$Start = splitlonger2$Start + ceiling(splitlonger2$length/2)
longer1 = longer2 = longer3 = rnacoord[which(rnacoord$length>20000),]
longer1$End = longer1$End - 2*ceiling(longer1$length/3) 
longer2$Start = longer2$End
longer2$End = longer2$End - ceiling(longer2$length/3)
longer3$Start = longer3$End - ceiling(longer3$length/3)
rnacoord = rbind(rnacoord, splitlonger1, splitlonger2, longer1, longer2, longer3)
rnacoord$length = abs(rnacoord$End - rnacoord$Start)
rnacoord = rnacoord[which(rnacoord$length<=10000),]
rnacoord$ID = paste0(rnacoord$Chr,":",rnacoord$Start,"-",rnacoord$End,":",rnacoord$Direction)
dim(rnacoord) # 2254 introns
write.table(rnacoord[,"ID"], file="./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/intron.coordinates.RBPmap.format.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)

###############################################################
######################### RNA binding #########################
###############################################################

### Check RNA binding protein motifs present in introns

# After running the code in split_RBPMap_results_byAsterisk.txt on server 4, scan in names and create data.frame of RNA binding proteins

setwd("cd /media/DATA/Amanda/RBPMap/intronRBPbyprotein")
ids = scan("../introns.txt", what = "character")
proteinIDs = proteins = list()
for (i in 1:length(ids)){
  proteinIDs[[i]] = scan(paste0(ids[i], "_ids.txt"), what="character")
}
proteinIDs = unlist(proteinIDs)
split = strsplit(proteinIDs, "_", fixed = T)

for (i in grep("No", proteinIDs, invert = T)){
  proteins[[i]] = read.table(proteinIDs[i], header = T, skip = 1)
  proteins[[i]][,"ID"] = proteinIDs[i]
  proteins[[i]][,"chromosome"] = split[[i]][1]
  proteins[[i]][,"interval"] = split[[i]][2]
  proteins[[i]][,"strand"] = split[[i]][3]
  proteins[[i]][,"proteinID"] = gsub("().txt","", split[[i]][6], fixed = T)
  proteins[[i]][,"intronID"] = paste0(split[[i]][1],":",split[[i]][2])
}
names(proteins) = proteinIDs
allrbp = do.call(rbind, proteins)

rbp = as.data.table(allrbp)
rbp = rbp[rbp[,.I[P.value == min(P.value)], by=ID]$V1]
rbp = rbp[rbp[,.I[Sequence_Position == min(Sequence_Position)], by=ID]$V1]
split.nos = strsplit(proteinIDs[grep("No", proteinIDs)], "_", fixed = T)
rbp = rbind(rbp, data.frame("Sequence_Position" = NA,"Genomic_Coordinate" = NA,"Motif" = NA,"K.mer" = NA,"Z.score" = NA,"P.value" = NA,
                            "ID" = proteinIDs[grep("No", proteinIDs)],
                            "chromosome" = unlist(lapply(split.nos, function(x) x[1])),
                            "interval" = unlist(lapply(split.nos, function(x) x[2])),          
                            "strand" = unlist(lapply(split.nos, function(x) x[3])),
                            "proteinID" = gsub(".txt","", unlist(lapply(split.nos, function(x) x[6])), fixed = T),
                            "intronID" = paste0(unlist(lapply(split.nos, function(x) x[1])),":",unlist(lapply(split.nos, function(x) x[2])))))
save(proteinIDs, proteins, allrbp, rbp, file="/media/DATA/Amanda/RBPMap/intronRBPbyprotein/RBP_objects.rda")

# load RBP results locally, and compare binding motifs by intron localization
load("./Dropbox/sorted_figures/new/github_controlled/intron_retention/data/intron_IR_comparisons/RBPMap/RBP_objects.rda")

### Do differentially retained introns by fraction have more RBP sites?
# Map RBP table to differentially retained intron list
rbpIntrons = lapply(introns, function(x) cbind(x, ))

### Do differentially retained introns by age have more RBP sites?


### What is the distribution of RBPs by different groups of introns?
length(unique(rbp$proteinID)) # 95 represented in entire list

### Any RBPs stick out?







  # more information can be found at http://rbpmap.technion.ac.il/