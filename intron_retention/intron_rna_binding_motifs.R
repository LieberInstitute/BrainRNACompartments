library(GenomicRanges)
library(ggplot2)
library(plyr)

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






  # more information can be found at http://rbpmap.technion.ac.il/
