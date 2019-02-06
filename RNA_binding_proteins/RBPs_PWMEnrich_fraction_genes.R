library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PWMEnrich)
library("VennDiagram")
library(clusterProfiler)
library(org.Hs.eg.db)


load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
# load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/retained.byAge.downsampled.rda")


### For RNA Binding protein analysis, create DNA strings for PWM analysis

genes = makeGRangesFromDataFrame(geneMap[!geneMap$Chr %in% c("GL000192.1", "GL000193.1", "GL000195.1", 
                                                             "GL000199.1", "GL000202.1", "GL000204.1", 
                                                             "GL000205.1", "GL000212.1", "GL000220.1", 
                                                             "GL000228.1", "GL000241.1"),], keep.extra.columns = T)
frac = lapply(sig[elementNROWS(sig)>0], function(x) makeGRangesFromDataFrame(geneMap[which(geneMap$gencodeID %in% x$geneID),], keep.extra.columns = T))
fracSeq = lapply(frac, function(x) getSeq(Hsapiens, x))
fracSeq = lapply(fracSeq, reverseComplement)
geneSeq = getSeq(Hsapiens, genes)
geneSeq = reverseComplement(geneSeq)
save(geneSeq, fracSeq, file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/makegenePWMbcgnd.rda")


## load RBP database from ATtRACT

rbpdb = read.table("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT/ATtRACT_db.txt", header=T, sep = "\t")
rbpdb = read.table("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT/ATtRACT_db.txt", header=T, sep = "\t")

rbpdb = rbpdb[rbpdb$Organism=="Homo_sapiens",]
dim(rbpdb) # 3256 RBP motifs included
ids = as.character(unique(rbpdb$Gene_name))
length(ids) # 160 RBPs reflected

## load motif databases

pwm = read.table("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT/pwm.txt", header = FALSE, 
                 sep = "\t", col.names = paste0("V",seq_len(4)), fill = TRUE)
pwm = read.table("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT/pwm.txt", header = FALSE, 
                 sep = "\t", col.names = paste0("V",seq_len(4)), fill = TRUE)

colnames(pwm) = c("A","C","G","T")
ast = grep(">", pwm$A, fixed = T)
sep = list()
for (i in 1:length(ast)) { sep[[i]] = pwm[c(ast[i]:(ast[i]+pwm[ast[i],"C"])),] }
names(sep) = lapply(sep, function(x) gsub(">","",as.character(x[1,1])))
sep = lapply(sep, function(x) x[2:nrow(x),])
for (i in 1:length(sep)) {
  tmp = sep[[i]]
  tmp$A = as.numeric(as.character(tmp$A))
  tmp$C = as.numeric(as.character(tmp$C))
  tmp$G = as.numeric(as.character(tmp$G))
  tmp$'T' = as.numeric(as.character(tmp$'T'))
  sep[[i]] = t(tmp)
  colnames(sep[[i]]) = NULL
}
sep = sep[names(sep) %in% rbpdb$Matrix_id]
head(sep)
pos = lapply(sep, function(x) round(x*10000))
pos = lapply(pos, function(m) apply (m, c (1, 2), function (x) { (as.integer(x)) }))

load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/makegenePWMbcgnd.rda")
nums = sample(1:length(geneSeq),1000,replace=F)

registerCoresPWMEnrich(4)
gene.BG = makeBackground(pos, bg.seq = geneSeq[nums], type = "logn", algorithm = "human", verbose=F)


## Calculate enrichment for all sequences together

res_frac = list()
for (i in 1:6){
  res_frac[[i]] = motifEnrichment(fracSeq[[i]], gene.BG, verbose=F)
  save(res_frac, file=paste0("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes",i,".rda"))
}
names(res_frac) = names(fracSeq)[1:6]

FracReport = lapply(res_frac, groupReport)

save(res_frac, FracReport, gene.BG, 
     file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.rda")

x = motifEnrichment(c(fracSeq[["both_retained"]],fracSeq[["Ad_retained"]]), gene.BG, verbose=F)

save(x, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.AllAd.Nuc.rda")

All.Prenatal.Nuc = motifEnrichment(c(fracSeq[["both_retained"]],fracSeq[["Fet_retained"]]), gene.BG, verbose=F)
save(All.Prenatal.Nuc, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.AllPren.Nuc.rda")

All.Adult.Cyt = motifEnrichment(c(fracSeq[["both_exported"]],fracSeq[["Ad_exported"]]), gene.BG, verbose=F)
save(All.Adult.Cyt, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.AllAd.Cyt.rda")

All.Prenatal.Cyt = motifEnrichment(c(fracSeq[["both_exported"]],fracSeq[["Fet_exported"]]), gene.BG, verbose=F)
save(All.Prenatal.Cyt, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.AllPren.Cyt.rda")

load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.rda")

res_frac = c(res_frac, list(All.Adult.Nuc = x, All.Prenatal.Nuc = All.Prenatal.Nuc, All.Adult.Cyt = All.Adult.Cyt, All.Prenatal.Cyt = All.Prenatal.Cyt))
FracReport = lapply(res_frac, groupReport)

save(res_frac, FracReport, gene.BG, 
     file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.rda")



### Explore results

#load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.rda")
load("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_fraction_genes.rda")

FracReport = lapply(FracReport, as.data.frame)

rbpdb = read.table("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT/ATtRACT_db.txt", header=T, sep = "\t")
#rbpdb = read.table("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_db.txt", header=T, sep = "\t")
FracReport = Map(cbind, FracReport, padj = lapply(FracReport, function(x) p.adjust(x$p.value, method = "fdr")),
                 lapply(FracReport, function(x) rbpdb[match(x$id,rbpdb$Matrix_id),]))
elementNROWS(FracReport)
elementNROWS(lapply(FracReport, function(x) x[x$padj<=0.01,]))
names(FracReport) = c("Nuclear:\nBoth", "Cytoplasmic:\nBoth", "Nuclear:\nPrenatal Only", 
                      "Nuclear:\nAdult Only","Cytoplasmic:\nPrenatal Only","Cytoplasmic:\nAdult Only",
                      "Nuclear\nin Adult", "Nuclear\nin Prenatal", "Cytoplasmic\nin Adult", "Cytoplasmic\nin Prenatal")


## filter to those expressed in brain

load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
#load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/rawCounts_combined_NucVSCyt_n23.rda")
rbps = geneMap[which(geneMap$ensemblID %in% unique(FracReport[[1]][,"Gene_id"])),"gencodeID"]
rbpCounts = geneCounts[which(rownames(geneCounts) %in% rbps),grep("polyA", colnames(geneCounts))]
dim(rbpCounts) # 158 RBPs measured

exprRBP = geneMap[which(geneMap$gencodeID %in% rownames(rbpCounts[which(apply(rbpCounts, 1, FUN=min)>5),])),"ensemblID"]
length(exprRBP) # 140 expressed

sigRBP = lapply(FracReport, function(x) x[which(x$padj<=0.01 & x$Gene_id %in% exprRBP),])
head(sigRBP[[1]])
elementNROWS(lapply(sigRBP, function(x) unique(x$Gene_id)))
# Nuclear:\nBoth          Cytoplasmic:\nBoth     Nuclear:\nPrenatal Only        Nuclear:\nAdult Only Cytoplasmic:\nPrenatal Only 
#             85                         114                          84                         116                          20 
# Cytoplasmic:\nAdult Only           Nuclear\nin Adult        Nuclear\nin Prenatal       Cytoplasmic\nin Adult    Cytoplasmic\nin Prenatal 
#                      123                         113                          86                         123                         115 

tbl = lapply(FracReport, function(x) x[which(x$Gene_id %in% exprRBP),which(colnames(x) %in% c("rank","id","raw.score","p.value","top.motif.prop",
                                                                 "padj","Gene_name","Gene_id","Motif","Len","Pubmed"))])
tbl = lapply(tbl, function(x) x[order(x$id),])
head(tbl[[1]])

df = data.frame("Motif ID" = tbl[[1]]$id, Motif = tbl[[1]]$Motif, "Motif Length" = tbl[[1]]$Len, "RBP Name" = tbl[[1]]$Gene_name, "Ensembl ID" = tbl[[1]]$Gene_id, 
                Pubmed = tbl[[1]]$Pubmed, "Cyt. In Adult Rank" = tbl[["Cytoplasmic\nin Adult"]]$rank, "Cyt. In Prenatal Rank" = tbl[["Cytoplasmic\nin Prenatal"]]$rank,
                "Nuc. In Adult Rank" = tbl[["Nuclear\nin Adult"]]$rank, "Nuc. In Prenatal Rank" = tbl[["Nuclear\nin Prenatal"]]$rank,
                "Cyt. In Adult Raw Score" = tbl[["Cytoplasmic\nin Adult"]]$raw.score, "Cyt. In Prenatal Raw Score" = tbl[["Cytoplasmic\nin Prenatal"]]$raw.score,
                "Nuc. In Adult Raw Score" = tbl[["Nuclear\nin Adult"]]$raw.score, "Nuc. In Prenatal Raw Score" = tbl[["Nuclear\nin Prenatal"]]$raw.score,
                "Cyt. In Adult Top Motif Prop." = tbl[["Cytoplasmic\nin Adult"]]$top.motif.prop, "Cyt. In Prenatal Top Motif Prop." = tbl[["Cytoplasmic\nin Prenatal"]]$top.motif.prop,
                "Nuc. In Adult Top Motif Prop." = tbl[["Nuclear\nin Adult"]]$top.motif.prop, "Nuc. In Prenatal Top Motif Prop." = tbl[["Nuclear\nin Prenatal"]]$top.motif.prop,
                "Cyt. In Adult P-value" = tbl[["Cytoplasmic\nin Adult"]]$p.value, "Cyt. In Prenatal P-value" = tbl[["Cytoplasmic\nin Prenatal"]]$p.value,
                "Nuc. In Adult P-value" = tbl[["Nuclear\nin Adult"]]$p.value, "Nuc. In Prenatal P-value" = tbl[["Nuclear\nin Prenatal"]]$p.value,
                "Cyt. In Adult FDR" = tbl[["Cytoplasmic\nin Adult"]]$padj, "Cyt. In Prenatal FDR" = tbl[["Cytoplasmic\nin Prenatal"]]$padj,
                "Nuc. In Adult FDR" = tbl[["Nuclear\nin Adult"]]$padj, "Nuc. In Prenatal FDR" = tbl[["Nuclear\nin Prenatal"]]$padj)
head(df)
write.csv(df, quote = F, file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/RBP_PWMEnrich_Results.csv")


## Venn diagrams of the overlap of RBPs in different groups of nuclear vs cytoplasmic genes

venn.diagram(lapply(sigRBP[grep("in", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.Alltogether.jpeg", 
             fill = c("#d95f02","#d95f02","#1b9e77","#1b9e77"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, 
             cat.cex=1.5, cat.dist=0.3, cat.fontfamily = "Arial", margin=.2)

venn.diagram(lapply(sigRBP[grep("Cytoplasmic:", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.cytoplasmic.jpeg", 
             fill = c("#1b9e77","#1b9e77","#1b9e77"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=.2)

venn.diagram(lapply(sigRBP[grep("Nuclear:", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.nuclear.jpeg", 
             fill = c("#d95f02","#d95f02", "#d95f02"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.35,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Both", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.both.jpeg", 
             fill = c("#d95f02","#1b9e77"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Adult Only", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.Adult.jpeg", 
             fill = c("#d95f02","#1b9e77"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Prenatal Only", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.Prenatal.jpeg", 
             fill = c("#d95f02","#1b9e77"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(list(Nuclear = do.call(rbind, sigRBP[grep("Nuclear:", names(sigRBP))]),
                         Cytoplasmic = do.call(rbind, sigRBP[grep("Cytoplasmic:", names(sigRBP))])), function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.all.jpeg", 
             fill = c("#d95f02","#1b9e77"), alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)



## Look up the uniquely enriched RBPs by fraction

head(sigRBP[[1]])
allRBPs = unique(as.character(FracReport[[1]][which(FracReport[[1]][,"Gene_id"] %in% exprRBP),"Gene_id"]))

frgenes = lapply(sigRBP, function(x) unique(as.character(x$Gene_id)))
names(frgenes)


rbps = list("Nuclear (In Adult)" = frgenes$"Nuclear\nin Adult"[-which(frgenes$"Nuclear\nin Adult" %in% frgenes$"Cytoplasmic\nin Adult")],
            "Nuclear (In Prenatal)" = frgenes$"Nuclear\nin Prenatal"[-which(frgenes$"Nuclear\nin Prenatal" %in% frgenes$"Cytoplasmic\nin Prenatal")],
            "Cytoplasmic (In Adult)" = frgenes$"Cytoplasmic\nin Adult"[-which(frgenes$"Cytoplasmic\nin Adult" %in% frgenes$"Nuclear\nin Adult")],
            "Cytoplasmic (In Prenatal)" = frgenes$"Cytoplasmic\nin Prenatal"[-which(frgenes$"Cytoplasmic\nin Prenatal" %in% frgenes$"Nuclear\nin Prenatal")])
rbps$"Nuclear (In Both)" = rbps$"Nuclear (In Adult)"[which(rbps$"Nuclear (In Adult)" %in% rbps$"Nuclear (In Prenatal)")]
rbps$"Nuclear (Adult-Specific)" = rbps$"Nuclear (In Adult)"[-which(rbps$"Nuclear (In Adult)" %in% rbps$"Nuclear (In Prenatal)")]
rbps$"Nuclear (Prenatal-Specific)" = rbps$"Nuclear (In Prenatal)"[-which(rbps$"Nuclear (In Prenatal)" %in% rbps$"Nuclear (In Adult)")]
rbps$"Cytoplasmic (In Both)" = rbps$"Cytoplasmic (In Adult)"[which(rbps$"Cytoplasmic (In Adult)" %in% rbps$"Cytoplasmic (In Prenatal)")]
rbps$"Cytoplasmic (Adult-Specific)" = rbps$"Cytoplasmic (In Adult)"[-which(rbps$"Cytoplasmic (In Adult)" %in% rbps$"Cytoplasmic (In Prenatal)")]
rbps$"Cytoplasmic (Prenatal-Specific)" = rbps$"Cytoplasmic (In Prenatal)"[-which(rbps$"Cytoplasmic (In Prenatal)" %in% rbps$"Cytoplasmic (In Adult)")]

elementNROWS(rbps)


FR = do.call(rbind, Map(cbind, FracReport, Group = as.list(names(FracReport))))
FR$gencodeID = geneMap[match(FR$Gene_id, geneMap$ensemblID),"gencodeID"]
head(FR)

exclFR = lapply(rbps, function(x) FR[which(FR$Gene_id %in% x),])
lapply(exclFR, head)

# unique RBPs to nuclear RNA

x = lapply(exclFR, function(x) unique(as.character(x$Gene_name))[order(unique(as.character(x$Gene_name)))])
x
#$`Nuclear (In Adult)`
# "CELF4"     "CELF5"     "HNRNPA1L2" "PTBP2"     "RBFOX1"    "RBFOX2"    "RBM24"     "RBM6"     

#$`Nuclear (In Prenatal)`
# "ANKHD1"  "CELF4"   "HNRNPA3" "HNRNPM"  "RBFOX2"  "RBM24"   "SAMD4A"  "SUPV3L1" "ZC3H10" 

#$`Nuclear (In Both)`
# "CELF4"  "RBFOX2" "RBM24" 

#$`Nuclear (Adult-Specific)`
# "CELF5"     "HNRNPA1L2" "PTBP2"     "RBFOX1"    "RBM6"     

#$`Nuclear (Prenatal-Specific)`
# "ANKHD1"  "HNRNPA3" "HNRNPM"  "SAMD4A"  "SUPV3L1" "ZC3H10" 

#$`Cytoplasmic (In Adult)`
# "AGO1"    "CPEB2"   "ELAVL2"  "ELAVL3"  "FXR1"    "HNRNPA0" "HNRNPAB" "MSI1"    "NUDT21"  "PABPC5"  "RBM28"   "RBM3"    "RBM41"   "RBM42"   "RBMS1"  
# "SNRNP70" "TUT1"    "ZCRB1"  

#$`Cytoplasmic (In Prenatal)`
# "AGO2"    "CMTR1"   "CPEB2"   "CPEB4"   "DDX19B"  "EIF4A3"  "ELAVL2"  "ELAVL3"  "ELAVL4"  "FXR1"    "HNRNPA0" "HNRNPD"  "IGHMBP2" "KHDRBS1" "KHDRBS2"
# "KHDRBS3" "MSI1"    "NUDT21"  "PABPC1"  "PABPC4"  "PABPC5"  "PPIE"    "RALY"    "RBM3"    "RBM41"   "RBM42"   "RBM6"    "RBMS1"   "RBMS3"   "RNASEL" 
# "SART3"   "SNRNP70" "SYNCRIP" "TIA1"    "TIAL1"   "TUT1"    "U2AF2"   "ZCRB1"  

#$`Cytoplasmic (In Both)`
# "CPEB2"   "ELAVL2"  "ELAVL3"  "FXR1"    "HNRNPA0" "MSI1"    "NUDT21"  "PABPC5"  "RBM3"    "RBM41"   "RBM42"   "RBMS1"   "SNRNP70" "TUT1"    "ZCRB1"  

#$`Cytoplasmic (Adult-Specific)`
# "AGO1"    "HNRNPAB" "RBM28"  

#$`Cytoplasmic (Prenatal-Specific)`
# "AGO2"    "CMTR1"   "CPEB4"   "DDX19B"  "EIF4A3"  "ELAVL4"  "HNRNPD"  "IGHMBP2" "KHDRBS1" "KHDRBS2" "KHDRBS3" "PABPC1"  "PABPC4"  "PPIE"    "RALY"   
# "RBM6"    "RBMS3"   "RNASEL"  "SART3"   "SYNCRIP" "TIA1"    "TIAL1"   "U2AF2"

x$`Nuclear (In Adult)`[which(x$`Nuclear (In Adult)` %in% sigRBP$"Cytoplasmic\nin Prenatal"$Gene_name)]
# "RBFOX1" "RBM6"

unlist(lapply(frgenes, function(y) y[which(y==geneMap[which(geneMap$Symbol=="AGO1"),"ensemblID"])]))
#Cytoplasmic\nin Adult 

unlist(lapply(frgenes, function(y) y[which(y==geneMap[which(geneMap$Symbol=="AGO2"),"ensemblID"])]))
# Nuclear\nin Adult    Cytoplasmic\nin Adult Cytoplasmic\nin Prenatal 


max(elementNROWS(lapply(exclFR, function(x) unique(as.character(x$Gene_name))[order(unique(as.character(x$Gene_name)))]))) # 38

l = lapply(exclFR, function(x) c(unique(as.character(x$Gene_name))[order(unique(as.character(x$Gene_name)))], 
                                 rep.int("NA", 38-length(unique(as.character(x$Gene_name))))))
write.csv(data.frame(do.call(cbind, l)), quote = F, file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/unique_RBPS_byFrac.csv")
write.csv(data.frame(unique(unlist(l))), quote = F, file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/unique_RBPS_byFrac_total.csv")


## What proportion of genes include the unique RBPs' motifs?

exclFR = lapply(exclFR, function(x) x[order(x$top.motif.prop, decreasing = T),
                                      which(colnames(x) %in% c("rank","id","raw.score","top.motif.prop","padj","Gene_name","Gene_id",
                                                                "Motif","Len","Matrix_id","Group","gencodeID"))]) 

N.inAd = data.table(exclFR$`Nuclear (In Adult)`[which(exclFR$`Nuclear (In Adult)`$Group %in% c("Nuclear\nin Adult","Cytoplasmic\nin Adult")),]
                    )[,list(Max.Prop = max(top.motif.prop)),by=c("Gene_name","Group")][order(Gene_name),,]
N.inAd.FC = data.frame(ID = N.inAd[grep("Nuc",N.inAd$Group),"Gene_name"], 
                       Fold.Diff = N.inAd[grep("Nuc",N.inAd$Group),"Max.Prop"]/N.inAd[grep("Cyt",N.inAd$Group),"Max.Prop"])
range(N.inAd[Group == "Nuclear\nin Adult","Max.Prop"])
# 0.01211827 0.03587009
mean(unlist(N.inAd[Group == "Nuclear\nin Adult","Max.Prop",]))
# 0.018238
N.inAd[which(N.inAd$Group=="Nuclear\nin Adult"),]
N.inAd.FC


C.inAd = data.table(exclFR$`Cytoplasmic (In Adult)`[which(exclFR$`Cytoplasmic (In Adult)`$Group %in% c("Nuclear\nin Adult","Cytoplasmic\nin Adult")),]
                    )[,list(Max.Prop = max(top.motif.prop)),by=c("Gene_name","Group")][order(Gene_name),,]
C.inAd.FC = data.frame(ID = C.inAd[grep("Nuc",C.inAd$Group),"Gene_name"], 
                       Fold.Diff = C.inAd[grep("Cyt",C.inAd$Group),"Max.Prop"]/C.inAd[grep("Nuc",C.inAd$Group),"Max.Prop"])
range(C.inAd[Group == "Cytoplasmic\nin Adult","Max.Prop"])
# 0.02753005 0.10314075
mean(unlist(C.inAd[Group == "Cytoplasmic\nin Adult","Max.Prop",]))
# 0.06729568
C.inAd[which(C.inAd$Group=="Cytoplasmic\nin Adult"),]
C.inAd.FC


N.inPren = data.table(exclFR$`Nuclear (In Prenatal)`[which(exclFR$`Nuclear (In Prenatal)`$Group %in% c("Nuclear\nin Prenatal","Cytoplasmic\nin Prenatal")),]
                      )[,list(Max.Prop = max(top.motif.prop)),by=c("Gene_name","Group")][order(Gene_name),,]
N.inPren.FC = data.frame(ID = N.inPren[grep("Nuc",N.inPren$Group),"Gene_name"], 
                         Fold.Diff = N.inPren[grep("Nuc",N.inPren$Group),"Max.Prop"]/N.inPren[grep("Cyt",N.inPren$Group),"Max.Prop"])
range(N.inPren[Group == "Nuclear\nin Prenatal","Max.Prop"])
# 0.0000000 0.3272059
mean(unlist(N.inPren[Group == "Nuclear\nin Prenatal","Max.Prop",]))
# 0.05718954
N.inPren[which(N.inPren$Group=="Nuclear\nin Prenatal"),]
N.inPren.FC


C.inPren = data.table(exclFR$`Cytoplasmic (In Prenatal)`[which(exclFR$`Cytoplasmic (In Prenatal)`$Group %in% c("Nuclear\nin Prenatal","Cytoplasmic\nin Prenatal")),]
                      )[,list(Max.Prop = max(top.motif.prop)),by=c("Gene_name","Group")][order(Gene_name),,]
C.inPren.FC = data.frame(ID = C.inPren[grep("Nuc",C.inPren$Group),"Gene_name"], 
                         Fold.Diff = C.inPren[grep("Nuc",C.inPren$Group),"Max.Prop"]/C.inPren[grep("Cyt",C.inPren$Group),"Max.Prop"])
range(C.inPren[Group == "Cytoplasmic\nin Prenatal","Max.Prop"])
# 0.04285714 0.27142857
mean(unlist(C.inPren[Group == "Cytoplasmic\nin Prenatal","Max.Prop",]))
# 0.1466165
C.inPren[which(C.inPren$Group=="Cytoplasmic\nin Prenatal"),]
C.inPren.FC


N.AdS = data.table(exclFR$`Nuclear (Adult-Specific)`[which(exclFR$`Nuclear (Adult-Specific)`$Group %in% c("Nuclear\nin Adult","Nuclear\nin Prenatal")),]
                   )[,list(Max.Prop = max(top.motif.prop)),by=c("Gene_name","Group")][order(Gene_name),,]
N.AdS.FC = data.frame(ID = N.AdS[grep("Adult", N.AdS$Group),"Gene_name"], 
                      Fold.Diff = N.AdS[grep("Adult", N.AdS$Group),"Max.Prop"]/N.AdS[grep("Prenatal", N.AdS$Group),"Max.Prop"])
range(N.AdS[Group == "Nuclear\nin Adult","Max.Prop"])
# 0.01211827 0.01938924
mean(unlist(N.AdS[Group == "Nuclear\nin Adult","Max.Prop",]))
# 0.01560834
N.AdS[which(N.AdS$Group=="Nuclear\nin Adult"),]
N.AdS.FC


N.PrenS = data.table(exclFR$`Nuclear (Prenatal-Specific)`[which(exclFR$`Nuclear (Prenatal-Specific)`$Group %in% c("Nuclear\nin Adult","Nuclear\nin Prenatal")),]
                     )[,list(Max.Prop = max(top.motif.prop)),by=c("Gene_name","Group")][order(Gene_name),,]
N.PrenS.FC = data.frame(ID = N.PrenS[grep("Adult", N.PrenS$Group),"Gene_name"], 
                        Fold.Diff = N.PrenS[grep("Prenatal", N.PrenS$Group),"Max.Prop"]/N.PrenS[grep("Adult", N.PrenS$Group),"Max.Prop"])
range(N.PrenS[Group == "Nuclear\nin Prenatal","Max.Prop"])
# 0.0000000 0.3272059
mean(unlist(N.PrenS[Group == "Nuclear\nin Prenatal","Max.Prop",]))
# 0.07414216
N.PrenS[which(N.PrenS$Group=="Nuclear\nin Prenatal"),]
N.PrenS.FC


## How about the difference statistic?

load("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_fraction_genes.rda")

vs = list(NvC.A = motifDiffEnrichment(res1 = res_frac[["All.Adult.Nuc"]], res2 = res_frac[["All.Adult.Cyt"]], pwms = gene.BG),
          NvC.P = motifDiffEnrichment(res1 = res_frac[["All.Prenatal.Nuc"]], res2 = res_frac[["All.Prenatal.Cyt"]], pwms = gene.BG),
          AvP.Nuc = motifDiffEnrichment(res1 = res_frac[["All.Adult.Nuc"]], res2 = res_frac[["All.Prenatal.Nuc"]], pwms = gene.BG),
          AvP.Cyt = motifDiffEnrichment(res1 = res_frac[["All.Adult.Cyt"]], res2 = res_frac[["All.Prenatal.Cyt"]], pwms = gene.BG))

pdf("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/change_in_RBP_Difference_Statistic.pdf",w=5.5, h=3.5)
ggplot(data.frame(vs$NvC.A), aes(group.bg)) + geom_histogram() + xlab("(Nuclear - Cytoplasmic)") + ylab("Count") +
  theme(title = element_text(size = 20), text = element_text(size = 20)) + ggtitle("Difference Statistic\nBy Fraction in Adult")
ggplot(data.frame(vs$NvC.P), aes(group.bg)) + geom_histogram() + xlab("(Nuclear - Cytoplasmic)") + ylab("Count") +
  theme(title = element_text(size = 20), text = element_text(size = 20)) + ggtitle("Difference Statistic\nBy Fraction in Prenatal")
ggplot(data.frame(vs$AvP.Nuc), aes(group.bg)) + geom_histogram() + xlab("(Adult - Prenatal)") + ylab("Count") +
  theme(title = element_text(size = 20), text = element_text(size = 20)) + ggtitle("Difference Statistic\nBy Age in Nuclear RNA")
ggplot(data.frame(vs$AvP.Cyt), aes(group.bg)) + geom_histogram() + xlab("(Adult - Prenatal)") + ylab("Count") +
  theme(title = element_text(size = 20), text = element_text(size = 20)) + ggtitle("Difference Statistic\nBy Age in Cytoplasmic RNA")
dev.off()

do.call(rbind,lapply(vs, function(x) data.frame(lower = range(x$group.bg[is.finite(x$group.bg)])[1], 
                                                higher = range(x$group.bg[is.finite(x$group.bg)])[2])))
#                lower       higher
#NvC.A   -4.364430e+16 3.041708e+16
#NvC.P   -4.723720e+05 2.878363e+16
#AvP.Nuc -1.467873e+03 3.041708e+16
#AvP.Cyt -1.071220e+05 5.083724e+16


## 1 standard deviation difference statistic (no background)

vs.1sd = lapply(vs, function(x) x$group.nobg[which(abs(x$group.nobg)>=(mean(x$group.nobg)+sd(x$group.nobg[is.finite(x$group.nobg)])))])
elementNROWS(vs.1sd)
vs.1sd = lapply(vs.1sd, function(x) cbind(rbpdb[match(names(x), rbpdb$Matrix_id),], names = names(x), diff = x))

lapply(vs.1sd, function(x) unique(as.character(x[which(x$diff<0),"Gene_name"])))
#$NvC.A
# "TARDBP" "CELF1" 
#$NvC.P
# "TARDBP" "CELF1" 
#$AvP.Nuc
# "HNRNPH2" "HNRNPH1" "SRSF6"  
#$AvP.Cyt
# "TARDBP" "CELF1"  "NOVA1"

lapply(vs.1sd, function(x) unique(as.character(x[which(x$diff>0),"Gene_name"])))
#$NvC.A
# "HNRNPH1"   "HNRNPA1"   "HNRNPA2B1" "HNRNPH2"   "HNRNPH3"   "SRSF6"     "HNRNPF"    "SRSF3"     "YBX1"     
#$NvC.P
# "HNRNPH1" "HNRNPH2" "HNRNPH3" "HNRNPF" 
#$AvP.Nuc
# "TARDBP" "CELF1" 
#$AvP.Cyt
# "ELAVL1" "TIA1"   "AGO2"   "ELAVL4"
