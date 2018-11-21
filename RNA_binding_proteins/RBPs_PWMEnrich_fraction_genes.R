library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PWMEnrich)
library("venn")
library("VennDiagram")


load("./Dropbox/sorted_figures/github_controlled/RNA_localization_and_age/data/retained.byAge.downsampled.rda")
load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/retained.byAge.downsampled.rda")


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
dim(rbpdb) # 3256 RBPs included
ids = as.character(unique(rbpdb$Gene_name))


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


## Calculate for disease gene sets

aej_sets = openxlsx::read.xlsx('./Dropbox/sorted_figures/github_controlled/Birnbaum_2013_AJP_Supplementary_table.xlsx')

AEJmap = geneMap[which(geneMap$Symbol %in% as.character(aej_sets$Gene.Symbol)),]
AEJmap = cbind(AEJmap, Gene.Set = aej_sets[match(AEJmap$Symbol, aej_sets$Gene.Symbol),"Gene.Set"])
AEJmap = split(AEJmap, AEJmap$Gene.Set)
elementNROWS(AEJmap)
AEJmap = lapply(AEJmap, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

DiseaseSeq = lapply(AEJmap, function(x) getSeq(Hsapiens, x))
DiseaseSeq = lapply(DiseaseSeq, reverseComplement)
save(DiseaseSeq, file = "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/DiseaseSeq.rda")

res_Disease = list()
for (i in 1:length(DiseaseSeq)){
  res_Disease[[i]] = motifEnrichment(DiseaseSeq[[i]], gene.BG, verbose=F)
  save(res_Disease, file=paste0("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_Disease_genes",i,".rda"))
}
names(res_Disease) = names(DiseaseSeq)

DiseaseReport = lapply(res_Disease, groupReport)

save(res_Disease, DiseaseReport, gene.BG, 
     file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_Disease_genes.rda")


## Explore results

load("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_fraction_genes.rda")

FracReport = lapply(FracReport, as.data.frame)

rbpdb = read.table("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT/ATtRACT_db.txt", header=T, sep = "\t")
FracReport = Map(cbind, FracReport, padj = lapply(FracReport, function(x) p.adjust(x$p.value, method = "fdr")),
                 lapply(FracReport, function(x) rbpdb[match(x$id,rbpdb$Matrix_id),]))

elementNROWS(FracReport)
elementNROWS(lapply(FracReport, function(x) x[x$padj<=0.001,]))
names(FracReport) = c("Nuclear:\nBoth", "Cytoplasmic:\nBoth", "Nuclear:\nPrenatal Only", 
                      "Nuclear:\nAdult Only","Cytoplasmic:\nPrenatal Only","Cytoplasmic:\nAdult Only")
sigRBP = lapply(FracReport, function(x) x[x$padj<=0.001,])
head(sigRBP[[1]])
elementNROWS(lapply(sigRBP, function(x) unique(x$Gene_id)))
venn(lapply(sigRBP, function(x) unique(x$Gene_id)))


venn.diagram(lapply(sigRBP[grep("Cytoplasmic", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.cytoplasmic.jpeg", 
             main="RBPs Enriched in Cytoplasmic Genes\n(LFC≥1, FDR<0.05)",
             fill = c("#1b9e77","#1b9e77","#1b9e77"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.2,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Nuclear", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.nuclear.jpeg", 
             main="RBPs Enriched in Nuclear Genes\n(LFC≥1, FDR<0.05)",
             fill = c("#d95f02","#d95f02", "#d95f02"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Both", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.both.jpeg", 
             main="RBPs Enriched in Nuclear and Cytoplasmic Genes\n(LFC≥1, FDR<0.05)",
             fill = c("#d95f02","#1b9e77"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Adult", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.Adult.jpeg", 
             main="RBPs Enriched in Nuclear and Cytoplasmic Genes\n(LFC≥1, FDR<0.05)",
             fill = c("#d95f02","#1b9e77"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigRBP[grep("Prenatal", names(sigRBP))], function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.Prenatal.jpeg", 
             main="RBPs Enriched in Nuclear and Cytoplasmic Genes\n(LFC≥1, FDR<0.05)",
             fill = c("#d95f02","#1b9e77"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(list(Cytoplasmic = do.call(rbind, sigRBP[grep("Cytoplasmic", names(sigRBP))]),
                         Nuclear = do.call(rbind, sigRBP[grep("Nuclear", names(sigRBP))])), 
                    function(x) unique(x$Gene_id)), 
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.byFrac.all.jpeg", 
             main="RBPs Enriched in Nuclear and Cytoplasmic Genes\n(LFC≥1, FDR<0.05)",
             fill = c("#d95f02","#1b9e77"),
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)


## Enrichment of Nuclear or Cytoplasmic RBP motifs in disease gene sets?

fisher.test()















