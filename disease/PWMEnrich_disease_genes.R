library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PWMEnrich)
library("VennDiagram")


load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_fraction_genes.rda")


## Calculate RBP motif enrichment for disease gene sets

registerCoresPWMEnrich(5)

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


pgc = read.delim("./Dropbox/sorted_figures/github_controlled/disease/data/tableS3_Nature2014_pgc2_loci.txt",as.is=TRUE)
pgc$chr = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[1]))
pgc$range = unlist(lapply(strsplit(pgc$Position..hg19.,":"), function(x) x[2]))
pgc$start = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[1])))
pgc$end = as.numeric(unlist(lapply(strsplit(pgc$range,"-"), function(x) x[2])))
pgcGR = GRanges(pgc$chr, IRanges(pgc$start,pgc$end))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep.extra.columns = T)
pgc2 = geneMap[queryHits(findOverlaps(geneMapGR, pgcGR)),]
pgc2 = makeGRangesFromDataFrame(pgc2, keep.extra.columns = T)

pgc2Seq = getSeq(Hsapiens, pgc2)
pgc2Seq = reverseComplement(pgc2Seq)
save(pgc2Seq, file="./Dropbox/sorted_figures/github_controlled/disease/data/pgc2Seq.rda")

load("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_Disease_genes.rda")

res_pgc2 = motifEnrichment(pgc2Seq, gene.BG, verbose=F)

res_Disease = c(res_Disease, list(PGC2 = res_pgc2))
DiseaseReport = c(DiseaseReport, list(PGC2 = groupReport(res_pgc2)))

save(res_Disease, DiseaseReport, gene.BG, 
     file="./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_Disease_genes.rda")


load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/DiseaseSeq.rda")
load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_Disease_genes.rda")
load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/pgc2Seq.rda")

Disease.Nuc.Both = motifEnrichment(c(DiseaseSeq[["ASD CNV"]],DiseaseSeq[["ASD DATABASE"]], DiseaseSeq[["SCZ CNV"]]), gene.BG, verbose=F)
save(Disease.Nuc.Both, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_disease_genes.Nuc.Both.rda")

Disease.Nuc.Ad = motifEnrichment(c(DiseaseSeq[["BPAD GWAS"]],DiseaseSeq[["SCZ SNV"]], pgc2Seq), gene.BG, verbose=F)
save(Disease.Nuc.Ad, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_disease_genes.Nuc.Ad.rda")

Disease.Nuc = motifEnrichment(c(DiseaseSeq[["ASD CNV"]],DiseaseSeq[["ASD DATABASE"]], DiseaseSeq[["SCZ CNV"]],
                                DiseaseSeq[["BPAD GWAS"]],DiseaseSeq[["SCZ SNV"]], pgc2Seq)
                              [unique(names(c(DiseaseSeq[["ASD CNV"]],DiseaseSeq[["ASD DATABASE"]], DiseaseSeq[["SCZ CNV"]],
                                              DiseaseSeq[["BPAD GWAS"]],DiseaseSeq[["SCZ SNV"]], pgc2Seq)))], gene.BG, verbose=F)
save(Disease.Nuc, file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_disease_genes.Nuc.rda")

res_Disease = c(res_Disease, list(Disease.Nuc = Disease.Nuc, Disease.Nuc.Both = Disease.Nuc.Both, Disease.Nuc.Ad = Disease.Nuc.Ad))


Disease.notNucAll = motifEnrichment(c(DiseaseSeq[["ID"]],DiseaseSeq[["NDD"]], DiseaseSeq[["Neurodegenerative"]],
                                DiseaseSeq[["SCZ Meta-analysis"]])
                                [unique(names(c(DiseaseSeq[["ID"]],DiseaseSeq[["NDD"]], DiseaseSeq[["Neurodegenerative"]],
                                              DiseaseSeq[["SCZ Meta-analysis"]])))], gene.BG, verbose=F)

Disease.notNucAd = motifEnrichment(c(DiseaseSeq[["ID"]],DiseaseSeq[["NDD"]], DiseaseSeq[["Neurodegenerative"]],
                                      DiseaseSeq[["SCZ Meta-analysis"]],DiseaseSeq[["ASD CNV"]],DiseaseSeq[["ASD DATABASE"]], DiseaseSeq[["SCZ CNV"]])
                                    [unique(names(c(DiseaseSeq[["ID"]],DiseaseSeq[["NDD"]], DiseaseSeq[["Neurodegenerative"]],
                                                    DiseaseSeq[["SCZ Meta-analysis"]],DiseaseSeq[["ASD CNV"]],DiseaseSeq[["ASD DATABASE"]], DiseaseSeq[["SCZ CNV"]])))],
                                   gene.BG, verbose=F)

Disease.notNucBoth = motifEnrichment(c(DiseaseSeq[["ID"]],DiseaseSeq[["NDD"]], DiseaseSeq[["Neurodegenerative"]],
                                      DiseaseSeq[["SCZ Meta-analysis"]], DiseaseSeq[["BPAD GWAS"]],DiseaseSeq[["SCZ SNV"]], pgc2Seq)
                                    [unique(names(c(DiseaseSeq[["ID"]],DiseaseSeq[["NDD"]], DiseaseSeq[["Neurodegenerative"]],
                                                    DiseaseSeq[["SCZ Meta-analysis"]],DiseaseSeq[["BPAD GWAS"]],DiseaseSeq[["SCZ SNV"]], pgc2Seq)))],
                                    gene.BG, verbose=F)


res_Disease = c(res_Disease, list(Disease.notNucAll = Disease.notNucAll, Disease.notNucAd = Disease.notNucAd, Disease.notNucBoth = Disease.notNucBoth))

DiseaseReport = lapply(res_Disease, groupReport)

save(res_Disease, DiseaseReport, gene.BG, 
     file="/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/ATtRACT_Disease_genes.rda")


## Explore disease results:

load("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_Disease_genes.rda")

DiseaseReport = lapply(DiseaseReport, as.data.frame)
rbpdb = read.table("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT/ATtRACT_db.txt", header=T, sep = "\t")
DiseaseReport = Map(cbind, DiseaseReport, padj = lapply(DiseaseReport, function(x) p.adjust(x$p.value, method = "fdr")),
                    lapply(DiseaseReport, function(x) rbpdb[match(x$id,rbpdb$Matrix_id),]))

elementNROWS(DiseaseReport)
elementNROWS(lapply(DiseaseReport, function(x) x[x$padj<=0.01,]))
names(DiseaseReport) = c("ASD\n(CNV)","ASD\n(Database)","BPAD\n(GWAS)","Intellectual\nDisability","Neuro-\ndevel.","Neuro-\ndegen.","SCZ\n(CNV)",
                         "SCZ\n(Meta\nanalysis)","SCZ\n(PGC1)","SCZ\n(SNV)","SCZ\n(PGC2)", "All Nuclear-\nEnriched","Nuclear-Enriched\nin Both","Nuclear-Enriched\nin Adult",
                         "Not Nuclear-\nEnriched","Not Nuclear-Enriched\nin Adult","Not Nuclear-Enriched\nin Both")


## filter to those expressed in brain

load("./Dropbox/sorted_figures/github_controlled/QC_section/data/rawCounts_combined_NucVSCyt_n23.rda")
#load("/dcl01/lieber/ajaffe/Amanda/BrainRNACompartments/rdas/rawCounts_combined_NucVSCyt_n23.rda")

rbps = geneMap[which(geneMap$ensemblID %in% unique(DiseaseReport[[1]][,"Gene_id"])),"gencodeID"]
rbpCounts = geneCounts[which(rownames(geneCounts) %in% rbps),grep("polyA", colnames(geneCounts))]
dim(rbpCounts) # 158 RBPs measured

exprRBP = geneMap[which(geneMap$gencodeID %in% rownames(rbpCounts[which(apply(rbpCounts, 1, FUN=min)>5),])),"ensemblID"]
length(exprRBP) # 140 expressed


sigdRBP = lapply(DiseaseReport[which(names(DiseaseReport)!="SCZ\n(PGC1)")], function(x) x[which(x$padj<=0.01 & x$Gene_id %in% exprRBP),])
head(sigdRBP[[1]])
elementNROWS(lapply(sigdRBP, function(x) unique(x$Gene_id)))
#                    ASD\n(CNV)                ASD\n(Database)                   BPAD\n(GWAS)       Intellectual\nDisability                 Neuro-\ndevel. 
#110                            132                            114                            122                            119 
#Neuro-\ndegen.                     SCZ\n(CNV)          SCZ\n(Meta\nanalysis)                     SCZ\n(SNV)                    SCZ\n(PGC2) 
#110                            106                            100                            128                            123 
#All Nuclear-\nEnriched      Nuclear-Enriched\nin Both     Nuclear-Enriched\nin Adult         Not Nuclear-\nEnriched Not Nuclear-Enriched\nin Adult 
#128                            128                            124                            129                            130 
#Not Nuclear-Enriched\nin Both 
#128 


for (i in 1:length(sigdRBP)) {
  venn.diagram(lapply(c(sigdRBP[which(names(sigdRBP)==names(sigdRBP)[i])],
                        list(Others = do.call(rbind, sigdRBP[which(names(sigdRBP)!=names(sigdRBP)[i])]))), function(x) unique(x$Gene_id)), 
               paste0("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.disease", i, ".vs.AllOthers.jpeg"), 
               main=paste0("RBPs Enriched in ",names(sigdRBP)[i], "\n(FDR<0.05)"),
               alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
               cat.fontfamily = "Arial", margin=0.2)
}

venn.diagram(lapply(sigdRBP[which(names(sigdRBP) %in% c("All Nuclear-\nEnriched","Nuclear-Enriched\nin Both","Nuclear-Enriched\nin Adult"))], 
                    function(x) unique(x$Gene_id)),
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.disease_combined.jpeg",
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)


venn.diagram(lapply(sigdRBP[which(names(sigdRBP) %in% c("Nuclear-Enriched\nin Both", "Not Nuclear-Enriched\nin Both"))], function(x) unique(x$Gene_id)),
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.disease_nucBoth.jpeg",
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigdRBP[which(names(sigdRBP) %in% c("Nuclear-Enriched\nin Adult", "Not Nuclear-Enriched\nin Adult"))], function(x) unique(x$Gene_id)),
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.disease_nucAdult.jpeg",
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)

venn.diagram(lapply(sigdRBP[which(names(sigdRBP) %in% c("All Nuclear-\nEnriched", "Not Nuclear-\nEnriched"))], function(x) unique(x$Gene_id)),
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.disease_nucAll.jpeg",
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.15,
             cat.fontfamily = "Arial", margin=0.2)
                                                        
venn.diagram(lapply(sigdRBP[which(names(sigdRBP) %in% c("Neuro-\ndegen.","Intellectual\nDisability","Neuro-\ndevel.","Nuclear-Enriched\nin Both"))], 
                    function(x) unique(x$Gene_id)),
             "./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/figures/venn.RBPs.disease_nuclear_vs_not.jpeg",
             alpha = 0.50, fontfamily = "Arial", fontface = "bold", cex=1.5, cat.cex=1.5, cat.dist=0.25,
             cat.fontfamily = "Arial", margin=0.2)



dgenes = lapply(sigdRBP, function(x) unique(as.character(x$Gene_id)))
names(dgenes)
head(dgenes[[1]])

Dis = list(NucOnly = dgenes$`All Nuclear-\nEnriched`[-which(dgenes$`All Nuclear-\nEnriched` %in% dgenes$`Not Nuclear-\nEnriched`)],
           notNucOnly = dgenes$`Not Nuclear-\nEnriched`[-which(dgenes$`Not Nuclear-\nEnriched` %in% dgenes$`All Nuclear-\nEnriched`)],
           Overlap = dgenes$`All Nuclear-\nEnriched`[which(dgenes$`All Nuclear-\nEnriched` %in% dgenes$`Not Nuclear-\nEnriched`)])
elementNROWS(Dis)

lapply(Dis[1:2], function(x) geneMap[which(geneMap$ensemblID %in% x),"Symbol"])
#$NucOnly
# "ZNF638" "RBM6"   "G3BP2" 

#$notNucOnly
# "HNRNPAB"   "RBM28"     "HNRNPA1L2" "YBX2"

# genes NOT enriched in adult nuclear enriched disease sets but in others
geneMap[which(geneMap$ensemblID %in% dgenes$`Not Nuclear-Enriched\nin Adult`[!(dgenes$`Not Nuclear-Enriched\nin Adult` %in% dgenes$`Nuclear-Enriched\nin Adult`)]),"Symbol"]
# "PTBP2"     "HNRNPAB"   "HNRNPA1L2" "CELF4"     "CELF5"     "RBM42"

# Gene enriched in "nuclear in both" disease genes but not others
geneMap[which(geneMap$ensemblID %in% dgenes$`Nuclear-Enriched\nin Both`[!(dgenes$`Nuclear-Enriched\nin Both` %in% dgenes$`Not Nuclear-Enriched\nin Both`)]),"Symbol"]
# "CELF5"

# Gene Not enriched in "nuclear in both" disease genes but others
geneMap[which(geneMap$ensemblID %in% dgenes$`Not Nuclear-Enriched\nin Both`[!(dgenes$`Not Nuclear-Enriched\nin Both` %in% dgenes$`Nuclear-Enriched\nin Both`)]),"Symbol"]
# "RBM28"


DisExcl = lapply(DiseaseReport[1:11], function(x) x[which(x$Gene_id %in% Dis$NucOnly),])
DisExcl = do.call(rbind, Map(cbind, DisExcl, GeneSet = as.list(names(DisExcl))))
rownames(DisExcl) = NULL
DisExcl[which(DisExcl$padj<=0.05),colnames(DisExcl) %in% c("top.motif.prop","padj","Gene_name","GeneSet")]
#top.motif.prop         padj Gene_name                  GeneSet
#1      0.03448276 4.497498e-07    ZNF638               ASD\n(CNV)
#2      0.03448276 7.088524e-05     G3BP2               ASD\n(CNV)
#4      0.01421801 2.284680e-08      RBM6          ASD\n(Database)
#5      0.02843602 6.389431e-07    ZNF638          ASD\n(Database)
#6      0.01895735 3.472428e-06     G3BP2          ASD\n(Database)
#7      0.04761905 1.798971e-03     G3BP2             BPAD\n(GWAS)
#10     0.03846154 3.506470e-02     G3BP2 Intellectual\nDisability
#13     0.03333333 1.074522e-03    ZNF638           Neuro-\ndevel.
#16     0.04878049 4.570470e-02    ZNF638           Neuro-\ndegen.
#19     0.02020202 1.048663e-02    ZNF638               SCZ\n(CNV)
#20     0.01010101 1.114916e-02     G3BP2               SCZ\n(CNV)
#22     0.00000000 1.695824e-02     G3BP2    SCZ\n(Meta\nanalysis)
#25     0.03000000 1.677077e-05     G3BP2              SCZ\n(PGC1)
#28     0.05102041 1.170330e-16     G3BP2               SCZ\n(SNV)
#29     0.03061224 1.507827e-07    ZNF638               SCZ\n(SNV)
#30     0.04081633 4.903027e-06      RBM6               SCZ\n(SNV)
#31     0.06771654 6.116565e-12     G3BP2              SCZ\n(PGC2)
#32     0.04724409 8.403527e-05      RBM6              SCZ\n(PGC2)
#33     0.04409449 4.920842e-04    ZNF638              SCZ\n(PGC2)


## Enrichment of Nuclear or Cytoplasmic RBP motifs in disease gene sets?

load("./Dropbox/sorted_figures/github_controlled/RNA_binding_proteins/data/ATtRACT_fraction_genes.rda")

FracReport = lapply(FracReport, as.data.frame)
FracReport = Map(cbind, FracReport, padj = lapply(FracReport, function(x) p.adjust(x$p.value, method = "fdr")),
                 lapply(FracReport, function(x) rbpdb[match(x$id,rbpdb$Matrix_id),]))
elementNROWS(FracReport)
names(FracReport) = c("Nuclear:\nBoth", "Cytoplasmic:\nBoth", "Nuclear:\nPrenatal Only", 
                      "Nuclear:\nAdult Only","Cytoplasmic:\nPrenatal Only","Cytoplasmic:\nAdult Only",
                      "Nuclear\nin Adult", "Nuclear\nin Prenatal", "Cytoplasmic\nin Adult", "Cytoplasmic\nin Prenatal")

sigRBP = lapply(FracReport, function(x) x[which(x$padj<=0.01 & x$Gene_id %in% exprRBP),])
elementNROWS(lapply(sigRBP, function(x) unique(x$Gene_id)))


allRBPs = unique(as.character(FracReport[[1]][which(FracReport[[1]][,"Gene_id"] %in% exprRBP),"Gene_id"]))
frgenes = lapply(sigRBP, function(x) unique(as.character(x$Gene_id)))

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

ct = list(adult = data.frame(DisNuc = c(length(Dis$NucOnly[which(Dis$NucOnly %in% rbps$`Nuclear (In Adult)`)]), 
                                        length(Dis$NucOnly[-which(Dis$NucOnly %in% rbps$`Nuclear (In Adult)`)])),
                             DisNotNuc = c(length(allRBPs[-which(allRBPs %in% Dis$NucOnly)][which(allRBPs[-which(allRBPs %in% Dis$NucOnly)] %in% rbps$`Nuclear (In Adult)`)]), 
                                           length(allRBPs[-which(allRBPs %in% rbps$`Nuclear (In Adult)`)][-which(allRBPs[-which(allRBPs %in% rbps$`Nuclear (In Adult)`)] %in% Dis$NucOnly)])),
                             row.names = c("FracNuc", "FracNotNuc")),
          prenatal = data.frame(DisNuc = c(length(Dis$NucOnly[which(Dis$NucOnly %in% rbps$`Nuclear (In Prenatal)`)]),
                                           length(Dis$NucOnly[!(Dis$NucOnly %in% rbps$`Nuclear (In Prenatal)`)])),
                                DisNotNuc = c(length(allRBPs[-which(allRBPs %in% Dis$NucOnly)][which(allRBPs[-which(allRBPs %in% Dis$NucOnly)] %in% rbps$`Nuclear (In Prenatal)`)]),
                                              length(allRBPs[-which(allRBPs %in% rbps$`Nuclear (In Prenatal)`)][-which(allRBPs[-which(allRBPs %in% rbps$`Nuclear (In Prenatal)`)] %in% Dis$NucOnly)])),
                                row.names = c("FracNuc", "FracNotNuc")),
          both = data.frame(DisNuc = c(length(Dis$NucOnly[which(Dis$NucOnly %in% rbps$`Nuclear (In Both)`)]),
                                       length(Dis$NucOnly[!(Dis$NucOnly %in% rbps$`Nuclear (In Both)`)])),
                            DisNotNuc = c(length(allRBPs[-which(allRBPs %in% Dis$NucOnly)][which(allRBPs[-which(allRBPs %in% Dis$NucOnly)] %in% rbps$`Nuclear (In Both)`)]),
                                          length(allRBPs[-which(allRBPs %in% rbps$`Nuclear (In Both)`)][-which(allRBPs[-which(allRBPs %in% rbps$`Nuclear (In Both)`)] %in% Dis$NucOnly)])),
                            row.names = c("FracNuc", "FracNotNuc")))

res = lapply(ct, fisher.test)
res = do.call(rbind, Map(cbind, Comparison = names(res), lapply(res, function(x) data.frame(OR = x$estimate, pval = x$p.value))))
res$FDR = p.adjust(res$pval, method = "fdr")
row.names(res) = NULL
write.csv(res, file = "./Dropbox/sorted_figures/github_controlled/disease/figures/disease_nuclearDEGs_RBP_overlap_fishertest.csv")
res
#  Comparison       OR      pval       FDR
#1      adult 8.932413 0.1629206 0.4887618
#2   prenatal 0.000000 1.0000000 1.0000000
#3       both 0.000000 1.0000000 1.0000000


## The one overlap is RBM6, the prevalence of the motif in disease sets and in fraction differences?

# disease genes
do.call(rbind, Map(cbind, lapply(sigdRBP, function(x) x[which(x$Gene_name=="RBM6"),colnames(x) %in% c("top.motif.prop","padj","Gene_name")])
                   [elementNROWS(lapply(sigdRBP, function(x) x[which(x$Gene_name=="RBM6"),]))>0]))
#                               top.motif.prop         padj Gene_name
#ASD\n(Database)                    0.01421801 2.284680e-08      RBM6
#SCZ\n(SNV)                         0.04081633 4.903027e-06      RBM6
#SCZ\n(PGC2)                        0.04724409 8.403527e-05      RBM6
#All Nuclear-\nEnriched             0.03491078 7.734558e-13      RBM6
#Nuclear-Enriched\nin Both          0.01511335 4.095709e-05      RBM6
#Nuclear-Enriched\nin Adult         0.04262295 2.585492e-09      RBM6
#Not Nuclear-Enriched\nin Adult     0.02245250 8.983884e-06      RBM6
#Not Nuclear-Enriched\nin Both      0.04231831 2.631014e-09      RBM6

# fraction genes
do.call(rbind, Map(cbind, lapply(sigRBP, function(x) x[which(x$Gene_name=="RBM6"),colnames(x) %in% c("top.motif.prop","padj","Gene_name")])
                   [elementNROWS(lapply(sigRBP, function(x) x[which(x$Gene_name=="RBM6"),]))>0])) 
#                         top.motif.prop         padj Gene_name
#Cytoplasmic:\nBoth           0.04918033 9.318325e-03      RBM6
#Nuclear:\nAdult Only         0.01305767 1.182742e-08      RBM6
#Nuclear\nin Adult            0.01211827 1.370404e-05      RBM6
#Cytoplasmic\nin Prenatal     0.05714286 6.666902e-03      RBM6